#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from hashlib import sha256
from typing import Any
import re
import pandas as pd
from Bio import SeqIO

"""
ARNm complet fusionné
=
5′UTR_A + CDS_A_partiel + CDS_B/UTR_B selon breakpoint + 3′UTR_B
"""

# =========================
# INTERFACE SNAKEMAKE
# =========================

# Entrées déclarées 
ARRIBA_TSV = Path(snakemake.input.arriba_fusion)
STAR_FUSION_TSV = Path(snakemake.input.star_fusion)  
GTF_FILE = Path(snakemake.input.gtf)
TRANSCRIPTOME_FASTA = Path(snakemake.input.transcriptome_fasta)
KALLISTO_ABUNDANCE = Path(snakemake.input.kallisto_abundance)

# Sorties déclarées 
OUT_ARRIBA_FASTA = Path(snakemake.output.arriba_transcripts)
OUT_STAR_FASTA = Path(snakemake.output.star_fusion_transcripts)
OUT_COMBINED_FASTA = Path(snakemake.output.fusions_transcripts)

OUT_ARRIBA_SUMMARY = Path(snakemake.output.arriba_summary)
OUT_STAR_SUMMARY = Path(snakemake.output.star_summary)
OUT_COMBINED_SUMMARY = Path(snakemake.output.summary)

OUT_STAR_CANDIDATES = Path(snakemake.output.star_candidates)

# Paramètres 
MIN_TPM = float(snakemake.params.get("min_tpm", 0.0))
ALLOW_INTRONIC_SNAP = bool(snakemake.params.get("allow_intronic_snap", True))
MAX_SNAP_DISTANCE = snakemake.params.get("max_snap_distance", None)


# =========================
# UTILITAIRES
# =========================

def strip_version(x):
    if x is None:
        return None
    x = str(x)
    if x == "." or x == "":
        return "."
    return x.split(".")[0]


def normalize_chrom(chrom):
    chrom = str(chrom)
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    return chrom


def parse_gtf_attributes(attr_text):
    attrs = {}
    for field in attr_text.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        parts = field.split(" ", 1)
        if len(parts) != 2:
            continue
        key, value = parts
        attrs[key] = value.strip().strip('"')
    return attrs


def parse_breakpoint(bp):
    """
    Exemple Arriba:
    12:50474167
    chr12:50474167
    """
    bp = str(bp)
    chrom, pos = bp.split(":")
    return normalize_chrom(chrom), int(pos)


def clean_gene_symbol(gene):
    """
    Arriba peut écrire:
    COL5A2(62657),KRT18P19(23294)
    On retourne une liste de symboles possibles.
    """
    if gene is None or gene == ".":
        return []

    genes = []
    for part in str(gene).split(","):
        part = part.strip()
        part = re.sub(r"\(.*?\)", "", part)
        if part and part != ".":
            genes.append(part)
    return genes


def wrap_fasta(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


# =========================
# LOAD KALLISTO
# =========================

def load_kallisto(path):
    df = pd.read_csv(path, sep="\t")

    if "target_id" not in df.columns or "tpm" not in [c.lower() for c in df.columns]:
        raise ValueError("Le fichier kallisto doit contenir target_id et TPM.")

    tpm_col = [c for c in df.columns if c.lower() == "tpm"][0]

    tpm = {}
    for _, row in df.iterrows():
        tx = strip_version(row["target_id"])
        tpm[tx] = float(row[tpm_col])

    return tpm


# =========================
# LOAD FASTA
# =========================

def load_transcriptome_fasta(path):
    fasta = {}

    for record in SeqIO.parse(path, "fasta"):
        raw_id = record.id
        tx_id = strip_version(raw_id.split("|")[0].split()[0])
        fasta[tx_id] = str(record.seq).upper()

    return fasta


# =========================
# PARSING GTF
# =========================

def load_gtf_transcripts(gtf_path):
    transcripts = {}
    gene_to_transcripts = {}

    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs_txt = fields

            if feature not in {"exon", "CDS"}:
                continue

            attrs = parse_gtf_attributes(attrs_txt)

            transcript_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id", ".")
            gene_name = attrs.get("gene_name", gene_id)
            transcript_biotype = attrs.get("transcript_biotype", attrs.get("transcript_type", "."))
            gene_biotype = attrs.get("gene_biotype", attrs.get("gene_type", "."))

            if transcript_id is None:
                continue

            tx = strip_version(transcript_id)
            chrom = normalize_chrom(chrom)

            if tx not in transcripts:
                transcripts[tx] = {
                    "transcript_id": tx,
                    "gene_id": strip_version(gene_id),
                    "gene_name": gene_name,
                    "gene_biotype": gene_biotype,
                    "transcript_biotype": transcript_biotype,
                    "chrom": chrom,
                    "strand": strand,
                    "exons_genomic": [],
                    "exons_tx_order": [],
                    "cds_genomic": [],
                    "cds_cdna": None,
                }

            if feature == "exon":
                transcripts[tx]["exons_genomic"].append((int(start), int(end)))

            if feature == "CDS":
                transcripts[tx]["cds_genomic"].append((int(start), int(end)))
            gene_to_transcripts.setdefault(gene_name, set()).add(tx)
            gene_to_transcripts.setdefault(strip_version(gene_id), set()).add(tx)

    for tx, info in transcripts.items():
        exons = sorted(info["exons_genomic"], key=lambda x: x[0])

        if info["strand"] == "+":
            info["exons_tx_order"] = exons
        else:
            info["exons_tx_order"] = sorted(exons, key=lambda x: x[0], reverse=True)
        cds_positions = []

        for cds_start, cds_end in info["cds_genomic"]:
            cds_positions.extend(
                interval_genomic_to_cdna(
                    info["exons_tx_order"],
                    cds_start,
                    cds_end,
                    info["strand"]
                )
            )

        if cds_positions:
            cds_min = min(x[0] for x in cds_positions)
            cds_max = max(x[1] for x in cds_positions)
            info["cds_cdna"] = (cds_min, cds_max)
        else:
            info["cds_cdna"] = None
    gene_to_transcripts = {g: sorted(list(v)) for g, v in gene_to_transcripts.items()}

    return transcripts, gene_to_transcripts

def interval_genomic_to_cdna(exons_tx_order, cds_start, cds_end, strand):
    positions = []

    cdna_offset = 0

    for exon_start, exon_end in exons_tx_order:
        exon_len = exon_end - exon_start + 1

        overlap_start = max(exon_start, cds_start)
        overlap_end = min(exon_end, cds_end)

        if overlap_start <= overlap_end:
            if strand == "+":
                cdna_start = cdna_offset + (overlap_start - exon_start)
                cdna_end = cdna_offset + (overlap_end - exon_start)
            else:
                cdna_start = cdna_offset + (exon_end - overlap_end)
                cdna_end = cdna_offset + (exon_end - overlap_start)

            positions.append((cdna_start, cdna_end))

        cdna_offset += exon_len

    return positions

def region_of_pos(pos, cds_cdna):
    if pos is None:
        return "intronic_or_unmapped"

    if cds_cdna is None:
        return "noncoding"

    cds_start, cds_end = cds_cdna

    if pos < cds_start:
        return "5UTR"

    if cds_start <= pos <= cds_end:
        return "CDS"

    return "3UTR"


def is_frame_preserved(bp1_pos, bp2_pos, cds1, cds2):
    if cds1 is None or cds2 is None:
        return False

    cds1_start, cds1_end = cds1
    cds2_start, cds2_end = cds2

    left_cds_len = bp1_pos - cds1_start + 1
    right_cds_offset = bp2_pos - cds2_start

    return (left_cds_len + right_cds_offset) % 3 == 0


def classify_fusion(bp1_pos, bp2_pos, cds1, cds2):
    r1 = region_of_pos(bp1_pos, cds1)
    r2 = region_of_pos(bp2_pos, cds2)

    if r1 == "CDS" and r2 == "CDS" and cds1 is not None and cds2 is not None:
        if is_frame_preserved(bp1_pos, bp2_pos, cds1, cds2):
            return "CDS_CDS_in_frame", r1, r2
        else:
            return "CDS_CDS_out_of_frame", r1, r2

    elif r1 == "5UTR" and r2 == "5UTR":
        return "promoter_swap_or_5UTR_fusion", r1, r2

    elif r1 == "3UTR" and r2 == "3UTR":
        return "3UTR_fusion", r1, r2

    elif r1 == "CDS" and r2 == "3UTR":
        return "CDS_to_3UTR", r1, r2

    elif r1 == "5UTR" and r2 == "CDS":
        return "5UTR_to_CDS", r1, r2

    else:
        return "complex_or_noncanonical", r1, r2
    
# =========================
# Choix du transcrit
# =========================

def transcript_exists(tx, transcripts, fasta):
    tx = strip_version(tx)
    return tx in transcripts and tx in fasta


def pick_transcript_from_gene(gene_symbols, gene_to_transcripts, transcripts, fasta, tpm_dict, chrom=None):
    candidates = []

    for gene in gene_symbols:
        for tx in gene_to_transcripts.get(gene, []):
            if tx not in transcripts:
                continue
            if tx not in fasta:
                continue
            if chrom is not None and transcripts[tx]["chrom"] != normalize_chrom(chrom):
                continue

            seq_len = len(fasta[tx])
            tpm = tpm_dict.get(tx, 0.0)

            candidates.append({
                "tx": tx,
                "tpm": tpm,
                "seq_len": seq_len,
                "gene_name": transcripts[tx]["gene_name"],
                "biotype": transcripts[tx]["transcript_biotype"]
            })

    if not candidates:
        return None, "no_candidate"

    candidates = sorted(
        candidates,
        key=lambda x: (
            x["tpm"] > MIN_TPM,
            x["tpm"],
            x["biotype"] == "protein_coding",
            x["seq_len"]
        ),
        reverse=True
    )

    return candidates[0]["tx"], f"picked_by_kallisto_tpm={candidates[0]['tpm']}"


def resolve_transcript(arriba_tx, gene_name, gene_id, chrom, gene_to_transcripts, transcripts, fasta, tpm_dict):
    arriba_tx = strip_version(arriba_tx)

    if transcript_exists(arriba_tx, transcripts, fasta):
        return arriba_tx, "arriba_transcript_id"

    gene_symbols = []

    for g in clean_gene_symbol(gene_name):
        gene_symbols.append(g)

    if gene_id is not None and gene_id != ".":
        gene_symbols.append(strip_version(gene_id))

    picked, reason = pick_transcript_from_gene(
        gene_symbols=gene_symbols,
        gene_to_transcripts=gene_to_transcripts,
        transcripts=transcripts,
        fasta=fasta,
        tpm_dict=tpm_dict,
        chrom=chrom
    )

    return picked, reason


# =========================
# GENOMIC to cDNA
# =========================

def genomic_to_cdna_indices(exons_tx_order, breakpoint, strand):
    """
    Retourne:
    - left_cut_index: index Python pour garder seq[:left_cut_index]
      donc inclut la base du breakpoint si le breakpoint est dans un exon.

    - right_start_index: index Python pour garder seq[right_start_index:]
      donc commence à la base du breakpoint si le breakpoint est dans un exon.

    Si intronique, retourne None.
    """
    cdna_offset = 0

    for start, end in exons_tx_order:
        exon_len = end - start + 1

        if start <= breakpoint <= end:
            if strand == "+":
                zero_based_in_exon = breakpoint - start
            else:
                zero_based_in_exon = end - breakpoint

            right_start_index = cdna_offset + zero_based_in_exon
            left_cut_index = right_start_index + 1

            return left_cut_index, right_start_index, "exonic"

        cdna_offset += exon_len

    return None, None, "not_exonic"


def snap_intronic_to_exon_boundary(exons_tx_order, breakpoint, side):
    """
    Si breakpoint intronique, approxime au bord exonique le plus proche.

    side = "left":
        pour le gène 5', on garde la partie amont transcriptomique.
        On coupe à la fin de l'exon précédent dans l'ordre transcript.

    side = "right":
        pour le gène 3', on garde la partie aval transcriptomique.
        On commence au début de l'exon suivant dans l'ordre transcript.

    Retourne:
    - left_cut_index
    - right_start_index
    - status
    """
    cdna_offset = 0
    exon_boundaries = []

    for start, end in exons_tx_order:
        exon_len = end - start + 1
        exon_boundaries.append({
            "start": start,
            "end": end,
            "cdna_start": cdna_offset,
            "cdna_end_exclusive": cdna_offset + exon_len
        })
        cdna_offset += exon_len

    compatible = []

    for b in exon_boundaries:
        start = b["start"]
        end = b["end"]

        dist_to_start = abs(breakpoint - start)
        dist_to_end = abs(breakpoint - end)

        compatible.append((dist_to_start, b["cdna_start"], b))
        compatible.append((dist_to_end, b["cdna_end_exclusive"], b))

    if not compatible:
        return None, None, "snap_failed_no_exon"

    compatible = sorted(compatible, key=lambda x: x[0])
    best_dist, best_cdna_boundary, best_exon = compatible[0]

    if side == "left":
        left_cut_index = best_cdna_boundary
        right_start_index = best_cdna_boundary
    else:
        left_cut_index = best_cdna_boundary
        right_start_index = best_cdna_boundary

    return left_cut_index, right_start_index, f"intronic_snapped_distance={best_dist}"


# =========================
# FUSION mRNA
# =========================

def build_fusion_mrna(row, tx1, tx2, transcripts, fasta):
    chrom1, bp1 = parse_breakpoint(row["breakpoint1"])
    chrom2, bp2 = parse_breakpoint(row["breakpoint2"])

    t1 = transcripts[tx1]
    t2 = transcripts[tx2]

    if normalize_chrom(chrom1) != t1["chrom"]:
        return None, {
            "status": "fail",
            "reason": f"chrom_mismatch_tx1_arriba={chrom1}_gtf={t1['chrom']}"
        }

    if normalize_chrom(chrom2) != t2["chrom"]:
        return None, {
            "status": "fail",
            "reason": f"chrom_mismatch_tx2_arriba={chrom2}_gtf={t2['chrom']}"
        }

    # 1. Mapping réel du breakpoint vers le cDNA
    cut1, start1, status1 = genomic_to_cdna_indices(
        t1["exons_tx_order"], bp1, t1["strand"]
    )

    cut2, start2, status2 = genomic_to_cdna_indices(
        t2["exons_tx_order"], bp2, t2["strand"]
    )

    # 2. Garder les vraies positions AVANT snap
    real_bp1_pos = cut1 - 1 if cut1 is not None else None
    real_bp2_pos = start2 if start2 is not None else None

    # 3. Snap seulement pour reconstruire la séquence si breakpoint intronique
    if cut1 is None and ALLOW_INTRONIC_SNAP:
        cut1, start1, status1 = snap_intronic_to_exon_boundary(
            t1["exons_tx_order"],
            bp1,
            side="left"
        )

    if start2 is None and ALLOW_INTRONIC_SNAP:
        cut2, start2, status2 = snap_intronic_to_exon_boundary(
            t2["exons_tx_order"],
            bp2,
            side="right"
        )

    # 4. Vérifier qu'on peut reconstruire la fusion
    if cut1 is None:
        return None, {
            "status": "fail",
            "reason": "breakpoint1_not_mappable_to_cdna",
            "bp1_status": status1,
            "bp2_status": status2
        }

    if start2 is None:
        return None, {
            "status": "fail",
            "reason": "breakpoint2_not_mappable_to_cdna",
            "bp1_status": status1,
            "bp2_status": status2
        }

    # 5. Classification avec les vraies positions, pas les positions snap
    fusion_type, region1, region2 = classify_fusion(
        bp1_pos=real_bp1_pos,
        bp2_pos=real_bp2_pos,
        cds1=t1["cds_cdna"],
        cds2=t2["cds_cdna"]
    )

    seq1 = fasta[tx1]
    seq2 = fasta[tx2]

    left_part = seq1[:cut1]
    right_part = seq2[start2:]

    fusion_seq = left_part + right_part

    diagnostics = {
        "status": "ok",
        "reason": "built",
        "bp1_status": status1,
        "bp2_status": status2,
        "tx1_cut_index": cut1,
        "tx2_start_index": start2,
        "real_bp1_pos": real_bp1_pos,
        "real_bp2_pos": real_bp2_pos,
        "left_len": len(left_part),
        "right_len": len(right_part),
        "fusion_len": len(fusion_seq),
        "fusion_type": fusion_type,
        "bp1_region": region1,
        "bp2_region": region2,
        "tx1_cds_cdna": t1["cds_cdna"],
        "tx2_cds_cdna": t2["cds_cdna"],
    }

    return fusion_seq, diagnostics

# =========================
# LOAD ARRIBA
# =========================

def load_arriba(path):
    with open(path, "r") as f:
        first_line = f.readline().rstrip("\n")

    if first_line.startswith("#"):
        header = first_line[1:].split("\t")
        df = pd.read_csv(path, sep="\t", comment="#", names=header)
    else:
        df = pd.read_csv(path, sep="\t")

    required = [
        "gene1", "gene2",
        "breakpoint1", "breakpoint2",
        "gene_id1", "gene_id2",
        "transcript_id1", "transcript_id2"
    ]

    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Colonnes Arriba manquantes: {missing}")

    return df



# =========================
# STAR-FUSION
# =========================

def parse_star_gene(value):
    """PVT1^ENSG00000249859.14 -> (PVT1, ENSG00000249859)."""
    text = str(value).strip()
    if "^" in text:
        symbol, gene_id = text.split("^", 1)
    else:
        symbol = text
        gene_id = text if text.startswith("ENSG") else "."
    return symbol.strip(), strip_version(gene_id)


def parse_star_breakpoint(value):
    """chr8:127989291:+ -> ('8', 127989291, '+')."""
    parts = str(value).strip().rsplit(":", 2)
    if len(parts) != 3:
        raise ValueError(f"Breakpoint STAR-Fusion invalide: {value!r}")
    chrom, pos, strand = parts
    if strand not in {"+", "-"}:
        raise ValueError(f"Brin STAR-Fusion invalide: {value!r}")
    return normalize_chrom(chrom), int(pos), strand


def load_star_fusion(path):
    df = pd.read_csv(path, sep="\t", dtype=str)
    required = [
        "#FusionName", "LeftGene", "LeftBreakpoint",
        "RightGene", "RightBreakpoint", "JunctionReadCount",
        "SpanningFragCount", "SpliceType", "LargeAnchorSupport", "FFPM",
    ]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"Colonnes STAR-Fusion manquantes: {missing}")
    return df


def star_boundary_candidates(info):
    boundaries = []
    cdna_offset = 0
    for exon_start, exon_end in info["exons_tx_order"]:
        exon_len = exon_end - exon_start + 1
        boundaries.extend([
            {"genomic": exon_start, "cdna": cdna_offset},
            {"genomic": exon_end, "cdna": cdna_offset + exon_len},
        ])
        cdna_offset += exon_len
    return boundaries


def evaluate_star_breakpoint(info, breakpoint):
    cut, start, status = genomic_to_cdna_indices(
        info["exons_tx_order"], breakpoint, info["strand"]
    )
    if cut is not None and start is not None:
        return {
            "compatible": True, "mapping": "exonic", "snap_distance": 0,
            "cut_index": cut, "start_index": start,
        }

    if not ALLOW_INTRONIC_SNAP:
        return {
            "compatible": False, "mapping": "not_exonic", "snap_distance": None,
            "cut_index": None, "start_index": None,
        }

    boundaries = star_boundary_candidates(info)
    if not boundaries:
        return {
            "compatible": False, "mapping": "no_exon", "snap_distance": None,
            "cut_index": None, "start_index": None,
        }

    best = min(boundaries, key=lambda b: abs(breakpoint - b["genomic"]))
    distance = abs(breakpoint - best["genomic"])
    if MAX_SNAP_DISTANCE is not None and distance > int(MAX_SNAP_DISTANCE):
        return {
            "compatible": False, "mapping": "snap_too_far", "snap_distance": distance,
            "cut_index": None, "start_index": None,
        }

    boundary = best["cdna"]
    return {
        "compatible": True, "mapping": "intronic_snap", "snap_distance": distance,
        "cut_index": boundary, "start_index": boundary,
    }


def choose_star_transcript(gene_name, gene_id, chrom, star_strand, breakpoint,
                           side, gene_to_transcripts, transcripts, fasta,
                           tpm_dict, fusion_index):
    candidate_ids = (
        set(gene_to_transcripts.get(gene_name, [])) |
        set(gene_to_transcripts.get(gene_id, []))
    )
    rows = []
    for tx in sorted(candidate_ids):
        info = transcripts.get(tx)
        if info is None or tx not in fasta:
            continue
        if info["chrom"] != normalize_chrom(chrom):
            continue

        mapping = evaluate_star_breakpoint(info, breakpoint)
        rows.append({
            "source": "STAR-Fusion",
            "fusion_index": fusion_index,
            "side": side,
            "gene_name": gene_name,
            "gene_id": gene_id,
            "breakpoint_chrom": normalize_chrom(chrom),
            "breakpoint_pos": breakpoint,
            "star_strand": star_strand,
            "tx": tx,
            "tx_strand": info["strand"],
            "strand_matches_star": info["strand"] == star_strand,
            "biotype": info["transcript_biotype"],
            "tpm": tpm_dict.get(tx, 0.0),
            "sequence_length": len(fasta[tx]),
            **mapping,
        })

    compatible = [row for row in rows if row["compatible"]]
    compatible.sort(
        key=lambda row: (
            row["mapping"] == "exonic",
            row["strand_matches_star"],
            row["tpm"] > MIN_TPM,
            row["tpm"],
            row["biotype"] == "protein_coding",
            -int(row["snap_distance"] or 0),
            row["sequence_length"],
        ),
        reverse=True,
    )

    compatible_ids = {id(row) for row in compatible}
    for rank, candidate in enumerate(compatible, start=1):
        candidate["compatible_rank"] = rank
        candidate["selected"] = rank == 1
    for candidate in rows:
        if id(candidate) not in compatible_ids:
            candidate["compatible_rank"] = "."
            candidate["selected"] = False

    if not compatible:
        return None, None, rows
    return compatible[0]["tx"], compatible[0], rows


def write_fasta(path, records):
    with open(path, "w") as handle:
        for header, sequence in records:
            handle.write(header + "\n")
            handle.write(wrap_fasta(sequence) + "\n")


def deduplicate_records(arriba_records, star_records):
    """Déduplique par séquence exacte et conserve toutes les provenances."""
    unique = {}
    order = []
    membership = []

    for source, records in (("Arriba", arriba_records), ("STAR-Fusion", star_records)):
        for header, sequence in records:
            digest = sha256(sequence.encode()).hexdigest()
            record_id = header[1:].split("|", 1)[0]
            if digest not in unique:
                unique[digest] = {
                    "sequence": sequence,
                    "ids": [record_id],
                    "sources": [source],
                    "headers": [header],
                }
                order.append(digest)
                duplicate_of = "."
            else:
                unique[digest]["ids"].append(record_id)
                if source not in unique[digest]["sources"]:
                    unique[digest]["sources"].append(source)
                unique[digest]["headers"].append(header)
                duplicate_of = unique[digest]["ids"][0]
            membership.append({
                "source": source,
                "record_id": record_id,
                "sequence_sha256": digest,
                "duplicate_of": duplicate_of,
            })

    combined = []
    for index, digest in enumerate(order, start=1):
        item = unique[digest]
        combined_id = f"fusion_unique_{index}"
        header = (
            f">{combined_id}|sources={';'.join(item['sources'])}|"
            f"original_ids={';'.join(item['ids'])}|sequence_sha256={digest}"
        )
        combined.append((header, item["sequence"]))

    return combined, membership

# =========================
# MAIN
# =========================

def validate_inputs():
    required_inputs = {
        "Arriba": ARRIBA_TSV,
        "STAR-Fusion": STAR_FUSION_TSV,
        "GTF": GTF_FILE,
        "transcriptome FASTA": TRANSCRIPTOME_FASTA,
        "Kallisto abundance": KALLISTO_ABUNDANCE,
    }
    missing = [f"{name}: {path}" for name, path in required_inputs.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Fichiers d'entrée absents:\n- " + "\n- ".join(missing))


def main():
    validate_inputs()
    outputs = (
        OUT_ARRIBA_FASTA, OUT_STAR_FASTA, OUT_COMBINED_FASTA,
        OUT_ARRIBA_SUMMARY, OUT_STAR_SUMMARY, OUT_COMBINED_SUMMARY,
        OUT_STAR_CANDIDATES,
    )
    for output in outputs:
        output.parent.mkdir(parents=True, exist_ok=True)

    print("Loading Kallisto...")
    tpm_dict = load_kallisto(KALLISTO_ABUNDANCE)
    print("Loading transcriptome FASTA...")
    fasta = load_transcriptome_fasta(TRANSCRIPTOME_FASTA)
    print("Loading GTF...")
    transcripts, gene_to_transcripts = load_gtf_transcripts(GTF_FILE)

    # -------------------------------------------------------------------------
    # 1. ARRIBA — logique existante conservée
    # -------------------------------------------------------------------------
    print("Processing Arriba...")
    arriba = load_arriba(ARRIBA_TSV)
    arriba_summary = []
    arriba_records = []

    for idx, row in arriba.iterrows():
        fusion_id = re.sub(
            r"[^A-Za-z0-9_.|-]+", "_",
            f"arriba_{idx + 1}_{row['gene1']}_{row['gene2']}"
        )
        chrom1, _ = parse_breakpoint(row["breakpoint1"])
        chrom2, _ = parse_breakpoint(row["breakpoint2"])

        tx1, reason_tx1 = resolve_transcript(
            row.get("transcript_id1", "."), row.get("gene1", "."),
            row.get("gene_id1", "."), chrom1, gene_to_transcripts,
            transcripts, fasta, tpm_dict,
        )
        tx2, reason_tx2 = resolve_transcript(
            row.get("transcript_id2", "."), row.get("gene2", "."),
            row.get("gene_id2", "."), chrom2, gene_to_transcripts,
            transcripts, fasta, tpm_dict,
        )

        summary = {
            "source": "Arriba", "fusion_index": idx + 1, "fusion_id": fusion_id,
            "gene1": row.get("gene1", "."), "gene2": row.get("gene2", "."),
            "breakpoint1": row.get("breakpoint1", "."),
            "breakpoint2": row.get("breakpoint2", "."),
            "site1": row.get("site1", "."), "site2": row.get("site2", "."),
            "type": row.get("type", "."), "confidence": row.get("confidence", "."),
            "reading_frame_arriba": row.get("reading_frame", "."),
            "tx1": tx1 or ".", "tx2": tx2 or ".",
            "tx1_choice_reason": reason_tx1, "tx2_choice_reason": reason_tx2,
            "tx1_tpm": tpm_dict.get(tx1, 0.0) if tx1 else 0.0,
            "tx2_tpm": tpm_dict.get(tx2, 0.0) if tx2 else 0.0,
        }

        if tx1 is None or tx2 is None:
            summary.update({"status": "fail", "reason": "missing_transcript_choice"})
            arriba_summary.append(summary)
            continue

        fusion_seq, diag = build_fusion_mrna(row, tx1, tx2, transcripts, fasta)
        summary.update(diag)
        arriba_summary.append(summary)
        if fusion_seq is not None:
            header = (
                f">{fusion_id}|source=Arriba|gene1={row.get('gene1', '.')}|"
                f"gene2={row.get('gene2', '.')}|tx1={tx1}|tx2={tx2}|"
                f"bp1={row.get('breakpoint1', '.')}|bp2={row.get('breakpoint2', '.')}|"
                f"tx1_tpm={tpm_dict.get(tx1, 0.0)}|tx2_tpm={tpm_dict.get(tx2, 0.0)}"
            )
            arriba_records.append((header, fusion_seq))

    # -------------------------------------------------------------------------
    # 2. STAR-FUSION — sélection breakpoint-compatible + GTF + Kallisto
    # -------------------------------------------------------------------------
    print("Processing STAR-Fusion...")
    star = load_star_fusion(STAR_FUSION_TSV)
    star_summary = []
    star_candidates = []
    star_records = []

    for idx, row in star.iterrows():
        fusion_index = idx + 1
        fusion_name = str(row["#FusionName"])
        fusion_id = re.sub(r"[^A-Za-z0-9_.|-]+", "_", f"starfusion_{fusion_index}_{fusion_name}")
        gene1, gene_id1 = parse_star_gene(row["LeftGene"])
        gene2, gene_id2 = parse_star_gene(row["RightGene"])
        chrom1, bp1, strand1 = parse_star_breakpoint(row["LeftBreakpoint"])
        chrom2, bp2, strand2 = parse_star_breakpoint(row["RightBreakpoint"])

        tx1, choice1, candidates1 = choose_star_transcript(
            gene1, gene_id1, chrom1, strand1, bp1, "left",
            gene_to_transcripts, transcripts, fasta, tpm_dict, fusion_index,
        )
        tx2, choice2, candidates2 = choose_star_transcript(
            gene2, gene_id2, chrom2, strand2, bp2, "right",
            gene_to_transcripts, transcripts, fasta, tpm_dict, fusion_index,
        )
        star_candidates.extend(candidates1)
        star_candidates.extend(candidates2)

        summary = {
            "source": "STAR-Fusion", "fusion_index": fusion_index,
            "fusion_id": fusion_id, "fusion_name": fusion_name,
            "gene1": gene1, "gene_id1": gene_id1,
            "breakpoint1": row["LeftBreakpoint"],
            "gene2": gene2, "gene_id2": gene_id2,
            "breakpoint2": row["RightBreakpoint"],
            "junction_reads": row["JunctionReadCount"],
            "spanning_fragments": row["SpanningFragCount"],
            "splice_type": row["SpliceType"],
            "large_anchor_support": row["LargeAnchorSupport"],
            "ffpm": row["FFPM"],
            "tx1": tx1 or ".", "tx2": tx2 or ".",
            "tx1_tpm": choice1["tpm"] if choice1 else 0.0,
            "tx2_tpm": choice2["tpm"] if choice2 else 0.0,
            "tx1_mapping": choice1["mapping"] if choice1 else ".",
            "tx2_mapping": choice2["mapping"] if choice2 else ".",
            "tx1_snap_distance": choice1["snap_distance"] if choice1 else ".",
            "tx2_snap_distance": choice2["snap_distance"] if choice2 else ".",
        }

        if tx1 is None or tx2 is None or choice1 is None or choice2 is None:
            summary.update({"status": "fail", "reason": "no_compatible_transcript_pair"})
            star_summary.append(summary)
            continue

        left_cut = int(choice1["cut_index"])
        right_start = int(choice2["start_index"])
        left_seq = fasta[tx1][:left_cut]
        right_seq = fasta[tx2][right_start:]
        fusion_seq = left_seq + right_seq
        digest = sha256(fusion_seq.encode()).hexdigest()
        summary.update({
            "status": "ok", "reason": "built",
            "left_cut_index": left_cut, "right_start_index": right_start,
            "left_len": len(left_seq), "right_len": len(right_seq),
            "fusion_len": len(fusion_seq), "sequence_sha256": digest,
        })
        star_summary.append(summary)
        header = (
            f">{fusion_id}|source=STAR-Fusion|fusion={fusion_name}|"
            f"gene1={gene1}|gene2={gene2}|tx1={tx1}|tx2={tx2}|"
            f"bp1={row['LeftBreakpoint']}|bp2={row['RightBreakpoint']}|"
            f"tx1_tpm={choice1['tpm']}|tx2_tpm={choice2['tpm']}|"
            f"map1={choice1['mapping']}|map2={choice2['mapping']}"
        )
        star_records.append((header, fusion_seq))

    # -------------------------------------------------------------------------
    # 3. SORTIES INDIVIDUELLES ET COMBINAISON DÉDUPLIQUÉE
    # -------------------------------------------------------------------------
    write_fasta(OUT_ARRIBA_FASTA, arriba_records)
    write_fasta(OUT_STAR_FASTA, star_records)
    pd.DataFrame(arriba_summary).to_csv(OUT_ARRIBA_SUMMARY, sep="\t", index=False)
    pd.DataFrame(star_summary).to_csv(OUT_STAR_SUMMARY, sep="\t", index=False)
    pd.DataFrame(star_candidates).to_csv(OUT_STAR_CANDIDATES, sep="\t", index=False)

    combined_records, membership = deduplicate_records(arriba_records, star_records)
    write_fasta(OUT_COMBINED_FASTA, combined_records)
    pd.DataFrame(membership).to_csv(OUT_COMBINED_SUMMARY, sep="\t", index=False)

    print("Done.")
    print(f"Arriba reconstructed: {len(arriba_records)}/{len(arriba)}")
    print(f"STAR-Fusion reconstructed: {len(star_records)}/{len(star)}")
    print(f"Combined before deduplication: {len(arriba_records) + len(star_records)}")
    print(f"Combined unique sequences: {len(combined_records)}")


if __name__ == "__main__":
    main()
