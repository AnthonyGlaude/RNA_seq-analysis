#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import re
import pandas as pd


# =========================
# CONFIGURATION
# =========================

GTF = Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/data/references/gtf/Homo_sapiens.GRCh38.113.gtf")

SAMPLES = {
    "SH-SY5Y": {
        "srr": "SRR4787043",
        "fasta": Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/final/SRR4787043/combined_transcripts.fa"),
        "log_dir": Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/logs/variants/SRR4787043"),
    },
    "SK-N-Be(2)": {
        "srr": "SRR4787045",
        "fasta": Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/final/SRR4787045/combined_transcripts.fa"),
        "log_dir": Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/logs/variants/SRR4787045"),
    },
}

OUTDIR = Path("transcript_reconstruction_summary")
OUTDIR.mkdir(exist_ok=True)


# =========================
# BASIC FUNCTIONS
# =========================

def extract_attr(attr_string, key):
    target = f'{key} "'
    if target not in attr_string:
        return None
    return attr_string.split(target, 1)[1].split('"', 1)[0]


def strip_version(tx_id):
    return tx_id.split(".")[0]


def parse_transcript_id(header):
    """
    Extrait le premier identifiant ENST trouvé dans le header FASTA.
    """
    m = re.search(r"(ENST[0-9]+(?:\.[0-9]+)?)", header)
    if not m:
        return None
    return m.group(1)


def is_mutated_header(header):
    """
    Détecte les entrées mutées dans le FASTA.
    À ajuster seulement si tes headers n'utilisent pas '_mut'.
    """
    h = header.lower()
    return "_mut" in h or "|mut" in h or "mut_" in h or "mutated" in h


def median_from_counter_values(counter):
    """
    Médiane des valeurs d'un Counter sans passer par pandas Series.
    Ici, counter = transcript_id -> nombre de versions mutées.
    """
    values = sorted(counter.values())
    n = len(values)

    if n == 0:
        return 0

    mid = n // 2

    if n % 2 == 1:
        return values[mid]

    return (values[mid - 1] + values[mid]) / 2


# =========================
# BIOTYPE GROUPING
# =========================

def group_biotype(biotype):
    """
    Regroupe les biotypes Ensembl en grandes classes lisibles pour les résultats.
    """

    if biotype is None or biotype == "NA":
        return "unknown"

    protein_coding_related = {
        "protein_coding",
        "protein_coding_CDS_not_defined",
        "nonsense_mediated_decay",
        "non_stop_decay",
    }

    non_coding_transcript = {
        "lncRNA",
        "processed_transcript",
        "retained_intron",
        "TEC",
        "macro_lncRNA",
        "bidirectional_promoter_lncRNA",
        "antisense",
        "sense_intronic",
        "sense_overlapping",
        "3prime_overlapping_ncRNA",
    }

    small_or_structural_rna = {
        "miRNA",
        "snoRNA",
        "snRNA",
        "scaRNA",
        "rRNA",
        "ribozyme",
        "misc_RNA",
        "sRNA",
        "vault_RNA",
        "Mt_rRNA",
        "Mt_tRNA",
    }

    if biotype in protein_coding_related:
        return "protein_coding_related"

    if biotype in non_coding_transcript:
        return "non_coding_transcript"

    if "pseudogene" in biotype:
        return "pseudogene"

    if biotype in small_or_structural_rna:
        return "small_or_structural_RNA"

    if biotype.startswith("IG_") or biotype.startswith("TR_"):
        return "IG_TR"

    return "other"


# =========================
# GTF LOADING
# =========================

def load_gtf_biotypes(gtf_path):
    """
    Retourne transcript_id -> biotype brut.
    Gère les transcript_id avec et sans version.
    """
    tx_to_biotype = {}

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            feature = parts[2]
            attrs = parts[8]

            if feature != "transcript":
                continue

            tx = extract_attr(attrs, "transcript_id")
            biotype = (
                extract_attr(attrs, "transcript_biotype")
                or extract_attr(attrs, "transcript_type")
                or extract_attr(attrs, "gene_biotype")
                or extract_attr(attrs, "gene_type")
                or "NA"
            )

            if tx:
                tx_to_biotype[tx] = biotype
                tx_to_biotype[strip_version(tx)] = biotype

    return tx_to_biotype


# =========================
# FASTA SUMMARY
# =========================

def summarize_fasta(fasta_path, tx_to_biotype):
    """
    Lit seulement les headers FASTA.
    Ne charge pas les séquences en mémoire.
    """

    total_records = 0
    ref_records = 0
    mut_records = 0
    no_enst_records = 0

    unique_ref_tx = set()
    unique_mut_tx = set()
    unique_all_tx = set()

    mut_versions_per_tx = Counter()

    # Biotypes comptés par entrée FASTA
    ref_biotypes_raw_entries = Counter()
    mut_biotypes_raw_entries = Counter()
    ref_biotypes_grouped_entries = Counter()
    mut_biotypes_grouped_entries = Counter()

    # Biotypes comptés par transcrit unique
    ref_tx_to_biotype = {}
    mut_tx_to_biotype = {}

    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue

            total_records += 1
            header = line[1:].strip()

            tx = parse_transcript_id(header)

            if tx is None:
                no_enst_records += 1
                continue

            tx_nover = strip_version(tx)
            unique_all_tx.add(tx_nover)

            raw_biotype = tx_to_biotype.get(tx) or tx_to_biotype.get(tx_nover) or "NA"
            grouped_biotype = group_biotype(raw_biotype)

            if is_mutated_header(header):
                mut_records += 1
                unique_mut_tx.add(tx_nover)
                mut_versions_per_tx[tx_nover] += 1

                mut_biotypes_raw_entries[raw_biotype] += 1
                mut_biotypes_grouped_entries[grouped_biotype] += 1
                mut_tx_to_biotype[tx_nover] = raw_biotype

            else:
                ref_records += 1
                unique_ref_tx.add(tx_nover)

                ref_biotypes_raw_entries[raw_biotype] += 1
                ref_biotypes_grouped_entries[grouped_biotype] += 1
                ref_tx_to_biotype[tx_nover] = raw_biotype

    vals = list(mut_versions_per_tx.values())

    summary = {
        "total_fasta_entries": total_records,
        "reference_entries": ref_records,
        "mutated_entries": mut_records,
        "entries_without_enst": no_enst_records,
        "unique_reference_transcripts": len(unique_ref_tx),
        "unique_transcripts_with_mutated_version": len(unique_mut_tx),
        "unique_transcripts_total": len(unique_all_tx),
        "mean_mutated_versions_per_transcript": round(sum(vals) / len(vals), 3) if vals else 0,
        "median_mutated_versions_per_transcript": median_from_counter_values(mut_versions_per_tx),
        "max_mutated_versions_per_transcript": max(vals) if vals else 0,
    }

    # Biotypes par transcrit unique
    ref_biotypes_raw_unique_tx = Counter(ref_tx_to_biotype.values())
    mut_biotypes_raw_unique_tx = Counter(mut_tx_to_biotype.values())

    ref_biotypes_grouped_unique_tx = Counter(
        group_biotype(biotype) for biotype in ref_tx_to_biotype.values()
    )
    mut_biotypes_grouped_unique_tx = Counter(
        group_biotype(biotype) for biotype in mut_tx_to_biotype.values()
    )

    biotype_counters = {
        "raw_reference_entries": ref_biotypes_raw_entries,
        "raw_mutated_entries": mut_biotypes_raw_entries,
        "grouped_reference_entries": ref_biotypes_grouped_entries,
        "grouped_mutated_entries": mut_biotypes_grouped_entries,
        "raw_reference_unique_transcripts": ref_biotypes_raw_unique_tx,
        "raw_mutated_unique_transcripts": mut_biotypes_raw_unique_tx,
        "grouped_reference_unique_transcripts": ref_biotypes_grouped_unique_tx,
        "grouped_mutated_unique_transcripts": mut_biotypes_grouped_unique_tx,
    }

    return summary, biotype_counters, mut_versions_per_tx


# =========================
# LOG SUMMARY
# =========================

def summarize_logs(log_dir):
    """
    Analyse tous les logs apply_variants_*.log.
    """

    log_files = sorted(log_dir.glob("apply_variants_*.log"))

    ignored_mutations = 0
    ignored_combos = 0

    tx_with_ignored_mutations = set()
    tx_with_ignored_combos = set()

    tx_re = re.compile(r"(ENST[0-9]+(?:\.[0-9]+)?)")

    for log_file in log_files:
        with open(log_file) as f:
            for line in f:
                if "Mutation ignorée" in line:
                    ignored_mutations += 1
                    m = tx_re.search(line)
                    if m:
                        tx_with_ignored_mutations.add(strip_version(m.group(1)))

                if "Combo ignoré" in line:
                    ignored_combos += 1
                    m = tx_re.search(line)
                    if m:
                        tx_with_ignored_combos.add(strip_version(m.group(1)))

    return {
        "n_log_files": len(log_files),
        "ignored_mutations": ignored_mutations,
        "ignored_combos": ignored_combos,
        "transcripts_with_ignored_mutations": len(tx_with_ignored_mutations),
        "transcripts_with_ignored_combos": len(tx_with_ignored_combos),
    }


# =========================
# OUTPUT TABLES
# =========================

def counter_to_df(counter, sample, counter_type):
    rows = []
    total = sum(counter.values())

    for category, n in counter.most_common():
        rows.append({
            "sample": sample,
            "counter_type": counter_type,
            "category": category,
            "n": n,
            "pct": (n / total * 100) if total else 0,
        })

    return pd.DataFrame(rows)


def mut_versions_to_df(counter, sample):
    """
    Distribution du nombre de versions mutées générées par transcrit.
    """
    dist = Counter(counter.values())

    rows = []

    for n_versions, n_transcripts in sorted(dist.items()):
        rows.append({
            "sample": sample,
            "mutated_versions_per_transcript": n_versions,
            "n_transcripts": n_transcripts,
        })

    return pd.DataFrame(rows)


# =========================
# MAIN
# =========================

def main():
    print("[1] Chargement des biotypes du GTF...")
    tx_to_biotype = load_gtf_biotypes(GTF)

    all_summary_rows = []
    all_biotype_dfs = []
    all_version_dfs = []

    for sample, info in SAMPLES.items():
        print(f"\n[2] Analyse de {sample} ({info['srr']})")

        fasta_summary, biotype_counters, mut_versions = summarize_fasta(
            info["fasta"],
            tx_to_biotype,
        )

        log_summary = summarize_logs(info["log_dir"])

        row = {
            "sample": sample,
            "srr": info["srr"],
            **fasta_summary,
            **log_summary,
        }

        all_summary_rows.append(row)

        for counter_type, counter in biotype_counters.items():
            all_biotype_dfs.append(counter_to_df(counter, sample, counter_type))

        all_version_dfs.append(mut_versions_to_df(mut_versions, sample))

    summary_df = pd.DataFrame(all_summary_rows)
    biotype_df = pd.concat(all_biotype_dfs, ignore_index=True)
    version_df = pd.concat(all_version_dfs, ignore_index=True)

    summary_path = OUTDIR / "summary_transcript_reconstruction.tsv"
    biotype_path = OUTDIR / "biotypes_transcript_reconstruction.tsv"
    version_path = OUTDIR / "mutated_versions_per_transcript.tsv"

    summary_df.to_csv(summary_path, sep="\t", index=False)
    biotype_df.to_csv(biotype_path, sep="\t", index=False)
    version_df.to_csv(version_path, sep="\t", index=False)

    print("\n=== SUMMARY ===")
    print(summary_df.to_string(index=False))

    print("\n=== GROUPED BIOTYPES - MUTATED UNIQUE TRANSCRIPTS ===")

    grouped_mut_unique = biotype_df[
        biotype_df["counter_type"] == "grouped_mutated_unique_transcripts"
    ]

    for sample in SAMPLES:
        print(f"\n{sample}")
        print(
            grouped_mut_unique[grouped_mut_unique["sample"] == sample]
            .sort_values("n", ascending=False)
            .to_string(index=False)
        )

    print("\n=== GROUPED BIOTYPES - MUTATED ENTRIES ===")

    grouped_mut_entries = biotype_df[
        biotype_df["counter_type"] == "grouped_mutated_entries"
    ]

    for sample in SAMPLES:
        print(f"\n{sample}")
        print(
            grouped_mut_entries[grouped_mut_entries["sample"] == sample]
            .sort_values("n", ascending=False)
            .to_string(index=False)
        )

    print("\nFichiers générés :")
    print(summary_path)
    print(biotype_path)
    print(version_path)


if __name__ == "__main__":
    main()