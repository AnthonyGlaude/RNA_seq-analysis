import os
import sys
import datetime
import pandas as pd
import numpy as np
from itertools import chain, combinations
import concurrent.futures
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import IntervalTree  
import copy
import re
import pickle
from multiprocessing import Process, Queue
from pathlib import Path
from queue import Empty


# ------------------------------------------------------------
# constante globale
# ------------------------------------------------------------
# Variables globales partagées par les workers
global_genome = None
global_gtf_exons = None
global_trees = None
#Variable pour le nombre de combinaison mutée/transcrit (_all_subsets)
MAX_COMBINATION_SIZE = 10
MAX_WORKER=3
#variable pour debug_loglog
DEBUG_MODE = False
# ------------------------------------------------------------
# 1) Configuration du log : redirige stdout & stderr vers un fichier
# ------------------------------------------------------------
output_path = Path(snakemake.output.mutated_output_fasta)
log_path = output_path.with_suffix(".log")  # ex: 1_to_9_mutated_transcripts.log
log_fh = open(log_path, "w", buffering=1)

def log(msg: str):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_fh.write(f"{ts} - {msg}\n")
    log_fh.flush()

# Redirige stdout & stderr vers le log
sys.stdout = log_fh
sys.stderr = log_fh

def debug_log(msg: str):
    if DEBUG_MODE:
        log(f"[DEBUG] {msg}")

# ------------------------------------------------------------
# 2) Fonctions du pipeline
# ------------------------------------------------------------

def read_vcf(vcf_file: str):
    """
    Génère les mutations lues depuis un VCF simple (pas bgzippé).
    Chaque mutation est retournée comme un dictionnaire.
    """
    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith('#'):
                continue
            chrom, pos, _, ref, alt, *_ = line.strip().split('\t')
            for alt_index, a in enumerate(alt.split(',')):
                yield {
                    'chromosome': chrom,
                    'position': int(pos),
                    'ref_nucleotide': ref,
                    'alt_nucleotide': a,
                    'alt_index': alt_index
                }


def _parse_cds_range(description):
    m = re.search(r'CDS=(\d+)-(\d+)', description)
    if m:
        return tuple(map(int, m.groups()))
    return None

# Génère toutes les combinaisons non vides
def _all_subsets(iterable, max_size=None):
    s = list(iterable)
    max_len = len(s) if max_size is None else min(max_size, len(s))
    for r in range(1, max_len + 1):
        for combo in combinations(s, r):
            pos_seen = set()
            skip = False
            for m in combo:
                if m['position'] in pos_seen:
                    skip = True
                    break
                pos_seen.add(m['position'])
            if not skip:
                yield combo

def _apply_mutations_to_transcript(transcript_id, seq, mapping, muts, strand):
    """
    Applique les mutations à une séquence de transcript.
    REF/ALT sont ajustés si le transcrit est sur le brin négatif.
    """
    muts_in = [m.copy() for m in muts if m['position'] in mapping]

    applied_mutations = []
    for m in muts_in:
        m['transcript_pos'] = mapping[m['position']]
    muts_in.sort(key=lambda x: x['transcript_pos'])

    new_seq = seq
    shift = 0

    for m in muts_in:
        p = m['transcript_pos'] + shift
        ref = m['ref_nucleotide']
        alt = m['alt_nucleotide']

        adj_ref, adj_alt = _adjust_for_strand(ref, alt, strand)
        debug_log(f"[DEBUG] Adj REF/ALT for strand {strand}: {ref}/{alt} -> {adj_ref}/{adj_alt}")

        # Calcul de l’indice de départ
        if strand == '-':
            start_idx = p - len(adj_ref) + 1
        else:
            start_idx = p

        # Tronquage si mutation dépasse la fin
        max_len = len(new_seq) - start_idx
        adj_ref_trimmed = adj_ref[:max_len]
        adj_alt_trimmed = adj_alt[:max_len]

        seq_check = new_seq[start_idx : start_idx + len(adj_ref_trimmed)]
        debug_log(f"[DEBUG] seq[{start_idx}:{start_idx + len(adj_ref_trimmed)}] = {seq_check} vs REF attendue = {adj_ref_trimmed} (transcript {transcript_id})")

        if seq_check != adj_ref_trimmed:
            log(f"[WARN] Mutation ignorée pour {transcript_id} à pos transcript {p+1}: attendu {adj_ref_trimmed}, trouvé {seq_check}")
            continue

        if len(adj_ref_trimmed) < len(adj_ref):
            log(f"[INFO] Troncage de REF/ALT pour {transcript_id} à pos transcript {p+1}: {adj_ref}/{adj_alt} -> {adj_ref_trimmed}/{adj_alt_trimmed}")

        debug_log(f"[DEBUG] Applying: pos={p}, ref={adj_ref_trimmed}, alt={adj_alt_trimmed}")
        debug_log(f"[DEBUG] Avant mutation  [{transcript_id}]: {new_seq[start_idx-5:start_idx+len(adj_ref_trimmed)+5]}")

        if strand == '-':
            new_seq = new_seq[:start_idx] + adj_alt_trimmed + new_seq[p + 1:]
        else:
            new_seq = new_seq[:start_idx] + adj_alt_trimmed + new_seq[start_idx + len(adj_ref_trimmed):]

        debug_log(f"[DEBUG] Après mutation  [{transcript_id}]: {new_seq[start_idx-5:start_idx+len(adj_alt_trimmed)+5]}")

        shift += len(adj_alt_trimmed) - len(adj_ref_trimmed)

        # Met à jour REF/ALT pour le suivi
        m['ref_nucleotide'] = adj_ref_trimmed
        m['alt_nucleotide'] = adj_alt_trimmed

        applied_mutations.append(m)

    return new_seq, applied_mutations

# Traite un transcript dans un worker
def _process_transcript_worker(transcript_id, mutations):
    global global_genome, global_gtf_exons
    log(f"Début du traitement pour {transcript_id} ({len(mutations)} mutations)")

    exons = global_gtf_exons[global_gtf_exons['transcript_id'] == transcript_id]
    seq, mapping = _print_transcript_exons(transcript_id)
    if seq is None:
        return []
    out = []
    # WT : on garde le header complet avec description originale
    wt_rec = SeqRecord(
        Seq(seq),
        id=transcript_id,
        description=global_genome[transcript_id].description
    )
    out.append(wt_rec)
    strand = exons.iloc[0]['strand']
    cds_range = _parse_cds_range(global_genome[transcript_id].description)
    if len(mutations) > MAX_COMBINATION_SIZE:
        if cds_range:
            _, cds_end = cds_range
            mutations = [m for m in mutations if m['position'] in mapping and mapping[m['position']] + 1 <= cds_end]
            log(f"[FILTER] {transcript_id} contient >15 mutations, filtrées sur CDS (pos ≤ {cds_end}) -> {len(mutations)} restantes.")
        else:
            mutations = sorted(
                [m for m in mutations if m['position'] in mapping],
                key=lambda m: mapping[m['position']]
            )[:15]
            log(f"[FILTER] {transcript_id} contient >15 mutations, CDS absent -> 15 premières positions gardées.")

    for combo in _all_subsets(mutations, max_size=MAX_COMBINATION_SIZE):
        if any(m['position'] not in mapping for m in combo):
            continue

        ms, applied = _apply_mutations_to_transcript(transcript_id, seq, mapping, combo, strand)

        if len(applied) != len(combo):
            log(f"[INFO] Combo ignoré pour {transcript_id} : demandé {len(combo)} mutations, appliquées {len(applied)}")
            continue

        signature = "_".join(f"{m['ref_nucleotide']}({m['alt_nucleotide']})" for m in applied)
        positions = ";".join(str(mapping[m['position']] + 1) for m in applied)
        new_id = f"{transcript_id}_mut_{signature}_pos={positions}"

        #  Ne pas inclure de description sur les transcrits mutés
        rec = SeqRecord(
            ms if isinstance(ms, Seq) else Seq(ms),
            id=new_id,
            description=""
        )
        out.append(rec)
        for m in applied:
            cpos = mapping[m['position']]
            ref_len = len(m['ref_nucleotide'])
            if strand == '-':
                start_idx = cpos - ref_len + 1
                ref_codon = seq[start_idx:cpos + 1]
            else:
                ref_codon = seq[cpos:cpos + ref_len]
            log(f"[CHECK] {transcript_id} | GENOME pos {m['position']} -> cDNA pos {cpos} | REF VCF: {m['ref_nucleotide']} vs REF_SEQ: {ref_codon}")
    return out

def _worker_loop(task_queue, result_queue, fasta_path, exon_parquet_path, trees_pkl_path):
    """
    Fonction exécutée par chaque worker. Initialise les données locales et traite les transcripts à la chaîne.
    """
    _init_worker(fasta_path, exon_parquet_path, trees_pkl_path)
    while True:
        item = task_queue.get()
        if item is None:
            break  # Fin du traitement
        transcript_id, mutations = item
        try:
            recs = _process_transcript_worker(transcript_id, mutations)
            start_put = time.time()
            log(f"[DEBUG] {transcript_id} - trying to put {len(recs)} recs in result_queue")
            result_queue.put(recs)
            end_put = time.time()
            log(f"[DEBUG] {transcript_id} - put done in {end_put - start_put:.2f} sec")

        except Exception as e:
            result_queue.put([])
            log(f"[ERROR WORKER] {transcript_id} - {e}")


def _run_with_queue(tx_muts, num_workers, fasta_path, exon_parquet_path, trees_pkl_path):
    task_queue = Queue()
    result_queue = Queue()

    # Enfile tous les transcripts à traiter
    for tx, muts in tx_muts.items():
        task_queue.put((tx, muts))
    for _ in range(num_workers):
        task_queue.put(None)  # Sentinelle pour arrêt

    # Lance les workers
    workers = []
    for _ in range(num_workers):
        p = Process(target=_worker_loop, args=(task_queue, result_queue, fasta_path, exon_parquet_path, trees_pkl_path))
        p.start()
        workers.append(p)
    total = len(tx_muts)
    received = 0

    # with tqdm(total=total, desc="Traitement des transcrits", file=sys.__stdout__) as pbar:
    while received < total:
        try:
            recs = result_queue.get(timeout=1)  # récupère dès que dispo
            received += 1
            for r in recs:
                yield r
            # pbar.update(1)
        except Empty:
            continue  # aucun résultat pour l'instant, on réessaie
    # Attendre la fin des workers
    for p in workers:
        p.join()


def _init_worker(fasta_path: str, exon_parquet_path: str, trees_pkl_path="trees.pkl"):
    """Initialise chaque worker : charge le génome, le DataFrame d'exons et construit ou charge les IntervalTrees"""
    global global_genome, global_gtf_exons, global_trees

    # 1) Genome en mémoire
    global_genome = {
        rec.id.split('.')[0]: rec
        for rec in SeqIO.parse(fasta_path, "fasta")
    }
    # 2) DataFrame d'exons
    global_gtf_exons = pd.read_parquet(exon_parquet_path)
    # 3) Chargement ou construction des IntervalTrees
    if os.path.exists(trees_pkl_path):
        with open(trees_pkl_path, "rb") as f:
            global_trees = pickle.load(f)
        log(f"[INIT] IntervalTrees chargés depuis {trees_pkl_path}")
    else:
        global_trees = {}
        for exon in global_gtf_exons.itertuples():
            chrom = exon.chromosome
            start, end = exon.start, exon.end
            if start > end:
                start, end = end, start
            end += 1  # intervaltree half-open
            tree = global_trees.setdefault(chrom, IntervalTree())
            tree.addi(start, end, exon.transcript_id)

        with open(trees_pkl_path, "wb") as f:
            pickle.dump(global_trees, f)
        log(f"[BUILD] IntervalTree construit et sauvegardé dans {trees_pkl_path}")
    # Infos debug_logour debug
    log(f"[INFO] IntervalTrees disponibles pour chromosomes : {list(global_trees.keys())}")
    chrom = '1'
    df1 = global_gtf_exons[global_gtf_exons["chromosome"] == chrom] #### chrom = '1' retiré
    debug_log(f"[DEBUG] Exons chr{chrom}: start min={df1['start'].min()}, max={df1['start'].max()}")


# Récupère et assemble les exons d'un transcript
def _print_transcript_exons(transcript_id):
    global global_genome, global_gtf_exons

    exons = global_gtf_exons[global_gtf_exons['transcript_id'] == transcript_id]
    if exons.empty or transcript_id not in global_genome:
        log(f"[ERROR] Transcript {transcript_id} non trouvé.")
        return None, None
    strand = exons.iloc[0]['strand']
    sorted_exons = exons.sort_values('start', ascending=(strand == '+'))
    seq = str(global_genome[transcript_id].seq)
    mapping = {}
    offset = 0
    all_entries = []
    for _, exon in sorted_exons.iterrows():
        start, end = exon.start, exon.end
        length = end - start + 1
        coords = (
            range(start, end + 1) if strand == '+' else
            range(end, start - 1, -1)
        )
        for i, gpos in enumerate(coords):
            cpos = offset + i
            mapping[gpos] = cpos
            base = seq[cpos] if cpos < len(seq) else '?'
            all_entries.append((gpos, cpos, base))
        offset += length

    # Log premier et dernier nucléotide
    log(f"Transcript: {transcript_id} (strand = {strand})")
    log(f"  FIRST:  GENOME {all_entries[0][0]} -> cDNA {all_entries[0][1]} -> BASE {all_entries[0][2]}")
    log(f"  LAST :  GENOME {all_entries[-1][0]} -> cDNA {all_entries[-1][1]} -> BASE {all_entries[-1][2]}")
    log(f"  TOTAL LENGTH: {len(seq)}")

    return seq, mapping

# Applique les mutations sur la séquence du transcript
def _adjust_for_strand(ref, alt, strand):
    if strand == "-":
        return str(Seq(ref).reverse_complement()), str(Seq(alt).reverse_complement())
    return ref, alt

def _generate_mutated_transcripts_parallel(
    fasta, exon_parquet, max_workers=MAX_WORKER, trees_pkl="trees.pkl",vcf_file=None
):
    log("2) Chargement des exons pour le main")
    _init_worker(fasta, exon_parquet, trees_pkl)
    for chrom, tree in global_trees.items():
        intervals = list(tree)
        if intervals:
            min_start = min(iv.begin for iv in intervals)
            max_end = max(iv.end for iv in intervals)
            debug_log(f"IntervalTree {chrom}: {len(intervals)} intervalles, start min={min_start}, end max={max_end}")
        else:
            log(f"IntervalTree {chrom}: AUCUN intervalle.")
    log(f"   -> {len(global_gtf_exons)} exons chargés.")
    tx_muts = {}
    for mut in read_vcf(vcf_file):
        chrom = mut['chromosome']
        pos   = mut['position']
        if chrom not in global_trees:
            continue
        for iv in global_trees[chrom][pos]:
            tx = iv.data
            tx_muts.setdefault(tx, []).append(mut)
    for tx, muts in tx_muts.items():
        log(f"--> {tx} a {len(muts)} mutations (positions uniques : {len({m['position'] for m in muts})})")
    total = len(tx_muts)
    log(f"3) Mapping construit -> {total} transcripts à traiter (≤15 mutations chacun).")
    log("4) Lancement des workers avec Queue…")

    return _run_with_queue(tx_muts, max_workers, fasta, exon_parquet, trees_pkl)


# Ecrit les transcrits mutés par batch pour économiser la mémoire
def save_mutated_transcripts_in_batches_parallel(
    fasta, exon_parquet, out_fasta,
    batch_size=5, max_workers=MAX_WORKER, trees_pkl="trees.pkl",vcf_file=None
):
    batch = []
    with open(out_fasta, "w") as out:
        for rec in _generate_mutated_transcripts_parallel(
            fasta, exon_parquet, max_workers, trees_pkl,vcf_file
        ):
            batch.append(rec)
            if len(batch) >= batch_size:
                SeqIO.write(batch, out, "fasta")
                out.flush(); os.fsync(out.fileno())
                batch.clear()
        if batch:
            SeqIO.write(batch, out, "fasta")
            out.flush()
        os.fsync(out.fileno())


def filter_transcripts(input_fasta, output_fasta):
    """
    Parcourt le FASTA d'entrée et regroupe les séquences par identifiant de transcript.
    Ne conserve que les groupes qui contiennent à la fois :
      - Un transcript WT (pas de "_mut_" dans l'ID)
      - Au moins un transcript muté (ID contenant "_mut_")
    """
    groups = {}
    total_records = 0
    for record in SeqIO.parse(input_fasta, "fasta"):
        total_records += 1
        m = re.search(r"(ENST\d+)", record.id)
        if not m:
            continue
        tid = m.group(1)
        groups.setdefault(tid, []).append(record)
    filtered = []
    kept = 0
    for tid, recs in groups.items():
        has_wt  = any("_mut_" not in r.id for r in recs)
        has_mut = any("_mut_" in r.id     for r in recs)
        if has_wt and has_mut:
            filtered.extend(recs)
            kept += 1
            
    SeqIO.write(filtered, output_fasta, "fasta")
    print(f" Total lus : {total_records}")
    print(f" Groupes gardés (WT+mut) : {kept}")
    print(f" Total écrits : {len(filtered)}")

# ------------------------------------------------------------
# 3) Fonction main()
# ------------------------------------------------------------
def main():
    log("=== DÉMARRAGE DU PIPELINE ===")
# ------------------------------------------------------------
#  Input Paths
# ------------------------------------------------------------

    vcf_file = snakemake.input.vcf
    transcripts_fasta = snakemake.input.transcripts_fasta
    mutated_output_fa = snakemake.output.mutated_output_fasta
    exon_parquet = snakemake.input.exon_parquet
    genome_pkl = snakemake.input.genome_pkl
    trees_pkl = snakemake.input.tree_pkl
    
    
#    vcf_file = "20QC_variant_tronque.vcf"
#    transcripts_fasta= "breast_cancer/pickle/transcriptome.fa"
#    mutated_output_fa  = "mutated_transcripts.fa"
#    exon_parquet = "exon_data.parquet"  
#    trees_pkl = "trees.pkl"

# ------------------------------------------------------------
# output Path + constante
# ------------------------------------------------------------
    final_fasta = "final_transcriptome.fa"
    
    
# ------------------------------------------------------------
# 3) logic :)
# ------------------------------------------------------------

    try:
#        df = read_vcf(vcf_file)
        save_mutated_transcripts_in_batches_parallel(
            fasta=transcripts_fasta,
           # vcf_df=df,
            exon_parquet=exon_parquet,
            out_fasta=mutated_output_fa,
            batch_size=25,
            max_workers=MAX_WORKER,
            trees_pkl=trees_pkl,
            vcf_file=vcf_file
        )
        filter_transcripts(mutated_output_fa, final_fasta)
        log(f"=== FASTA filtré écrit dans {final_fasta} ===")
        log("=== PIPELINE TERMINÉ AVEC SUCCÈS ===")
    except Exception as e:
        log(f"!!! ERREUR CRITIQUE : {e}")
    finally:
        log_fh.close()

if __name__ == "__main__":
    main()
