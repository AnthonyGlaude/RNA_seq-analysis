import pandas as pd
import numpy as np
import pickle
from Bio import SeqIO
from intervaltree import IntervalTree


def extract_exons(gtf_file):
    exons = []
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'exon':
                chromosome = fields[0]
                start, end = int(fields[3]), int(fields[4])
                strand = fields[6]
                attributes = {
                    kv.split(" ")[0]: kv.split(" ")[1].strip('";')
                    for kv in fields[8].split("; ") if kv
                }
                gene_name = attributes.get("gene_name", attributes.get("gene_id", None))
                transcript_id = attributes.get("transcript_id", None)
                gene_biotype = attributes.get("gene_biotype", None)
                exons.append((chromosome, start, end, strand, gene_name, transcript_id, gene_biotype))
    exons_df = pd.DataFrame(
        exons,
        columns=['chromosome', 'start', 'end', 'strand', 'gene_name', 'transcript_id', 'gene_biotype']
    )
    exons_df['exon_number'] = (
        exons_df.sort_values(['transcript_id', 'start'])
                .groupby('transcript_id')
                .cumcount()
    )
    exons_df[['start', 'end']] = exons_df[['start', 'end']].astype(np.uint32)
    exons_df['exon_interval'] = pd.IntervalIndex.from_arrays(exons_df['start'], exons_df['end'], closed="both")
    return exons_df


def save_genome_dict(fasta_path, output_path):
    genome_dict = {
        rec.id.split('.')[0]: str(rec.seq)
        for rec in SeqIO.parse(fasta_path, "fasta")
    }
    with open(output_path, "wb") as f:
        pickle.dump(genome_dict, f)


def save_intervaltrees(exons_df, output_path):
    tree_dict = {}
    for exon in exons_df.itertuples():
        tree = tree_dict.setdefault(exon.chromosome, IntervalTree())
        start, end = min(exon.start, exon.end), max(exon.start, exon.end) + 1
        tree.addi(start, end, exon.transcript_id)
    with open(output_path, "wb") as f:
        pickle.dump(tree_dict, f)


def main():
    gtf_file = snakemake.input.gtf
    fasta_file = snakemake.input.transcripts_fasta
    output_parquet = snakemake.output.parquet
    output_genome_pkl = snakemake.output.genome_pkl
    output_tree_pkl = snakemake.output.tree_pkl

    exons_df = extract_exons(gtf_file)
    exons_df.to_parquet(output_parquet, index=False)

    save_genome_dict(fasta_file, output_genome_pkl)
    save_intervaltrees(exons_df, output_tree_pkl)


if __name__ == "__main__":
    main()