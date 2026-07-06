import sys
import pandas as pd
from Bio import SeqIO

fasta_in = sys.argv[1]
gtf_parquet = sys.argv[2]
output_dir = sys.argv[3]

# Charger exon dataframe
exon_df = pd.read_parquet(gtf_parquet)

# Créer un mapping transcript_id → chromosome
transcript_chrom = exon_df.drop_duplicates("transcript_id")[["transcript_id", "chromosome"]]
transcript_chrom = transcript_chrom.rename(columns={"chromosome": "chrom"})

# Définir les groupes de chromosomes
GROUPS = {
    "1_to_4": [str(c) for c in range(1, 5)],
    "5_to_9": [str(c) for c in range(5, 10)],
    "10_to_16": [str(c) for c in range(10, 17)],
}
excluded = set(GROUPS["1_to_4"] + GROUPS["5_to_9"] + GROUPS["10_to_16"])
all_chroms = set(transcript_chrom["chrom"].unique())
GROUPS["rest"] = sorted(all_chroms - excluded)

# Créer un dict transcript_id → group
tx_to_group = {}
for group, chroms in GROUPS.items():
    tx_subset = transcript_chrom[transcript_chrom["chrom"].isin(chroms)]
    for tx in tx_subset["transcript_id"]:
        tx_to_group[tx] = group

# Lire le FASTA et répartir directement dans le bon groupe
records_by_group = {group: [] for group in GROUPS}
for record in SeqIO.parse(fasta_in, "fasta"):
    tx_id = record.id.split()[0]
    group = tx_to_group.get(tx_id)
    if group:
        records_by_group[group].append(record)

# Écrire les fichiers
for group, records in records_by_group.items():
    output_path = f"{output_dir}/{group}_transcripts.fa"
    with open(output_path, "w") as out_f:
        SeqIO.write(records, out_f, "fasta")