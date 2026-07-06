
import pandas as pd
from gtfparse import read_gtf
import os

tx2gene = snakemake.input.tx2gene
gtf = snakemake.input.gtf
outfile = snakemake.output.tpm
log_file_path = snakemake.log[0]

import pandas as pd
from Bio import SeqIO

# Charger le fichier VCF
vcf_file =  snakemake.output.call_variants
vcf_data = pd.read_csv(vcf_file, sep='\t', comment='#', header=None, 
                       names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])

# Charger les séquences de transcrits à partir du fichier FASTA
fasta_file =  snakemake.output.build_transcriptome
transcripts = {rec.id: rec.seq for rec in SeqIO.parse(fasta_file, "fasta")}

# Parcourir les variants et appliquer les modifications aux séquences
for idx, variant in vcf_data.iterrows():
    chrom = variant['CHROM']
    pos = int(variant['POS']) - 1  # L'index dans Python commence à 0
    ref = variant['REF']
    alt = variant['ALT']

    # Identifier le transcrit concerné par le variant
    # Ici, il faut une logique qui associe chaque variant à un transcrit
    # Supposons que vous ayez une fonction 'get_transcript_from_position' qui fait cette association
    transcript_id = get_transcript_from_position(chrom, pos)
    if transcript_id:
        transcript_seq = transcripts[transcript_id]

        # Vérifier que la base à la position correspond à REF
        if transcript_seq[pos] == ref:
            # Créer la séquence modifiée avec ALT
            modified_seq = transcript_seq[:pos] + alt + transcript_seq[pos + 1:]

            # Stocker la nouvelle séquence avec un nom modifié
            new_transcript_id = f"{transcript_id}_variant_chr{chrom}_{pos+1}"
            transcripts[new_transcript_id] = modified_seq
        else:
            print(f"Mismatch: expected {ref} at position {pos}, but found {transcript_seq[pos]}")

# Écrire les séquences modifiées dans un nouveau fichier FASTA
with open("modified_transcripts.fasta", "w") as output_fasta:
    for transcript_id, seq in transcripts.items():
        output_fasta.write(f">{transcript_id}\n{str(seq)}\n")
