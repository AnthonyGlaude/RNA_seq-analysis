
import re
from Bio import SeqIO

def filter_transcripts(input_fasta, output_fasta):
    """
    Parcourt le FASTA d'entrée et regroupe les séquences par identifiant de transcript.
    Ne conserve que les groupes qui contiennent à la fois :
      - Un transcrit de référence (wild-type, c'est-à-dire dont le header ne contient pas "_mut_")
      - Au moins un transcrit muté (header contenant "_mut_")
    """
    groups = {}
    # Parcours du FASTA d'entrée
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Extraction de l'identifiant transcript (ex : ENST00000372748)
        match = re.search(r"(ENST\d+)", record.id)
        if not match:
            continue  # On ignore les séquences dont le header n'a pas le format attendu
        transcript_id = match.group(1)
        groups.setdefault(transcript_id, []).append(record)
    
    filtered_records = []
    # Pour chaque groupe, on vérifie qu'il contient à la fois une séquence wild-type et une séquence mutée.
    for transcript_id, records in groups.items():
        has_wildtype = any("_mut_" not in rec.id for rec in records)
        has_mutated = any("_mut_" in rec.id for rec in records)
        if has_wildtype and has_mutated:
            filtered_records.extend(records)
    
    # Écriture du FASTA filtré dans le fichier de sortie
    SeqIO.write(filtered_records, output_fasta, "fasta")

def main():
    # Récupération des chemins depuis snakemake
    combined_fasta = snakemake.input.combined_fasta 
    filtered_fasta = snakemake.output.filtered_fasta 

    filter_transcripts(combined_fasta, filtered_fasta)

if __name__ == "__main__":
    main()