import pickle
import os
import sys
from collections import defaultdict

def build_exon_dict(gtf_file):
    exon_dict = defaultdict(list) 
    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'exon':
                chromosome = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8].split(';')
                gene_name = None
                transcript_id = None

                for attribute in attributes:
                    if 'gene_name' in attribute:
                        gene_name = attribute.split(' ')[-1].replace('"', '')
                    elif 'transcript_id' in attribute:
                        transcript_id = attribute.split(' ')[-1].replace('"', '')

                exon_dict[transcript_id].append({
                    'chromosome': chromosome,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_name': gene_name
                })
    return exon_dict

def dict_to_pickle(exon_dict, output_file):
    with open(output_file, "wb") as f:
        pickle.dump(exon_dict, f)
    
def main():
    gtf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/data/references/gtf/homo_sapiens.gtf"
    output_file = "../results/exon_dict.pkl"  # Définissez le nom du fichier de sortie

    exon_dict = build_exon_dict(gtf_file)
    dict_to_pickle(exon_dict, output_file)  # Passez le fichier de sortie ici

if __name__ == "__main__":
    main()
