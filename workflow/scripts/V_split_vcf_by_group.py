
import os
import gzip

vcf_file = snakemake.input.vcf
output_path = snakemake.output.vcf_split
chromosomes = set(snakemake.params.chromosomes)
include_rest = snakemake.params.get("include_rest", False)

# Créer le dossier de sortie si besoin
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# Gestion .vcf vs .vcf.gz
def open_vcf(path):
    return gzip.open(path, 'rt') if path.endswith(".gz") else open(path, 'r')

with open_vcf(vcf_file) as vcf_in, open(output_path, 'w') as vcf_out:
    for line in vcf_in:
        if line.startswith('#'):
            vcf_out.write(line)
        else:
            chrom = line.split('\t')[0]
            if include_rest:
                if chrom not in chromosomes:
                    vcf_out.write(line)
            else:
                if chrom in chromosomes:
                    vcf_out.write(line)