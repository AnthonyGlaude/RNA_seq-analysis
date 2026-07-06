import os
import pandas as pd
import matplotlib.pyplot as plt

# Chemins des fichiers
gtf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/data/references/gtf/homo_sapiens.gtf"
vcf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/results/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant.vcf"
gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, low_memory=False)
vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None, low_memory=False)

gtf_df.columns = ['CHROM', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NEW_COLUMN']

# Fonction pour obtenir les attributs comme 'gene_id' et 'transcript_id'
def extract_attributes(attribute_str):
    attributes = {}
    for attr in attribute_str.split(';'):
        if attr.strip():
            key, value = attr.split(' ', 1)
            key = key.strip().lower()
            value = value.strip().strip('"')
            attributes[key] = value
    return attributes

def get_feature_for_variant(row):
    chrom = row['#CHROM']
    pos = row['POS']

    # Filtrer le DataFrame GTF pour trouver les fonctionnalités correspondantes
    matching_gtf = gtf_df[(gtf_df['CHROM'] == chrom) & (gtf_df['start'] <= pos) & (gtf_df['end'] >= pos)]

    if not matching_gtf.empty:
        feature_counts = matching_gtf['feature'].value_counts()
        
        # Liste des fonctionnalités prioritaires
        priority_features = ['transcript', '5-utr', '3-utr', 'exon', 'gene']

        # Vérifier les fonctionnalités en fonction de l'ordre de priorité
        for feature in priority_features:
            if feature in feature_counts.index:
                return feature  # Retourne la première fonctionnalité prioritaire rencontrée

        # Gérer les fonctionnalités non prioritaires
        all_features = feature_counts.index.tolist()
        other_features = [feature for feature in all_features if feature not in priority_features]

        if other_features:
            return other_features  # Retourne les fonctionnalités autres trouvées

        # Si aucune fonctionnalité prioritaire ou autre n'est trouvée, retourner la plus fréquente
        return feature_counts.idxmax() if not feature_counts.empty else "Aucune fonctionnalité"

    return "Région Intragénique"  # Si aucune correspondance

def main():
    gtf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/data/references/gtf/homo_sapiens.gtf"
    global gtf_df
    gtf_df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
    gtf_df.columns = ['CHROM', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    
    variants_folder = "/mnt/c/Users/Antho/Documents/nouveau_projet/workflow/results/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/"
    total_counts = {"transcript": 0, "exon": 0, "gene": 0, "5-utr": 0, "3-utr": 0, "Région Intragénique": 0}  # Dictionnaire pour compter les fonctionnalités

    for patient_folder in os.listdir(variants_folder):
        patient_path = os.path.join(variants_folder, patient_folder)
        vcf_file = os.path.join(patient_path, "20QC_variant.vcf")

        if os.path.isfile(vcf_file):
            vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
            vcf_df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NEW_COLUMN']

            for index, row in vcf_df.iterrows():
                feature = get_feature_for_variant(row)
                if isinstance(feature, list):
                    # Compter chaque fonctionnalité non prioritaire
                    for f in feature:
                        total_counts[f] = total_counts.get(f, 0) + 1  # Utiliser get pour éviter une KeyError
                else:
                    total_counts[feature] += 1  # Incrémentation des comptages

    # Créer le pie chart
    plt.figure(figsize=(8, 6), dpi=600)
    plt.pie(total_counts.values(), labels=total_counts.keys(), autopct='%1.1f%%', startangle=90)
    plt.ylabel('')
    plt.title(f"Distribution des variants par localisation intragénique (N={len(os.listdir(variants_folder))})")
    plt.savefig("piechart_high_res.png", dpi=600)
    plt.show()

    # Afficher les comptages totaux
    print(total_counts)

if __name__ == "__main__":
    main()
