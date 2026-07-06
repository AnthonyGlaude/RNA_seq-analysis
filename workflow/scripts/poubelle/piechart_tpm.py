import pandas as pd
import matplotlib.pyplot as plt

# Charger le fichier VCF avec la fonction
def load_vcf_with_pandas(vcf_file):
    """
    Charge les données VCF à partir du fichier en utilisant pandas
    """
    vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
    vcf_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [f'Sample{i+1}' for i in range(vcf_df.shape[1] - 9)]
    return vcf_df

# Charger les données VCF
vcf_file = "F:/breast_cancer/workflow/results/fraction/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant.vcf"  # Remplacer par ton fichier VCF
vcf_df = load_vcf_with_pandas(vcf_file)

# Exemple : Supposons que les expressions de chaque variant sont dans les colonnes Sample1, Sample2, etc.
# Et que nous avons une fonction ou une règle pour calculer les TPM.
# Calculer les TPM (remplace par ta propre méthode si nécessaire)
# Ici, on fait simplement une moyenne des valeurs des échantillons comme exemple
vcf_df['TPM'] = vcf_df.iloc[:, 9:].mean(axis=1)  # Moyenne des valeurs des échantillons pour chaque variant

# Filtrer les variants avec TPM > 10
filtered_vcf_df = vcf_df[vcf_df['TPM'] > 10]

# Segmenter en fonction des catégories (Feature) - On suppose qu'il existe une colonne 'Feature' ou que tu as besoin de la créer
# Exemple de valeurs possibles : 'transcript', 'exon', 'intergénique', etc.
# Pour l'exemple, imaginons une colonne 'Feature' déjà définie dans ton dataframe
# Sinon, il faut l'extraire à partir du fichier GTF comme tu l'as fait précédemment

# Appliquer un processus pour segmenter les données (par exemple, transcrit, exon, intergénique)
# Ici, je vais supposer que tu as une colonne 'Feature' dans le VCF, mais si tu n'as pas cette information, 
# il faudra la rajouter en fonction de l'alignement avec le fichier GTF (comme fait dans ton code précédent).

# Exemple de création d'une colonne 'Feature' fictive pour démonstration (à adapter à ton cas)
filtered_vcf_df['Feature'] = 'transcript'  # Ceci est un exemple. Remplace-le par ta propre logique

# Comptage des catégories d'expression
feature_counts = filtered_vcf_df['Feature'].value_counts()

# Liste de couleurs spécifiques pour chaque catégorie
colors = ['#548235', '#A3D08E', '#FF6347', '#FFD700', '#98FB98']  # Exemples de couleurs

# Tracer le pie chart
plt.figure(figsize=(8, 6), dpi=600)
feature_counts.plot(kind='pie', autopct='%1.1f%%', startangle=90, cmap='Set3', colors=colors, legend=False)

plt.ylabel('')
plt.title("Distribution des variants avec TPM > 10 par catégorie d'expression")
plt.show()
plt.savefig("piechart_tpm_filtered.png", dpi=600)

print(feature_counts)  # Affiche les comptages des différentes catégories
