import pandas as pd
import matplotlib.pyplot as plt
import io

# Charger le fichier VCF (ignorer les lignes de commentaires commençant par ##)
vcf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/results/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant.vcf"
with open(vcf_file, 'r') as f:
    lines = [line for line in f if not line.startswith('##')]

# Charger le VCF dans un DataFrame (ici on conserve uniquement la colonne 'POS')
df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')

# Renommer les colonnes pour plus de clarté
df.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

# Charger le fichier GTF
gtf_file = "/mnt/c/Users/Antho/Documents/breast_cancer/breast_cancer/workflow/data/references/gtf/homo_sapiens.gtf"
gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)
gtf.columns = ['CHROM', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# Extraire les annotations du champ INFO du fichier VCF
annotations = []
for info in df["INFO"]:
    ann_field = [field for field in info.split(';') if field.startswith('ANN=')]
    if ann_field:
        ann_value = ann_field[0].split('=')[1]
        annotations.append(ann_value)
    else:
        annotations.append("")

# Extraire les types de mutations
mutation_types = []
for ann in annotations:
    if ann:
        mutation = ann.split('|')[1]  # Le type de région est dans la 2e colonne après '|'
        mutation_types.append(mutation)
    else:
        mutation_types.append("Unknown")

# Ajouter les types de mutations au DataFrame
df["Mutation_Type"] = mutation_types
pd.set_option('display.max_rows', None)  # Afficher toutes les lignes
pd.set_option('display.max_columns', None)  # Afficher toutes les colonnes
pd.set_option('display.width', None)  # Pas de troncature de largeur
pd.set_option('display.max_colwidth', None)  # Pas de troncature de largeur des colonnes

# Filtrer les gènes dans le fichier GTF
genes = gtf[gtf['feature'] == 'gene']

# Associer les annotations GTF aux variants VCF (en fonction de la position)
# Ici, on suppose qu'on veut ajouter le nom du gène pour chaque variante
def get_gene_name(row):
    chrom = row['CHROM']
    pos = row['POS']
    
    # Filtrer les gènes correspondant à la même position (ou dans une plage proche)
    gene_info = genes[(genes['CHROM'] == chrom) & (genes['start'] <= pos) & (genes['end'] >= pos)]
    
    # Si un gène est trouvé, on retourne son nom, sinon "Unknown"
    if not gene_info.empty:
        attributes = gene_info['attribute'].iloc[0]
        gene_name = [x for x in attributes.split(';') if x.startswith('gene_name')]
        return gene_name[0].split('=')[1] if gene_name else "Unknown"
    return "Unknown"

# Ajouter une colonne 'Gene' avec les noms de gènes associés aux variantes
df['Gene'] = df.apply(get_gene_name, axis=1)

# Compter les occurrences des types de mutation
mutation_counts = df["Mutation_Type"].value_counts()

# Créer un graphique en camembert
plt.figure(figsize=(8, 8))
mutation_counts.plot.pie(autopct='%1.1f%%', startangle=90, colormap='tab10')
plt.title("Répartition des mutations par type de région")
plt.ylabel("")  # Supprimer l'étiquette de l'axe Y

# Afficher le graphique
plt.show()
