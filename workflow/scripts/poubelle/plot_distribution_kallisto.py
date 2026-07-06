import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def load_vcf_with_pandas(vcf_file):
    """
    Charge les données VCF à partir du fichier en utilisant pandas
    """
    vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
    vcf_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + [f'Sample{i+1}' for i in range(vcf_df.shape[1] - 9)]
    return vcf_df

def extract_variant_type(info_column):
    """
    Extrait le type de variant à partir de la colonne 'INFO' du fichier VCF.
    """
    info_str = info_column
    type_info = [entry.split('=')[1] for entry in info_str.split(';') if entry.startswith('TYPE=')]
    return type_info[0] if type_info else 'Autre'

def plot_variant_types_with_legend(variants_df, output_plot_path):
    """
    Crée un pie chart des différents types de variants détectés avec une légende.
    """
    variants_df['variant_type'] = variants_df['INFO'].apply(extract_variant_type) 
    variant_types = variants_df['variant_type'].value_counts()
    fig, ax = plt.subplots(figsize=(7, 7)) 
    wedges, _ = ax.pie(variant_types, 
                       startangle=90, 
                       colors=sns.color_palette("Set1", len(variant_types)),
                       wedgeprops={'edgecolor': 'black', 'linewidth': 0.5})
    legend_labels = [f"{variant_types.index[i]}: {variant_types.iloc[i]} ({variant_types.iloc[i] / variant_types.sum() * 100:.1f}%)" 
                     for i in range(len(variant_types))]
    ax.legend(wedges, legend_labels, title="Types de Variants", loc="center left", bbox_to_anchor=(1, 0, 0.5, 1),
              fontsize=12, frameon=False) 
    ax.set_title("Types de Variants Détectés", fontsize=16)
    plt.savefig(output_plot_path, dpi=720)  
    plt.close()

def plot_abundance_distribution(abundance_data, output_plot_path):
    """
    Crée un graphique de distribution des abondances.
    """
    plt.figure(figsize=(10, 6))
    sns.histplot(abundance_data, kde=True, color='blue', bins=30)
    plt.title("Distribution des Abondances des Variants", fontsize=16)
    plt.xlabel("Abondance", fontsize=12)
    plt.ylabel("Fréquence", fontsize=12)
    plt.savefig(output_plot_path, dpi=720) 
    plt.close()

def plot_variant_types_piechart(variants_df, output_plot_path):
    """
    Crée un pie chart des types de variants détectés, incluant les labels avec ajustement.
    """
    variants_df['variant_type'] = variants_df['INFO'].apply(extract_variant_type)
    variant_types = variants_df['variant_type'].value_counts()
    fig, ax = plt.subplots(figsize=(8, 8))  
    wedges, texts, autotexts = ax.pie(variant_types, 
                                      labels=variant_types.index,
                                      autopct='%1.1f%%',
                                      startangle=90, 
                                      colors=sns.color_palette("Set1", len(variant_types)),
                                      wedgeprops={'edgecolor': 'black', 'linewidth': 0.5})
    
    for text in texts:
        text.set_fontsize(10)
        text.set_horizontalalignment('center')
    for autotext in autotexts:
        autotext.set_fontsize(10)
        autotext.set_horizontalalignment('center')
    
    ax.set_title("Types de Variants Détectés", fontsize=16)
    ax.legend(wedges, [f"{variant_types.index[i]}: {variant_types.iloc[i]} ({variant_types.iloc[i] / variant_types.sum() * 100:.1f}%)" 
                       for i in range(len(variant_types))],
              title="Types de Variants", loc="center left", bbox_to_anchor=(1.05, 0.5),
              fontsize=12, frameon=False)  
    plt.tight_layout()
    plt.savefig(output_plot_path, dpi=720) 
    plt.close()

def main():
    #  Path pour le fichier VCF a utiliser
    variants_file = "F:/breast_cancer/workflow/results/fraction/variants/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/20QC_variant.vcf"
    
    # Charger les données VCF dans un DataFrame
    vcf_df = load_vcf_with_pandas(variants_file)
    
    # Path de sortie 
    base_output_path = "F:/breast_cancer/workflow/results/fraction/dge/kallisto/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/"
    output_plot_path_variant_types_legend = f"{base_output_path}variant_types_with_legend.png"
    output_plot_path_abundance = f"{base_output_path}abundance_distribution.png"
    output_plot_path_variant_types_piechart = f"{base_output_path}variant_types_piechart.png"
    
    # Générer le pie chart des types de variants avec légende
    plot_variant_types_with_legend(vcf_df, output_plot_path_variant_types_legend)

    # Simuler des données d'abondance pour le graphique de distribution
    abundance_data = np.random.uniform(0, 100, size=1000)
    
    # Générer le graphique de distribution des abondances
    plot_abundance_distribution(abundance_data, output_plot_path_abundance)

    # Générer le pie chart des types de variants
    plot_variant_types_piechart(vcf_df, output_plot_path_variant_types_piechart)

    # Suivi troubleshooting et savoir ou sont enregistré les graphiques
    print(f"Graphiques générés et enregistrés sous :")
    print(f"- Types de variants avec légende : {output_plot_path_variant_types_legend}")
    print(f"- Distribution des abondances : {output_plot_path_abundance}")
    print(f"- Types de variants (pie chart détaillé) : {output_plot_path_variant_types_piechart}")

if __name__ == "__main__":
    main()
