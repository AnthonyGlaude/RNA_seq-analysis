import pandas as pd

def main():
    # Chemins des fichiers
    kallisto_file = "F:/breast_cancer/workflow/results/dge/kallisto/171992_SIMG0590_T_totalRNA_sarcoma_43378_S9_L002/abundance_filtered.tsv"
    gtf_file = "F:/breast_cancer/workflow/data/references/gtf/homo_sapiens.gtf"
    output_bedgraph = "F:/breast_cancer/workflow/final_data/top_gene_abundance.bedgraph"
    igv_session_file = "F:/breast_cancer/workflow/final_data/session.igv"

    # Charger les résultats de Kallisto
    kallisto_data = pd.read_csv(kallisto_file, sep='\t')

    # Identifier le transcript avec la plus grande abondance
    top_transcript = kallisto_data.loc[kallisto_data['tpm'].idxmax()]
    transcript_id = top_transcript['target_id']
    abundance = top_transcript['tpm']

    print(f"Transcript le plus abondant : {transcript_id} avec une TPM de {abundance}")

    # Extraire les coordonnées du transcript à partir du GTF
    coordinates = []

    with open(gtf_file, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):  # Ignorer les commentaires
                continue
            fields = line.strip().split('\t')
            if fields[2] == "transcript":  # Filtrer pour les transcrits
                attributes = fields[8]
                if f'transcript_id "{transcript_id}"' in attributes:
                    chrom = fields[0]
                    start = int(fields[3]) - 1  # Convertir en 0-based pour le format BedGraph
                    end = int(fields[4])
                    strand = fields[6]  # La chaîne n'est pas utilisée dans le BedGraph
                    coordinates.append((chrom, start, end, strand))
                    print(chrom, 'sssssssssssssssssss')
                    break

    # Vérifier si les coordonnées ont été trouvées
    if not coordinates:
        print(f"Aucune coordonnée trouvée pour le transcript {transcript_id} dans le GTF.")
    else:
        # Générer le fichier BedGraph
        with open(output_bedgraph, 'w') as bg:
            bg.write("track type=bedGraph name=\"Top Gene Abundance\"\n")
            for chrom, start, end, _ in coordinates:  # _ pour ignorer le strand
                bg.write(f"{chrom}\t{start}\t{end}\t{abundance}\n")
        print(f"Fichier BedGraph créé : {output_bedgraph}")

        # Générer le fichier de session IGV
        genome = "hg38"  # Changez cela si vous utilisez un autre génome
        with open(igv_session_file, 'w') as igv:
            igv.write(f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<session genome="{genome}" hasGeneTrack="true" hasSequenceTrack="true" version="8">
    <Resources>
        <Resource path="{output_bedgraph}"/>
    </Resources>
    <Panel>
        <Track id="{output_bedgraph}" name="Top Gene Abundance" color="0,0,255" height="50"/>
    </Panel>
</session>
""")
            # Ajouter la ligne goto correctement formatée
            chrom = coordinates[0][0]  # Obtenir le chromosome
            start = coordinates[0][1] + 1  # Convertir en 1-based pour IGV
            end = coordinates[0][2]
            print(chrom, "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee")
            # Vérifier et corriger le nom du chromosome
            if chrom.startswith('chr'):
                igv.write(f"goto {chrom}:{start}-{end}\n")  # Écrire la ligne goto correctement
            else:
                igv.write(f"goto chr{chrom}:{start}-{end}\n")  # Ajouter 'chr' si nécessaire

            igv.write("sortByFile\n")  # Ajouter la commande pour trier par fichier

        print(f"Fichier IGV créé : {igv_session_file}")

if __name__ == "__main__":
    main()
