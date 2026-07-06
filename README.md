# RNA-seq analysis workflow

Ce dépôt contient un pipeline Snakemake pour l’analyse de données RNA-seq. Le pipeline automatise les étapes principales du traitement des lectures brutes jusqu’à la génération de transcrits personnalisés.

Le workflow inclut le contrôle qualité, le nettoyage des lectures, l’alignement, la quantification transcriptomique, l’appel de variants, la détection de fusions géniques et l’application des variants aux séquences transcriptomiques.

```text
FASTQ → QC → trimming → alignement → quantification → variants/fusions → transcrits personnalisés
```

---

## Résumé du pipeline

Le pipeline effectue les étapes suivantes :

1. Téléchargement des références génomiques
2. Contrôle qualité des lectures RNA-seq avec FastQC
3. Nettoyage des lectures avec Trim Galore
4. Contrôle qualité après trimming
5. Alignement des lectures avec STAR
6. Quantification transcriptomique avec kallisto
7. Appel de variants
8. Détection de fusions géniques
9. Application des variants aux transcrits
10. Génération de transcrits personnalisés

---

## Structure du dépôt

```text
AGlaude_2024/
│
├── config/
│   └── config.json
│
├── profile_local/
│   └── config.yaml
│
├── profile_slurm/
│   ├── config.yaml
│   ├── cluster.yaml
│   └── slurmSubmit.py
│
└── workflow/
    ├── Snakefile
    ├── rules/
    ├── scripts/
    ├── envs/
    ├── data/
    ├── results/
    └── logs/
```

Le fichier principal du pipeline est :

```text
workflow/Snakefile
```

La configuration générale est chargée depuis :

```text
config/config.json
```

Les règles Snakemake sont séparées dans le dossier :

```text
workflow/rules/
```

Les scripts utilisés par certaines règles sont dans :

```text
workflow/scripts/
```

Les environnements conda associés aux différentes étapes sont dans :

```text
workflow/envs/
```

---

## Configuration

Le pipeline utilise le fichier suivant :

```text
config/config.json
```

Ce fichier contient notamment les chemins vers les données RNA-seq et les ressources utilisées par le pipeline.

Dans le `Snakefile`, les échantillons sont détectés automatiquement à partir du répertoire défini par :

```python
config["path"]["RNAseq_dir"]
```

Chaque sous-dossier de ce répertoire est considéré comme un échantillon à analyser.

---

## Exécution du pipeline

Le pipeline doit être lancé à partir du dossier `workflow/`.

```bash
cd workflow
```

L’environnement conda utilisé doit contenir Snakemake.

```bash
conda activate snakemake
```

---

## Profils Snakemake

Deux profils sont disponibles pour exécuter le pipeline :

```text
profile_local/
profile_slurm/
```

### Exécution locale

Le profil local permet de lancer le pipeline directement, sans soumission SLURM.

```bash
snakemake --profile ../profile_local/
```

Il est aussi utilisé pour les étapes nécessitant un accès Internet, comme le téléchargement des références.

Si les nœuds de calcul du cluster n’ont pas accès à Internet, les références doivent d’abord être téléchargées localement avec :

```bash
snakemake download_genome --profile ../profile_local/
```

### Exécution sur SLURM

Une fois les références téléchargées, le pipeline complet peut être lancé sur le cluster avec le profil SLURM.

```bash
snakemake --profile ../profile_slurm/
```

Cette commande exécute la cible par défaut définie dans `rule all`.

---

## Règles principales

Le `Snakefile` principal inclut les fichiers de règles suivants :

```text
rules/I_download_genome.smk
rules/II_preprocess.smk
rules/III_star.smk
rules/III_kallisto.smk
rules/IV_calling_variants.smk
rules/IV_calling_genefusion.smk
rules/V_variant_split_apply_merge.smk
rules/IV_DRPs.smk
```

Ces fichiers contiennent les différentes étapes du workflow.

---

## Sorties principales

Les sorties principales sont générées dans le dossier :

```text
workflow/results/
```

Exemples de fichiers générés :

```text
results/STAR/{id}/Aligned.sortedByCoord.out.bam
results/STAR/{id}/Chimeric.out.junction
results/kallisto/{id}/abundance.tsv
results/variants/{id}/20QC_variant.vcf
results/star_fusion/{id}/star-fusion.fusion_predictions.abridged.tsv
results/final/{id}/combined_transcripts.fa
```

Dans ces chemins, `{id}` correspond au nom de l’échantillon.

---

## Fichiers non suivis par Git

Les dossiers contenant les données, résultats, logs et fichiers temporaires ne sont généralement pas suivis par Git.

Exemples :

```text
workflow/.snakemake/
workflow/data/
workflow/logs/
workflow/reference/
workflow/results/
```

Ces dossiers sont générés ou remplis pendant l’exécution du pipeline.

---

## Exemple d’utilisation

Télécharger les références localement :

```bash
cd workflow
snakemake download_genome --profile ../profile_local/
```

Lancer le pipeline complet sur SLURM :

```bash
cd workflow
snakemake --profile ../profile_slurm/
```

Lancer le pipeline localement :

```bash
cd workflow
snakemake --profile ../profile_local/
```
