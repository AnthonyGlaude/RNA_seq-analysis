#!/bin/bash
#SBATCH --job-name=mutation_histogram  # Nom du job mis à jour
#SBATCH --output=mutation_histogram_output.txt  # Fichier de sortie
#SBATCH --error=mutation_histogram_error.txt    # Fichier d'erreurs
#SBATCH --time=00:30:00                         # Temps limite
#SBATCH --ntasks=1                              # Nombre de tâches
#SBATCH --cpus-per-task=1                       # Nombre de cœurs par tâche
#SBATCH --mem=8G                               # Mémoire RAM demandée

# Chargement de l'environnement Conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate python_env  # Remplace par ton env exact si besoin

# Exécution du nouveau script Python
python distribution_mutation_by_gene.py
