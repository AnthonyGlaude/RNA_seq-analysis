import pandas as pd
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
import os

# Constantes
VENN_DIR = "venn_sets"
BACKGROUND_FILE = os.path.join(VENN_DIR, "background_genes.txt")
GAF = os.path.join(VENN_DIR, "data", "goa_human.gaf")
OBO = os.path.join(VENN_DIR, "data", "go-basic.obo")

# Lire les gènes depuis un fichier
def lire_genes(fichier):
    with open(fichier) as f:
        return [l.strip().upper() for l in f if l.strip()]

# Enrichissement GO
def gene_enrichment(study, background, gaf_file, obo_file, namespace='BP'):
    dag = GODag(obo_file)
    gene2gos = {}
    ns_code = {'BP': 'P', 'MF': 'F', 'CC': 'C'}[namespace]

    with open(gaf_file) as f:
        for line in f:
            if line.startswith('!'):
                continue
            parts = line.strip().split('\t')
            gene_symbol = parts[2].strip().upper()
            go_term = parts[4]
            ns = parts[8]
            if ns == ns_code and gene_symbol in background:
                gene2gos.setdefault(gene_symbol, set()).add(go_term)

    if not gene2gos:
        raise RuntimeError("Aucune annotation GO trouvée pour les gènes du background.")

    enricher = GOEnrichmentStudy(
        background,
        gene2gos,
        dag,
        propagate_counts=False,
        alpha=0.05,
        methods=['fdr_bh']
    )

    results = enricher.run_study(study)

    # Résultats significatifs
    df_sig = pd.DataFrame([{
        'GO': r.GO,
        'Name': r.name,
        'NS': r.NS,
        'p_uncorrected': r.p_uncorrected,
        'p_fdr_bh': getattr(r, 'p_fdr_bh'),
        'Study_Count': r.study_count,
        'Pop_Count': r.pop_count
    } for r in results if getattr(r, 'p_fdr_bh', 1) < 0.05])

    # Tous les résultats
    df_all = pd.DataFrame([{
        'GO': r.GO,
        'Name': r.name,
        'NS': r.NS,
        'p_uncorrected': r.p_uncorrected,
        'p_fdr_bh': getattr(r, 'p_fdr_bh'),
        'Study_Count': r.study_count,
        'Pop_Count': r.pop_count
    } for r in results])

    return (
        df_sig.sort_values("p_fdr_bh") if not df_sig.empty else df_sig,
        df_all.sort_values("p_fdr_bh") if not df_all.empty else df_all
    )

# Lire le background
background = set(lire_genes(BACKGROUND_FILE))

# Boucle sur tous les fichiers .txt sauf background
for fname in os.listdir(VENN_DIR):
    if not fname.endswith(".txt") or fname == "background_genes.txt":
        continue

    study_file = os.path.join(VENN_DIR, fname)
    prefix = fname.replace(".txt", "")
    study = set(lire_genes(study_file))

    print(f"\n Analyse de {prefix} ({len(study)} gènes)")
    output_dir = os.path.join(VENN_DIR, "results", prefix)
    os.makedirs(output_dir, exist_ok=True)

    for ns in ["BP", "MF", "CC"]:
        try:
            print(f" Namespace {ns}")
            df_sig, df_all = gene_enrichment(study, background, GAF, OBO, namespace=ns)
            print(f" {len(df_sig)} termes significatifs trouvés pour {ns}")
            
            # Sauvegardes
            df_sig.to_csv(os.path.join(output_dir, f"{prefix}_GO_{ns}_significant.csv"), index=False)
            df_all.to_csv(os.path.join(output_dir, f"{prefix}_GO_{ns}_all.csv"), index=False)

            # Affichage console : top 10
            print(df_all[['GO', 'Name', 'p_fdr_bh']].head(10))

        except RuntimeError as e:
            print(f" {ns} ignoré : {str(e)}")