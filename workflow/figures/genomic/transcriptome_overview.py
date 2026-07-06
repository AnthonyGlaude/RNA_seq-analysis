#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# =========================
# Paths
# =========================
kallisto_files = {
    "SH_SY5Y": "/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/kallisto/SRR4787043/abundance.tsv",
    "SK_N_BE2": "/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/kallisto/SRR4787045/abundance.tsv",
    "SK_N_BE2_C": "/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/kallisto/SRR4787046/abundance.tsv",
}

outdir = Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/figures/genomic/transcriptome_overview")
outdir.mkdir(parents=True, exist_ok=True)


blue_white_orange = LinearSegmentedColormap.from_list(
    "blue_white_orange",
    ["#2166AC", "#FFFFFF", "#F46D43"],
    N=256
)

# =========================
# Functions
# =========================
def load_kallisto_tpm(tsv_path: str, sample_name: str) -> pd.DataFrame:
    """
    Load one kallisto abundance.tsv and return a dataframe with:
    - target_id
    - TPM column renamed to sample_name
    """
    df = pd.read_csv(tsv_path, sep="\t")

    required_cols = {"target_id", "tpm"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"{tsv_path} is missing columns: {missing}")

    df = df[["target_id", "tpm"]].copy()
    df = df.rename(columns={"tpm": sample_name})
    return df


def build_tpm_matrix(files_dict: dict) -> pd.DataFrame:
    """
    Merge TPM columns from several kallisto files by target_id.
    """
    merged = None

    for sample_name, filepath in files_dict.items():
        df = load_kallisto_tpm(filepath, sample_name)

        if merged is None:
            merged = df
        else:
            merged = merged.merge(df, on="target_id", how="outer")

    merged = merged.fillna(0)
    merged = merged.set_index("target_id")
    return merged


# =========================
# Main
# =========================
def main():
    print("Loading kallisto TPM files...")
    tpm_matrix = build_tpm_matrix(kallisto_files)

    # Save raw TPM matrix
    raw_tpm_out = outdir / "kallisto_tpm_matrix.tsv"
    tpm_matrix.to_csv(raw_tpm_out, sep="\t")

    # log2(TPM + 1)
    log_tpm = np.log2(tpm_matrix + 1)

    log_tpm_out = outdir / "kallisto_log2_tpm_plus1_matrix.tsv"
    log_tpm.to_csv(log_tpm_out, sep="\t")

    # Correlation between samples
    corr = log_tpm.corr(method="pearson")
    corr_out = outdir / "sample_correlation_matrix.tsv"
    corr.to_csv(corr_out, sep="\t")

    print("\nCorrelation matrix:")
    print(corr)

    # Plot heatmap with matplotlib only
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        corr.values,
        aspect="auto",
        cmap=blue_white_orange,
        vmin=0.85,   # important pour corrélation
        vmax=1
    )
    ax.set_xticks(range(len(corr.columns)))
    ax.set_yticks(range(len(corr.index)))
    ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticklabels(corr.index)

    # Write correlation values inside cells
    for i in range(corr.shape[0]):
        for j in range(corr.shape[1]):
            ax.text(
                j, i,
                f"{corr.iloc[i, j]:.3f}",
                ha="center", va="center"
            )

    ax.set_title("Cell-line correlation\nlog2(TPM + 1)")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Pearson correlation")

    plt.tight_layout()

    png_out = outdir / "heatmap_sample_correlation_kallisto.png"
    pdf_out = outdir / "heatmap_sample_correlation_kallisto.pdf"

    plt.savefig(png_out, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_out, bbox_inches="tight")
    plt.close()

    print(f"\nFiles written to: {outdir}")
    print(f"- {raw_tpm_out}")
    print(f"- {log_tpm_out}")
    print(f"- {corr_out}")
    print(f"- {png_out}")
    print(f"- {pdf_out}")


if __name__ == "__main__":
    main()