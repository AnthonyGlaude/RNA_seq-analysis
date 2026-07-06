#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# =========================
# Inputs
# =========================
kallisto_files = {
    "SH_SY5Y": "/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/kallisto/SRR4787043/abundance.tsv",
    "SK_N_BE2": "/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/kallisto/SRR4787045/abundance.tsv",
    "SK_N_BE2_C": "/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/kallisto/SRR4787046/abundance.tsv",
}

gtf_path = "/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/data/references/gtf/Homo_sapiens.GRCh38.113.gtf"

outdir = Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/figures/genomic/transcriptome_overview")
outdir.mkdir(parents=True, exist_ok=True)


# =========================
# Optional custom colormap
# =========================
blue_white_orange = LinearSegmentedColormap.from_list(
    "blue_white_orange",
    ["#2166AC", "#FFFFFF", "#F46D43"],
    N=256
)


# =========================
# Helpers
# =========================
def parse_gtf_attributes(attr_string: str) -> dict:
    """
    Parse the GTF attributes column into a dictionary.
    Example:
    gene_id "ENSG..."; transcript_id "ENST..."; gene_name "MYCN";
    """
    attrs = {}
    for match in re.finditer(r'(\S+)\s+"([^"]+)"', attr_string):
        key, value = match.groups()
        attrs[key] = value
    return attrs


def get_transcripts_for_gene(gtf_file: str, gene_name: str) -> list[str]:
    """
    Extract transcript_id values for a given gene_name from a GTF.
    Restricts to feature == 'transcript' when possible.
    """
    transcripts = set()

    with open(gtf_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            attrs = parse_gtf_attributes(fields[8])

            if attrs.get("gene_name") == gene_name:
                tx_id = attrs.get("transcript_id")
                if tx_id is None:
                    continue

                # Prefer transcript rows, but if absent in some GTF contexts,
                # this still works because we deduplicate with a set.
                if feature == "transcript":
                    transcripts.add(tx_id)

    # Fallback in case transcript rows are absent or unusual
    if not transcripts:
        with open(gtf_file, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue

                attrs = parse_gtf_attributes(fields[8])
                if attrs.get("gene_name") == gene_name:
                    tx_id = attrs.get("transcript_id")
                    if tx_id is not None:
                        transcripts.add(tx_id)

    return sorted(transcripts)


def clean_kallisto_target_id(series: pd.Series) -> pd.Series:
    """
    Remove version suffix if present:
    ENST00000377970.8 -> ENST00000377970
    """
    return series.astype(str).str.replace(r"\.\d+$", "", regex=True)


# =========================
# Main
# =========================
def main():
    gene_of_interest = "MYCN"

    # 1) Get MYCN transcripts from GTF
    mycn_transcripts = get_transcripts_for_gene(gtf_path, gene_of_interest)

    if not mycn_transcripts:
        raise RuntimeError(f"No transcripts found for gene {gene_of_interest} in {gtf_path}")

    print(f"[INFO] {gene_of_interest} transcripts found in GTF: {len(mycn_transcripts)}")
    for tx in mycn_transcripts:
        print("   ", tx)

    # Save transcript list
    tx_out = outdir / f"{gene_of_interest}_transcripts_from_gtf.txt"
    with open(tx_out, "w", encoding="utf-8") as out:
        for tx in mycn_transcripts:
            out.write(tx + "\n")

    # 2) Extract MYCN TPM from each sample
    results = []

    for sample, path in kallisto_files.items():
        df = pd.read_csv(path, sep="\t")

        if "target_id" not in df.columns or "tpm" not in df.columns:
            raise ValueError(f"{path} must contain 'target_id' and 'tpm' columns")

        df["target_id_clean"] = clean_kallisto_target_id(df["target_id"])
        df_mycn = df[df["target_id_clean"].isin(mycn_transcripts)].copy()

        total_tpm = df_mycn["tpm"].sum()
        log2_tpm = np.log2(total_tpm + 1)

        results.append({
            "sample": sample,
            "MYCN_TPM_sum": total_tpm,
            "MYCN_log2_TPM_plus1": log2_tpm,
            "n_MYCN_transcripts_detected": (df_mycn["tpm"] > 0).sum(),
            "n_MYCN_transcripts_in_gtf": len(mycn_transcripts),
        })

        # Save per-sample MYCN transcript table
        df_mycn[["target_id", "target_id_clean", "tpm", "est_counts"]].sort_values(
            by="tpm", ascending=False
        ).to_csv(
            outdir / f"{sample}_{gene_of_interest}_transcripts.tsv",
            sep="\t",
            index=False
        )

    plot_df = pd.DataFrame(results).sort_values("MYCN_log2_TPM_plus1", ascending=False)
    plot_df.to_csv(outdir / f"{gene_of_interest}_expression_summary.tsv", sep="\t", index=False)

    print("\n[INFO] Summary:")
    print(plot_df)

    # 3) Plot
    fig, ax = plt.subplots(figsize=(6, 5))
    bars = ax.bar(
        plot_df["sample"],
        plot_df["MYCN_log2_TPM_plus1"]
    )

    # Optional: color bars with your palette endpoints / midpoint feel
    bar_colors = ["#F46D43", "#E6E6E6", "#2166AC"]
    for bar, color in zip(bars, bar_colors[:len(bars)]):
        bar.set_color(color)

    ax.set_ylabel("log2(sum TPM + 1)")
    ax.set_title("MYCN expression across samples")
    plt.xticks(rotation=45, ha="right")

    # Add values above bars
    for bar, value in zip(bars, plot_df["MYCN_log2_TPM_plus1"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{value:.2f}",
            ha="center",
            va="bottom"
        )

    plt.tight_layout()
    plt.savefig(outdir / f"{gene_of_interest}_expression_barplot.png", dpi=300, bbox_inches="tight")
    plt.savefig(outdir / f"{gene_of_interest}_expression_barplot.pdf", bbox_inches="tight")
    plt.savefig(outdir / f"{gene_of_interest}_expression_barplot.svg", bbox_inches="tight")
    plt.close()

    print(f"\n[OK] Outputs written to: {outdir}")


if __name__ == "__main__":
    main()