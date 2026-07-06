#!/usr/bin/env python3

from pathlib import Path
import subprocess
import shutil
import pandas as pd
import matplotlib.pyplot as plt

GTF = Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/data/references/gtf/Homo_sapiens.GRCh38.113.gtf")

VCFS = [
    Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/variants/SRR4787043/20QC_variant.vcf"),
    Path("/home/glaudea/scratch/glaudea/AGlaude_2024/workflow/results/variants/SRR4787046/20QC_variant.vcf"),
]

OUT_ROOT = Path("mut_hist_transcript_out")

MAX_X = 10
GREEN = "#7CFC90"


def run(cmd):
    subprocess.run(cmd, shell=True, check=True)


def check_dependencies():
    if shutil.which("bedtools") is None:
        raise RuntimeError(
            "bedtools introuvable. Fais: module load bedtools "
            "ou conda install -c bioconda -c conda-forge bedtools"
        )


def extract_attr(attr_string, key):
    target = f'{key} "'
    if target not in attr_string:
        return None
    return attr_string.split(target, 1)[1].split('"', 1)[0]


def gtf_to_bed(feature, outbed):
    with open(GTF) as f, open(outbed, "w") as out:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _, feat, start, end, _, _, _, attrs = parts

            if feat != feature:
                continue

            tx = extract_attr(attrs, "transcript_id")
            if tx is None:
                continue

            start0 = int(start) - 1
            end = int(end)

            out.write(f"{chrom}\t{start0}\t{end}\t{tx}\n")


def vcf_to_bed(vcf, outbed):
    with open(vcf) as f, open(outbed, "w") as out:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alts = parts[4].split(",")

            start0 = pos - 1
            end = start0 + len(ref)

            for alt in alts:
                if alt and alt != ".":
                    out.write(f"{chrom}\t{start0}\t{end}\t{ref}\t{alt}\n")


def intersect_and_count(outdir, label, bed):
    hits = outdir / f"hits_{label}.tsv"
    mut_per_tx = outdir / f"mut_per_transcript_{label}.txt"
    dist = outdir / f"dist_{label}.tsv"

    run(f"bedtools intersect -a {outdir/'variants.bed'} -b {bed} -wa -wb > {hits}")

    if not hits.exists() or hits.stat().st_size == 0:
        mut_per_tx.write_text("")
        dist.write_text("")
        return

    df = pd.read_csv(hits, sep="\t", header=None)

    if df.empty:
        mut_per_tx.write_text("")
        dist.write_text("")
        return

    # colonnes:
    # 0 chrom variant
    # 1 start variant
    # 2 end variant
    # 3 ref
    # 4 alt
    # 5 chrom annotation
    # 6 start annotation
    # 7 end annotation
    # 8 transcript_id
    df = df[[0, 1, 2, 3, 4, 8]].drop_duplicates()
    df.columns = ["chrom", "start", "end", "ref", "alt", "transcript_id"]

    counts = (
        df.groupby("transcript_id")
        .size()
        .reset_index(name="n_mut")
        .sort_values("n_mut", ascending=False)
    )

    counts[["n_mut", "transcript_id"]].to_csv(
        mut_per_tx, sep="\t", index=False, header=False
    )

    d = counts["n_mut"].value_counts().sort_index().reset_index()
    d.columns = ["k_mut", "n_transcripts"]
    d.to_csv(dist, sep="\t", index=False, header=False)


def capture_summary(outdir, labels=("transcript", "exon", "cds")):
    rows = []

    for label in labels:
        p = outdir / f"dist_{label}.tsv"

        if not p.exists() or p.stat().st_size == 0:
            continue

        d = pd.read_csv(p, sep="\t", header=None, names=["k_mut", "n_transcripts"])

        total_tx = d["n_transcripts"].sum()
        captured_tx = d.loc[d["k_mut"] <= MAX_X, "n_transcripts"].sum()

        total_mut_assignments = (d["k_mut"] * d["n_transcripts"]).sum()
        captured_mut_assignments = (
            d.loc[d["k_mut"] <= MAX_X, "k_mut"]
            * d.loc[d["k_mut"] <= MAX_X, "n_transcripts"]
        ).sum()

        rows.append({
            "level": label,
            "threshold": MAX_X,
            "captured_transcripts": captured_tx,
            "total_mutated_transcripts": total_tx,
            "pct_transcripts_captured": captured_tx / total_tx * 100,
            "captured_mutation_assignments": captured_mut_assignments,
            "total_mutation_assignments": total_mut_assignments,
            "pct_mutation_assignments_captured": captured_mut_assignments / total_mut_assignments * 100,
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(outdir / "capture_summary_leq10.tsv", sep="\t", index=False)

    print("\nCapture summary:")
    print(summary.to_string(index=False))


def plot_mutation_distribution_2x2():
    out_png = Path("figure_supp_distribution_mutations_exon_cds.png")

    panels = [
        ("A", "SRR4787043", "exon", "SH-SY5Y : régions exoniques"),
        ("B", "SRR4787043", "cds",  "SH-SY5Y : CDS"),
        ("C", "SRR4787046", "exon", "SK-N-Be2 : régions exoniques"),
        ("D", "SRR4787046", "cds",  "SK-N-Be2 : CDS"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 9.5), sharex=True, sharey=False)
    axes = axes.flatten()

    for ax, (panel, srr, level, title) in zip(axes, panels):
        path = OUT_ROOT / srr / f"dist_{level}.tsv"

        if not path.exists() or path.stat().st_size == 0:
            ax.text(
                0.5,
                0.5,
                f"Fichier manquant:\n{path}",
                ha="center",
                va="center",
                fontsize=12
            )
            ax.axis("off")
            continue

        d = pd.read_csv(path, sep="\t", header=None, names=["k_mut", "n_transcripts"])
        d = d.sort_values("k_mut")

        d_zoom = d[d["k_mut"] <= MAX_X].copy()

        total_tx = d["n_transcripts"].sum()
        n_leq10 = d.loc[d["k_mut"] <= MAX_X, "n_transcripts"].sum()
        n_gt10 = d.loc[d["k_mut"] > MAX_X, "n_transcripts"].sum()
        pct_leq10 = n_leq10 / total_tx * 100 if total_tx > 0 else 0

        ax.bar(
            d_zoom["k_mut"],
            d_zoom["n_transcripts"],
            width=0.82,
            color=GREEN,
            edgecolor="black",
            linewidth=0.35,
            alpha=0.95,
        )

        ax.set_title(
            f"{panel}. {title}",
            loc="left",
            fontsize=16,
            fontweight="bold",
            pad=10,
        )

        ax.text(
            0.97,
            0.88,
            f"≤10 mutations : {pct_leq10:.1f} %\n>10 mutations : {n_gt10:,}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=14,
            fontweight="bold",
            linespacing=1.4,
            bbox=dict(
                boxstyle="round,pad=0.55",
                facecolor="white",
                edgecolor="gray",
                linewidth=0.9,
                alpha=0.97,
            ),
        )

        ax.set_xlim(0.5, 10.8)
        ax.set_xticks(range(1, MAX_X + 1))

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.tick_params(axis="both", labelsize=13)

    for ax in axes[2:]:
        ax.set_xlabel("Nombre de mutations par transcrit", fontsize=15)

    for ax in axes[::2]:
        ax.set_ylabel("Nombre de transcrits", fontsize=15)

    plt.tight_layout(pad=2.0)

    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    print(f"\n[OK] Figure écrite: {out_png}")
    print(f"[OK] Figure PDF écrite: {out_png.with_suffix('.pdf')}")
    print(f"[OK] Figure SVG écrite: {out_png.with_suffix('.svg')}")


def main():
    check_dependencies()

    OUT_ROOT.mkdir(exist_ok=True)

    refdir = OUT_ROOT / "_refbeds"
    refdir.mkdir(exist_ok=True)

    print("[0] Build transcript-level BEDs from GTF")
    gtf_to_bed("transcript", refdir / "transcripts.bed")
    gtf_to_bed("exon", refdir / "exons.bed")
    gtf_to_bed("CDS", refdir / "cds.bed")

    for vcf in VCFS:
        srr = vcf.parent.name
        outdir = OUT_ROOT / srr
        outdir.mkdir(exist_ok=True)

        print(f"\n[1] Processing {srr}")

        vcf_to_bed(vcf, outdir / "variants.bed")

        intersect_and_count(outdir, "transcript", refdir / "transcripts.bed")
        intersect_and_count(outdir, "exon", refdir / "exons.bed")
        intersect_and_count(outdir, "cds", refdir / "cds.bed")

        capture_summary(outdir)

    plot_mutation_distribution_2x2()

    print(f"\n[DONE] Outputs in: {OUT_ROOT}/<SRR>/")


if __name__ == "__main__":
    main()