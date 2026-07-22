
import sys
from pathlib import Path
from Bio import SeqIO


"""
Crée un transcriptome personnalisé pour un patient en combinant :
- tous les transcrits originaux (de transcripts_fasta)
- tous les transcrits mutés générés à partir des différents morceaux de VCF

Ajout :
- retire les doublons de transcrits en se basant sur record.id
- conserve la première occurrence rencontrée
- traite les fichiers un par un pour éviter de charger les FASTA complets en mémoire
"""


def save_combined_transcripts(
    transcripts_fasta,
    mutated_output_fastas,
    combined_output_fasta,
    report_output=None,
):
    seen_ids = set()

    stats = []

    Path(combined_output_fasta).parent.mkdir(parents=True, exist_ok=True)

    with open(combined_output_fasta, "w") as combined_out:

        # 1. Ajouter le transcriptome original en premier
        source_name = "transcripts_fasta"
        total = 0
        written = 0
        skipped = 0

        for record in SeqIO.parse(transcripts_fasta, "fasta"):
            total += 1
            record_id = record.id

            if record_id in seen_ids:
                skipped += 1
                continue

            seen_ids.add(record_id)
            SeqIO.write(record, combined_out, "fasta")
            written += 1

        stats.append((source_name, transcripts_fasta, total, written, skipped))

        # 2. Ajouter ensuite les FASTA mutés, en retirant les IDs déjà vus
        for i, mutated_fasta in enumerate(mutated_output_fastas, start=1):
            source_name = f"mutated_output_fasta{i}"
            total = 0
            written = 0
            skipped = 0

            for record in SeqIO.parse(mutated_fasta, "fasta"):
                total += 1
                record_id = record.id

                if record_id in seen_ids:
                    skipped += 1
                    continue

                seen_ids.add(record_id)
                SeqIO.write(record, combined_out, "fasta")
                written += 1

            stats.append((source_name, mutated_fasta, total, written, skipped))

    # 3. Écrire un petit rapport optionnel
    if report_output is not None:
        Path(report_output).parent.mkdir(parents=True, exist_ok=True)

        with open(report_output, "w") as report:
            report.write(
                "source\tinput_fasta\ttotal_records\twritten_records\tskipped_duplicates\n"
            )

            for source_name, fasta_path, total, written, skipped in stats:
                report.write(
                    f"{source_name}\t{fasta_path}\t{total}\t{written}\t{skipped}\n"
                )

    print("[INFO] Combined transcriptome created.")
    print(f"[INFO] Output: {combined_output_fasta}")
    print(f"[INFO] Unique transcript IDs written: {len(seen_ids)}")

    for source_name, fasta_path, total, written, skipped in stats:
        print(
            f"[INFO] {source_name}: "
            f"total={total}, written={written}, skipped_duplicates={skipped}"
        )


try:
    transcripts_fasta = snakemake.input.transcripts_fasta

    mutated_output_fastas = [
        snakemake.input.mutated_output_fasta1,
        snakemake.input.mutated_output_fasta2,
        snakemake.input.mutated_output_fasta3,
        snakemake.input.mutated_output_fasta4,
        snakemake.input.fusions_transcript
    ]

    combined_output_fasta = snakemake.output.combined_transcripts

    try:
        report_output = snakemake.output.report
    except AttributeError:
        report_output = None

except NameError:
    if len(sys.argv) < 4:
        sys.exit(
            "Usage: combine_transcripts.py "
            "<transcripts_fasta> <combined_output_fasta> <mutated_output_fastas...>"
        )

    transcripts_fasta = sys.argv[1]
    combined_output_fasta = sys.argv[2]
    mutated_output_fastas = sys.argv[3:]
    report_output = None


save_combined_transcripts(
    transcripts_fasta=transcripts_fasta,
    mutated_output_fastas=mutated_output_fastas,
    combined_output_fasta=combined_output_fasta,
    report_output=report_output,
)
