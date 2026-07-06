# Règle pour la quantification des transcrits avec Kallisto

rule build_transcriptome:  # Sert de transcriptome de référence avec tous les ID transcrits
    """UTRs and exons are extracted using gffread but introns are not included."""
    input:
        genome = rules.download_human_genome.output.genome,
        gtf = rules.download_human_gtf.output.gtf
    output:
        transcriptome = config["path"]["transcriptome"]
    conda:
        "../envs/gffread.yml"
    message:
        "Build a reference transcriptome using gffread."
    log:
        "logs/kallisto/build_transcriptome.log"
    shell:
        """
        gffread {input.gtf} -g {input.genome} -w {output}
        """


rule kallisto_index:  # L'index qui sert à quantifier
    input:
        rules.build_transcriptome.output
    output:
        "data/references/kallisto.idx"
    threads:
        32
    conda:
        "../envs/kallisto.yml"
    log:
        "logs/kallisto/index.log"
    message:
        "Builds an index from the FASTA file."
    shell:
        "kallisto index "
        "--index={output} "
        "{input} "
        "&> {log}"


rule kallisto_quant:  # Donne l'abondance des transcrits chez l'individu
    input:
        idx = rules.kallisto_index.output,
        fq1 = rules.trim_reads.output.gal_trim1,
        fq2 = rules.trim_reads.output.gal_trim2
    output:
        abundance_tsv = "results/kallisto/{id}/abundance.tsv"
    params:
        bootstrap = "50",
        outdir = "results/kallisto/{id}/"
    threads:
        32
    conda:
        "../envs/kallisto.yml"
    log:
        "logs/kallisto/{id}/log/{id}.log"
    message:
        "Perform pseudoalignment and quantify transcript abundance for {wildcards.id}."
    shell:
        "kallisto quant "
        "--bias "
        "--index={input.idx} "
        "--output-dir={params.outdir} "
        "--bootstrap-samples={params.bootstrap} "
        "--threads={threads} "
        "{input.fq1} {input.fq2} "
        "&> {log}"


rule filter_abundance:  # Ne contient que les transcrits présents dans l'ARNseq (donc retire les 0)
    input:
        abundance = rules.kallisto_quant.output.abundance_tsv
    output:
        filtered_abundance = "results/kallisto/{id}/abundance_filtered.tsv",  # Contient seulement les gènes quantifiés
        transcript_ids = "results/kallisto/{id}/transcript_ids.txt"  # Contient seulement l'ID des gènes quantifiés
    log:
        "logs/kallisto/{id}/log/filter_{id}.log"
    shell:
        "awk '$4 > 0' {input.abundance} > {output.filtered_abundance} && "
        "cut -f1 {output.filtered_abundance} > {output.transcript_ids}"


rule build_filtered_transcriptome:
    input:
        gtf = rules.download_human_gtf.output.gtf,
        genome = rules.download_human_genome.output.genome,
        transcript_ids = rules.filter_abundance.output.transcript_ids
    output:
        filtered_gtf = "results/{id}/filtered_transcripts.gtf", #utile?
        transcriptome_final_custom = "results/{id}/transcriptome_final_custom.fa"
    conda:
        "../envs/gffread.yml"
    message:
        "Build a filtered reference transcriptome using gffread."
    log:
        "logs/kallisto/{id}/log/build_filtered_transcriptome_{id}.log"
    shell:
        """
        # Filtrer le GTF pour ne garder que les transcrits désirés
        grep -F -f {input.transcript_ids} {input.gtf} > {output.filtered_gtf}
        
        # Modifier le GTF pour remplacer 'geneID' par 'gene_name'
        awk '{{ gsub("geneID", "gene_name"); print }}' {output.filtered_gtf} > {output.filtered_gtf}_modified.gtf
        
        # Générer le transcriptome à partir du GTF filtré et modifié
        gffread {output.filtered_gtf}_modified.gtf -g {input.genome} -w {output.transcriptome_final_custom}
        """
