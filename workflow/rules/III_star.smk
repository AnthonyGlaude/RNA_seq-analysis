rule star_index:
    """ Generates the genome index for STAR """
    input:
        fasta = rules.download_human_genome.output.genome,
        gtf = rules.download_human_gtf.output.gtf #GRCh38.p113 reference
    output:
        chrNameLength = "data/references/star_index/chrNameLength.txt"
    params:
        dir = config['path']['star_index']  
    log:
        "logs/index.log"
    conda:
        "../envs/star.yml"
    threads:
        32
    shell:
        """
        mkdir -p {params.dir} && \
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params.dir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99 \
        &> {log}
        """

rule star_alignreads:
    """ Generates a bam file using STAR
        unimapping stric
    """
    input:
        idx = rules.star_index.output,
        fq1 = rules.trim_reads.output.gal_trim1,
        fq2 = rules.trim_reads.output.gal_trim2
    output:
        bam = "results/STAR/{id}/Aligned.sortedByCoord.out.bam",
        bam_logs = "results/STAR/{id}/Log.final.out" 
    params:
        index = config['path']['star_index'],
        output_dir = "results/STAR/{id}/"
    log:
        "logs/{id}/alignreads.log"
    threads:
        32
    conda:
        "../envs/star.yml"
    shell:
        """
        mkdir -p {params.output_dir} && \
        STAR --runMode alignReads \
            --genomeDir {params.index} \
            --readFilesIn {input.fq1} {input.fq2} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outReadsUnmapped Fastx \
            --outFilterType BySJout \
            --outStd Log \
            --outSAMunmapped None \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.output_dir} \
            --outFilterScoreMinOverLread 0.66 \
            --outFilterMatchNminOverLread 0.66 \
            --outFilterMultimapNmax 1 \
            --winAnchorMultimapNmax 1 \
            --limitBAMsortRAM 120000000000 \
            --outTmpDir /tmp/{wildcards.id} \
            &> {log}
        """

rule build_star_index_ctat:
    input:
        fasta = "data/references/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa",
        gtf   = "data/references/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
    output:
        sa = "data/references/star_ctat_index/SA"
    params:
        dir = "data/references/star_ctat_index"
    threads: 32
    conda:
        "../envs/star.yml"
    log:
        "logs/build_star_index_ctat.log"
    shell:
        """
        mkdir -p {params.dir} && \
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params.dir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99 \
        &> {log}
        """

rule star_chimeric:
    input:
        fq1 = rules.trim_reads.output.gal_trim1,
        fq2 = rules.trim_reads.output.gal_trim2,
        sa  = rules.build_star_index_ctat.output.sa
    output:
        bam = "results/STAR_fusion/{id}/Aligned.sortedByCoord.out.bam",
        junction = "results/STAR_fusion/{id}/Chimeric.out.junction",
        sj = "results/STAR_fusion/{id}/SJ.out.tab"
    params:
        index = "data/references/star_ctat_index",
        outdir = "results/STAR_fusion/{id}/"
    log:
        "logs/{id}/star_chimeric.log"
    threads: 32
    conda:
        "../envs/star.yml"
    shell:
        """
        mkdir -p {params.outdir}
        STAR \
            --runMode alignReads \
            --genomeDir {params.index} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {params.outdir} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --outSAMunmapped Within \
            --outReadsUnmapped None \
            --twopassMode Basic \
            --chimOutType Junctions WithinBAM HardClip \
            --chimOutJunctionFormat 1 \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --chimScoreDropMax 30 \
            --chimScoreSeparation 1 \
            --chimSegmentReadGapMax 3 \
            --chimScoreJunctionNonGTAG 0 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --alignIntronMax 1000000 \
            --outFilterMultimapNmax 50 \
            --limitBAMsortRAM 120000000000 \
            --outTmpDir /tmp/{wildcards.id}_chim \
            &> {log}
        """