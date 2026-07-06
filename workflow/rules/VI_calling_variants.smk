rule index_genome:
    input:
        genome=rules.download_human_genome.output.genome
    output:
        fai="data/references/genome_fa/homo_sapiens_genome.fa.fai"
    conda:
        "../envs/freebayes.yml"
    shell:
        "samtools faidx {input.genome}"

rule create_targets:
    input:
        fai=rules.index_genome.output.fai
    output:
        targets="results/variants/{id}/targets.txt"
    conda:
        "../envs/freebayes.yml"
    shell:
        "scripts/fasta_generate_regions.py {input.fai} 1000000 > {output.targets}"

rule index_bam:
    input:
        bam=rules.star_alignreads.output.bam
    output:
        bai="results/STAR/{id}/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/freebayes.yml"
    shell:
        "samtools index {input.bam}"

rule call_variants:
    input:
        bam = rules.star_alignreads.output.bam,
        bai = rules.index_bam.output.bai,
        genome = rules.download_human_genome.output.genome,
        fai = rules.index_genome.output.fai,
        targets = rules.create_targets.output.targets
    output:
        vcf = "results/variants/{id}/variants.vcf"
    params:
        min_alternate_count = 10,            #Exige au moins 10 lectures supportant l'allele alternatif (reduit faux positif)
        min_coverage = 10,                  # couverture min de 10reads a  une position pour etre positif
        min_base_quality = 20,              # ignore Q phred <20 (≥ 1% erreur) lors des appels
        min_mapping_quality = 30,           # lecture bien alignée MAPQ = 30 (1/1000 erreur)
        min_alternate_fraction = 0.2,       # 20 % des lectures doivent soutenir l’allèle alternatif
    conda:
        "../envs/freebayes.yml"
    log:
        "logs/freebayes_{id}/freebayes.log"
    threads: 32
    shell:
        """
        freebayes-parallel {input.targets} {threads} \
            -f {input.genome} {input.bam} \
            --min-alternate-count {params.min_alternate_count} \
            --min-coverage {params.min_coverage} \
            --min-base-quality {params.min_base_quality} \
            --min-mapping-quality {params.min_mapping_quality} \
            --min-alternate-fraction {params.min_alternate_fraction} \
            > {output.vcf} 2> {log}
        """

rule snpeff_annotation:
    input:
        vcf = rules.call_variants.output.vcf,
        db = rules.download_snpeff_db.output.done
    output:
        vcf = "results/variants/{id}/variants.snpeff.vcf"
    params:
        db = "GRCh38.99",
        datadir = "data/references/snpeff",
        snpeff = config["path"]["snpeff"]
    log:
        "logs/snpeff_{id}.log"
    shell:
        """
        {params.snpeff} -dataDir {params.datadir} -nodownload -v {params.db} {input.vcf} > {output.vcf} 2> {log}
        """
        
rule filter_variants:
    """
    Filtration des variants FreeBayes :
    - QUAL <20 : moins de 99% de confiance (Phred)
    - FORMAT/DP <5 : profondeur trop faible pour un appel fiable
    - AO <10 : au moins 10 lectures supportant l’allèle alternatif
    """
    input:
        vcf = rules.call_variants.output.vcf
    output:
        vcf_filtered = "results/variants/{id}/20QC_variant.vcf"
    conda:
        "../envs/bcftools.yml"
    log:
        "logs/freebayes_{id}/filter_variants.log"
    shell:
        """
        bcftools filter -s LowQual -e 'QUAL<20 || INFO/DP<5 || INFO/AO<10' {input.vcf} |
        bcftools view -f PASS > {output.vcf_filtered} 2> {log}
        """


#        "bcftools filter -s LowQual -e 'QUAL<20' {input.vcf} -o {output.vcf_filtered} > {log} 2>&1"
