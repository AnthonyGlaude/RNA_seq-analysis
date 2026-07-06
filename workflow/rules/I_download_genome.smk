rule download_human_gtf:
    """ Download gtf of human genome from Ensembl """
    output:
        gtf = 'data/references/gtf/Homo_sapiens.GRCh38.113.gtf'
    params:
        link = config['download']['human_gtf']
    shell:
        """
        mkdir -p data/references/gtf &&
        wget -O temp_gtf.gz {params.link} &&
        gunzip temp_gtf.gz &&
        mv temp_gtf {output.gtf} &&  # Renommer le fichier temporaire en gtf
        rm -f temp_gtf.gz
        """


rule download_human_genome:
    """Download the reference genome (fasta file) of human from ENSEMBL ftp servers."""
    output:
        genome = 'data/references/genome_fa/homo_sapiens_genome.fa'
    params:
        link = config['download']['human_genome_fa']
    shell:
        """
        mkdir -p data/references/genome_fa && 
        wget -O temp_genome.gz {params.link} &&
        gunzip temp_genome.gz && 
        mv temp_genome {output.genome} &&
        rm -f temp_genome.gz
        """

rule install_snpeff:
    output:
        jar  = "data/references/snpeff/snpEff/snpEff.jar",
        done = "data/references/snpeff/snpEff/.install_done"
    params:
        url = "https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip",
        outdir = "data/references/snpeff"
    log:
        "logs/install_snpeff.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}
        wget -O {params.outdir}/snpEff_latest_core.zip {params.url} > {log} 2>&1
        unzip -o {params.outdir}/snpEff_latest_core.zip -d {params.outdir} >> {log} 2>&1
        touch {output.done}
        """
        
rule download_snpeff_db:
    input:
        jar = rules.install_snpeff.output.jar
    output:
        done2 = "data/references/snpeff/snpEff/data/GRCh38.115/.download_done"
    params:
        db = "GRCh38.115",
        config = "data/references/snpeff/snpEff/snpEff.config"
    log:
        "logs/snpeff_download.log"
    shell:
        r"""
        set -euo pipefail

        module load java/21.0.1

        java -jar {input.jar} download \
            -c {params.config} \
            -configOption data.dir=data \
            -v {params.db} > {log} 2>&1

        test -d data/references/snpeff/snpEff/data/{params.db}
        touch {output.done2}
        """

rule prepare_arriba_database:
    output:
        blacklist = "data/references/arriba/blacklist_hg38_GRCh38_v2.5.0.tsv.gz",
        known     = "data/references/arriba/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz",
        domains   = "data/references/arriba/protein_domains_hg38_GRCh38_v2.5.0.gff3",
        done3      = "data/references/arriba/.download_done"
    params:
        url    = "https://github.com/suhrig/arriba/releases/download/v2.5.0/arriba_v2.5.0.tar.gz",
        outdir = "data/references/arriba",
    log:
        "logs/arriba_database.log"
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.outdir}
        wget -O {params.outdir}/arriba_v2.5.0.tar.gz {params.url} > {log} 2>&1

        tar -xzf {params.outdir}/arriba_v2.5.0.tar.gz -C {params.outdir} >> {log} 2>&1

        cp {params.outdir}/arriba_v2.5.0/database/blacklist_hg38_GRCh38_v2.5.0.tsv.gz {output.blacklist}
        cp {params.outdir}/arriba_v2.5.0/database/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz {output.known}
        cp {params.outdir}/arriba_v2.5.0/database/protein_domains_hg38_GRCh38_v2.5.0.gff3 {output.domains}

        touch {output.done3}
        """

rule download_ctat_lib:
    output:
        tar = "data/references/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz"
    params:
        url = "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play.tar.gz"
    log:
        "logs/download_ctat_lib.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p data/references/ctat
        wget -O {output.tar} {params.url} > {log} 2>&1
        """

rule extract_ctat_lib:
    input:
        tar = rules.download_ctat_lib.output.tar
    output:
        done5 = touch("data/references/ctat/GRCh38_gencode_v44_CTAT_lib_Oct292023/.extract_done")
    log:
        "logs/extract_ctat_lib.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p data/references/ctat
        tar -xzf {input.tar} -C data/references/ctat > {log} 2>&1
        """