rule get_genome:
    output:
        "resources/genome.fasta"
    log:
        "logs/ref/get_genome/all.log"
    benchmark:
        "benchmarks/ref/get_genome/all.benchmark"
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-sequence"


checkpoint genome_faidx:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.fasta.fai"
    log:
        "logs/ref/genome_faidx/all.log"
    cache: True
    wrapper:
        "0.74.0/bio/samtools/faidx"


rule genome_dict:
    input:
        "resources/genome.fasta"
    output:
        "resources/genome.dict"
    log:
        "logs/ref/genome_dict/all.log"
    benchmark:
        "benchmarks/ref/genome_dict/all.benchmark"
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


rule get_known_variation:
    input:
        # use fai to annotate contig lengths for GATK BQSR
        fai="resources/genome.fasta.fai"
    output:
        vcf="resources/variation.vcf.gz"
    log:
        "logs/ref/get_known_variation/all.log"
    benchmark:
        "benchmarks/ref/get_known_variation/all.benchmark"
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        type="all"
    cache: True
    wrapper:
        "0.74.0/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz"
    output:
        "resources/variation.noiupac.vcf.gz"
    log:
        "logs/ref/remove_iupac_codes/all.log"
    benchmark:
        "benchmarks/ref/remove_iupac_codes/all.benchmark"
    conda:
        "../envs/rbt.yaml"
    cache: True
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule tabix_vcf:
    input:
        "{file_to_index}.vcf.gz"
    output:
        "{file_to_index}.vcf.gz.tbi"
    log:
        "logs/ref/tabix_vcf/{file_to_index}.log"
    benchmark:
        "benchmarks/ref/tabix_vcf/{file_to_index}.benchmark"
    params:
        "-p vcf",
    wrapper:
        "0.74.0/bio/tabix"


rule bwa_index:
    input:
        "resources/genome.fasta"
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/ref/bwa_index/all.log"
    benchmark:
        "benchmarks/ref/bwa_index/all.benchmark"
    threads: 20
    resources:
        mem_mb=369000
    cache: True
    wrapper:
        "0.74.0/bio/bwa/index"


rule get_vep_cache:
    output:
        directory("resources/vep/cache")
    log:
        "logs/ref/get_vep_cache/all.log"
    benchmark:
        "benchmarks/ref/get_vep_cache/all.benchmark"
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"]
    wrapper:
        "0.74.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins")
    log:
        "logs/ref/get_vep_plugins/all.log"
    benchmark:
        "benchmarks/ref/get_vep_plugins/all.benchmark"
    params:
        release=config["ref"]["release"]
    wrapper:
        "0.74.0/bio/vep/plugins"
