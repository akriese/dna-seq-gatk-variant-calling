if "restrict-regions" in config["processing"]:

    rule compose_regions:
        input:
            config["processing"]["restrict-regions"],
        output:
            "results/called/{contig}.regions.bed",
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
        regions=(
            "results/called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else []
        ),
    output:
        gvcf="results/called/gatk/{sample}.{contig}.g.vcf.gz",
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log",
    threads: 16
    resources:
        mem_mb=10000
    params:
        extra=get_call_variants_params,
    benchmark:
        "benchmarks/results/gatk/haplotypecaller/{sample}_{contig}.benchmark"
    wrapper:
        "0.59.0/bio/gatk/haplotypecaller"


rule freebayes_snv_caller:
    input:
        samples = get_all_sample_bams(),
        ref = "resources/genome.fasta",
        # optional BED file specifying chromosomal regions on which freebayes
        # should run, e.g. all regions that show coverage
        regions = config["processing"].get("restrict-regions")
    output:
        "results/genotyped/freebayes/all.vcf",
    log:
        "logs/results/freebayes/all.log",
    benchmark:
        "benchmarks/results/freebayes/all.benchmark",
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 50
    resources:
        mem_mb = 100000,
    wrapper:
        "v1.3.2/bio/freebayes"


rule combine_calls:
    input:
        ref="resources/genome.fasta",
        gvcfs=expand(
            "results/called/gatk/{sample}.{{contig}}.g.vcf.gz", sample=samples.index
        ),
    output:
        gvcf="results/called/gatk/all.{contig}.g.vcf.gz",
    log:
        "logs/gatk/combinegvcfs.gatk_{contig}.log",
    resources:
        mem_mb=10000
    benchmark:
        "benchmarks/results/gatk/combinegvcfs/gatk_{contig}.benchmark"
    wrapper:
        "0.74.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref="resources/genome.fasta",
        gvcf="results/called/gatk/all.{contig}.g.vcf.gz",
    output:
        vcf=temp("results/genotyped/gatk/all.{contig}.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
    log:
        "logs/gatk/genotypegvcfs.gatk_{contig}.log",
    resources:
        mem_mb=10000
    benchmark:
        "benchmarks/results/gatk/genotypegvcfs/gatk_{contig}.benchmark"
    wrapper:
        "0.74.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand(
            "results/genotyped/gatk/all.{contig}.vcf.gz", contig=get_contigs()
        ),
    output:
        vcf="results/genotyped/gatk/all.vcf.gz",
    log:
        "logs/picard/merge-genotyped_gatk.log",
    resources:
        mem_mb=10000
    benchmark:
        "benchmarks/results/picard/mergevcfs_gatk.benchmark"
    wrapper:
        "0.74.0/bio/picard/mergevcfs"

rule vcf_gzip:
    input:
        "{prefix}.vcf"
    output:
        "{prefix}.vcf.gz"
    threads: 8
    log:
        "logs/bgzip/{prefix}.benchmark"
    shell:
        """
        (bgzip -f --threads {threads} {input} > {output}) &> {log}
        """

