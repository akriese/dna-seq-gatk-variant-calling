if "restrict-regions" in config["processing"]:

    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "results/called/{contig}.regions.bed"
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"


rule gatk_snv_caller:
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
        )
    output:
        gvcf="results/called/gatk/{sample}.{contig}.g.vcf.gz"
    log:
        "logs/calling/gatk_snv_caller/{sample}/{contig}.log"
    benchmark:
        "benchmarks/calling/gatk_snv_caller/{sample}/{contig}.benchmark"
    threads: 5
    resources:
        mem_mb = 20000, # sometimes hits 7g, give it 20
    params:
        extra=get_call_variants_params
    wrapper:
        "v1.10.0/bio/gatk/haplotypecaller"


rule freebayes_snv_caller:
    input:
        samples = get_all_sample_bams(),
        ref = "resources/genome.fasta",
        # optional BED file specifying chromosomal regions on which freebayes
        # should run, e.g. all regions that show coverage
        regions = config["processing"].get("restrict-regions", [])
    output:
        "results/genotyped/freebayes/all.vcf"
    log:
        "logs/calling/freebayes_snv_caller/all.log"
    benchmark:
        "benchmarks/calling/freebayes_snv_caller/all.benchmark"
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: lambda wc: min(100, len(get_all_sample_bams())*30)  # TODO with many samples, split this job up into multiple
    resources:
        mem_mb = 500000,
    wrapper:
        "v1.3.2/bio/freebayes"


rule combine_calls:
    input:
        ref="resources/genome.fasta",
        gvcfs=expand(
            "results/called/gatk/{sample}.{{contig}}.g.vcf.gz", sample=samples.index
        ),
    output:
        gvcf="results/called/gatk/all.{contig}.g.vcf.gz"
    log:
        "logs/calling/combine_calls/{contig}.log"
    benchmark:
        "benchmarks/calling/combine_calls/{contig}.benchmark"
    resources:
        mem_mb=10000
    wrapper:
        "0.74.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref="resources/genome.fasta",
        gvcf="results/called/gatk/all.{contig}.g.vcf.gz"
    output:
        vcf=temp("results/genotyped/gatk/all.{contig}.vcf.gz")
    log:
        "logs/calling/genotype_variants/{contig}.log"
    benchmark:
        "benchmarks/calling/genotype_variants/{contig}.benchmark"
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    resources:
        mem_mb=10000
    wrapper:
        "0.74.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand(
            "results/genotyped/gatk/all.{contig}.vcf.gz", contig=get_contigs()
        )
    output:
        vcf="results/genotyped/gatk/all.vcf.gz"
    log:
        "logs/calling/merge_variants/all.log"
    benchmark:
        "benchmarks/calling/merge_variants/all.benchmark"
    resources:
        mem_mb=10000
    wrapper:
        "0.74.0/bio/picard/mergevcfs"

rule vcf_gzip:
    input:
        "results/genotyped/{prefix}.vcf"
    output:
        "results/genotyped/{prefix}.vcf.gz"
    log:
        "logs/calling/vcf_gzip/{prefix}.log"
    benchmark:
        "benchmarks/calling/vcf_gzip/{prefix}.benchmark"
    threads: 8
    shell:
        """
        (bgzip -f --threads {threads} {input} > {output}) &> {log}
        """

