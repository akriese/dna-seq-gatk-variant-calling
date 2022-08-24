rule select_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/genotyped/{caller}/all.vcf.gz",
        tbi="results/genotyped/{caller}/all.vcf.gz.tbi",
    output:
        vcf=temp("results/filtered/all.{vartype}_{caller}.vcf.gz"),
    params:
        extra=get_vartype_arg,
    log:
        "logs/gatk/selectvariants/{vartype}_{caller}.log",
    wrapper:
        "0.59.0/bio/gatk/selectvariants"


rule hard_filter_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/filtered/all.{vartype}_{caller}.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}_{caller}.hardfiltered.vcf.gz"),
    params:
        filters=get_filter,
    log:
        "logs/gatk/variantfiltration/{vartype}_{caller}.log",
    wrapper:
        "0.74.0/bio/gatk/variantfiltration"


rule recalibrate_calls:
    input:
        vcf="results/filtered/all.{vartype}_{caller}.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}_{caller}.recalibrated.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"],
    log:
        "logs/gatk/variantrecalibrator/{vartype}_{caller}.log",
    threads: 40
    resources:
        mem_mb = 50000
    benchmark:
        "benchmarks/results/gatk/variantrecalibrator/{vartype}_{caller}.benchmark"
    wrapper:
        "0.74.0/bio/gatk/variantrecalibrator"


rule merge_calls:
    input:
        vcfs=expand(
            "results/filtered/all.{vartype}_{{caller}}.{filtertype}.vcf.gz",
            vartype=["snvs", "indels"],
            filtertype="recalibrated"
            if config["filtering"]["vqsr"]
            else "hardfiltered",
        ),
    output:
        vcf="results/filtered/all.{caller}.vcf.gz",
    log:
        "logs/picard/merge-filtered.{caller}.log",
    benchmark:
        "benchmarks/results/picard/mergevcfs.{caller}.benchmark"
    wrapper:
        "0.74.0/bio/picard/mergevcfs"
