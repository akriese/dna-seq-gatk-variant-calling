rule select_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/genotyped/{caller}/all.vcf.gz",
        tbi="results/genotyped/{caller}/all.vcf.gz.tbi"
    output:
        vcf=temp("results/filtered/all.{vartype}_{caller}.vcf.gz")
    log:
        "logs/filtering/select_calls/{caller}/{vartype}.log"
    benchmark:
        "benchmarks/filtering/select_calls/{caller}/{vartype}.benchmark"
    params:
        extra=get_vartype_arg
    wrapper:
        "0.59.0/bio/gatk/selectvariants"


rule hard_filter_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/filtered/all.{vartype}_{caller}.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}_{caller}.hardfiltered.vcf.gz")
    log:
        "logs/filtering/hard_filter_calls/{caller}/{vartype}.log"
    benchmark:
        "benchmarks/filtering/hard_filter_calls/{caller}/{vartype}.benchmark"
    params:
        filters=get_filter
    wrapper:
        "0.74.0/bio/gatk/variantfiltration"


rule recalibrate_calls:
    input:
        vcf="results/filtered/all.{vartype}_{caller}.vcf.gz"
    output:
        vcf=temp("results/filtered/all.{vartype}_{caller}.recalibrated.vcf.gz")
    log:
        "logs/filtering/recalibrate_calls/{caller}/{vartype}.log"
    benchmark:
        "benchmarks/filtering/recalibrate_calls/{caller}/{vartype}.benchmark"
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]
    threads: 40
    resources:
        mem_mb = 50000
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
        vcf="results/filtered/all.{caller}.vcf.gz"
    log:
        "logs/filtering/merge_calls/{caller}.log"
    benchmark:
        "benchmarks/filtering/merge_calls/{caller}.benchmark"
    resources:
        mem_mb = 400000
    wrapper:
        "0.74.0/bio/picard/mergevcfs"

rule merge_technologies:
    input:
        ref="resources/genome.fasta",
        vcfs = ["results/filtered/all.gatk.vcf.gz",
                "results/filtered/all.freebayes.vcf.gz"]
    output:
        vcf="results/filtered/all.combined.vcf.gz"
    log:
        "logs/filtering/merge_technologies/all.log"
    benchmark:
        "benchmarks/filtering/merge_technologies/all.benchmark"
    params:
        extra="",
    resources:
        mem_mb=10000
    wrapper:
        "v1.12.0/bio/picard/mergevcfs"
