rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        temp("results/trimmed/{sample}-{unit}.fastq.gz"),
    params:
        **config["params"]["trimmomatic"]["se"],
        extra="",
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    benchmark:
        "benchmarks/results/trim/se_{sample}_{unit}.benchmark"
    wrapper:
        "0.74.0/bio/trimmomatic/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp("results/trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp("results/trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp("results/trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog="results/trimmed/{sample}-{unit}.trimlog.txt",
    params:
        **config["params"]["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    benchmark:
        "benchmarks/results/trim/pe_{sample}_{unit}.benchmark"
    wrapper:
        "0.74.0/bio/trimmomatic/pe"


"""
rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        "results/mapped/{sample}-{unit}.sorted.bam",
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: 48
    resources:
        mem_mb = 20000
    benchmark:
        "benchmarks/results/bwa_mem/{sample}_{unit}.benchmark"
    wrapper:
        "0.74.0/bio/bwa/mem"


"""
rule make_link_to_mapped_reads:
    input:
        "resources/bam_files/{sample}-{unit}.sorted.bam"
    output:
        "results/mapped/{sample}-{unit}.sorted.bam"
    shell:
        " ln -s $(readlink -e {input}) {output} "

"""
"""
rule mark_duplicates:
    input:
        "results/mapped/{sample}-{unit}.sorted.bam",
    output:
        bam="results/dedup/{sample}-{unit}.bam",
        metrics="results/qc/dedup/{sample}-{unit}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}-{unit}.log",
    params:
        config["params"]["picard"]["MarkDuplicates"],
    threads: 30
    resources:
        mem_mb = 40000,
    benchmark:
        "benchmarks/results/picard/markduplicates/{sample}_{unit}.benchmark"
    wrapper:
        "0.74.0/bio/picard/markduplicates"

"""

rule sambamba_markdup:
    input:
        "results/mapped/{sample}-{unit}.sorted.bam"
    output:
        "results/dedup/{sample}-{unit}.bam"
    params:
        extra=("--overflow-list-size=600000 " +
            "--hash-table-size=600000 " +
            "-p " +
            ("-r" if config["params"]["picard"]["MarkDuplicates"] == "REMOVE_DUPLICATES=true" else ""))
    log:
        "logs/sambamba-markdup/{sample}-{unit}.log"
    benchmark:
        "benchmarks/results/sambamba/markduplicates/{sample}_{unit}.benchmark"
    threads: 16
    wrapper:
        "v1.8.0/bio/sambamba/markdup"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        known_idx="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table="results/recal/{sample}-{unit}.grp",
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    threads: 20
    resources:
        mem_mb=40000,
    benchmark:
        "benchmarks/results/gatk/bqsr/{sample}_{unit}.benchmark"
    wrapper:
        "0.74.0/bio/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        recal_table="results/recal/{sample}-{unit}.grp",
    output:
        bam="results/recal/{sample}-{unit}.bam",
    log:
        "logs/gatk/apply-bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param(),
    threads: 20
    resources:
        mem_mb=40000,
    benchmark:
        "benchmarks/results/gatk/apply-bqsr/{sample}_{unit}.benchmark"
    wrapper:
        "0.74.0/bio/gatk/applybqsr"


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools/index/{prefix}.log",
    benchmark:
        "benchmarks/results/samtools/index_{prefix}.benchmark"
    wrapper:
        "0.74.0/bio/samtools/index"
