GiB = int(1024**3)
MiB = int(1024**2)

localrules: link_or_map, make_link_to_mapped_reads

rule trim_reads_se:
    input:
        unpack(get_fastq)
    output:
        temp("results/trimmed/{sample}-{unit}.fastq.gz")
    log:
        "logs/mapping/trim_reads_se/{sample}-{unit}.log"
    benchmark:
        "benchmarks/mapping/trim_reads_se/{sample}-{unit}.benchmark"
    params:
        **config["params"]["trimmomatic"]["se"],
        extra=""
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
        trimlog="results/trimmed/{sample}-{unit}.trimlog.txt"
    log:
        "logs/mapping/trim_reads_pe/{sample}-{unit}.log"
    benchmark:
        "benchmarks/mapping/trim_reads_pe/{sample}-{unit}.benchmark"
    params:
        **config["params"]["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {}".format(output.trimlog)
    wrapper:
        "0.74.0/bio/trimmomatic/pe"

rule map_reads_from_fastq:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        "results/mapped/{sample}-{unit}.sorted.bam"
    log:
        "logs/mapping/map_reads_from_fastq/{sample}-{unit}.log"
    benchmark:
        "benchmarks/mapping/map_reads_from_fastq/{sample}-{unit}.benchmark"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 48
    resources:
        # mem_mb = 35000, # sometimes hits 30g
        time_min = 1000, # ~16h should suffice
    wrapper:
        "0.74.0/bio/bwa/mem"

rule make_link_to_mapped_reads:
    input:
        get_unit_bam
    output:
        "results/mapped/{sample}-{unit}_link.sorted.bam"
    shell:
        " ln -s $(readlink -e {input}) {output} "

# this might be a bit too hacky...
rule link_or_map:
    input:
        link_or_mapped_output
    output:
        rules.map_reads_from_fastq.output
    shell:
        """ [[ {input} == {output} ]] || mv {input} {output} """


rule mark_duplicates:
    input:
        bams="results/mapped/{sample}-{unit}.sorted.bam",
    output:
        bam="results/dedup/{sample}-{unit}.bam",
        metrics="results/qc/dedup/{sample}-{unit}.metrics.txt"
    log:
        "logs/mapping/mark_duplicates/{sample}-{unit}.log"
    benchmark:
        "benchmarks/mapping/mark_duplicates/{sample}-{unit}.benchmark"
    params:
        extra = config["params"]["picard"]["MarkDuplicates"]
    resources:
        time_min = 960, # usually takes 3-5h, make it 16
        # mem_mb = 50000
        mem_mb = lambda _wc, input: (input.size//MiB) * 3, # 3x input size
        tmpdir_gb = lambda _wc, input: (input.size//GiB) * 2 #
    wrapper:
        "v1.14.0/bio/picard/markduplicates"

"""

rule sambamba_markdup:
    input:
        "results/mapped/{sample}-{unit}.sorted.bam"
    output:
        "results/dedup/{sample}-{unit}.bam"
    log:
        "logs/mapping/sambamba_markdup/{sample}-{unit}.log"
    benchmark:
        "benchmarks/mapping/sambamba_markdup/{sample}-{unit}.benchmark"
    params:
        extra=("--overflow-list-size=600000 " +
            "--hash-table-size=600000 " +
            "-p " +
            ("-r" if config["params"]["picard"]["MarkDuplicates"] == "REMOVE_DUPLICATES=true" else ""))
    threads: 16
    resources:
        mem_mb = 10000, # 10g, usually at 3g
        time_min = 300, # 5h should suffice
        tmpdir_gb = 100 # 100g cluster tmp dir size
    wrapper:
        "v1.8.0/bio/sambamba/markdup"
"""


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        known="resources/variation.noiupac.vcf.gz",
        known_idx="resources/variation.noiupac.vcf.gz.tbi"
    output:
        recal_table="results/recal/{sample}-{unit}.grp"
    log:
        "logs/mapping/recalibrate_base_qualities/{sample}-{unit}.log"
    benchmark:
        "benchmarks/mapping/recalibrate_base_qualities/{sample}-{unit}.benchmark"
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    resources:
        mem_mb = 40000, # usually takes up to 15g
        time_min = 1440 # usually takes about 6h, make it 24h
    wrapper:
        "0.74.0/bio/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        recal_table="results/recal/{sample}-{unit}.grp"
    output:
        bam="results/recal/{sample}-{unit}.bam"
    log:
        "logs/mapping/apply_base_quality_recalibration/{sample}-{unit}.log"
    benchmark:
        "benchmarks/mapping/apply_base_quality_recalibration/{sample}-{unit}.benchmark"
    params:
        extra=get_regions_param()
    resources:
        time_min = 960 # usually takes up to 4h, give it 16
    #     mem_mb=40000,
    wrapper:
        "0.74.0/bio/gatk/applybqsr"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "logs/mapping/samtools_index/{prefix}.log"
    benchmark:
        "benchmarks/mapping/samtools_index/{prefix}.benchmark"
    wrapper:
        "0.74.0/bio/samtools/index"
