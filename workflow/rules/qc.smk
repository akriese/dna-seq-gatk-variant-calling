rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="results/qc/fastqc/{sample}-{unit}.html",
        zip="results/qc/fastqc/{sample}-{unit}.zip",
    log:
        "logs/qc/fastqc/{sample}-{unit}.log"
    benchmark:
        "benchmarks/qc/fastqc/{sample}-{unit}.benchmark"
    wrapper:
        "0.74.0/bio/fastqc"


rule samtools_stats:
    input:
        "results/recal/{sample}-{unit}.bam"
    output:
        "results/qc/samtools-stats/{sample}-{unit}.txt"
    log:
        "logs/qc/samtools_stats/{sample}-{unit}.log"
    benchmark:
        "benchmarks/qc/samtools_stats/{sample}-{unit}.benchmark"
    wrapper:
        "0.74.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(
            [
                "results/qc/samtools-stats/{u.sample}-{u.unit}.txt",
                "results/qc/fastqc/{u.sample}-{u.unit}.zip",
                "results/qc/dedup/{u.sample}-{u.unit}.metrics.txt",
            ],
            u=units.itertuples()
        )
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control"
        )
    log:
        "logs/qc/multiqc/all.log"
    benchmark:
        "benchmarks/qc/multiqc/all.benchmark"
    wrapper:
        "0.74.0/bio/multiqc"
