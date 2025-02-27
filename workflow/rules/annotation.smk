rule annotate_variants:
    input:
        calls="results/filtered/all.combined.vcf.gz",
        cache="resources/vep/cache",
        plugins="resources/vep/plugins"
    output:
        calls=report(
            "results/annotated/all.vcf.gz",
            caption="../report/vcf.rst",
            category="Calls"
        ),
        stats=report(
            "results/stats/all.stats.html",
            caption="../report/stats.rst",
            category="Calls"
        )
    log:
        "logs/annotation/annotate_variants/all.log"
    benchmark:
        "benchmarks/annotation/annotate_variants/all.benchmark"
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=config["params"]["vep"]["plugins"],
        extra=config["params"]["vep"]["extra"]
    threads: 4
    resources:
        time_min = lambda w: max(960, len(get_all_sample_bams())*180) # give at least 16h, otherwise #samples*3h
    wrapper:
        "0.74.0/bio/vep/annotate"

