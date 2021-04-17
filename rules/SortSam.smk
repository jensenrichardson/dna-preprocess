rule SortSam:
    input:
        "02-mapping/{sample}/{readgroup}.aligned.sam"
    output:
        temp("03-Sorted/{sample}/{readgroup}.sorted.bam"),
    log:
        "03-Sorted/{sample}/{readgroup}.log"
    resources:
        cores=16,
	runtime=lambda wildcards, attempt: 30 * attempt
    shell:
        "gatk SortSam "
        "-I {input} "
        "-O {output} "
        "-SO coordinate "
        "&> {log}"
