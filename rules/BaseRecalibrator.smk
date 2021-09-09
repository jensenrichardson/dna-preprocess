configfile: "config.yaml"

rule BaseRecalibration:
    input:
        bam="04-Markdup/{sample}.markdup.bam",
        ref=config["ref_gen"]
    output:
        table="05-BaseRecalibrator/{sample}.table",
    params:
        known_sites="--known-sites " + " --known-sites ".join(config["known_sites"])
    log:
        "05-BaseRecalibrator/{sample}.log"
    resources:
        cores=16,
	runtime=lambda wildcards, attempt: 60 * attempt + 240 * ( attempt - 1 )
    shell:
        "OMP_NUM_THREADS={resources.cores} "
        "gatk BaseRecalibrator "
        "-R {input.ref} "
        "-I {input.bam} "
	"{params.known_sites} "
        "-O {output.table} "
        "&> {log}"

