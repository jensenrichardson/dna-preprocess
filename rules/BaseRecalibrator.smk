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
	runtime=lambda wildcards, attempt: 30 * attempt
    shell:
        "gatk BaseRecalibrator "
        "-R {input.ref} "
        "-I {input.bam} "
	"{params.known_sites} "
        "-O {output.table} "
        "&> {log}"

