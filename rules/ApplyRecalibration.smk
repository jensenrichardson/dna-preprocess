configfile: "config.yaml"

rule ApplyCalibration:
    input:
        bam="04-Markdup/{sample}.markdup.bam",
        ref=config["ref_gen"],
        table="05-BaseRecalibrator/{sample}.table"
    output:
        bam=protected("06-ApplyRecalibration/{sample}.recalibrated.bam"),
        bai="06-ApplyRecalibration/{sample}.recalibrated.bai"
    log:
        "06-ApplyRecalibration/{sample}.log"
    container: 
        "broadinstitute/gatk"
    resources:
        cores=16,
	runtime=lambda wildcards, attempt: 15 + 30 * attempt + 60 * (attempt - 1)
    shell:
        "gatk ApplyBQSR "
        "-R {input.ref} "
        "-I {input.bam} "
	"--bqsr-recal-file {input.table} "
        "-O {output.bam} "
        "&> {log}"
