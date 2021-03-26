rule MarkDuplicates:
    input:
        "02-mapping/{sample}/{sample}.Aligned.out.bam"
    output:
        bam="03-markdup/{sample}.markdup.bam",
        bai="03-markdup/{sample}.markdup.bam.bai",
        sbi="03-markdup/{sample}.markdup.bam.sbi"
    log:
        "03-markdup/{sample}.log"
    resources:
        cores=16,
	runtime=45
    shell:
        "gatk MarkDuplicates "
        "-I {input} "
        "-O {output.bam} "
        "&> {log}"
