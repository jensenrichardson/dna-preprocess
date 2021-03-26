import pandas as pd
import ast
configfile: "config.yaml"
samples = pd.read_table(config["samples_tsv"], converters={"files": ast.literal_eval}).set_index("sample_name", drop=False)
wildcard_constraints:
    sample ="|".join(samples.index.tolist())
#print(expand("02-mapping/ALMC1/{readgroup}.sam", readgroup=list(samples.loc["ALMC1", "files"].keys())))
rule MarkDuplicates:
    input:
        lambda wildcards: expand("03-Sorted/{sample}/{readgroup}.sorted.bam", sample=wildcards.sample, readgroup=list(samples.loc[wildcards.sample, "files"].keys()))
    output:
        bam="04-Markdup/{sample}.markdup.bam",
        bai="04-Markdup/{sample}.markdup.bam.bai",
        sbi="04-Markdup/{sample}.markdup.bam.sbi"
    log:
        "03-markdup/{sample}.log"
    params:
        lambda wildcards, input: '-I ' + ' -I '.join(input)
    resources:
        cores=16,
	    runtime=45
    shell:
        "gatk MarkDuplicates "
        "{params}"
        "-O {output.bam} "
        "&> {log}"
