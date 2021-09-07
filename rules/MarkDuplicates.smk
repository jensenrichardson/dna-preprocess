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
        bam=temp("04-Markdup/{sample}.markdup.bam"),
	metrics="04-Markdup/{sample}.metrics"
    log:
        "04-Markdup/{sample}.log"
    params:
        lambda wildcards, input: '-I ' + ' -I '.join(input)
    resources:
        cores=16,
	runtime=lambda wildcards, attempt: 30 + 60 * attempt
    shell:
        "gatk MarkDuplicates "
        "{params} "
        "-O {output.bam} "
	"-M {output.metrics} "
        "&> {log}"
