import pandas as pd
import ast
configfile: "config.yaml"
samples = pd.read_table(config["samples_tsv"], converters={"files": ast.literal_eval}).set_index("sample_name", drop=False)
wildcard_constraints:
    sample ="|".join(samples.index.tolist())
rule bwa_map:
    input:
        fastq=lambda wildcards: samples.loc[wildcards.sample, "files"][wildcards.readgroup][:-1]
    params:
        command=lambda wildcards: samples.loc[wildcards.sample, "files"][wildcards.readgroup][-1]
    output:
        sam="02-mapping/{sample}/{readgroup}.aligned.sam"
    log:
        "02-mapping/{sample}/{readgroup}.log",
    resources:
        runtime=120,
        cores=48
    shell:
        "bwa-mem2 mem "
        "-t {resources.cores} "
        "-o {output.sam} "
        "{params.command} "
        "{input.fastq} "
        "&> {log}"
