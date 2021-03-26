import pandas as pd
import ast
configfile: "config.yaml"
samples = pd.read_table(config["samples_tsv"], converters={"files": ast.literal_eval}).set_index("sample_name", drop=False)
wildcard_constraints:
    sample ="|".join(samples.index.tolist())
#print(expand("02-mapping/ALMC1/{readgroup}.sam", readgroup=list(samples.loc["ALMC1", "files"].keys())))
rule SortSam:
    input:
        "02-mapping/{sample}/{readgroup}.aligned.sam"
    output:
        "03-Sorted/{sample}/{readgroup}.sorted.bam",
    log:
        "03-Sorted/{sample}/{readgroup}.log"
    resources:
        cores=16,
	    runtime=45
    shell:
        "gatk SortSam "
        "-I {input} "
        "-O {output} "
        "-SO coordinate "
        "&> {log}"
