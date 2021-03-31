import pandas as pd
import ast
configfile: "config.yaml"
samples = pd.read_table(config["samples_tsv"], converters={"files": ast.literal_eval}).set_index("sample_name", drop=False)
wildcard_constraints:
    sample ="|".join(samples.index.tolist())

rule all:
	#input: expand("02-mapping/{sample}/{sample}.Aligned.out.bam", sample=samples.to_dict('index'))
    input: expand("06-ApplyRecalibration/{sample}.recalibrated.bam", sample=samples.to_dict('index'))
	#input: expand("04-SplitCigar/{sample}.splitcigar.bam", sample=samples.to_dict('index'))
	#input: expand("06-ApplyRecalibration/{sample}.recalibrated.bam", sample=samples.to_dict('index'))
	#input: expand("07-HaplotypeCaller/{sample}.hapcall.vcf", sample=samples.to_dict('index'))

include: "rules/Bwa2.smk"
include: "rules/MarkDuplicates.smk"
include: "rules/SortSam.smk"
include: "rules/BaseRecalibrator.smk"
include: "rules/ApplyRecalibration.smk"
