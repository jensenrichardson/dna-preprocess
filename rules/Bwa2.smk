import pandas as pd
import ast
configfile: "config.yaml"
samples = pd.read_table(config["samples_tsv"], converters={"files": ast.literal_eval}).set_index("sample_name", drop=False)
wildcard_constraints:
    sample ="|".join(samples.index.tolist())

rule bwa_map:
    input:
        lambda wildcards: samples.loc[wildcards.sample, "files"],
        genome=config["star_genome"]
    params:
        command=lambda wildcards: samples.loc[wildcards.sample, "command"]
    output:
        bam="02-mapping/{sample}/{sample}.Aligned.out.bam"
    log:
        fin="02-mapping/{sample}/{sample}.Log.final.out",
        full="02-mapping/{sample}/{sample}.Log.out",
        prog="02-mapping/{sample}/{sample}.Log.progress.out",
        sjo="02-mapping/{sample}/{sample}.SJ.out.tab",
        sji="02-mapping/{sample}/{sample}._STARgenome/sjdbInfo.txt",
        sjl="02-mapping/{sample}/{sample}._STARgenome/sjdbList.out.tab",
        pass1="02-mapping/{sample}/{sample}._STARpass1/Log.final.out",
        pass1s="02-mapping/{sample}/{sample}._STARpass1/SJ.out.tab"
    resources:
        runtime=120,
	cores=48
    shell:
        "STAR "
        "--runThreadN {resources.cores} " 
        "--genomeDir {input.genome} "
        "{params.command} " 
        "--outSAMtype BAM Unsorted " 
        "--twopassMode Basic "
        "--outFileNamePrefix ./02-mapping/{wildcards.sample}/{wildcards.sample}. &> /dev/null"
