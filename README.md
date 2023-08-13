# Snakemake workflow: dna-preprocess
For preprocessing dna fastq files according to the gatk preprocessing pipeline using Snakemake.

## Quick Start
This quick start guide assumes you know how to configue a snakemake workflow for deployment on a cluster.
1. Fill out the required fields in `config.yaml`. Ensure that your references has bwa index (`bwa index <reference>`)and gatk indices (`gatk CreateSequenceDictionary -R <reference>` and `samtools faidx <reference>`).
2. Go over the resources required for each rule in the `rules` directory. Change these as required.
3. Make sure that you have the dependencies managed properly. Dependencies for this workflow include:
    * Bwa Mem 2 (you can use bwa mem just change the command from `bwa-mem2 mem` to just `bwa mem` in the file `rules/Bwa2.smk`).
        * Please read [the documentation](https://github.com/bwa-mem2/bwa-mem2) for bwa mem 2.
        It has significantly different system and *especially memory* requirements.
    * Gatk 4
4. Use the `parse_samples.py` script to parse your samples of the form `SAMPLE_READGROUP...R1...<extension>`.
The script should be used as `./parse_samples.py -b <fastq-directory>`.
5. Preview your run with `snakemake -npr`. Please look at at least one of each command to make sure that your samples were processed.
6. Queue your jobs with `snakemake -c1` or `snakemake --profile <profile>` (recommended).

## Detailed(ish) Instructions
### A note on `parse_samples.py`
This is the most important part of the whole repository.
Using this script saves hours of formatting and config file headaches and it is highly recommended that you use it, not only for this workflow but in other places.
It has several different output formats.
It can parse samples of the format that are usable in any snakemake workflow using
```
import pandas as pd
import ast
configfile: "config.yaml"
samples = pd.read_table(config["samples_tsv"], converters={"files": ast.literal_eval}).set_index("sample_name", drop=False)
wildcard_constraints:
    sample ="|".join(samples.index.tolist())
```
See `rules/Bwa2.smk` for examples on how to leverage the properties of the DataFrame in your workflows.

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in `config.yaml` to configure the workflow execution.

Execute the `parse_samples.py` script:

    `./parse_samples.py -b <fastq-directory>`
    
to parse your samples of the form `SAMPLE_READGROUP...R1...<extension>`. Essentially the only requirement is that the sample and readgroup are separated by an underscore and R1 and R2 are differentiated with an 'R1' and an 'R2' after the readgroup.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally (not recommended) via

    snakemake --cores $N

using `$N` cores.

Conda integration is not provided because I have found it does not play nicely with clusters.

It is recomended that you use singularity to do any sort of actual work or testing on an actual cluster and I have provided docker containers for each rule. The containers do not provide any of the files requested by the `config.yaml` and without those files (and indices of each of them) the workflow will not run.

    snakemake --use-singularity

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.
An example (using some trivial test data) can be seen [here](https://cdn.rawgit.com/snakemake-workflows/rna-seq-kallisto-sleuth/master/.test/report.html).

### Step 6: Commit changes

Whenever you change something, don't forget to commit the changes back to your github copy of the repository:

    git commit -a
    git push

### Step 7: Obtain updates from upstream

Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: `git remote add -f upstream git@github.com:snakemake-workflows/neoant-pred.git` or `git remote add -f upstream https://github.com/snakemake-workflows/neoant-pred.git` if you do not have setup ssh keys.
2. Update the upstream version: `git fetch upstream`.
3. Create a diff with the current version: `git diff HEAD upstream/master workflow > upstream-changes.diff`.
4. Investigate the changes: `vim upstream-changes.diff`.
5. Apply the modified diff via: `git apply upstream-changes.diff`.
6. Carefully check whether you need to update the config files: `git diff HEAD upstream/master config`. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.


### Step 8: Contribute back

In case you have also changed or added steps, please consider contributing them back to the original repository:

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the original repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., `cp -r workflow path/to/fork`. Make sure to **not** accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary.
4. Commit and push your changes to your fork.
5. Create a [pull request](https://help.github.com/en/articles/creating-a-pull-request) against the original repository.


## Authors

* Jensen Richardson (@jensenrichardson)