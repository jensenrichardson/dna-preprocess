#!/usr/bin/env python3
# This script is responsible for parsing a directory of
# fastq samples to produce a tsv that is easy to manage
# with snakemake or nextflow.
from os import listdir
from os.path import isfile, join
import os
from pathlib import Path
# for regex
import re
# for arguments
import argparse
import pandas as pd
import glob


# Defines class Sample with name string and readgroups list
class Sample:
    """Sample class is for representing a single genetic samples with all readgroups.

    Attributes:
        name (string): Name of the sample
        readgroups (list): List of the readgroups associated with the sample. List of Readgroup objects
        command (string): Optional command associated with sample. Used with STAR mapper. 
    """
    def __init__(self, name=""):
        self.name = name
        self.readgroups = []
        self.command = []

    def __repr__(self):
        rep = f'Sample(name: {self.name}, readgroups: {self.readgroups}, command: {self.command})'
        return rep

    def addReadgroup(self, rg):
        self.readgroups.append(
            rg) if rg not in self.readgroups else self.readgroups


class Readgroup:
    """Defines Readgroup Object to represent readgroup

    Attributes:
        rg (string): Name of readgroup
        r1 (file): Read 1 file of readgroup
        r2 (file): Read 2 file of readgroup
        command (string): Readgroup string for readgroup (bwa). Includes illumina platform

    Returns:
        _type_: _description_
    """
    # name of readgroup
    def __init__(self, rg=""):
        self.rg = rg
        self._r1 = ""
        self._r2 = ""
        self.command = ""

    # Equality
    def _is_valid_operand(self, other):
        return hasattr(other, "rg")

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.rg == other.rg

    # better printing
    def __repr__(self):
        rep = f'Readgroup(rg: {self.rg}, r1: {self.r1}, r2: {self.r2}, command: {self.command})'
        return rep

    # getters and setters for r1 and r2 that ensure that the file exists
    @property
    def r1(self):
        return self._r1

    @r1.setter
    def r1(self, r1):
        if not os.path.isfile(r1):
            raise OSError(f'{r1} is not a file')
        self._r1 = r1

    @property
    def r2(self):
        return self._r2

    @r2.setter
    def r2(self, r2):
        if not os.path.isfile(r2):
            raise OSError(f'{r2} is not a file')
        self._r2 = r2


# Used for ensuring directory provided as argument is a valid directory
def dir_path(string: str):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)


def main(args):
    verb = args.verbose
    directory = args.fastq_dir
    directory = Path(directory).resolve()
    fastq_files = [f for f in listdir(directory) if isfile(join(directory, f))]
    samples = []
    if not args.input_samples:
        samples = get_samples(fastq_files)
    elif args.input_samples:
        samples_df = pd.read_table(args.input_samples)
        samples_list = samples_df[samples_df.columns[0]]
        samples = [Sample(sample) for sample in samples_list]
    # gets list of fastq files
    # sets up samples with the format:
    # Sample(
    #   Name: str
    #   Readgroups: [
    #       Readgroup(
    #           rg: str
    #           r1: str
    #           r2: str
    # )])
    samples = get_readgroups(samples, fastq_files, directory, verb)
    if verb > 1:
        for s in samples:
            print(f'{s.name} has {s.readgroups} readgroups')
    if verb > 0:
        for s in samples:
            print(f'{s.name} has {len(s.readgroups)} readgroups')
    if args.star:
        for s in samples:
            s.command = getStarCommand(s)
    if args.bwa or args.nextflow:
        for s in samples:
            s = getBwaCommand(s)
    # Constructs a dictionary with the following format:
    # "sample": ["sample name", ["rg1"...], "star command"]
    sample_dict = constructDict(samples, args.star, args.bwa, args.nextflow)
    if args.output:
        directory = directory.parent
        tsv_file = os.path.join(directory, args.output)
    else:
        directory = directory.parent
        tsv_file = os.path.join(directory, "samples.tsv")
    print_tsv(sample_dict, tsv_file, args.star, args.bwa, args.nextflow)
    print(f'Formed {tsv_file} of {len(samples)} samples.')


def get_samples(fastq_files: list):
    """Gets incomplete samples from list of fastq files.

    Just extracts the sample name. Does not format the readgroups or ensure that files exists.
    Sample names are assumed to have the form "SAMPLE_..." (i.e. are followed by an underscore).

    Args:
        fastq_files (list): List of fastq files

    Returns:
        list: List of incomplete Sample objects
    """
    # Establishes list of samples
    samples = []
    for file in fastq_files:
        # Assumes that sample name will be followed by _ and then finds the name
        sample_name = re.search('^.+?(?=_)', file).group(0)
        # Initializes a sample with that name
        sample = Sample(sample_name)
        # Makes a regex looking for that sample name
        name_match = re.compile(f"^{sample.name}")
        # Adds that sample to the list of samples if there isn't already one with that name initialized
        samples if any(name_match.match(s.name)
                       for s in samples) else samples.append(sample)
    # sorts the list of samples alphabetically (had to use a lambda function because I used Sample Class)
    samples.sort(key=lambda s: s.name)
    return samples


# Establishes readgroups
def get_readgroups(samples: list, fastq_files: list, fastq_dir: Path, verb: int):
    """Produces a list of samples with readgroups

    Args:
        samples (list): List of incomplete sample objects missing readgroup information
        fastq_files (list): List of files to form readgroups from
        fastq_dir (directory): Directory of fastq files
        verb (int): Amount of debug verbosity

    Raises:
        Exception: Unable to find r1 or r2 file associated with 

    Returns:
        list: List of complete sample objects with defined name and readgroups.
    """
    # Iterates through samples
    final_samples = []
    for sample in samples:
        # each sample must have at least one R1 and R2 file to be valid
        match_r1 = re.compile(f"^{sample.name}.+?(?=R1)")
        match_r2 = re.compile(f"^{sample.name}.+?(?=R2)")
        # Takes the first read, up to but not including the R1, as long as there is a fastq file
        # for both R1 and R2 for the sample
        readgroups = [match_r1.match(file).group(0) for file in fastq_files if
                      match_r1.match(file) and any(match_r2.match(file) for file in fastq_files)]
        # Makes sure that there were some readgroups found for each sample
        if not readgroups and verb > 0:
            print(f"{sample.name} has no readgroups. {readgroups}")
        for rg in readgroups:
            # Ensures that each readgroup has 2 reads files
            reads = glob.glob(f'{os.path.join(fastq_dir, rg)}*')
            if len(reads) != 2:
                if verb > 0:
                    print(f'{reads} does not have a mate')
            # Separates out the reads into r1 and r2
            elif len(reads) == 2:
                # Defines read 1 and read 2
                r1 = glob.glob(f'{os.path.join(fastq_dir, rg)}*R1*')
                r2 = glob.glob(f'{os.path.join(fastq_dir, rg)}*R2*')
                if not r1 or not r2:
                    print(r1, r2)
                    raise Exception(
                        f'Either R1: {r1} or R2: {r2} was not found. Make sure your samples are formatted correctly. \n A good example is one names "SAMPLE_READGROUP_R1.fastq".')
                # Creates readgroup object
                final_rg = Readgroup()
                final_rg.rg = rg
                final_rg.r1 = r1[0]
                final_rg.r2 = r2[0]
                # Adds the readgroup to the sample's list
                sample.addReadgroup(final_rg)
        if len(sample.readgroups) != 0:
            final_samples.append(sample)
    return final_samples


def getStarCommand(sample: Sample):
    """Returns a string with a partial STAR command for that sample

    Args:
        sample (Sample): Sample you're lookign for

    Returns:
        string: Command string to read files into STAR mapper
    """
    if len(sample.readgroups) == 1:
        return f'--readFilesIn {sample.readgroups[0].r1} {sample.readgroups[0].r2} --outSAMattributes All ' \
               f'--outSAMattrRGline ID:{sample.readgroups[0].rg} PL:ILLUMINA SM:{sample.name} '
    elif len(sample.readgroups) > 1:
        r1s = []
        r2s = []
        rgs = []
        for readg in sample.readgroups:
            r1s.append(readg.r1)
            r2s.append(readg.r2)
            rgs.append(readg.rg)
        samattr = f'--outSAMattributes All --outSAMattrRGline ID:{rgs[0]} PL:ILLUMINA SM:{sample.name} '
        for rg in rgs[1:]:
            samattr = samattr + f', ID:{rg} PL:ILLUMINA SM:{sample.name} '
        return f'--readFilesIn {",".join(r1s)} {",".join(r2s)} {samattr}'


def getBwaCommand(sample: Sample):
    """Formats a readgroup attribute flag for bwa

    Args:
        sample (Sample): A sample object

    Returns:
        Sample: modified sample object with readgroups with commands
    """
    for rg in sample.readgroups:
        # need two back slashed on the tabs for snakemake
        command = f' -R \'@RG\\tID:{rg.rg}\\tPL:ILLUMINA\\tSM:{sample.name}\' '
        rg.command = command
    return sample


def constructDict(samples: list, star: bool, bwa: bool, nextflow: bool):
    """Constructs a dictionary for each output type

    Args:
        samples (List[Sample]): List of fully formatted samples
        star (bool): Format for star output. Outputs sample [files] command tsv.
        bwa (bool): Format for bwa output. Outputs sample Readgroup(r1, r2, command) tsv.
        nextflow (bool): Format for nextflow output. Outputs Sample Readgroup r1 r2 bwa_flag_command tsv.

    Returns:
        dict: Dictionary of the described format
    """
    sample_dict = {}
    if nextflow:
        sample_dict["nextflow_table"] = []
    for s in samples:
        files = [rg.r1 for rg in s.readgroups] + [rg.r2 for rg in s.readgroups]
        if star:
            sample_dict[s.name] = [s.name, files, s.command]
        elif bwa:
            groups = {}
            for readgroup in s.readgroups:
                name = readgroup.rg
                file1 = readgroup.r1
                file2 = readgroup.r2
                command = readgroup.command
                groups[name] = [file1, file2, command]
            sample_dict[s.name] = [s.name, groups]
        elif nextflow:
            for rg in s.readgroups:
                sample_dict["nextflow_table"].append(
                    [s.name, rg.rg, rg.r1, rg.r2, rg.command])
        else:
            sample_dict[s.name] = [s.name, files]
    return sample_dict


def print_tsv(sample_d: dict, tsv_file: str, star: bool, bwa: bool, nextflow: bool):
    """Outputs the tsv of the sample dict

    Args:
        sample_d (dict): Dictionary output by @constructDict
        tsv_file (str): Path to desired tsv file
        star (bool): Format for STAR output.
        bwa (bool): Format for bwa output.
        nextflow (bool): Format for nextflow output.
    """
    if star:
        data = pd.DataFrame.from_dict(sample_d, orient='index', columns=[
                                      'sample_name', 'files', 'command'])
    elif bwa:
        data = pd.DataFrame.from_dict(sample_d, orient='index', columns=[
                                      'sample_name', 'files'])
    elif nextflow:
        data = pd.DataFrame(sample_d["nextflow_table"], columns=[
                            'sample_name', 'readgroup', 'r1', 'r2', 'bwa_read_group_string'])
    else:
        data = pd.DataFrame.from_dict(sample_d, orient='index', columns=[
                                      'sample_name', 'files'])
    data.to_csv(tsv_file, sep='\t', index=False)


if __name__ == "__main__":
    # Parses arugments
    # Only required argument is for the sample directory
    parser = argparse.ArgumentParser(
        description="Create sample tsv for RNAseq mapping with STAR")
    parser.add_argument("fastq_dir", type=dir_path,
                        help="Directory containing fastq files")
    parser.add_argument("-o", "--output", type=str, help="Name of output tsv")
    parser.add_argument("-v", "--verbose", action="count",
                        default=0, help="Enable debug output")
    parser.add_argument("-s", "--star", action="store_true",
                        help="Form commands for mapping with STAR.")
    parser.add_argument("-b", "--bwa", action="store_true",
                        help="Form commands for mapping with bwa (mem 2).")
    parser.add_argument("-n", "--nextflow", action="store_true",
                        help="Format for nextflow. Creates table with sample_name, readgroup, r1, r2, bwa_read_group_string as a tsv.")
    parser.add_argument("-i", "--input-samples", type=str, help="Provide a tsv of sample names to look for in the "
                                                                "fastq directory.")
    args = parser.parse_args()
    main(args)
