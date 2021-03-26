#!/usr/bin/env python3
# This script is responsible for parsing a directory of
# fastq samples to produce a yaml that is easy to manage
# with snakemake.
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
    def __init__(self, name=""):
        self.name = name
        self.readgroups = []
        self.command = ""

    def __repr__(self):
        rep = f'Sample(name: {self.name}, readgroups: {self.readgroups}, command: {self.command})'
        return rep

    def addReadgroup(self, rg):
        self.readgroups.append(rg) if rg not in self.readgroups else self.readgroups


class Readgroup:
    # name of readgroup
    def __init__(self, rg=""):
        self.rg = rg
        self._r1 = ""
        self._r2 = ""

    # Equality
    def _is_valid_operand(self, other):
        return hasattr(other, "rg")

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.rg == other.rg

    # better printing
    def __repr__(self):
        rep = f'Readgroup(rg: {self.rg}, r1: {self.r1}, r2: {self.r2})'
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
def dir_path(string):
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
            s.command = getCommand(s)
    # Constructs a dictionary with the following format:
    # "sample": ["sample name", ["rg1"...], "star command"]
    sample_dict = constructDict(samples, args.star)
    if args.output:
        tsv_file = args.output
    else:
        directory = directory.parent
        tsv_file = os.path.join(directory, "samples.tsv")
    print_tsv(sample_dict, tsv_file, args.star)
    print(f'Formed {tsv_file}/samples.tsv of {len(samples)} samples.')


def get_samples(fastq_files):
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
        samples if any(name_match.match(s.name) for s in samples) else samples.append(sample)
    # sorts the list of samples alphabetically (had to use a lambda function because I used Sample Class)
    samples.sort(key=lambda s: s.name)
    return samples


# Establishes readgroups
def get_readgroups(samples, fastq_files, fastq_dir, verb):
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
                    raise Exception(f'Either R1: {r1} or R2: {r2} was not caught. File bug report')
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


def getCommand(sample):
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


def constructDict(samples, star):
    sample_dict = {}
    for s in samples:
        rgs = [rg.r1 for rg in s.readgroups] + [rg.r2 for rg in s.readgroups]
        if star:
            sample_dict[s.name] = [s.name, rgs, s.command]
        else:
            sample_dict[s.name] = [s.name, rgs]
    return sample_dict


def print_tsv(sample_d, tsv_file, star):
    if star:
        data = pd.DataFrame.from_dict(sample_d, orient='index', columns=['sample_name', 'files', 'command'])
    else:
        data = pd.DataFrame.from_dict(sample_d, orient='index', columns=['sample_name', 'files'])
    data.to_csv(tsv_file, sep='\t', index=False)



if __name__ == "__main__":
    # Parses arugments
    # Only argument is for the directory
    parser = argparse.ArgumentParser(description="Create sample tsv for RNAseq mapping with STAR")
    parser.add_argument("fastq_dir", type=dir_path, help="Directory containing fastq files")
    parser.add_argument("-o", "--output", type=str, help="Name of output tsv")
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Enable debug output")
    parser.add_argument("-s", "--star", action="store_true", help="Form commands for mapping with STAR.")
    parser.add_argument("-i", "--input-samples", type=str, help="Provide a tsv of sample names to look for in the "
                                                                "fastq directory.")
    args = parser.parse_args()
    main(args)
    exit(0)
