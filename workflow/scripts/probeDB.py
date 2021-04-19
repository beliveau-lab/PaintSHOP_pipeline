#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 16:11:26 2019

@author: hershe

A script to generate a database of refseq annotations with a row for each
probe that overlaps the region.

"""

import os
import timeit
import pybedtools
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd

# configure file paths
PROBE_TSVS = snakemake.input.probes
ANNOTATIONS = snakemake.input.annotations
OUTPUT_PATH = snakemake.output.tsv

def main():

    start_time = timeit.default_timer()

    # creates a list of the all of the probes from each
    # probe file, each represented as a single string
    probes = []
    for f in PROBE_TSVS:
        with open(f) as probe_file:
            for line in probe_file:
                probes.append(line.strip())

    # creates a bedtool object with the entire probe set for the assembly
    probe_bedtool = pybedtools.BedTool(probes)

    # creates a bedtool object for the refseq annotations
    # (snakemake input needs to be refseq BED file)
    refseq_bedtool = pybedtools.BedTool(ANNOTATIONS)

    # create output name for raw intersect result
    intersect_result = OUTPUT_PATH.split('.')[0] + '_raw_intersect.bed'

    # wo is the same as -wo from BedTools, and writes the original
    # A and B entries plus the number of base pairs overlap of the
    # two features. f is the same as -f in BedTools, which is the
    # minimum overlap percentage accepted. 1 means 100% overlap
    probe_result = probe_bedtool.intersect(refseq_bedtool, wo=True, f=1,
                                           output = intersect_result)
    # remove since I parse output file
    del probe_result

    # read in the bed file of the results of intersection
    df = pd.read_csv(intersect_result, sep='\t', header=None)
    df = df.rename(columns={
        0: 'probe_chrom',
        1: 'probe_start',
        2: 'probe_stop',
        3: 'probe_seq',
        4: 'probe_tm',
        5: 'probe_on_target',
        6: 'probe_off_target',
        7: 'probe_repeat',
        8: 'probe_max_kmer',
        9: 'probe_prob',
        10: 'probe_strand',
        11: 'ref_chrom',
        12: 'ref_start',
        13: 'ref_stop',
        14: 'ref_refseq',
        15: 'ref_score',
        16: 'ref_strand',
    })
    df.columns = [*df.columns[:-1], 'overlap']

    # reverse complement probe sequences as needed based on reference strand
    df['probe_seq'] = df.apply(check_polarity, axis=1)
    df['probe_strand'] = df.apply(get_strand, axis=1)

    # drop unnecessary columns
    drop_cols = [
        'ref_chrom',
        'ref_start',
        'ref_stop',
        'ref_score',
        'ref_strand',
        'overlap',
    ]
    df = df.drop(drop_cols, axis=1).sort_values(['probe_chrom', 'probe_start'])

    # save data frame to disk
    df.to_csv(OUTPUT_PATH,
                  sep = '\t',
                  index = False,
                  index_label = False,
                  header = False,
                  float_format='%.3f')

    # remove temp file
    os.remove(intersect_result)

    # success
    print('Program took %f seconds' % (timeit.default_timer() - start_time))


def check_polarity(row):
    """"Checks polarity and flips probe sequence if reference is +."""
    if row['ref_strand'] == '+':
        return str(Seq(row['probe_seq'], IUPAC.unambiguous_dna).reverse_complement())
    else:
        return str(row['probe_seq'])


def get_strand(row):
    """"Checks polarity and returns strand as '+' or '-'."""
    if row['ref_strand'] == '+':
        return '-'
    else:
        return '+'


if __name__ == "__main__":
    main()
