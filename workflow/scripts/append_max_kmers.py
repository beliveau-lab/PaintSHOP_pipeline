#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:01:32 2019

@author: hershe

A script to append the max kmer count for a given probe to the BED row of the
probe. Written for use in my OligoServer/PaintSHOP Snakemake prediction 
pipeline.

"""

import pandas as pd

# read in the probes ('probes/hg38_chrM_DNA-FISH_probes.bed')
probes = pd.read_csv(snakemake.input[0],
                     sep = '\t',
                     header = None)

probes.columns = ['chrom',
                  'start',
                  'stop',
                  'parent',
                  'Tm',
                  'on_target_score',
                  'off_target_score',
                  'repeat']

# read in max kmer counts ('max_kmer_counts/hg38_chrM_DNA-FISH_kmer_max.txt')
with open(snakemake.input[1]) as file:
    max_kmers = [line.strip() for line in file]

# append counts
probes['max_kmers'] = max_kmers

# add hard-coded + strand (blockparse output)
probes['strand'] = '+'

# sort probes on chromosomal start coordinate
probes = probes.sort_values('start')

# save to disk
probes.to_csv(snakemake.output[0],
              sep = '\t',
              index = False,
              index_label = False,
              header = False)
