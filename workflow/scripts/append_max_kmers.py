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

# load probes with metadata
df = pd.read_csv(snakemake.input[0], sep='\t', header=None)
df.columns = [
    'chrom',
    'start',
    'stop',
    'parent',
    'Tm',
    'on_target_score',
    'off_target_score',
    'repeat',
    'prob'
]

# read in max kmer counts
with open(snakemake.input[1], 'r') as file:
    max_kmers = [line.strip() for line in file]

# append counts
df['max_kmers'] = max_kmers

# add hard-coded + strand (blockparse output)
df['strand'] = '+'

# sort probes on chromosomal start coordinate
df = df.sort_values('start')

# save to disk
df.to_csv(snakemake.output[0], sep='\t', index=False, header=False)
