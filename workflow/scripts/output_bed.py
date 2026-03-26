#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 10:08:45 2019

@author: hershe
"""

import pandas as pd
from Bio.SeqUtils import MeltingTemp as mt

# function to calculate Tm for each probe
def probeTm(row):
    """Calculates the Tm of a given sequence."""
    tmval = float(('%0.2f' \
                   % mt.Tm_NN(row['parent'], Na=390, dnac1=25, dnac2=25)))
    fcorrected = ('%0.2f' % mt.chem_correction(tmval, fmd=50))
    return fcorrected

# read in data
df = pd.read_csv(snakemake.input[0])

# compute on target score as the duplex probability for the on target probe
on_target_scores = df[df['on_target'] == 1] \
                        .groupby(['probe_ID', 'parent'])['duplex_pred'] \
                        .sum().reset_index(name='on_target_score')

# compute off target score as the sum of all duplex probabilities 
# of all off target regions
off_target_scores = df[df['on_target'] == 0] \
                        .groupby('probe_ID')['duplex_pred'] \
                        .sum().reset_index(name='off_target_score')

# create a data set with both scores by a join
scores_df = on_target_scores.set_index('probe_ID') \
                        .join(off_target_scores.set_index('probe_ID')) \
                        .fillna(0)

scores_df = scores_df.reset_index()

# parse info for BED from probe_ID (format: "chrom:start-stop|repeat")
id_and_repeat = scores_df['probe_ID'].str.split('|', n=1)
coords = id_and_repeat.str[0]  # "chrom:start-stop"
scores_df['repeat'] = id_and_repeat.str[1]
scores_df['chrom'] = coords.str.split(':').str[0]
start_stop = coords.str.split(':').str[1].str.split('-', n=1)
scores_df['start'] = start_stop.str[0]
scores_df['stop'] = start_stop.str[1]

# remove unnecesary column
scores_df = scores_df.drop('probe_ID', axis = 1)

# calculate melting temps
scores_df['Tm'] = scores_df.apply(probeTm, axis = 1)

# rearrange data to BED format
scores_df = scores_df[['chrom',
                      'start',
                      'stop',
                      'parent',
                      'Tm',
                      'on_target_score',
                      'off_target_score',
                      'repeat']]

# save to disk
scores_df.to_csv(snakemake.output[0],
                 sep = '\t',
                 index = False,
                 index_label = False,
                 header = False)