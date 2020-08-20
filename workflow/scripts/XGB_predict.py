#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:41:16 2019

@author: hershe

Script to generate predictions of the duplexing probability of probe
candidates and their Bowtie2 alignments.

Does preprocessing necessary to generate features for the model, loads
the pickled model, and predicts the probability for each row in the
data frame.

Written to work with the results of parse_pairwise.py in a snakemake pipeline.

"""

import pandas as pd
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pickle

# read in data
# in snakemake pipeline first argument should be snakemake.input[0]
df = pd.read_csv(snakemake.input[0])

# get reverse complement of derived to mirror format of training data
def reverse_comp(seq):
    forward_seq = Seq(seq, generic_dna)
    return str(forward_seq.reverse_complement())

df['derived'] = df['derived'].apply(reverse_comp)

# calculate GC content of both sequences
df['parent_gc'] = df['parent'].apply(GC)
df['derived_gc'] = df['derived'].apply(GC)

# length of each seq
df['parent_len'] = df['parent'].apply(len)
df['derived_len'] = df['derived'].apply(len)

# count dinucleotides for each sequence for each row
dinucleotides = ['AA','AT','AG','AC',
                 'TA','TT','TG','TC',
                 'GA','GT','GG','GC',
                 'CA','CT','CG','CC']

for pair in dinucleotides:
  parent_name = 'parent_%s' % (pair)
  derived_name = 'derived_%s' % (pair)
  
  df[parent_name] = df['parent'].str.count(pair)
  df[derived_name] = df['derived'].str.count(pair)

# rearrange in to the order for model
features = df
features = features.drop(['parent', 'derived', 'probe_ID',
                          'align_chr', 'align_start'],
                         axis=1)

features = features[['bowtie', 'parent_gc', 'derived_gc',
               'parent_len', 'derived_len', 'parent_AA', 'parent_AT', 'parent_AG',
               'parent_AC', 'parent_TA', 'parent_TT', 'parent_TG', 'parent_TC',
               'parent_GA', 'parent_GT', 'parent_GG', 'parent_GC', 'parent_CA',
               'parent_CT', 'parent_CG', 'parent_CC', 'derived_AA', 'derived_AT',
               'derived_AG', 'derived_AC', 'derived_TA', 'derived_TT', 'derived_TG',
               'derived_TC', 'derived_GA', 'derived_GT', 'derived_GG', 'derived_GC',
               'derived_CA', 'derived_CT', 'derived_CG', 'derived_CC']]

# load in pickled model
model = pickle.load(open(snakemake.input[1], "rb"))

result = model.predict(features.values)

# function to set neg values to 0, and values greater than 100 to 100
def correct_range(val):
    if val < 0:
        return 0
    if val > 100:
        return 100;
    else:
        return val

for i in range(0, len(result), 1):
    result[i] = correct_range(result[i])

df['duplex_pred'] = result

# save data frame to disk
df.to_csv(snakemake.output[0], index_label=False)

