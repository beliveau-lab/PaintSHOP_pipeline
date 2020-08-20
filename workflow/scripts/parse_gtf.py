
import os
import pickle

import numpy as np
import pandas as pd

# configure file paths
CHROM_NAMES = snakemake.input.chrom_names
GTF_PATH = snakemake.input.gtf
BED_PATH = snakemake.output.bed
CHROM_DF_DIR = snakemake.output.chrom_df_dir

# ensure target directory exists
os.makedirs(CHROM_DF_DIR, exist_ok=True)

def main():

    # load the input annotation file
    df = load_gtf(GTF_PATH)

    # filter annotations to only exons on canonical chromosomes
    df = parse_and_filter(df)

    # split gtf dataframe into individual chromosome files
    chroms = df['seqid'].unique()
    for chrom in chroms:

        # filter gtf to only the current chromosome
        chrom_df = df[df['seqid'] == chrom]

        # save filtered gtf data to disk for this chromosome
        df_path = os.path.join(CHROM_DF_DIR, f'{chrom}_filtered_gtf.dat')
        pickle.dump(chrom_df, open(df_path, 'wb'))

    # save exon-resolved bed file
    print('Saving isoform-resolved BED annotation file...')
    bed_df = df[[
        'seqid', 
        'start', 
        'end', 
        'transcript_id', 
        'score', 
        'strand',
        'transcript_version',
        'gene_id',
    ]]
    bed_df.to_csv(BED_PATH, sep='\t', index=False, header=None)


def load_gtf(file_path):

    # read input gtf file and rename columns
    df = pd.read_csv(GTF_PATH, sep='\t', header=None)
    df.rename(columns={
        0: 'seqid',
        1: 'source',
        2: 'type',
        3: 'start',
        4: 'end',
        5: 'score',
        6: 'strand',
        7: 'phase',
        8: 'attributes',
    }, inplace=True)

    # success
    return(df)


def parse_and_filter(df):

    # load canonical chromosome names
    chrom_names = get_chrom_names()

    # filter annotations to canonical chromosomes only
    df = df[df['seqid'].isin(chrom_names)]

    # filter annotations to exon records only
    df = df[df['type'] == 'exon']

    # parse attributes string into a dict
    attr_data = df['attributes'].apply(parse_attributes)

    # extract attributes as new columns
    df['gene_id'] = attr_data.apply(lambda x: x['gene_id'])
    df['transcript_version'] = attr_data.apply(lambda x: x['transcript_id'])

    df['transcript_id'] = df['transcript_version'].apply(lambda x: x.split('.')[0])

    # df['exon_id'] = attr_data.apply(lambda x: x['exon_id'])

    # remove unnecessary column
    df.drop(['attributes'], axis=1, inplace=True)

    # success
    return(df)


def get_chrom_names():

    # load chromosome names from filtered fasta
    with open(CHROM_NAMES, 'r') as infile:
        text = infile.read().strip(' \n')
    chrom_names = text.split('\n')

    # success
    return(chrom_names)


def parse_attributes(data):

    # parse attributes
    split_data = data.replace('"', '').strip(' ;').split(';')
    attr_data = {}
    for field in split_data:
        split_field = field.strip(' ').split(' ')
        key = split_field[0]
        val = ' '.join(split_field[1:])
        attr_data[key] = val

    # success        
    return(attr_data)


if __name__ == '__main__':
    main()
