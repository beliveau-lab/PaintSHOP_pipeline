
import pickle

import numpy as np
import pandas as pd

# configure file paths
PICKLED_DF = snakemake.input.df
OUTPUT_BED = snakemake.output.bed

def main():

    # load the input annotation file
    df = pickle.load(open(PICKLED_DF, 'rb'))

    # flatten isoforms to shared segments only
    flat_df = flatten_isoforms(df)

    # save isoform-flattened bed file
    flat_df.to_csv(OUTPUT_BED, sep='\t', index=False, header=None)


def flatten_isoforms(df):

    # group exons by gene
    grouped_df = df.groupby('gene_id').aggregate({
        'seqid': 'first',
        'start': lambda x: tuple(x),
        'end': lambda x: tuple(x),
        'score': 'first',
        'strand': lambda x: tuple(x),
        'transcript_version': lambda x: tuple(x),
    }).reset_index(drop=False)

    # perform isoform flattening
    exon_lists = grouped_df.apply(get_flattened_exons, axis=1).reset_index(drop=True)
    merged_exon_list = [exon for exon_list in exon_lists.values for exon in exon_list]
    flat_df = pd.DataFrame(merged_exon_list)

    # success
    return(flat_df)


def get_flattened_exons(row):
    
    # convert input coordinates to numpy arrays
    starts = np.array(row.start, dtype=int)
    ends = np.array(row.end, dtype=int)
    
    # calculate the total span of all exons
    num_exons = starts.size
    min_start = np.min(starts)
    max_end = np.max(ends)

    # zero shift coordinates
    starts_norm = starts - min_start
    ends_norm = ends - min_start
    max_end_norm = np.max(ends_norm)
    
    # calculate exon coverage
    coverage = np.zeros(max_end - min_start, dtype=int)
    for i in range(num_exons):
        coverage[starts_norm[i]:ends_norm[i]] += 1
    max_coverage = np.max(coverage)

    # collapse isoforms as needed
    output_starts = starts
    output_ends = ends
    if max_coverage > 1:
    
        # find start and end coordinates of maximally shared segments
        mask = (coverage == max_coverage).astype(int)
        shared_starts = np.where(np.append([mask[0]], np.diff(mask)) == 1)[0]
        shared_ends = np.where(np.diff(mask) == -1)[0] + 1
        if shared_ends.size < shared_starts.size:
            shared_ends = np.append(shared_ends, max_end_norm)

        # update output coordinates if shared segments were found
        if shared_starts.size > 0:
            
            # undo previous zero shift
            output_starts = shared_starts + min_start
            output_ends = shared_ends + min_start

    # format output annotation data for shared segments
    segments = []
    for i in range(output_starts.size):
        segment_data = (
            row.seqid,
            output_starts[i],
            output_ends[i],
            row.gene_id,
            row.score,
            row.strand[i],
            max_coverage,
        )
        segments.append(segment_data)

    # success        
    return(segments)


if __name__ == '__main__':
    main()
