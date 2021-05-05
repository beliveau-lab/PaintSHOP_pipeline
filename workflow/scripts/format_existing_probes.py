
import os
import pandas as pd

# configure file paths
PROBE_FILE = snakemake.config['existing_probes']
FASTQ_PATH = snakemake.output[0]
CHROM = os.path.basename(FASTQ_PATH).replace('.fastq', '')

def main():

    # load existing probes
    df = pd.read_csv(PROBE_FILE, sep='\t', header=None)

    # truncate df to first four columns and name columns
    df = df[df.columns[0:4]]
    df.columns = ['chrom', 'start', 'stop', 'seq']

    # filter probes to the current chromosome only
    df = df[df['chrom'] == CHROM].reset_index(drop=True)

    # create fastq formatted probe file with the current chromosome's probes
    with open(FASTQ_PATH, 'w') as outfile:

        # write a fastq entry for each probe
        for i in range(len(df)):

            # extract this probe's data
            start = df['start'][i]
            stop = df['stop'][i]
            seq = df['seq'][i]

            # write probe to fastq file
            fastq_str = f'@{CHROM}:{start}-{stop}|0\n' # assume no-repeat
            fastq_str += f'{seq}\n'
            fastq_str += '+\n'
            fastq_str += '~' * len(seq) + '\n'
            outfile.write(fastq_str)


if __name__ == '__main__':
    main()
