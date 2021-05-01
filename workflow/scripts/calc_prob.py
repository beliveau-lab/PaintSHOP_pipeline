
import numpy as np
import pandas as pd
import nupack

# configure nupack model
NUPACK_MODEL = nupack.Model(
    material=snakemake.config.get('nupack_material', 'dna'),
    ensemble=snakemake.config.get('nupack_ensemble', 'stacking'),
    celsius=snakemake.config.get('nupack_celsius', 74.5),
    sodium=snakemake.config.get('nupack_sodium', 0.39),
    magnesium=snakemake.config.get('nupack_magnesium', 0.0),
)

def main():

    # load probe data
    df = pd.read_csv(snakemake.input[0], sep='\t', header=None)
    df.columns = [
        'chrom',
        'start',
        'stop',
        'seq',
        'Tm',
        'on_target_score',
        'off_target_score',
        'repeat'
    ]

    # compute prob value for each sequence
    df['prob'] = df['seq'].apply(prob)

    # write outfile with prob value added
    df.to_csv(snakemake.output[0], sep='\t', header=None, index=False)


def prob(seq):
    '''Calculates the probability that the input sequence has no secondary structure.'''
    
    # invoke nupack's 'prob' command using input seq and a linear structure
    probability = nupack.structure_probability(strands=[seq], structure='.' * len(seq), model=NUPACK_MODEL)
    
    # success
    return(probability)


if __name__ == '__main__':
    main()
