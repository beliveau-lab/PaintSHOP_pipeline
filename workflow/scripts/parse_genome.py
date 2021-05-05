
from Bio import SeqIO

# configure file paths
INPUT_FA = snakemake.input.fasta
OUTPUT_FA = snakemake.output.multi_fasta
CHROM_NAMES = snakemake.output.chrom_names

# chromosome records to be filtered out
REMOVED_UNK = snakemake.output.removed_unk
REMOVED_EXCL = snakemake.output.removed_excl
REMOVED_RAND = snakemake.output.removed_rand
REMOVED_ALT = snakemake.output.removed_alt
REMOVED_HAP = snakemake.output.removed_hap
REMOVED_FIX = snakemake.output.removed_fix

# optionally manually filter out chromosome records
EXCL_CHROMS = []
if 'exclude_chroms' in snakemake.config:
    EXCL_CHROMS = snakemake.config['exclude_chroms']

def main():

    # load raw genomic fasta
    print('~~~~~~~~~~ Processing chromosomes ~~~~~~~~~~')
    recs = SeqIO.parse(snakemake.input[0], 'fasta')

    # store names of excluded chromosomes
    removed_unk = []
    removed_excl = []
    removed_rand = []
    removed_alt = []
    removed_hap = []
    removed_fix = []
    chrom_names = []

    # filter out alt seq records
    filtered = []
    num_recs = 0
    for seq in recs:

        # track total sequence records
        num_recs += 1
        
        # truncate fasta description to seqid only
        seq.description = seq.id

        # no associated chromosome, include in bowtie2 index
        if 'Un_' in seq.id:
            removed_unk.append(seq.id)  # keep a record of this chromosome
            filtered.append(seq)        # add to multifasta for bowtie2 index

        # user-excluded chrom, include in bowtie2 index but skip probe design
        elif seq.id in EXCL_CHROMS:
            removed_excl.append(seq.id) # keep a record of this chromosome
            filtered.append(seq)        # add to multifasta for bowtie2 index

        # un-placed on chromosome, include in bowtie2 index
        elif '_random' in seq.id:
            removed_rand.append(seq.id) # keep a record of this chromosome
            filtered.append(seq)        # add to multifasta for bowtie2 index

        # novel sequence, exclude entirely
        elif '_alt' in seq.id:
            removed_alt.append(seq.id)  # keep a record of this chromosome

        # alternate haplotype, exclude entirely
        elif '_hap' in seq.id:
            removed_hap.append(seq.id)  # keep a record of this chromosome

        # fix sequence, exclude entirely
        elif '_fix' in seq.id:
            removed_fix.append(seq.id)  # keep a record of this chromosome

        # presumed to be a canonical chromosome, include in all steps
        else:
            print(seq.id)
            chrom_names.append(seq.id)  # keep a record of this chromosome
            filtered.append(seq)        # add to multifasta for bowtie2 index

    # save filtered multi-fasta genome
    print('Writing filtered seq records to new multi-fasta file...')
    with open(OUTPUT_FA, 'w') as outfile:
        SeqIO.write(filtered, outfile, 'fasta')
        num_filtered = len(filtered)
        print(f'Saving {num_filtered} of {num_recs} total chromosomes.')

    # write list of presumptive canonical chromosomes
    print('Writing chromosome name files...')
    with open(CHROM_NAMES, 'w') as outfile:
        print(f'found {len(chrom_names)} presumptive canonical chromosomes')
        outfile.write('\n'.join(chrom_names))

    # write list of filtered chromosomes
    with open(REMOVED_UNK, 'w') as outfile:
        print(f'found {len(removed_unk)} Un_ chromosomes')
        outfile.write('\n'.join(removed_unk))
    with open(REMOVED_RAND, 'w') as outfile:
        print(f'found {len(removed_rand)} _random chromosomes')
        outfile.write('\n'.join(removed_rand))
    with open(REMOVED_EXCL, 'w') as outfile:
        print(f'found {len(removed_excl)} user-excluded chromosomes')
        outfile.write('\n'.join(removed_excl))
    with open(REMOVED_ALT, 'w') as outfile:
        print(f'found {len(removed_alt)} _alt chromosomes')
        outfile.write('\n'.join(removed_alt))
    with open(REMOVED_HAP, 'w') as outfile:
        print(f'found {len(removed_hap)} _hap chromosomes')
        outfile.write('\n'.join(removed_hap))
    with open(REMOVED_FIX, 'w') as outfile:
        print(f'found {len(removed_fix)} _fix chromosomes')
        outfile.write('\n'.join(removed_fix))

    # success
    print('DONE!')

if __name__ == '__main__':

    import time
    start = time.time()
    main()
    print('Script took {} seconds'.format(time.time() - start))
