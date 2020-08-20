
import os

from Bio import SeqIO

# configure file paths
INPUT_FA = snakemake.input.fasta
CHROM_NAMES = snakemake.input.chrom_names
CHROM_FASTA_DIR = snakemake.output.chrom_fasta_dir

# create output directory
os.makedirs(CHROM_FASTA_DIR, exist_ok=True)
	
def main():

	# load raw genomic fasta
	print('~~~~~~~~~~ Splitting genome into chromosome fastas ~~~~~~~~~~')
	recs = SeqIO.parse(INPUT_FA, 'fasta')
	

	# load filtered list of chromosome names
	chrom_names = get_chrom_names()

	# save individual chromosome files
	for seq in recs:

		# skip non-canonical chromosomes
		if seq.id not in chrom_names:
			print('Skipping non-canonical chromosome: ', seq.id)
			continue

		# save this chromosome
		file_path = os.path.join(CHROM_FASTA_DIR, f'{seq.id}.fa')
		with open(file_path, 'w') as outfile:
			print(f'Writing {seq.id}')
			SeqIO.write(seq, outfile, 'fasta')



def get_chrom_names():

    # load chromosome names from filtered fasta
    with open(CHROM_NAMES, 'r') as infile:
        text = infile.read().strip(' \n')
    chrom_names = text.split('\n')

    # success
    return(chrom_names)


if __name__ == '__main__':

	import time
	start = time.time()
	main()
	print('Script took {} seconds'.format(time.time() - start))
