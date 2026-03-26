"""Mine candidate probes from a chromosome FASTA using OligoMiner2."""

from oligominer import mine_fasta

probes = mine_fasta(
    snakemake.input[0],
    min_length=snakemake.params.min_length,
    max_length=snakemake.params.max_length,
    min_tm=snakemake.params.min_tm,
    max_tm=snakemake.params.max_tm,
    overlap=False,
)

# Write fastq with 1-based coords and "|0" repeat suffix to match downstream
# pipeline expectations (bowtie2/SAM uses 1-based positions for on-target matching,
# blockParse wrote "|0" for unmasked, "|1" for repeat-masked probes)
with open(snakemake.output[0], 'w') as f:
    for seq_id, start, stop, probe_seq, tm in probes:
        f.write(f"@{seq_id}:{start + 1}-{stop}|0\n{probe_seq}\n+\n{'~' * len(probe_seq)}\n")
