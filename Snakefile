# load config file
configfile: "config.yaml"

#glob all mm9 DNA files
SAMPLES, = glob_wildcards(config["fastq_directory"] + "{sample}.fastq")

# a pseudo-rule that collects the target files
rule all:
    input:
        expand("final_probes/{sample}_probes.bed", sample=SAMPLES)

rule align:
    input:
        config["fastq_directory"] + "{sample}.fastq"
    output:
        "alignments/{sample}.bam"
    shell:
        "bowtie2 -x {config[bowtie_index]} -U {input} --very-sensitive-local -k 100 | "
        "samtools view -bS - > {output}"

rule sam_to_pairwise:
    input:
        "alignments/{sample}.bam"
    output:
        "pairwise_alignments/{sample}_pairwise.out"
    shell:
        "samtools view {input} | sam2pairwise > {output}"

rule extract_alignment_scores:
    input:
        "alignments/{sample}.bam"
    output:
        "alignment_scores/{sample}_AS.txt"
    shell:
        "samtools view {input} | awk '{{print $12}}' > {output}"

rule parse_pairwise:
    input:
        "pairwise_alignments/{sample}_pairwise.out",
        "alignment_scores/{sample}_AS.txt"
    output:
        "data_frames/{sample}_alignments.csv"
    script:
        "scripts/parse_pairwise.py"

rule predict_duplex:
    input:
        "data_frames/{sample}_alignments.csv",
        expand("pickled_models/{param}_all_fixed_xgb.pickle.dat", param=config["model_temp"])
    output:
        "predictions/{sample}_predictions.csv"
    script:
        "scripts/XGB_predict.py"

rule score:
    input:
        "predictions/{sample}_predictions.csv"
    output:
        temp("probes/{sample}_probes.bed")
    script:
        "scripts/output_bed.py"

rule max_kmer:
    input:
        "probes/{sample}_probes.bed",
        config["jf_file"]
    output:
        "max_kmer_counts/{sample}_kmer_max.txt"
    script:
        "scripts/kmer_frequency.py"

rule final_probes:
    input:
        "probes/{sample}_probes.bed",
        "max_kmer_counts/{sample}_kmer_max.txt"
    output:
        "final_probes/{sample}_probes.bed"
    script:
        "scripts/append_max_kmers.py"

