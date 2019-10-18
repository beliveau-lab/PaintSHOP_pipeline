<div align="center">
  <img src="images/PaintSHOP-logo.png" width="200">
</div>

# PaintSHOP Overview [![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.4-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

PaintSHOP is a technology that enables the interactive design of modern oligonucleotide FISH experiments at genome and transcriptome-scale.

PaintSHOP is comprised of two components:
1. A scalable machine learning pipeline for probe specifity prediction
2. An interactive Shiny web application for probe design

This repository contains the code for the machine learning pipeline.

We provide this open source software without any warranty under the [MIT license](https://opensource.org/licenses/MIT).

A manuscript describing this work is in preparation.

# Description of Pipeline

The pipeline is built to take a set of candidate oligonucleotide (oligo) fluorescent *in-situ hybridization* (FISH) probes, and generate quantitative predictions about the likelihood of the candidates hybridizing with sequences other than their intended target in the genome. Off-target hybridization results in background noise, and makes FISH signal harder to interpret.

Here is a schematic overview of the pipeline:

<div align="center">
  <img src="images/HOP-schematic.png">
</div>

<div><br></div>

Detailed descriptions of the pipeline are provided in the PaintSHOP publication. Briefly, the pipeline works by aligning all probe candidates against the genome they are designed for using the [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) NGS aligner with very sensitive parameters. Pairwise alignments are reconstructed using [sam2pairwise](https://github.com/mlafave/sam2pairwise). A gradient boosting regression model implemented with [XGBoost](https://xgboost.readthedocs.io/en/latest/#) is used to predict the probability of the candidate hybridizing at each alignment site based on a thermodynamic partition function. The scores are aggregrated into an on-target and off-target score for each probe in the set.

# Installation & Use

In order to use the pipeline, you first need to download [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Next, you can run:

```
$ git clone https://github.com/elliothershberg/PaintSHOP_pipeline.git && cd PaintSHOP_pipeline/
$ conda env create -f envs/environment.yaml
$ conda activate PaintSHOP_pipeline
```

to clone this repository and create a conda environment with the dependencies required.

The last installation step is to install [sam2pairwise](https://github.com/mlafave/sam2pairwise) following the instructions provided.

You will also need to either create or download a Bowtie2 index and a [Jellyfish](http://www.genome.umd.edu/jellyfish.html) file. Information for downloading pre-computed Bowtie2 indeces or building new indeces can be found [here](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). The [OligoMiner](https://github.com/brianbeliveau/OligoMiner) repository also provides detailed information on creating Bowtie2 indeces and Jellyfish files.

This pipeline is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html), and distributed according to [best practices](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html). If you are familiar with Snakemake or similar workflow tools, you can edit the config.yaml file to include your information, and run the pipeline. A quick example of the workflow can be found [here](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html). If you are new to Snakemake, the [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) is a great place to get started to learn more.

The information needed by the config file are the 4 following things:
1. The path to the directory with your probe file(s) in fastq format.
2. The path to your Bowtie2 index.
3. The temperature at which which you want to predict probe hybridization (37, 42, 47, 52, 60)
4. The path to your Jellyfish file

# Pipeline Output

The pipeline will return a number of output folders. The most important are:

| Folder        | Description                                                       |
|---------------|-------------------------------------------------------------------|
| alignments/   | BAM files with Bowtie2 alignments                                 |
| predictions/  | data frames with the ML predictions used to generate final scores |
| final_probes/ | the final probe set with all k-mer counts and scores              |

The columns of the files in the ```final_probes/``` folder are:

| 0          | 1     | 2    | 3        | 4   | 5                       | 6                           | 7               | 8         |
|------------|-------|------|----------|-----|------------------------ |-----------------------------|-----------------|-----------| 
| chromosome | start | stop | sequence | Tm  | on-target score (0-100) | off-target score (0-10,000) | repeat (0 or 1) | max k-mer | 

The on-target score is the hybridization prediction at the target site of the probe.

The off-target score is the sum of the hybridization predictions for all off-target sites (up to 100) for the probe.

The repeat column is 0 if the sequence contains no repeat-masked bases, and 1 if it does. 

The max k-mer column is the max count out of all k-mers in the probe sequence.

# Utilities

Several addition utilities are provided for using the pipeline and the results. The pipeline is written to take fastq files with the following information in each description line:

```@chr:start-stop|0``` for a probe containing NO repeat-masked bases in its sequences, or ```@chr:start-stop|1``` for a probe with repeat-masked bases.

Since most probe sets are designed against hard-masked genomes (with no soft-masked bases) and are often stored in BED format, we provided the ```bedToFastq_repeat.py``` script to generate a fastq file with all |0 flags from a BED file of probes for use in the pipeline. This is a modified version of a utility from the [OligoMiner](https://github.com/brianbeliveau/OligoMiner) suite.

We also have provided a modified version of the ```blockParse.py``` script from OligoMiner called ```blockParse_unmasked.py``` that can mine candidate probes from soft-masked genome assemblies, and store repeat information in the format that is downstream compatible with this pipeline. We have also provided a separate conda environment in the envs/ folder that can be used to run this script to mine for probes. Thanks to [Conor Camplisson](https://github.com/conorcamplisson/oligo_miner_env) for putting that together.

Another useful utility is the script used to generate the backend database files for PaintSHOP based on intersections with UCSC RefSeq annotation sets. This database can be used to quickly retrieve probes that target any RefSeq annotation.

An example use is:

```
python probeDB.py -p final_probes/ -a hg38_refseq.bed -o hg38_refseq_db.tsv
```
# Questions

If you have questions or issues, feel free to either open an issue, or [contact me](mailto:eahershberg@gmail.com) directly.
