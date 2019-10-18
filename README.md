<div align="center">
  <img src="images/PaintSHOP-logo.png" width="200">
</div>

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.4-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

# PaintSHOP Overview

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

# Pipeline Output









