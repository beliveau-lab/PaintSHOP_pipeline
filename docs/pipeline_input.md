# Pipeline Input

PaintSHOP pipeline input specification

### Overview

The PaintSHOP pipeline takes two input files: 

1. A fasta file containing the genome sequence

2. A GTF/GFF file containing annotations

The paths to these files are specified in the [config.yml](../example_run/config.yml). Sample files are included in the provided [example genome assembly](../example_run/data/). 

For example, to run this pipeline on the hg38 assembly, the [hg38.fa.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) and [genes/hg38.refGene.gtf.gz](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz) files from [this page](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) were downloaded and de-compressed. 

### Config file

The [config.yml](../example_run/config.yml) file sets all the pipeline parameters for a particular run. A new config file should be created to accompany each pipeline run.

**Required parameters:**

* `assembly` is the name of the genome assembly e.g. `hg38`
* `genome_fasta` is the path to the genome sequence file e.g. `hg38.fa`
* `annotation_file` is the path to the annotation file e.g. `hg38.refGene.gtf`

**Optional parameters:**

* `bowtie2_threads` is the number of threads for bowtie2 to use during index building and alignment e.g. `4`
* `model_temp` is the temperature at which probe hybridization is predicted with XGBoost (e.g. `37`, `42`, `47`, `52`, or `60`)
* `nupack_<param>` settings configure the [NUPACK 4 model](https://piercelab-caltech.github.io/nupack-docs/model/#model-specification) used for secondary structure predictions
* `existing_probes` allows users to provide already designed DNA-FISH probes to the pipeline and then run downstream steps, including thermodynamic, specifity, and secondary structure analysis, as well as intersection with the genome annotations file to yield RNA probe sets.
* `exclude_chroms` allows users to manually exclude particular fasta records from probe design. This is useful in cases where no probes are available for a particular fasta record.

### Genome annotation files

#### GTF-formatted genome annotation file

A GTF-formatted annotation file is required to run the pipeline in order to generate the RNA probe sets. The pipeline can be run on other assemblies found [here](https://hgdownload.soe.ucsc.edu/downloads.html) and at a variety of sources. Currently, the pipeline supports annotations in the following GTF format:

```
chr1  refGene  transcript  11874  14409  .  +  .  gene_id "DDX11L1"; transcript_id "NR_046018"; gene_name "DDX11L1";
chr1  refGene  exon        11874  12227  .  +  .  gene_id "DDX11L1"; transcript_id "NR_046018"; exon_number "1"; exon_id "NR_046018.1"; gene_name "DDX11L1";
chr1  refGene  exon        12613  12721  .  +  .  gene_id "DDX11L1"; transcript_id "NR_046018"; exon_number "2"; exon_id "NR_046018.2"; gene_name "DDX11L1";
chr1  refGene  exon        13221  14409  .  +  .  gene_id "DDX11L1"; transcript_id "NR_046018"; exon_number "3"; exon_id "NR_046018.3"; gene_name "DDX11L1";
chr1  refGene  transcript  14362  29370  .  -  .  gene_id "WASH7P"; transcript_id "NR_024540"; gene_name "WASH7P";
chr1  refGene  exon        14362  14829  .  -  .  gene_id "WASH7P"; transcript_id "NR_024540"; exon_number "1"; exon_id "NR_024540.1"; gene_name "WASH7P";
chr1  refGene  exon        14970  15038  .  -  .  gene_id "WASH7P"; transcript_id "NR_024540"; exon_number "2"; exon_id "NR_024540.2"; gene_name "WASH7P";
chr1  refGene  exon        15796  15947  .  -  .  gene_id "WASH7P"; transcript_id "NR_024540"; exon_number "3"; exon_id "NR_024540.3"; gene_name "WASH7P";
```

Only `exon` features (the third column) will be used in the pipeline. In the `attribute` section (the last column), `gene_id` and `exon_id` are extracted and used in the pipeline. 

#### Converting GFF files to GTF format

In some cases, the only available annotation files are in GFF format:

```
Chr1  TAIR10  chromosome      1     30427671  .  .  .  ID=Chr1;Name=Chr1
Chr1  TAIR10  gene            3631  5899      .  +  .  ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
Chr1  TAIR10  exon            3631  3913      .  +  .  Parent=AT1G01010.1
Chr1  TAIR10  exon            3996  4276      .  +  .  Parent=AT1G01010.1
Chr1  TAIR10  CDS             3996  4276      .  +  2  Parent=AT1G01010.1,AT1G01010.1-Protein;
```

These files are readily converted to GTF format using the [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread) command-line utility: 

```
$ gffread input_file.gff -T -o output_file.gtf
```

**NOTE:** for convenience, `gffread` is included in the provided `paintshop_snakemake` conda environment.
After following the installation instructions, activate the conda environment using `conda activate paintshop_snakemake` 
prior to running the above command.

### A note on sequence IDs

Chromosomes must have identical names in the two input files. For example, if a fasta file contains a `chr4` sequence, then annotations for this sequence must have `chr4` as their `seqid` in the input GTF file, and `chrIV` or `4` would be invalid values. 

#### Inspecting seqid values

To inspect unique `seqid` values in the input fasta file:

```
$ cat example_run/data/example.fa | grep ">" | cut -c 2- | sortc
GTFhrI
chrII
chrIII
chrIV
```

To inspect unique `seqid` values in the input GTF file:

```
$ awk '{a[$1]++} END {for (b in a) {print b}}' example_run/data/example.gtf | sort
chrI
chrII
chrIII
chrIV
```

#### Modifying seqid values

If the `seqid` values do not match, then the names in at least one file must be changed such that they do match.

Changing the names in both files is demonstrated below. If you only need to change one file, you can 
call just one of the `rename_fasta_ids` / `rename_gtf_ids` functions and exclude the unused file paths. This code outputs new files leaving input files intact:

```python
from io import StringIO
import pandas as pd
from Bio import SeqIO


# old name: new name
RENAME_DICT = {
    'chrI': 'chr1',
    'chrII': 'chr2',
    'chrIII': 'chr3',
    'chrIV': 'chr4',
}

# configure file paths
INPUT_FA = 'example_run/data/example.fa'
OUTPUT_FA = 'example_run/data/example_renamed.fa'
INPUT_GTF = 'example_run/data/example.gtf'
OUTPUT_GTF = 'example_run/data/example_renamed.gtf'


# fasta helper function
def rename_fasta_ids():
    # load input fasta
    recs = SeqIO.parse(INPUT_FA, 'fasta')

    # rename fasta records according to lookup table
    output_seqs = []
    for seq in recs:
        seq.id = RENAME_DICT[seq.id]
        seq.description = seq.id
        output_seqs.append(seq)   

    # write output fasta
    with open(OUTPUT_FA, 'w') as outfile:
        SeqIO.write(output_seqs, outfile, 'fasta')


# gtf helper function
def rename_gtf_ids():
    
    # load input gtf file skipping commented lines
    csv_data = ''
    with open(INPUT_GTF, 'r') as infile:
        for line in infile.readlines():
            if not line.startswith('#'):
                csv_data += line
    
    # construct a pandas dataframe from gtf data tsv string
    df = pd.read_csv(StringIO(csv_data), sep='\t', header=None)
    
    # rename chromosome ids
    df[0] = df[0].apply(lambda x: RENAME_DICT[x])
    
    # save renamed GTF file
    df.to_csv(OUTPUT_GTF, sep='\t', header=None, index=False)


# create output fasta file with new seqids
rename_fasta_ids()

# create output gtf file with new seqids
rename_gtf_ids()
```

**NOTE:** the modules required to run the above code are included in the provided `paintshop_snakemake` conda environment.
After following the installation instructions, activate the conda environment using `conda activate paintshop_snakemake` 
prior to running the above code.
