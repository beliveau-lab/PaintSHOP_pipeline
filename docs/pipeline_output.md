# Pipeline Output

PaintSHOP pipeline output specification.

### Overview

Pipeline outputs are organized into three directories:

| Folder        | Description                                                       |
|---------------|-------------------------------------------------------------------|
| [01_reference_files/](../example_run/expected_pipeline_output/01_reference_files)   | files that may be of use in other pipelines or analyses |
| [02_intermediate_files/](../example_run/expected_pipeline_output/02_intermediate_files)  | large intermediate files, useful when debugging, but otherwise disposable |
| [03_output_files/](../example_run/expected_pipeline_output/03_output_files) | all DNA/RNA probe sets as both .tsv files and zip archives |

All unwanted files can be safely deleted once the pipeline is run. 

**NOTE:** To minimize disk usage, it may be desirable to keep only the three zip archives found in [03_output_files/04_zip_archives/](../example_run/expected_pipeline_output/03_output_files/04_zip_archives), as the intermediate files are very large (e.g. the pipeline will generate over 200 GB of intermediate files for hg38).

### Probes

#### File locations

| Item        | Location                                                       |
|---------------|-------------------------------------------------------------------|
| DNA FISH probes in .tsv format | [03_output_files/01_dna_probes/](../example_run/expected_pipeline_output/03_output_files/01_dna_probes) | 
| Isoform-resolved RNA FISH probes in .tsv format | [03_output_files/02_rna_probes_all/](../example_run/expected_pipeline_output/03_output_files/02_rna_probes_all) |
| Isoform-flattened RNA FISH probes in .tsv format | [03_output_files/03_rna_probes_iso/](../example_run/expected_pipeline_output/03_output_files/03_rna_probes_iso) |
| All DNA/RNA FISH probe sets as compressed .zip archives | [03_output_files/04_zip_archives/](../example_run/expected_pipeline_output/03_output_files/04_zip_archives) |

#### DNA probe sets

The columns in DNA probe files are: 

| # | Column | Description |
|---|--------|-------------|
| 0 | chromosome | the `seqid` that the probe targets |
| 1 | start | genomic start coordinate |
| 2 | stop | genomic stop coordinate |
| 3 | sequence | DNA sequence of the probe |
| 4 | Tm | Melting temperature of the probe DNA sequence |
| 5 | on-target score (0-100) | hybridization prediction at the target site of the probe |
| 6 | off-target score (0-10,000) | sum of the hybridization predictions for all off-target sites (up to 100) for the probe |
| 7 | repeat (0 or 1) | 0 if the sequence contains no repeat-masked bases, and 1 if it does  |
| 8 | prob | the probability that the probe has no secondary structure |
| 9 | max k-mer | max count out of all k-mers in the probe sequence |
| 10 | strand | probe strand, `+` or `-` |

#### Isoform-resolved RNA probe sets

The columns in isoform-resolved RNA probe sets are:

| # | Column | Description |
|---|--------|-------------|
| 0 | chromosome | the `seqid` that the probe targets |
| 1 | start | genomic start coordinate |
| 2 | stop | genomic stop coordinate |
| 3 | sequence | DNA sequence of the probe |
| 4 | Tm | Melting temperature of the probe DNA sequence |
| 5 | on-target score (0-100) | hybridization prediction at the target site of the probe |
| 6 | off-target score (0-10,000) | sum of the hybridization predictions for all off-target sites (up to 100) for the probe |
| 7 | repeat (0 or 1) | 0 if the sequence contains no repeat-masked bases, and 1 if it does  |
| 8 | prob | the probability that the probe has no secondary structure |
| 9 | max k-mer | max count out of all k-mers in the probe sequence |
| 10 | strand | probe strand, `+` or `-` |
| 11 | refseq | The transcript ID that the probe targets, stripped of any version suffixes e.g. NM_001180043 |
| 12 | transcript_id | The unmodified transcript ID of the transcript that the probe targets e.g. NM_001180043.1 |
| 13 | gene_id | The gene ID of the gene whose transcript the probe targets e.g. PAU8 |

#### Isoform-flattened RNA probe sets

The columns in isoform-flattened RNA probe sets are:

| # | Column | Description |
|---|--------|-------------|
| 0 | chromosome | the `seqid` that the probe targets |
| 1 | start | genomic start coordinate |
| 2 | stop | genomic stop coordinate |
| 3 | sequence | DNA sequence of the probe |
| 4 | Tm | Melting temperature of the probe DNA sequence |
| 5 | on-target score (0-100) | hybridization prediction at the target site of the probe |
| 6 | off-target score (0-10,000) | sum of the hybridization predictions for all off-target sites (up to 100) for the probe |
| 7 | repeat (0 or 1) | 0 if the sequence contains no repeat-masked bases, and 1 if it does  |
| 8 | prob | the probability that the probe has no secondary structure |
| 9 | max k-mer | max count out of all k-mers in the probe sequence |
| 10 | strand | probe strand, `+` or `-` |
| 11 | gene_id | The gene ID of the gene whose transcript the probe targets e.g. PAU8 |
| 12 | transcripts |  The number of isoforms that this probe targets |

### Reporting

An HTML report with diagnostics and detailed pipeline information can by generated with the following command:

```
$ snakemake --snakefile path/to/Snakefile --configfile path/to/config.yml --report pipeline_output/report.html
```

An example report is available [here](https://paintshop-bucket.s3.amazonaws.com/static/report.html). For a visualization of the pipeline DAG structure, see: [pipeline.pdf](../example_run/expected_pipeline_output/pipeline.pdf) or [pipeline.svg](../example_run/expected_pipeline_output/pipeline.svg)
