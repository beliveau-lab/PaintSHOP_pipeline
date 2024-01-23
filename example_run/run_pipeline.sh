
# configure file paths
CONFIG_FILE='config.yml'
SNAKE_FILE='../workflow/Snakefile'
CONDA_ENVS='../shared_conda_envs'

# activate conda environment
source activate paintshop_snakemake

# run the pipeline
snakemake --configfile config.yml --snakefile $SNAKE_FILE \
    --use-conda --conda-prefix $CONDA_ENVS --cores \
    --restart-times 3 --conda-frontend mamba

# export PDF and svg visualizations of the DAG structure of pipeline steps
echo -e "Exporting pipeline DAG to svg and pdf..."
snakemake --configfile config.yml --snakefile $SNAKE_FILE --dag > dag.dot
dot -Tpdf dag.dot > pipeline_output/pipeline.pdf
dot -Tsvg dag.dot > pipeline_output/pipeline.svg
rm dag.dot

echo -e "Generating pipeline HTML report..."
snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE --report pipeline_output/report.html --conda-frontend mamba

# success
echo -e "\nDONE!\n"
