# scnprep_scanpy

## Create working environment

```
# using pip
pip install -r requirements.txt

# using Conda
conda create --name <env_name> --file requirements.txt
conda activate <env_name>
```

## Run pipeline

```
snakemake -j 4
```