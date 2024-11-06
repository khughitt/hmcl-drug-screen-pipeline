# HMCL NCATS MIPE 4.0 Drug Screen Dataset Pipeline

## Overview

This repository contains code for a [Snakemake](https://snakemake.readthedocs.io/) pipeline used to
perform basic pre-processing and normalization for a large-scale human multiple myeloma cell line
("HMCL") drug screen dataset generated at the [National Center for Advancing Translational Sciences
(NCATS)](https://ncats.nih.gov/).

## Data availability

Raw + processed versions of the HMCL drugscreen data is available on Zenodo:

- [https://zenodo.org/records/13910207](https://zenodo.org/records/13910207)

## Usage

If you wish to run the data processing pipeline yourself, start by cloning the repo:

```
git clone https://github.com/khughitt/hmcl-drug-screen-pipeline
```

Next, download and install [conda](https://docs.conda.io/en/latest/) (or a variant such as
[mamba](https://mamba.readthedocs.io), which is used below) if it is not already present on your
system.

Once you have `conda` or `mamba` installed, create a new conda environment with the necessary
requiments:

```
mamba create -c conda-forge --file requirements.txt -n "hmcl"
```

If you are using `conda`, simply replace "mamba" with "conda".

Next, to configure the pipeline parameters, including where output should be stored, edit the config
file `config/config.yml` and adjust the settings as desired.

The default settings can be used to reproduce the tables, figures, and intermediate datasets
described in the manuscript and available on Zenodo.

Finally, to run the pipeline, activate the conda environment and call `snakemake`, specifying the
desired number of threads, e.g.:

```
mamba activate hmcl
snakemake -j2
```

## Related data

In addition to the drug screen data described here, other omics data including RNA-Seq, CNV, and
mutation data is available for the same same lines from the [Keats lab](https://sites.google.com/site/jonathankeatslab/data-repository).
