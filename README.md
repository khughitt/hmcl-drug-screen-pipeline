# HMCL NCATS MIPE 4.0 Drug Screen Dataset Pipeline

## Overview

This repository contains code for a [Snakemake](https://snakemake.readthedocs.io/) pipeline used to
perform basic pre-processing and normalization for a large-scale human multiple myeloma cell line
("HMCL") drug screen dataset generated at the [National Center for Advancing Translational Sciences
(NCATS)](https://ncats.nih.gov/).

## Data availability..

## Paper..

## Usage

If you wish to run the data processing pipeline yourself, the first step is to clone this repository 
and initialize a conda environment with the required dependencies:

```
git clone xx
```

Download and install [conda](https://docs.conda.io/en/latest/) (or a variant such as
[mamba](https://mamba.readthedocs.io), which is used below) if it is not already present on your
system.

TODO: add option specifying how one can use exported environment file / specific versions..

```
mamba create -c conda-forge --file requirements.txt -n "hmcl"
```

## Related data

In addition to the drug screen data described here, other omics data including RNA-Seq, CNV, and
mutation data is available for the same same lines from the [Keats lab](https://sites.google.com/site/jonathankeatslab/data-repository).
