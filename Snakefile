"""
HMCL NCATS MIPE 4.0 Drug Screen Dataset Pipeline
"""
import pathlib
import pandas as pd

# output directory
out_dir = pathlib.Path("/data/proj/hmcl/drug-screen-manuscript")

out_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

# get lists of plate ids, cell lines, and drugs
dat = pd.read_csv("data/raw.tsv.gz", sep="\t")

cell_lines = set(dat.cell_line)
plate_ids = set(dat.plate)

rule normalize_plates:
    input:
        out_dir.joinpath("drug_plates/combined/raw.tsv"),
        out_dir.joinpath("metadata/plate-metadata.tsv"),
    params:
        plate_ids=plate_ids
    output:
        out_dir.joinpath("drug_plates/combined/normed.tsv"),
        expand(out_dir.joinpath("drug_plates/{plate}/01-raw.tsv"), plate=plate_ids),
        expand(out_dir.joinpath("drug_plates/{plate}/02-normed.tsv"), plate=plate_ids),
    script:
        "scripts/normalize_plates.R"

rule create_combined_plate_matrix:
    input:
        "data/raw.tsv.gz"
    output:
        out_dir.joinpath("drug_plates/combined/raw.tsv"),
        out_dir.joinpath("metadata/plate-metadata.tsv"),
        out_dir.joinpath("metadata/drug-groups.tsv"),
    script:
        "scripts/create_combined_plate_matrix.R"
