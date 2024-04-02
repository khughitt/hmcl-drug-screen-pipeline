"""
HMCL NCATS MIPE 4.0 Drug Screen Dataset Pipeline
"""
import os
import pathlib
import pandas as pd

# output directory
out_dir = pathlib.Path("/data/proj/hmcl/drug-screen-manuscript")

out_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

# get lists of plate ids, cell lines, and drugs
dat = pd.read_csv("data/raw.tsv.gz", sep="\t")

cell_lines = set(dat.cell_line)
plate_ids = set(dat.plate)

rule visualize_plates:
    input:
        out_dir.joinpath("drug_plates/raw.tsv"),
        out_dir.joinpath("drug_plates/normed.tsv"),
        out_dir.joinpath("drug_plates/background_adjusted.tsv"),
        out_dir.joinpath("drug_plates/background.tsv"),
        out_dir.joinpath("metadata/plate-metadata.tsv")
    output:
        indiv=expand(os.path.join(out_dir, "fig/plates/{plate_id}.png"), plate_id=plate_ids),
        mean=out_dir.joinpath("fig/mean_plate.png"),
        median=out_dir.joinpath("fig/median_plate.png"),
        background=out_dir.joinpath("fig/background_plate.png")
    params:
        plate_ids=plate_ids
    script:
        "scripts/visualize_plates.R"

rule background_adjustment:
    input:
        out_dir.joinpath("drug_plates/normed.tsv"),
    output:
        out_dir.joinpath("drug_plates/background.tsv"),
        out_dir.joinpath("drug_plates/background_adjusted.tsv"),
    script:
        "scripts/background_adjustment.R"

rule normalize_plates:
    input:
        out_dir.joinpath("drug_plates/raw.tsv"),
    output:
        out_dir.joinpath("drug_plates/normed.tsv"),
    script:
        "scripts/normalize_plates.R"

rule create_plate_matrices:
    input:
        "data/raw.tsv.gz"
    output:
        out_dir.joinpath("drug_plates/raw.tsv"),
        out_dir.joinpath("drug_plates/concentrations.tsv"),
        out_dir.joinpath("metadata/plate-metadata.tsv"),
        out_dir.joinpath("metadata/drug-groups.tsv"),
    params:
        plate_ids=plate_ids
    script:
        "scripts/create_plate_matrices.R"
