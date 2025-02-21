"""
HMCL NCATS MIPE 4.0 Drug Screen Dataset Pipeline
https://github.com/khughitt/hmcl-drug-screen-pipeline
"""
import pathlib
import pandas as pd

configfile: "config/config.yml"

# output directory
out_dir = pathlib.Path(config["output_dir"])
out_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

# create a list of plate ids
dat = pd.read_csv("data/raw.tsv.gz", sep="\t")

plate_ids = sorted(list(set(dat.plate)))

mask = dat.cell_line.isin(config["outlier_cells"])
filtered_plate_ids = sorted(list(set(dat.plate[~mask])))

rule all:
    input:
        out_dir.joinpath("manuscript/figures/fig3a.tiff"),
        out_dir.joinpath("manuscript/figures/fig3b.tiff"),
        out_dir.joinpath("plates/images/background_plate.png"),
        out_dir.joinpath("datapackage.yml")

rule package_results:
    input:
        # 1. raw data
        out_dir.joinpath("raw/data.tsv"),

        # 2-3. manuscript tables
        out_dir.joinpath("manuscript/tables/table1.tsv"),
        out_dir.joinpath("manuscript/tables/table2.tsv"),

        # 4-10. plate-level data
        out_dir.joinpath("plates/raw.tsv"),
        out_dir.joinpath("plates/raw_filtered.tsv"),
        out_dir.joinpath("plates/normed.tsv"),
        out_dir.joinpath("plates/background_adjusted.tsv"),
        out_dir.joinpath("plates/background.tsv"),
        out_dir.joinpath("plates/concentrations.tsv"),
        out_dir.joinpath("plates/metadata.tsv"),

        # 11-14. drug-level data
        out_dir.joinpath("drugs/ac50.tsv"),
        out_dir.joinpath("drugs/drug_curves.tsv"),
        out_dir.joinpath("drugs/drug_indices.tsv"),
        out_dir.joinpath("drugs/metadata.tsv"),

        # 14. mutation data
        out_dir.joinpath("mutations/hmcl-predicted-mutations.tsv")
    output:
        out_dir.joinpath("datapackage.yml")
    script:
        "scripts/package_results.py"

rule create_figures:
    input:
        out_dir.joinpath("drugs/drug_curves.tsv"),
        out_dir.joinpath("drugs/ac50.tsv")
    output:
        out_dir.joinpath("manuscript/figures/fig3a.tiff"),
        out_dir.joinpath("manuscript/figures/fig3b.tiff")
    script:
        "scripts/create_figures.R"

rule visualize_plates:
    input:
        out_dir.joinpath("plates/raw_filtered.tsv"),
        out_dir.joinpath("plates/normed.tsv"),
        out_dir.joinpath("plates/background_adjusted.tsv"),
        out_dir.joinpath("plates/background.tsv"),
        out_dir.joinpath("plates/metadata.tsv")
    output:
        indiv=expand(out_dir.joinpath("plates/images/indiv/{plate_id}.jpg"), plate_id=filtered_plate_ids),
        mean=out_dir.joinpath("plates/images/mean_plate.png"),
        median=out_dir.joinpath("plates/images/median_plate.png"),
        background=out_dir.joinpath("plates/images/background_plate.png")
    params:
        plate_ids=filtered_plate_ids
    script:
        "scripts/visualize_plates.R"

rule create_cell_metadata:
    input:
        "data/cell-metadata.tsv",
    output:
        out_dir.joinpath("manuscript/tables/table1.tsv")
    script:
        "scripts/create_cell_metadata.R"

rule compute_ac50_summary:
  input:
      out_dir.joinpath("drugs/ac50.tsv")
  output:
      out_dir.joinpath("manuscript/tables/table2.tsv"),
  script:
      "scripts/compute_ac50_summary.R"

rule create_ac50_matrix:
    input:
        out_dir.joinpath("drugs/drug_curves.tsv"),
    output:
        out_dir.joinpath("drugs/ac50.tsv")
    script:
        "scripts/create_ac50_matrix.R"

rule fit_dose_response_curves:
    input:
        out_dir.joinpath("plates/background_adjusted.tsv"),
        out_dir.joinpath("drugs/drug_indices.tsv"),
    output:
        out_dir.joinpath("drugs/drug_curves.tsv"),
    script:
        "scripts/fit_dose_response_curves.R"

rule background_adjustment:
    input:
        out_dir.joinpath("plates/normed.tsv"),
    output:
        out_dir.joinpath("plates/background.tsv"),
        out_dir.joinpath("plates/background_adjusted.tsv"),
    script:
        "scripts/background_adjustment.R"

rule normalize_plates:
    input:
        out_dir.joinpath("plates/raw_filtered.tsv"),
        out_dir.joinpath("plates/metadata.tsv"),
    output:
        out_dir.joinpath("plates/normed.tsv"),
    script:
        "scripts/normalize_plates.R"

rule filter_outlier_cells:
    input:
        out_dir.joinpath("plates/raw.tsv"),
        out_dir.joinpath("plates/metadata.tsv"),
    output:
        out_dir.joinpath("plates/raw_filtered.tsv"),
    script:
        "scripts/filter_outlier_cells.R"

rule create_plate_matrices:
    input:
        out_dir.joinpath("raw/data.tsv")
    output:
        out_dir.joinpath("plates/raw.tsv"),
        out_dir.joinpath("plates/concentrations.tsv"),
        out_dir.joinpath("plates/metadata.tsv"),
        out_dir.joinpath("drugs/drug_indices.tsv"),
    params:
        plate_ids=plate_ids
    script:
        "scripts/create_plate_matrices.R"

rule copy_mutation_data:
    input:
        "data/hmcl-predicted-mutations.tsv"
    output:
        out_dir.joinpath("mutations/hmcl_predicted_mutations.tsv")
    shell:
        "cp {input} {output}"

rule copy_drug_metadata:
    input:
        "data/drug-metadata.tsv"
    output:
        out_dir.joinpath("drugs/metadata.tsv")
    shell:
        "cp {input} {output}"

rule copy_raw_data:
    input:
        "data/raw.tsv.gz"
    output:
        out_dir.joinpath("raw/data.tsv")
    shell:
        "cp {input} {output}.gz && gunzip {output}.gz"

# vi:ft=snakemake
