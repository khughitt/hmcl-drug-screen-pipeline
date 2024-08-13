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

cell_lines = sorted(list(set(dat.cell_line)))
plate_ids = sorted(list(set(dat.plate)))

# exclude outlier cell lines downstream plots
outliers = ["KMS21BM_JCRB", "Karpas417_ECACC"]
cell_lines = [x for x in cell_lines if x not in outliers]

# filepaths for cell line-specific drug cluster plots
cell_filepaths = [os.path.join(out_dir, "fig", "drugs", f"drug_curves_by_cluster_{cell}.png") for cell in cell_lines]

# drug response fields
num_conc = 11
response_fields = ["ac50", "lac50"] + [f"dose_{i}" for i in range(num_conc)]

rule all:
    input:
        out_dir.joinpath("similarity/cells.tsv"),
        out_dir.joinpath("similarity/drugs.tsv"),
        out_dir.joinpath("clusters/cells.tsv"),
        out_dir.joinpath("clusters/drugs.tsv"),
        out_dir.joinpath("enrichment/drug_cluster_annotation_enrichment.tsv"),
        out_dir.joinpath("fig/drugs/drug_umap_clusters.png"),
        out_dir.joinpath("fig/cells/cell_average_viability.png"),
        out_dir.joinpath("xlsx/drug_ac50.xslx"),
        out_dir.joinpath("xlsx/drug_clusters.xslx"),
        out_dir.joinpath("xlsx/cell_clusters.xslx")

rule visualize_drug_clusters:
    input:
        out_dir.joinpath("drug_curves/drug_curves_filtered.tsv"),
        out_dir.joinpath("projections/similarity/drugs_umap.tsv"),
        out_dir.joinpath("clusters/drugs.tsv"),
        out_dir.joinpath("clusters/cells.tsv"),
        "data/drug_metadata.tsv"
    output:
        os.path.join(out_dir, "fig/drugs/drug_umap_clusters.png"),
        os.path.join(out_dir, "fig/drugs/drug_umap_super_clusters.png"),
        os.path.join(out_dir, "fig/drugs/drug_cluster_mean_ac50.png"),
        os.path.join(out_dir, "fig/drugs/drug_curves_by_cluster_all_cells.png"),
        cell_filepaths
    script:
        "scripts/visualize_drug_clusters.R"

rule visualize_average_cell_response:
    input:
        out_dir.joinpath("drug_curves/drug_curves.tsv"),
    output:
        os.path.join(out_dir, "fig/cells/cell_average_viability.png")
    script:
        "scripts/visualize_average_cell_response.R"

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

rule create_result_tables:
    input:
        out_dir.joinpath("drug_curves/drug_curves.tsv"),
        out_dir.joinpath("clusters/cells.tsv"),
        out_dir.joinpath("clusters/drugs.tsv"),
        "data/drug_metadata.tsv"
    output:
        out_dir.joinpath("xlsx/drug_ac50.xslx"),
        out_dir.joinpath("xlsx/cell_clusters.xslx"),
        out_dir.joinpath("xlsx/drug_clusters.xslx")
    script:
        "scripts/create_result_tables.R"

rule quantify_drug_cluster_annotation_enrichment:
    input:
        out_dir.joinpath("similarity/drugs.tsv"),
        out_dir.joinpath("clusters/drugs.tsv"),
        "data/drug_metadata.tsv"
    output:
        out_dir.joinpath("enrichment/drug_cluster_annotation_enrichment.tsv")
    script:
        "scripts/quantify_drug_cluster_annotation_enrichment.R"

rule cluster_cells:
    input:
        out_dir.joinpath("similarity/cells.tsv"),
    output:
        out_dir.joinpath("clusters/cells.tsv"),
    script:
        "scripts/cluster_cells.R"

rule cluster_drugs:
    input:
        out_dir.joinpath("similarity/drugs.tsv"),
    output:
        out_dir.joinpath("clusters/drugs.tsv"),
    script:
        "scripts/cluster_drugs.R"

rule reduce_similarity_matrix_dimensions:
    input:
        out_dir.joinpath("similarity/cells.tsv"),
        out_dir.joinpath("similarity/drugs.tsv"),
    output:
        out_dir.joinpath("projections/similarity/cells_pca.tsv"),
        out_dir.joinpath("projections/similarity/cells_pca_var.txt"),
        out_dir.joinpath("projections/similarity/cells_umap.tsv"),
        out_dir.joinpath("projections/similarity/drugs_pca.tsv"),
        out_dir.joinpath("projections/similarity/drugs_pca_var.txt"),
        out_dir.joinpath("projections/similarity/drugs_umap.tsv"),
    script:
        "scripts/reduce_dimensions.R"

rule compute_drug_similarity:
    input:
        out_dir.joinpath("combined_viability_matrices/drugs.tsv"),
    output:
        out_dir.joinpath("similarity/drugs.tsv"),
    script:
        "scripts/compute_drug_similarity.R"

rule compute_cell_similarity:
    input:
        out_dir.joinpath("combined_viability_matrices/cells.tsv"),
    output:
        out_dir.joinpath("similarity/cells.tsv"),
    script:
        "scripts/compute_cell_similarity.R"

rule reduce_processed_data_dimensions:
    input:
        out_dir.joinpath("combined_viability_matrices/cells.tsv"),
        out_dir.joinpath("combined_viability_matrices/drugs.tsv"),
    output:
        out_dir.joinpath("projections/data/cells_pca.tsv"),
        out_dir.joinpath("projections/data/cells_pca_var.txt"),
        out_dir.joinpath("projections/data/cells_umap.tsv"),
        out_dir.joinpath("projections/data/drugs_pca.tsv"),
        out_dir.joinpath("projections/data/drugs_pca_var.txt"),
        out_dir.joinpath("projections/data/drugs_umap.tsv"),
    script:
        "scripts/reduce_dimensions.R"

rule create_combined_response_matrices:
    input:
        expand(os.path.join(out_dir, "drug_response_matrices/{response}.tsv"), response=response_fields)
    output:
        out_dir.joinpath("combined_viability_matrices/cells.tsv"),
        out_dir.joinpath("combined_viability_matrices/drugs.tsv"),
    script:
        "scripts/create_combined_response_matrices.R"

rule create_drug_response_matrices:
    input:
        out_dir.joinpath("drug_curves/drug_curves_filtered.tsv"),
    output:
        os.path.join(out_dir, "drug_response_matrices/{response}.tsv")
    script:
        "scripts/create_dose_response_matrices.R"

rule filter_outlier_cells:
    input:
        out_dir.joinpath("drug_curves/drug_curves.tsv"),
    output:
        out_dir.joinpath("drug_curves/drug_curves_filtered.tsv"),
    script:
        "scripts/filter_outlier_cells.R"

rule fit_dose_response_curves:
    input:
        out_dir.joinpath("drug_plates/background_adjusted.tsv"),
        out_dir.joinpath("metadata/drug-indices.tsv"),
    output:
        out_dir.joinpath("drug_curves/drug_curves.tsv"),
    script:
        "scripts/fit_dose_response_curves.R"

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
        out_dir.joinpath("metadata/plate-metadata.tsv"),
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
        out_dir.joinpath("metadata/drug-indices.tsv"),
        out_dir.joinpath("metadata/drug-groups.tsv"),
    params:
        plate_ids=plate_ids
    script:
        "scripts/create_plate_matrices.R"
