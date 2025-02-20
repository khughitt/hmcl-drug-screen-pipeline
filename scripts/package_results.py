"""
Creates a data package with the main result datasets + metadata

Schemas are excluded for matrix results since they don't provide a lot of useful information and
would clutter to the resulting datapackage.yml.
"""
import datetime
import os
from frictionless import describe, Package

snek = snakemake

# abstract from the manuscript
abstract="Multiple myeloma, a hematopoietic malignancy of terminally differentiated B cells, is the second most common hematological malignancy after leukemia. While patients have benefited from numerous advances in treatment in recent years resulting in significant increases to average survival time following diagnosis, myeloma remains incurable and relapse is common. To help identify novel therapeutic agents with efficacy against the disease and to search for biomarkers associated with differential response to treatment, a large-scale pharmacological screen was performed with 1,912 small molecule compounds tested at 11 doses for 47 human myeloma cell lines (HMCL). Raw and processed versions of the drug screen dataset are provided, as well as supportive information including drug and cell line metadata and high-level characterization of the most salient features of each. The dataset is publicly available at Zenodo and the workflow code used for data processing and generation of supporting figures and tables can be found at https://github.com/khughitt/hmcl-drug-screen-pipeline."

# helper function to create data resource objects
def create_resource(path:str, name:str, title:str, include_schema=True):
    res = describe(path, stats=True)
    res.name = name
    res.title = title

    if not include_schema:
        res.schema = None

    return res

# Create dataset Resource objects
resources = []

# 1-6. drug plates
res = create_resource(snek.input[0], "plates_raw", "Drug plate matrix (raw)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[1], "plates_raw_filtered", "Drug plate matrix (raw / filtered)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[2], "plates_normed", "Drug plate matrix (normed)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[3], "plates_bgadj", "Drug plate matrix (background-adjusted)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[4], "plate_background", "Drug plate matrix (background)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[5], "plate_concentrations", "Drug plate matrix (concentrations)", include_schema=False)
resources.append(res)

# 7. drug curves
res = create_resource(snek.input[6], "drug_curves", "Drug curves")
resources.append(res)

# 8-20. cell x drug matrices (AC-50, LAC-50, single doses)
res = create_resource(snek.input[7], "ac50_matrix", "AC-50 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[8], "lac50_matrix", "LAC-50 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[9], "dose_0_matrix", "Dose 0 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[10], "dose_1_matrix", "Dose 1 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[11], "dose_2_matrix", "Dose 2 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[12], "dose_3_matrix", "Dose 3 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[13], "dose_4_matrix", "Dose 4 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[14], "dose_5_matrix", "Dose 5 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[15], "dose_6_matrix", "Dose 6 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[16], "dose_7_matrix", "Dose 7 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[17], "dose_8_matrix", "Dose 8 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[18], "dose_9_matrix", "Dose 9 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[19], "dose_10_matrix", "Dose 10 Matrix", include_schema=False)
resources.append(res)

# 21-22. similarity matrices
res = create_resource(snek.input[20], "cell_similarity", "Cell Similarity Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[21], "drug_similarity", "Drug Similarity Matrix", include_schema=False)
resources.append(res)

# 23-26. similarity matrix PCA/UMAP projections
res = create_resource(snek.input[22], "cell_pca", "Cell Similarity Matrix (PCA)")
resources.append(res)

res = create_resource(snek.input[23], "cell_umap", "Cell Similarity Matrix (UMAP)")
resources.append(res)

res = create_resource(snek.input[24], "drug_pca", "Drug Similarity Matrix (PCA)")
resources.append(res)

res = create_resource(snek.input[25], "drug_umap", "Drug Similarity Matrix (UMAP)")
resources.append(res)

# 27-28. clusters
res = create_resource(snek.input[26], "cell_clusters", "Cell Clusters")
resources.append(res)

res = create_resource(snek.input[27], "drug_clusters", "Drug Clusters")
resources.append(res)

# 29. cell average viability curves
res = create_resource(snek.input[28], "average_cell_viability", "Average Cell Viability Curves")
resources.append(res)

# 30-31. combined viability matrices
res = create_resource(snek.input[29], "combined_viability_cells", "All viability measurements by cell", include_schema=False)
resources.append(res)

res = create_resource(snek.input[30], "combined_viability_drugs", "All viability measurements by drug", include_schema=False)
resources.append(res)

# 32-35. metadata
res = create_resource(snek.input[31], "cell_metadata", "Cell line metadata")
res.sources = [{
    "title": "Keats Lab Data Repository",
    "path": "https://www.keatslab.org/data-repository"
}]
resources.append(res)

res = create_resource(snek.input[32], "drug_indices", "Drug, cell plate and dose for each plate and position")
resources.append(res)

res = create_resource(snek.input[33], "drug_metadata", "Drug metadata")
resources.append(res)

res = create_resource(snek.input[34], "plate_metadata", "Plate metadata")
resources.append(res)

res = create_resource(snek.input[35], "mutation_data", "Predicted mutations")
res.sources = [{
    "title": "Keats Lab Data Repository",
    "path": "https://www.keatslab.org/data-repository"
}]

resources.append(res)

res = create_resource(snek.input[36], "raw_data", "Raw drugscreen data")
resources.append(res)

# convert absolute resource paths to relative paths
for resource in resources:
    resource.path = os.path.relpath(resource.path, snek.config["output_dir"])

# create data package
pkg = Package(
    resources=resources,
    name="hmcl_drug_screen",
    title="Human Myeloma Cell Line (HMCL) NCATS MIPE 4.0 Drug Screen Dataset",
    description="ABSTRACT: " + abstract,
    homepage="https://zenodo.org/records/14867030",
    created=datetime.datetime.now().isoformat(),
    licenses=[{
        "name": "CC-BY-4.0",
        "path": "https://creativecommons.org/licenses/by/4.0/",
        "title": "Creative Commons Attribution 4.0"
    }],
    contributors=[{
        "title": "V. Keith Hughitt",
        "email": "keith.hughitt@nih.gov",
        "roles": ["creator"]
    },
    {
        "title": "John Simmons",
        "roles": ["creator"]
    },
    {
        "title": "Beverly A. Mock",
        "email": "mockb@mail.nih.gov",
        "roles": ["contact"]
    }]
)

pkg.to_yaml(snek.output[0])
