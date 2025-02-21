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

# 1. raw data
res = create_resource(snek.input[0], "raw_data", "Raw drugscreen data")
resources.append(res)

# 2-3. manuscript tables
res = create_resource(snek.input[1], "cell_metadata", "Table 1) Cell line metadata")
res.sources = [{
    "title": "Keats Lab Data Repository",
    "path": "https://www.keatslab.org/data-repository"
}]
resources.append(res)

res = create_resource(snek.input[2], "ac50_summary", "Table 2) AC-50 summary")
resources.append(res)

# 4-10. drug plates
res = create_resource(snek.input[3], "plates_raw", "Drug plate matrix (raw)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[4], "plates_raw_filtered", "Drug plate matrix (raw / filtered)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[5], "plates_normed", "Drug plate matrix (normed)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[6], "plates_bgadj", "Drug plate matrix (background-adjusted)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[7], "plate_background", "Drug plate matrix (background)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[8], "plate_concentrations", "Drug plate matrix (concentrations)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[9], "plate_metadata", "Plate metadata")
resources.append(res)

# 11-14. drug-level data
res = create_resource(snek.input[10], "ac50_matrix", "AC-50 Matrix", include_schema=False)
resources.append(res)

res = create_resource(snek.input[11], "drug_curves", "Drug curves")
resources.append(res)

res = create_resource(snek.input[12], "drug_indices", "Drug, cell plate and dose for each plate and position")
resources.append(res)

res = create_resource(snek.input[13], "drug_metadata", "Drug metadata")
resources.append(res)

res = create_resource(snek.input[14], "mutation_data", "Predicted mutations")
res.sources = [{
    "title": "Keats Lab Data Repository",
    "path": "https://www.keatslab.org/data-repository"
}]

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
    homepage="https://zenodo.org/records/14902712",
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
