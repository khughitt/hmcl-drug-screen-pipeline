"""
Creates a data package with the main result datasets + metadata

Schemas are excluded for matrix results since they don't provide a lot of useful information and
would clutter to the resulting datapackage.yml.
"""
import os
from frictionless import describe, Package

snek = snakemake

# helper function to create data resource objects
def create_resource(path:str, name:str, title:str, include_schema=True):
    res = describe(path, stats=True)
    res.name = "plates_raw"
    res.title = "Drug plate matrix (raw)"

    if not include_schema:
        res.schema = None

    return res

# Create dataset Resource objects
resources = []

# 1-3. drug plates
res = create_resource(snek.input[0], "plates_raw", "Drug plate matrix (raw)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[1], "plates_normed", "Drug plate matrix (normed)", include_schema=False)
resources.append(res)

res = create_resource(snek.input[2], "plates_bgadj", "Drug plate matrix (background-adjusted)", include_schema=False)
resources.append(res)

# 4. plate background
res = create_resource(snek.input[3], "plate_background", "Drug plate matrix (background)", include_schema=False)
resources.append(res)

# 5. drug concentrations
res = create_resource(snek.input[4], "plate_concentrations", "Drug plate matrix (concentrations)", include_schema=False)
resources.append(res)

# 6. drug curves
res = create_resource(snek.input[5], "drug_curves", "Drug curves")
resources.append(res)

# 7. drug response matrix (AC-50)
res = create_resource(snek.input[6], "ac50_matrix", "AC-50 Matrix", include_schema=False)
resources.append(res)

# 8. drug pca
res = create_resource(snek.input[7], "drug_pca", "Drug PCA")
resources.append(res)

# 9. drug similarity matrix
res = create_resource(snek.input[8], "drug_similarity", "Drug Similarity Matrix", include_schema=False)
resources.append(res)

# 10. drug similarity matrix (umap projection)
res = create_resource(snek.input[9], "drug_umap", "Drug Similarity Matrix (UMAP)")
resources.append(res)

# 11. drug clusters
res = create_resource(snek.input[10], "drug_clusters", "Drug Clusters")
resources.append(res)

# 12. drug enrichment results
res = create_resource(snek.input[11], "drug_enrichment", "Drug Annotation Enrichment")
resources.append(res)

# 13. cell pca
res = create_resource(snek.input[12], "cell_pca", "Cell PCA")
resources.append(res)

# 14. cell similarity matrix

# 15. cell similarity matrix (umap projection)

# 16. cell clusters

# 17. cell average dose response curves

# 18. metadata

# x. figures?

# convert absolute resource paths to relative paths
base_dir = "/data/proj/hmcl/drug-screen-manuscript"

for resource in resources:
    resource.path = os.path.relpath(resource.path, base_dir)

# create data package
pkg = Package(
    resources=resources,
    name="hmcl_drug_screen",
    title="Human Myeloma Cell Line (HMCL) NCATS MIPE 4.0 Drug Screen Dataset"
)

pkg.to_yaml(snek.output[0])
