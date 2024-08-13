"""
Creates a data package with the main result datasets + metadata

Schemas are excluded for matrix results since they don't provide a lot of useful information and
would clutter to the resulting datapackage.yml.
"""
import os
from frictionless import describe, Package

snek = snakemake

# Create dataset Resource objects
resources = []

# 1. drug plates
res = describe(snek.input[0], stats=True)
res.schema = None
res.name = "plates_raw"
res.title = "Drug plate matrix (raw)"
resources.append(res)

res = describe(snek.input[1], stats=True)
res.schema = None
res.name = "plates_normed"
res.title = "Drug plate matrix (normed)"
resources.append(res)

res = describe(snek.input[2], stats=True)
res.schema = None
res.name = "plates_bgadj"
res.title = "Drug plate matrix (background-adjusted)"
resources.append(res)

res = describe(snek.input[3], stats=True)
res.schema = None
res.name = "plates_background"
res.title = "Drug plate matrix (background)"
resources.append(res)

res = describe(snek.input[4], stats=True)
res.schema = None
res.name = "plates_concentrations"
res.title = "Drug plate matrix (concentrations)"
resources.append(res)

# 2. drug curves

# 3. drug response matrices

# 4. drug pca/umap data

# 5. drug similarity matrix

# 6. drug clusters

# 7. drug enrichment results

# 8. cell pca/umap data

# 9. cell similarity matrix

# 10. cell clusters

# 11. cell average dose response curves

# 12. metadata

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
