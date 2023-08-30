"""
HMCL NCATS MIPE 4.0 Drug Screen Dataset Pipeline
"""
import pathlib

# output directory
out_dir = pathlib.Path("output/")
out_dir.mkdir(mode=0o755, parents=True, exist_ok=True)

rule create_drug_plate_matrix:
    input:
        "data/raw.tsv.gz"
    output:
        out_dir.joinpath("drug_plates/raw.tsv.gz"),
        out_dir.joinpath("metadata/plate-metadata.tsv"),
        out_dir.joinpath("metadata/drug-groups.tsv"),
    script:
        "scripts/create_drug_plate_matrix.R"
