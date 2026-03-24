"""Parse ichorCNA .params.txt and extract tumor fraction + ploidy."""
import sys

params_file = snakemake.input.params_txt
output_file = snakemake.output.tsv
library = snakemake.wildcards.library
window = snakemake.wildcards.window
pon = snakemake.wildcards.pon
preset = snakemake.wildcards.preset

tumor_fraction = "NA"
ploidy = "NA"

with open(params_file) as f:
    for line in f:
        if line.startswith("Tumor Fraction:"):
            tumor_fraction = line.strip().split("\t")[1]
        elif line.startswith("Ploidy:"):
            ploidy = line.strip().split("\t")[1]

with open(output_file, "w") as f:
    f.write("library_id\twindow\tpon\tpreset\tmethod\ttumor_fraction\tploidy\n")
    f.write(f"{library}\t{window}\t{pon}\t{preset}\tichorCNA\t{tumor_fraction}\t{ploidy}\n")
