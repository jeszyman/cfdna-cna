import pandas as pd

container: config["container"]

hg38_inputs = pd.read_table("test/inputs/hg38_libraries.tsv").set_index("library", drop = False)

HG38_LIBS = list(hg38_inputs.library.unique())
HG38_LIBS_NORMAL = list(hg38_inputs.library[(hg38_inputs["cohort"] == "healthy")].unique())
HG38_LIBS_NORMAL = list(hg38_inputs.library[(hg38_inputs["cohort"] == "mpnst")].unique())

rule all:
    input:
        expand(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.bam", library_id = HG38_LIBS, frag_distro = ["90_150"]),
        expand(config["wig_dir"] + "/{library_id}_frag{frag_distro}.wig", library_id = HG38_LIBS, frag_distro = ["90_150"]),
        expand(config["ichor_hg38_dir"] + "/{library_id}_frag{frag_distro}.cna.seg", library_id = HG38_LIBS, frag_distro = ["90_150"]),

include: "cfdna_wgs_cna.smk"
