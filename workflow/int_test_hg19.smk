import pandas as pd

container: config["container"]

libraries = pd.read_table("test/inputs/hg19_libraries.tsv").set_index("library", drop = False)

LIBRARY_IDS = list(libraries.library.unique())
LIBRARY_IDS_NORMAL = list(libraries.library[(libraries["cohort"] == "healthy")].unique())
LIBRARY_IDS_DISEASE = list(libraries.library[(libraries["cohort"] == "mpnst")].unique())

rule all:
    input:
        expand(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.bam", library_id = LIBRARY_IDS, frag_distro = ["90_150"]),
        expand(config["wig_dir"] + "/{library_id}_frag{frag_distro}.wig", library_id = LIBRARY_IDS, frag_distro = ["90_150"]),
        #expand(config["ichor_hg19_dir"] + "/{library_id}_frag{frag_distro}.cna.seg", library_id = LIBRARY_IDS, frag_distro = ["90_150"]),

include: "cfdna_wgs_cna.smk"
