container: config["container"]

IDS = ["lib001_hg19", "lib002_hg19"]

rule all:
    input:
        expand(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.bam", library_id = IDS, frag_distro = ["90_150"]),
        expand(config["wig_dir"] + "/{library_id}_frag{frag_distro}.wig", library_id = IDS, frag_distro = ["90_150"]),
        expand(config["ichor_dir"] + "/{library_id}_frag{frag_distro}.cna.seg", library_id = ["lib001_hg19", "lib002_hg19"], frag_distro = ["90_150"]),

include: "cfdna_wgs_cna.smk"
