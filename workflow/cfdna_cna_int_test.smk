container: config["container"]

IDS, = glob_wildcards(config["cna_bam_dir"] + "/{id}.bam")

rule all:
    input:
        expand(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.bam", library_id = IDS, frag_distro = ["90_150"]),

include: "cfdna_wgs_cna.smk"
