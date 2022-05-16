container: config["container"]

HG19_IDS = ["lib001_hg19", "lib002_hg19"]
HG38_IDS = ["lib003_hg38", "lib004_hg38"]

rule all:
    input:
        expand(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.bam", library_id = HG19_IDS, frag_distro = ["90_150"]),
        expand(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.bam", library_id = HG38_IDS, frag_distro = ["90_150"]),	
        expand(config["wig_dir"] + "/{library_id}_frag{frag_distro}.wig", library_id = HG19_IDS, frag_distro = ["90_150"]),
        expand(config["wig_dir"] + "/{library_id}_frag{frag_distro}.wig", library_id = HG38_IDS, frag_distro = ["90_150"]),	
        expand(config["ichor_hg19_dir"] + "/{library_id}_frag{frag_distro}.cna.seg", library_id = HG19_IDS, frag_distro = ["90_150"]),
        expand(config["ichor_hg38_dir"] + "/{library_id}_frag{frag_distro}.cna.seg", library_id = HG38_IDS, frag_distro = ["90_150"]),	

include: "cfdna_wgs_cna.smk"
