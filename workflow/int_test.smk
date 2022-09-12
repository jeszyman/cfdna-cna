#container: config["container"]

LIBRARIES = ["lib003",
	     "lib004",
	     "lib005",
	     "lib006"]

NORMAL_LIBRARIES = ["lib003", "lib004"]

cfdna_wgs_container = config["container"]["cfdna_wgs"]
cfdna_wgs_scripts = config["dir"]["cfdna_wgs_repo"] + "/workflow/scripts"

cfdna_wgs_cna_bam_inputs   = config["dir"]["data"] + "/bam/filt"
cfdna_wgs_cna_bam_fragfilt = config["dir"]["data"] + "/bam/frag"

wig = config["dir"]["data"] + "/wig"
ichor = config["dir"]["data"] + "/ichor"
cfdna_wgs_logs = config["dir"]["data"] + "logs/cfdna_wgs"
ichor_nopon = config["dir"]["data"] + "/ichor_nopon"

rule all:
    input:
	# Fragment-filtered bam
        #expand(cfdna_wgs_cna_bam_fragfilt + "/{library}_frag{frag_distro}.bam", library = LIBRARIES, frag_distro = ["90_150"]),
	# Wiggle
        #expand(wig + "/{library}_frag{frag_distro}.wig", library = LIBRARIES, frag_distro = ["90_150"]),
        # PoN list
        #wig + "/normal.txt",
        # PoN
        #wig + "/pon_median.rds",
	# Ichor
        expand(ichor + "/{library}_frag{frag_distro}.cna.seg", library = LIBRARIES, frag_distro = ["90_150"]),
        expand(ichor_nopon + "/{library}_frag{frag_distro}.cna.seg", library = LIBRARIES, frag_distro = ["90_150"]),

include: "cfdna_wgs_cna.smk"
