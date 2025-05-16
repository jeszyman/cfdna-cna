# import pandas as pd
# import os
# from tabulate import tabulate
# import sys


# ichor_bam_dir = config["ichor-bam-dir"]
# ichor_wig_dir = config["ichor-wig-dir"]
# ichor_out_main_dir = config["ichor-main-out-dir"]

# library_id = "NH_39_L1"


# rule all:
#     input:
#         expand(f"{ichor_out_main_dir}/{{library_id}}/{{library_id}}.cna.seg",
#                library_id=["NH_39_L1"]),


# Will follow symlinks
rule ichor_index_bam_check:
    input:
        bam = ancient(f"{ichor_bam_dir}/{{library_id}}.bam"),
    output:
        bai = f"{ichor_bam_dir}/{{library_id}}.bam.bai",
    shell:
        """
        samtools index -@ 8 {input.bam} {output.bai}
        """

rule make_wig:
    input:
        bam = f"{ichor_bam_dir}/{{library_id}}.bam",
        bai = f"{ichor_bam_dir}/{{library_id}}.bam.bai",
    output:
        wig = f"{ichor_wig_dir}/{{library_id}}.wig",
    params:
        window = "1000000",
        quality = 20,
    shell:
        """
        mkdir -p "{ichor_wig_dir}"
        readCounter \
        --window {params.window} \
        --quality {params.quality} \
	--chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \
        {input} > {output}
        """

rule run_ichor:
    input:
        wig = f"{ichor_wig_dir}/{{library_id}}.wig"
    output:
        f"{ichor_out_main_dir}/{{library_id}}/{{library_id}}.cna.seg",
        f"{ichor_out_main_dir}/{{library_id}}/{{library_id}}.params.txt",
    params:
        ichor_out_main_dir = ichor_out_main_dir,
        ichor_repo = ichor_repo,
    shell:
        """
        Rscript {params.ichor_repo}/scripts/runIchorCNA.R \
        --id {wildcards.library_id} \
        --WIG {input.wig} \
        --normal "c(0.95, 0.99, 0.995, 0.999)" \
        --genomeBuild "hg38" \
        --ploidy "c(2)" \
        --gcWig {params.ichor_repo}/inst/extdata/gc_hg38_1000kb.wig \
        --mapWig {params.ichor_repo}/inst/extdata/map_hg38_1000kb.wig \
        --centromere {params.ichor_repo}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
        --normalPanel {params.ichor_repo}/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
        --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" \
        --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
        --scStates "c()" --txnE 0.9999 --txnStrength 10000 --outDir {params.ichor_out_main_dir}/{wildcards.library_id}
         """

rule extract_tumor_fractions:
    input:
        expand(f"{ichor_out_main_dir}/{{library_id}}/{{library_id}}.params.txt", library_id=emseq_library_ids)
    output:
        f"{ichor_out_main_dir}/ichor_tumor_fractions.tsv"
    run:
        with open(output[0], "w") as out:
            out.write("library\ttf\n")
            for f in input:
                sample = f.split("/")[-1].replace(".params.txt", "")
                with open(f) as fh:
                    lines = fh.readlines()
                    if len(lines) >= 2:
                        tf = lines[1].split()[1]
                        out.write(f"{sample}\t{tf}\n")
