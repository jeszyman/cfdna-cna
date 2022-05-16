rule frag_filt:
    input:
        config["cna_bam_dir"] + "/{library_id}.bam"
    params:
        out_dir = config["frag_bam_dir"]
    output:
        nohead = temp(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.nohead"),
        onlyhead = temp(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.onlyhead"),
        final = config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.bam",
    shell:
        """
        frag_min=$(echo {wildcards.frag_distro} | sed -e "s/_.*$//g")
        frag_max=$(echo {wildcards.frag_distro} | sed -e "s/^.*_//g")
        workflow/scripts/frag_filt.sh {input} \
                                      {output.nohead} \
                                      $frag_min \
                                      $frag_max \
                                      {config[threads]} \
                                      {output.onlyhead} \
                                      {output.final}
        """

rule bam_to_wig:
    input: config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.bam",
    output: config["wig_dir"] + "/{library_id}_frag{frag_distro}.wig",
    params:
        chrs = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"
    shell:
        """
        /opt/hmmcopy_utils/bin/readCounter --window 1000000 --quality 20 \
        --chromosome {params.chrs} \
        {input} > {output}
        """

rule ichor:
    input:
        config["wig_dir"] + "/{library_id}_frag{frag_distro}.wig",
    output:
        config["ichor_dir"] + "/{library_id}_frag{frag_distro}.cna.seg",
    shell:
        """
        Rscript /opt/ichorCNA/scripts/runIchorCNA.R \
         --id {wildcards.library_id}_frag{wildcards.frag_distro} \
         --WIG {input} \
         --gcWig /opt/ichorCNA/inst/extdata/gc_hg19_1000kb.wig \
         --normal "c(0.95, 0.99)" \
         --ploidy "c(2)" \
         --maxCN 3 \
         --estimateScPrevalence FALSE \
         --scStates "c()" \
         --outDir {config[ichor_dir]}
        """
