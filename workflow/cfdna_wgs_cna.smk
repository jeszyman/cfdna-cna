rule frag_filt:
    input:
        config["data_dir"] + "/bam/{library_id}.bam"
    params:
        out_dir = config["data_dir"] + "/frag_bam"
    output:
        nohead = temp(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.nohead"),
        onlyhead = temp(config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.onlyhead"),
        final = config["frag_bam_dir"] + "/{library_id}_frag{frag_distro}.bam",
    shell:
        """
        frag_min=$(echo {wildcards.frag_distro} | sed -e "s/_.*$//g")
        frag_max=$(echo {wildcards.frag_distro} | sed -e "s/^.*_//g")
        {config[cfdna_cna_repo]}/workflow/scripts/frag_filt.sh {input} \
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

rule pon_txt:
    input:  expand(config["data_dir"] + "/wig/{library_id}_frag{frag_distro}.wig", library_id = H),
    output: config["data_dir"] + "/wig/normals.txt"
    shell:
        """
        {config[cfdna_cna_script_dir]}/pon_txt.sh \
        {input} {output}
        """
