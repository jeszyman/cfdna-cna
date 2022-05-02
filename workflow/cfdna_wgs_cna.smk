rule frag_filt:
    input:
        config["bam_dir"] + "/{library_id}.bam"
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
