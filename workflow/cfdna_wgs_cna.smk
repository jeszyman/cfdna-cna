rule frag_filt:
    input:
        cfdna_wgs_cna_bam_inputs + "/{library}.bam",
    params:
        script = cfdna_wgs_scripts + "/frag_filt.sh",
    output:
        nohead =   temp(cfdna_wgs_cna_bam_fragfilt + "/{library}_frag{frag_distro}.nohead"),
        onlyhead = temp(cfdna_wgs_cna_bam_fragfilt + "/{library}_frag{frag_distro}.onlyhead"),
        final =         cfdna_wgs_cna_bam_fragfilt + "/{library}_frag{frag_distro}.bam",
    container:
        cfdna_wgs_container
    shell:
        """
        frag_min=$(echo {wildcards.frag_distro} | sed -e "s/_.*$//g")
        frag_max=$(echo {wildcards.frag_distro} | sed -e "s/^.*_//g")
        {params.script} \
        {input} \
        {output.nohead} \
        $frag_min \
        $frag_max \
        {config[threads]} \
        {output.onlyhead} \
        {output.final}
        """

rule bam_to_wig:
    input:
        cfdna_wgs_cna_bam_fragfilt + "/{library}_frag{frag_distro}.bam",
    output:
        wig + "/{library}_frag{frag_distro}.wig",
    params:
        chrs = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"
    container:
        cfdna_wgs_container,
    shell:
        """
        /opt/hmmcopy_utils/bin/readCounter --window 1000000 --quality 20 \
        --chromosome {params.chrs} \
        {input} > {output}
        """

# Make ichorCNA panel of normals from healthy samples
rule pon_list:
    input:
        expand(wig + "/{library}_frag{frag_distro}.wig", library = NORMAL_LIBRARIES, frag_distro = ["90_150"]),
    output:
        wig + "/normal.txt",
    log:
        cfdna_wgs_logs + "/pon.log",
    container:
        cfdna_wgs_container,
    shell:
        """
        input_string=$(echo "{input}" | tr " " "\n")
        if [ -f {output} ]; then rm {output}; fi
        echo -e "${{input_string}}" >> {output}
        """

# Make ichorCNA panel of normals from healthy samples
rule pon:
    input:
        wig + "/normal.txt",
    params:
        script = cfdna_wgs_scripts + "/pon.sh",
	outdir = wig
    output:
        wig + "/pon_median.rds"
    log:
        cfdna_wgs_logs + "/pon.log",
    container:
        cfdna_wgs_container,
    shell:
        """
        {params.script} \
        {input} \
        {params.outdir} &> {log}
        """

rule ichor:
    input:
        wig = wig + "/{library}_frag{frag_distro}.wig",
	pon = wig + "/pon_median.rds",
    output:
        ichor + "/{library}_frag{frag_distro}.cna.seg",
    params:
        script = cfdna_wgs_scripts + "/MOD_runIchorCNA.R",
        out_dir = ichor,
    container:
        cfdna_wgs_container,
    shell:
        """
        Rscript {params.script} \
         --id {wildcards.library}_frag{wildcards.frag_distro} \
         --WIG {input.wig} \
         --gcWig /opt/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
         --mapWig /opt/ichorCNA/inst/extdata/map_hg38_1000kb.wig \
         --centromere /opt/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
         --normal "c(0.95, 0.99, 0.995, 0.999)" \
         --normalPanel {input.pon} \
         --ploidy "c(2)" \
         --maxCN 3 \
         --estimateScPrevalence FALSE \
         --scStates "c()" \
         --outDir {params.out_dir}
        """

rule ichor_nopon:
    input:
        wig = wig + "/{library}_frag{frag_distro}.wig",
	pon = wig + "/pon_median.rds",
    output:
        ichor_nopon + "/{library}_frag{frag_distro}.cna.seg",
    params:
        script = cfdna_wgs_scripts + "/MOD_runIchorCNA.R",
        out_dir = ichor_nopon,
    container:
        cfdna_wgs_container,
    shell:
        """
        Rscript {params.script} \
         --id {wildcards.library}_frag{wildcards.frag_distro} \
         --WIG {input.wig} \
         --gcWig /opt/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
         --mapWig /opt/ichorCNA/inst/extdata/map_hg38_1000kb.wig \
         --centromere /opt/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
         --normal "c(0.95, 0.99, 0.995, 0.999)" \
         --ploidy "c(2)" \
         --maxCN 3 \
         --estimateScPrevalence FALSE \
         --scStates "c()" \
         --outDir {params.out_dir}
        """
