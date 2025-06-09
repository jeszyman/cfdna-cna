# Will follow symlinks
rule ichor_index_bam_check:
    conda:
        "../config/cfdna-cna-conda-env.yaml"
    input:
        bam = ancient(f"{data_dir}/analysis/cfdna_cna/bams/{{library_id}}.bam"),
    output:
        bai = f"{data_dir}/analysis/cfdna_cna/bams/{{library_id}}.bam.bai",
    shell:
        """
        samtools index -@ 8 {input.bam} {output.bai}
        """

rule make_wig:
    conda:
        "../config/cfdna-cna-conda-env.yaml"
    input:
        bam = f"{data_dir}/analysis/cfdna_cna/bams/{{library_id}}.bam",
        bai = f"{data_dir}/analysis/cfdna_cna/bams/{{library_id}}.bam.bai",
    output:
        wig = f"{data_dir}/analysis/cfdna_cna/wigs/{{library_id}}.wig",
    params:
        window = "1000000",
        quality = 20,
        ichor_wig_dir = f"{data_dir}/analysis/cfdna_cna/wigs",
    shell:
        """
        mkdir -p "{params.ichor_wig_dir}"
        readCounter \
        --window {params.window} \
        --quality {params.quality} \
	--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
        {input} > {output}
        """

rule run_ichor:
    conda:
        "../config/cfdna-cna-conda-env.yaml"
    input:
        wig = f"{data_dir}/analysis/cfdna_cna/wigs/{{library_id}}.wig"
    output:
        f"{data_dir}/analysis/cfdna_cna/ichor/{{library_id}}/{{library_id}}.cna.seg",
    params:
        ichor_out_main_dir = f"{data_dir}/analysis/cfdna_cna/ichor",
        ichor_repo = ichor_repo
    shell:
        """
        mkdir -p $(dirname {output})
        Rscript {params.ichor_repo}/scripts/runIchorCNA.R \
            --id {wildcards.library_id} \
            --WIG {input.wig} \
            --normal "c(0.95, 0.99, 0.995, 0.999)" \
            --genomeBuild hg38 \
            --ploidy "c(2)" \
            --gcWig {params.ichor_repo}/inst/extdata/gc_hg38_1000kb.wig \
            --mapWig {params.ichor_repo}/inst/extdata/map_hg38_1000kb.wig \
            --centromere {params.ichor_repo}/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
            --normalPanel {params.ichor_repo}/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
            --includeHOMD FALSE \
            --chrs "c(1:22)" \
            --chrTrain "c(1:22)" \
            --estimateNormal TRUE \
            --estimatePloidy TRUE \
            --estimateScPrevalence TRUE \
            --scStates "c()" \
            --txnE 0.9999 \
            --txnStrength 10000 \
            --outDir {params.ichor_out_main_dir}/{wildcards.library_id} \
            --libdir {params.ichor_repo}
        """

rule extract_tumor_fractions:
    input:
        expand(f"{data_dir}/analysis/cfdna_cna/ichor/{{library_id}}/{{library_id}}.params.txt", library_id=cfdna_cna_library_ids)
    output:
        f"{data_dir}/analysis/cfdna_cna/ichor/ichor_tumor_fractions.tsv"
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
