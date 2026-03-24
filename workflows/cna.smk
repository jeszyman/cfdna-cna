# =============================================================================
# cna.smk — cfDNA Copy Number Analysis Pipeline Module
# =============================================================================
#
# Modular snakefile for cfDNA CNA analysis with ichorCNA and Fragle.
#
# Intended to be included via "include:" from a wrapper (e.g., test.smk).
# Expects the following to be defined before inclusion:
#   - resolve_config_paths() applied to config
#   - D_CNA, D_INPUTS, D_LOGS, D_BENCHMARK, ENV_ICHOR, ENV_FRAGLE
#   - LIBRARY_IDS, TUMOR_IDS, NORMAL_IDS, WINDOWS, PON_CONFIGS, PRESETS
#   - samples DataFrame
# =============================================================================

shell.prefix("set -e; export PYTHONNOUSERSITE=1; ")

# ---------------------------------------------------------------------------
# Helper: look up BAM path from sample table
# ---------------------------------------------------------------------------
def get_bam(wildcards):
    return samples.loc[samples.library_id == wildcards.library, "bam_path"].values[0]

# ---------------------------------------------------------------------------
# Rule 1a: Fragment filter (sized windows)
# ---------------------------------------------------------------------------
rule cna_fragment_filter:
    """Filter BAM by fragment length (TLEN) for a specific window."""
    input:
        bam=get_bam,
    output:
        bam=f"{D_CNA}/filt_bam/{{library}}.{{window}}.bam",
        bai=f"{D_CNA}/filt_bam/{{library}}.{{window}}.bam.bai",
    params:
        frag_min=lambda w: w.window.split("_")[0],
        frag_max=lambda w: w.window.split("_")[1],
    log:
        f"{D_LOGS}/{{library}}.{{window}}_fragment_filter.log"
    benchmark:
        f"{D_BENCHMARK}/{{library}}.{{window}}_fragment_filter.tsv"
    threads: 4
    conda: ENV_ICHOR
    wildcard_constraints:
        window=r"\d+_\d+"
    message: "Fragment filtering {wildcards.library} to {wildcards.window}"
    shell:
        """
        samtools view -h -@ {threads} \
          -e 'abs(tlen) >= {params.frag_min} && abs(tlen) <= {params.frag_max}' \
          {input.bam} | \
          samtools sort -@ {threads} -o {output.bam} 2> {log}
        samtools index {output.bam} 2>> {log}
        """

# ---------------------------------------------------------------------------
# Rule 1b: Fragment filter (full window — symlink)
# ---------------------------------------------------------------------------
rule cna_fragment_filter_full:
    """Symlink BAM for full (unfiltered) window."""
    input:
        bam=get_bam,
        bai=lambda w: get_bam(w) + ".bai",
    output:
        bam=f"{D_CNA}/filt_bam/{{library}}.full.bam",
        bai=f"{D_CNA}/filt_bam/{{library}}.full.bam.bai",
    log:
        f"{D_LOGS}/{{library}}.full_fragment_filter.log"
    benchmark:
        f"{D_BENCHMARK}/{{library}}.full_fragment_filter.tsv"
    message: "Symlinking full BAM for {wildcards.library}"
    shell:
        """
        ln -sf $(readlink -f {input.bam}) {output.bam} 2> {log}
        ln -sf $(readlink -f {input.bai}) {output.bai} 2>> {log}
        """

ruleorder: cna_fragment_filter_full > cna_fragment_filter

# ---------------------------------------------------------------------------
# Rule 2: BAM to WIG (readCounter)
# ---------------------------------------------------------------------------
rule cna_bam_to_wig:
    """Convert filtered BAM to WIG using readCounter (HMMcopy)."""
    input:
        bam=f"{D_CNA}/filt_bam/{{library}}.{{window}}.bam",
        bai=f"{D_CNA}/filt_bam/{{library}}.{{window}}.bam.bai",
    output:
        wig=f"{D_CNA}/wig/{{library}}.{{window}}.wig",
    params:
        chromosomes=config["readcounter-chromosomes"],
    log:
        f"{D_LOGS}/{{library}}.{{window}}_bam_to_wig.log"
    benchmark:
        f"{D_BENCHMARK}/{{library}}.{{window}}_bam_to_wig.tsv"
    threads: 1
    conda: ENV_ICHOR
    message: "readCounter: {wildcards.library}.{wildcards.window}"
    shell:
        """
        readCounter \
          --window 1000000 \
          --quality 20 \
          --chromosome "{params.chromosomes}" \
          {input.bam} > {output.wig} 2> {log}
        """

# ---------------------------------------------------------------------------
# Rule 3: PoN WIG list
# ---------------------------------------------------------------------------
rule cna_pon_wig_list:
    """Create text file listing normal-sample WIG paths for PoN construction."""
    input:
        wigs=expand(f"{D_CNA}/wig/{{normal}}.{{{{window}}}}.wig", normal=NORMAL_IDS),
    output:
        f"{D_CNA}/pon/{{pon}}.{{window}}_wiglist.txt",
    log:
        f"{D_LOGS}/{{pon}}.{{window}}_pon_wig_list.log"
    benchmark:
        f"{D_BENCHMARK}/{{pon}}.{{window}}_pon_wig_list.tsv"
    wildcard_constraints:
        pon="|".join([k for k, v in PON_CONFIGS.items() if v.get("normal-samples")])
    message: "Building WIG list for PoN {wildcards.pon}.{wildcards.window}"
    run:
        with open(output[0], "w") as f:
            for wig in input.wigs:
                f.write(wig + "\n")

# ---------------------------------------------------------------------------
# Rule 4: Build PoN
# ---------------------------------------------------------------------------
rule cna_build_pon:
    """Build panel of normals from normal sample WIGs."""
    input:
        wiglist=f"{D_CNA}/pon/{{pon}}.{{window}}_wiglist.txt",
    output:
        rds=f"{D_CNA}/pon/{{pon}}.{{window}}_median.rds",
    params:
        ichor_repo=config["ichor-repo"],
        gc_wig=config["ichor-ref"]["gc-wig"],
        map_wig=config["ichor-ref"]["map-wig"],
        centromere=config["ichor-ref"]["centromere"],
        chrs=config["ichor-chromosomes"],
        chrs_norm=config["ichor-chrs-normalize"],
    log:
        f"{D_LOGS}/{{pon}}.{{window}}_build_pon.log"
    benchmark:
        f"{D_BENCHMARK}/{{pon}}.{{window}}_build_pon.tsv"
    threads: 1
    conda: ENV_ICHOR
    message: "Building PoN {wildcards.pon} for window {wildcards.window}"
    shell:
        """
        Rscript {params.ichor_repo}/scripts/createPanelOfNormals.R \
          --filelist {input.wiglist} \
          --gcWig {params.gc_wig} \
          --mapWig {params.map_wig} \
          --centromere {params.centromere} \
          --chrs "{params.chrs}" \
          --chrNormalize "{params.chrs_norm}" \
          --outfile {output.rds} \
          --libdir {params.ichor_repo} \
          2> {log}
        """

# ---------------------------------------------------------------------------
# PoN input resolution
# ---------------------------------------------------------------------------
def get_pon_input(wildcards):
    """Return PoN RDS path as a list (empty for nopon)."""
    pon_cfg = PON_CONFIGS[wildcards.pon]
    if pon_cfg.get("prebuilt-path"):
        return [pon_cfg["prebuilt-path"]]
    elif pon_cfg.get("normal-samples"):
        return [f"{D_CNA}/pon/{wildcards.pon}.{wildcards.window}_median.rds"]
    else:  # nopon
        return []

# ---------------------------------------------------------------------------
# Rule 5: Run ichorCNA
# ---------------------------------------------------------------------------
rule cna_run_ichor:
    """Run ichorCNA with preset parameters and optional PoN."""
    input:
        wig=f"{D_CNA}/wig/{{library}}.{{window}}.wig",
        pon=get_pon_input,
    output:
        seg=f"{D_CNA}/{{library}}.{{window}}.{{pon}}.{{preset}}/{{library}}.cna.seg",
        params_txt=f"{D_CNA}/{{library}}.{{window}}.{{pon}}.{{preset}}/{{library}}.params.txt",
        corrected=f"{D_CNA}/{{library}}.{{window}}.{{pon}}.{{preset}}/{{library}}.correctedDepth.txt",
    params:
        ichor_repo=config["ichor-repo"],
        gc_wig=config["ichor-ref"]["gc-wig"],
        map_wig=config["ichor-ref"]["map-wig"],
        centromere=config["ichor-ref"]["centromere"],
        genome_build=config["genome-build"],
        genome_style=config["genome-style"],
        chrs=config["ichor-chromosomes"],
        chrs_norm=config["ichor-chrs-normalize"],
        outdir=lambda w: f"{D_CNA}/{w.library}.{w.window}.{w.pon}.{w.preset}",
        pon_flag=lambda w, input: f"--normalPanel {input.pon[0]}" if input.pon else "",
        normal=lambda w: PRESETS[w.preset]["normal"],
        ploidy=lambda w: PRESETS[w.preset]["ploidy"],
        maxCN=lambda w: PRESETS[w.preset]["maxCN"],
        includeHOMD=lambda w: PRESETS[w.preset]["includeHOMD"],
        txnE=lambda w: PRESETS[w.preset]["txnE"],
        txnStrength=lambda w: PRESETS[w.preset]["txnStrength"],
        estimateNormal=lambda w: PRESETS[w.preset]["estimateNormal"],
        estimatePloidy=lambda w: PRESETS[w.preset]["estimatePloidy"],
        estimateScPrevalence=lambda w: PRESETS[w.preset]["estimateScPrevalence"],
        scStates=lambda w: PRESETS[w.preset]["scStates"],
    log:
        f"{D_LOGS}/{{library}}.{{window}}.{{pon}}.{{preset}}_ichor.log"
    benchmark:
        f"{D_BENCHMARK}/{{library}}.{{window}}.{{pon}}.{{preset}}_ichor.tsv"
    threads: 1
    conda: ENV_ICHOR
    message: "ichorCNA: {wildcards.library}.{wildcards.window}.{wildcards.pon}.{wildcards.preset}"
    shell:
        """
        Rscript {params.ichor_repo}/scripts/runIchorCNA.R \
          --libdir {params.ichor_repo} \
          --id {wildcards.library} \
          --WIG {input.wig} \
          --gcWig {params.gc_wig} \
          --mapWig {params.map_wig} \
          --centromere {params.centromere} \
          --genomeBuild {params.genome_build} \
          --genomeStyle {params.genome_style} \
          --chrs "{params.chrs}" \
          --chrNormalize "{params.chrs_norm}" \
          --normal "{params.normal}" \
          --ploidy "{params.ploidy}" \
          --maxCN {params.maxCN} \
          --includeHOMD {params.includeHOMD} \
          --txnE {params.txnE} \
          --txnStrength {params.txnStrength} \
          --estimateNormal {params.estimateNormal} \
          --estimatePloidy {params.estimatePloidy} \
          --estimateScPrevalence {params.estimateScPrevalence} \
          --scStates "{params.scStates}" \
          {params.pon_flag} \
          --outDir {params.outdir}/ \
          2> {log}
        """

# ---------------------------------------------------------------------------
# Rule 6: Run Fragle
# ---------------------------------------------------------------------------
rule cna_run_fragle:
    """Run Fragle ctDNA quantification on filtered BAM."""
    input:
        bam=f"{D_CNA}/filt_bam/{{library}}.{{window}}.bam",
        bai=f"{D_CNA}/filt_bam/{{library}}.{{window}}.bam.bai",
    output:
        csv=f"{D_CNA}/fragle/{{library}}.{{window}}/Fragle.csv",
    params:
        fragle_repo=config["fragle"]["repo-path"],
        mode=config["fragle"]["mode"],
        genome_build=config["genome-build"],
        outdir=lambda w: f"{D_CNA}/fragle/{w.library}.{w.window}",
    log:
        f"{D_LOGS}/{{library}}.{{window}}_fragle.log"
    benchmark:
        f"{D_BENCHMARK}/{{library}}.{{window}}_fragle.tsv"
    threads: config["fragle"]["threads"]
    conda: ENV_FRAGLE
    message: "Fragle: {wildcards.library}.{wildcards.window}"
    shell:
        """
        python {params.fragle_repo}/main.py \
          --input {input.bam} \
          --output {params.outdir} \
          --mode {params.mode} \
          --genome_build {params.genome_build} \
          --cpu {threads} \
          --threads {threads} \
          2> {log}
        """

# ---------------------------------------------------------------------------
# Rule 7: Extract ichorCNA tumor fraction
# ---------------------------------------------------------------------------
rule cna_extract_ichor_tf:
    """Parse ichorCNA params.txt to extract tumor fraction."""
    input:
        params_txt=f"{D_CNA}/{{library}}.{{window}}.{{pon}}.{{preset}}/{{library}}.params.txt",
    output:
        tsv=f"{D_CNA}/{{library}}.{{window}}.{{pon}}.{{preset}}/{{library}}.tf.tsv",
    log:
        f"{D_LOGS}/{{library}}.{{window}}.{{pon}}.{{preset}}_extract_tf.log"
    benchmark:
        f"{D_BENCHMARK}/{{library}}.{{window}}.{{pon}}.{{preset}}_extract_tf.tsv"
    conda: ENV_ICHOR
    message: "Extract TF: {wildcards.library}.{wildcards.window}.{wildcards.pon}.{wildcards.preset}"
    script:
        "../scripts/extract_ichor_tf.py"

# ---------------------------------------------------------------------------
# Rule 8: Aggregate all TF results
# ---------------------------------------------------------------------------
rule cna_aggregate_results:
    """Merge ichorCNA and Fragle TF estimates into comparison table."""
    input:
        ichor_tfs=expand(
            f"{D_CNA}/{{lib}}.{{win}}.{{pon}}.{{preset}}/{{lib}}.tf.tsv",
            lib=TUMOR_IDS, win=WINDOWS,
            pon=PON_CONFIGS.keys(), preset=PRESETS.keys()
        ),
        fragle_csvs=expand(
            f"{D_CNA}/fragle/{{lib}}.{{win}}/Fragle.csv",
            lib=TUMOR_IDS, win=WINDOWS
        ),
    output:
        tsv=f"{D_CNA}/summary/tf_comparison.tsv",
    log:
        f"{D_LOGS}/aggregate_results.log"
    benchmark:
        f"{D_BENCHMARK}/aggregate_results.tsv"
    conda: ENV_ICHOR
    message: "Aggregating TF comparison table"
    script:
        "../scripts/aggregate_tf.R"

# ---------------------------------------------------------------------------
# Rule 9a: TF comparison plot
# ---------------------------------------------------------------------------
rule cna_plot_tf_comparison:
    input:
        tsv=f"{D_CNA}/summary/tf_comparison.tsv",
    output:
        pdf=f"{D_CNA}/plots/tf_comparison.pdf",
        png=f"{D_CNA}/plots/tf_comparison.png",
    log:
        f"{D_LOGS}/plot_tf_comparison.log"
    benchmark:
        f"{D_BENCHMARK}/plot_tf_comparison.tsv"
    conda: ENV_ICHOR
    message: "Plotting TF comparison"
    script:
        "../scripts/plot_tf_comparison.R"

# ---------------------------------------------------------------------------
# Rule 9b: LogR heatmap (one per pon/preset combo)
# ---------------------------------------------------------------------------
rule cna_plot_logr_heatmap:
    input:
        depth_files=expand(
            f"{D_CNA}/{{lib}}.{{win}}.{{{{pon}}}}.{{{{preset}}}}/{{lib}}.correctedDepth.txt",
            lib=TUMOR_IDS, win=WINDOWS,
        ),
    output:
        pdf=f"{D_CNA}/plots/logr_heatmap.{{pon}}.{{preset}}.pdf",
        png=f"{D_CNA}/plots/logr_heatmap.{{pon}}.{{preset}}.png",
    log:
        f"{D_LOGS}/plot_logr_heatmap.{{pon}}.{{preset}}.log"
    benchmark:
        f"{D_BENCHMARK}/plot_logr_heatmap.{{pon}}.{{preset}}.tsv"
    conda: ENV_ICHOR
    message: "Plotting logR heatmap for {wildcards.pon}.{wildcards.preset}"
    script:
        "../scripts/plot_logr_heatmap.R"

# ---------------------------------------------------------------------------
# Rule 9c: Window concordance
# ---------------------------------------------------------------------------
rule cna_plot_window_concordance:
    input:
        tsv=f"{D_CNA}/summary/tf_comparison.tsv",
    output:
        pdf=f"{D_CNA}/plots/window_concordance.pdf",
        png=f"{D_CNA}/plots/window_concordance.png",
    log:
        f"{D_LOGS}/plot_window_concordance.log"
    benchmark:
        f"{D_BENCHMARK}/plot_window_concordance.tsv"
    conda: ENV_ICHOR
    message: "Plotting window concordance"
    script:
        "../scripts/plot_window_concordance.R"

# ---------------------------------------------------------------------------
# Rule 9d: Fragle summary
# ---------------------------------------------------------------------------
rule cna_plot_fragle_summary:
    input:
        fragle_csvs=expand(
            f"{D_CNA}/fragle/{{lib}}.{{win}}/Fragle.csv",
            lib=TUMOR_IDS, win=WINDOWS
        ),
    output:
        pdf=f"{D_CNA}/plots/fragle_summary.pdf",
        png=f"{D_CNA}/plots/fragle_summary.png",
    log:
        f"{D_LOGS}/plot_fragle_summary.log"
    benchmark:
        f"{D_BENCHMARK}/plot_fragle_summary.tsv"
    conda: ENV_ICHOR
    message: "Plotting Fragle summary"
    script:
        "../scripts/plot_fragle_summary.R"
