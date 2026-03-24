# =============================================================================
# test.smk — Test wrapper for cfdna-cna pipeline
# =============================================================================
#
# Includes cna.smk and defines rule all for integration testing.
# Uses chr22-scoped test data.
# =============================================================================

import os
import pandas as pd

# ---------------------------------------------------------------------------
# Load YAML Configuration
# ---------------------------------------------------------------------------

configfile: "config/test.yaml"

def resolve_config_paths(config_dict):
    """Recursively expand ${HOME} and other env vars in config values."""
    for k, v in config_dict.items():
        if isinstance(v, str):
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)
        elif isinstance(v, list):
            config_dict[k] = [
                os.path.expandvars(os.path.expanduser(i)) if isinstance(i, str) else i
                for i in v
            ]

resolve_config_paths(config)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

D_CNA       = config["directories"]["cfdna-cna"]
D_INPUTS    = config["directories"]["inputs"]
D_LOGS      = config["directories"]["logs"]
D_BENCHMARK = config["directories"]["benchmark"]
ENV_ICHOR   = config["conda"]["ichor"]
ENV_FRAGLE  = config["conda"]["fragle"]

# ---------------------------------------------------------------------------
# Load Tabular Configuration
# ---------------------------------------------------------------------------

samples = pd.read_csv(config["sample-tsv-path"], sep="\t")
samples["bam_path"] = samples["bam_path"].apply(os.path.expandvars)
LIBRARY_IDS = samples["library_id"].tolist()
TUMOR_IDS = samples[samples["sample_type"] == "tumor"]["library_id"].tolist()
NORMAL_IDS = samples[samples["sample_type"] == "normal"]["library_id"].tolist()

# ---------------------------------------------------------------------------
# Window / PoN / Preset lists
# ---------------------------------------------------------------------------

WINDOWS = list(dict.fromkeys(config.get("windows", []) + ["full"]))
PON_CONFIGS = config["pon"]
PRESETS = config["ichor-presets"]

# ---------------------------------------------------------------------------
# Include module
# ---------------------------------------------------------------------------

include: "cna.smk"

# ---------------------------------------------------------------------------
# Rule all
# ---------------------------------------------------------------------------

rule all:
    input:
        # ichorCNA seg files
        expand(
            "{d}/{lib}.{win}.{pon}.{preset}/{lib}.cna.seg",
            d=D_CNA, lib=TUMOR_IDS, win=WINDOWS,
            pon=PON_CONFIGS.keys(), preset=PRESETS.keys()
        ),
        # Fragle CSVs
        expand(
            "{d}/fragle/{lib}.{win}/Fragle.csv",
            d=D_CNA, lib=TUMOR_IDS, win=WINDOWS
        ),
        # Aggregate TF table
        f"{D_CNA}/summary/tf_comparison.tsv",
        # Visualization outputs
        f"{D_CNA}/plots/tf_comparison.pdf",
        expand(
            "{d}/plots/logr_heatmap.{pon}.{preset}.pdf",
            d=D_CNA, pon=PON_CONFIGS.keys(), preset=PRESETS.keys()
        ),
        f"{D_CNA}/plots/window_concordance.pdf",
        f"{D_CNA}/plots/fragle_summary.pdf",
