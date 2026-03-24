# cfdna-cna

Modular bioinformatics pipeline for cfdna-cna processing.

## Quick reference

- **Conda env:** `cfdna-cna`
- **Dry-run:** `conda run -n basecamp snakemake -s workflows/test.smk --configfile config/test.yaml --dry-run`
- **Full run:** `conda run -n basecamp snakemake -s workflows/test.smk --configfile config/test.yaml --use-conda --cores 4`
- **Test data:** `tests/full/inputs/`
- **Org source:** `cfdna-cna.org` (edit here, tangle to generate workflow files)

## Development workflow

1. Edit `cfdna-cna.org` (all code lives in org-babel blocks)
2. Tangle to generate workflow/config/script files
3. Dry-run: verify DAG resolves
4. Run on test data
5. Commit
