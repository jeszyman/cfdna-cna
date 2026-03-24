# cfdna-cna

cfDNA copy number analysis biopipe using ichorCNA and Fragle.

## Quick Start

```bash
# Dry-run
conda run -n basecamp snakemake -s workflows/test.smk --use-conda --conda-frontend conda --dry-run

# Execute
conda run -n basecamp snakemake -s workflows/test.smk --use-conda --conda-frontend conda -j4
```

## Architecture

- Input: quality-filtered BAMs (not fragment-filtered)
- Processing: {library} x {window} x {pon} x {preset} matrix
- `full` window always included
- Two conda envs: ichor (R) and fragle (Python)
- Conditional PoN: prebuilt RDS or auto-generated from normal samples

## Key Conventions

- Always use `--conda-frontend conda` (mamba prefix bug)
- ichorCNA invoked via `--libdir` from `~/repos/ichorCNA-patched`
- Fragle cloned into `resources/FRAGLE/`
- Test data is chr22-only (DAG testing, not biological validation)
- Config keys use hyphens, not underscores

## Development workflow

1. Edit `cfdna-cna.org` (all code lives in org-babel blocks)
2. Tangle to generate workflow/config/script files
3. Dry-run: verify DAG resolves
4. Run on test data
5. Commit
