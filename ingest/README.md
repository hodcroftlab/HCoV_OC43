# Ingest Directory Overview
- Snakefile: Snakemake workflow file to run the ingest process.
- bin/: Directory containing scripts used in the ingest workflow.
- generate_from_genbank.py: Script to download and parse GenBank files into required formats.
- config/: Configuration files for the ingest process.
- data/: Directory containing input data for the ingest process.
- source-data/: Directory for annotations and geo-location rules.
- vendored/: Directory for vendored scripts utilized in the ingest.
- workflow/: Directory containing rules for the ingest workflow.


## Prepare Reference Files
Nextclade is used to align sequences to the reference and assign clades. You’ll need to provide reference files in the correct format (e.g. reference.fasta, annotation.gff3).

To generate these files from the GenBank record:

## Verify Config Settings
Open ingest/config/config.yaml and confirm that ncbi_taxon_id is set correctly.

## Run the generate_from_genbank.py script:

```bash
python3 bin/generate_from_genbank.py --reference "AY391777" --output-dir data/references/
```

## You may be prompted to specify CDS annotations.

Enter [product] or [0] to auto-select proteins, or leave blank for manual selection.

Note: For some pathogens, there may be no shared fields. Leave input blank and manually assign CDS names if necessary.

## Update Attributes
Ensure pathogen.json is up to date. You can consult the [Nextclade pathogen config documentation](https://docs.nextstrain.org/projects/nextclade/en/stable/user/input-files/05-pathogen-config.html) for additional optional attributes.

## Run the Ingest Workflow
Make scripts executable if needed:

```bash
chmod +x ./vendored/* ./bin/*
``` 

Then run the ingest pipeline from within the ingest directory:

```bash
snakemake --cores 9 all
```
This will download OC43 sequences and metadata, generate data/sequences.fasta and data/metadata.tsv, and prepare them for downstream analysis.

## Updating Vendored Scripts

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage vendored scripts in `vendored`.

### Steps to Update Vendored Scripts

1. Install `git subrepo` by following the [installation guide](https://github.com/ingydotnet/git-subrepo#installation).
2. Pull the latest changes from the central ingest repository by following the instructions in [`vendored/README.md`](vendored/README.md#vendoring).
