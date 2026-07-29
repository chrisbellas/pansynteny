# Pangenome explorer

A local, interactive gene-order viewer for a Panaroo pangenome. Search any
gene cluster by name or annotation, then build an anchor-gene gene-order
chart on demand: the anchor +/- a flanking window, sampled across carrier
genomes, RFE-important flanking genes colored, per-row gutter markers
showing whether each colored gene is visible in the window / present
elsewhere / absent, a genome-count toggle, and client-side metadata filter
boxes with a cohort-wide "% of \<value\> genomes carry this gene" stat.

Dataset-agnostic: a genome is viewable if it has a column in Panaroo's
`gene_presence_absence.csv` and a matching Prokka GFF -- nothing else is
required. No geNomad dependency.

## Requirements

Python 3, standard library only (no pip installs).

## Setup

Every external data path is an absolute path, given either as a CLI flag
or in a config file -- nothing is inferred from where you put this folder.

1. Copy `pangenome_viewer.config.example` to `pangenome_viewer.config` in
   this same directory (gitignored, since it's deployment-specific) and
   fill in absolute paths for your dataset. Only `panaroo_csv` and
   `prokka_dir` are required; every other line is optional and the
   corresponding feature just turns itself off if left blank -- see the
   comments in the example file and the docstring at the top of
   `pangenome_viewer.py` for exactly what each one does.
2. Run it:
   ```
   python3 pangenome_viewer.py
   ```
   (or override/skip the config file entirely with CLI flags, e.g.
   `python3 pangenome_viewer.py --panaroo-csv /abs/path.csv --prokka-dir /abs/dir`)
3. If this machine has no GUI (e.g. a remote server), reach it from your
   laptop via an SSH tunnel:
   ```
   ssh -L 8765:localhost:8765 <this-host>
   ```
   then open http://localhost:8765/ there.

## `pgv_lookup.py`

A separate, small CLI for instant single-cluster/single-genome spot-checks
(byte-offset index over the Panaroo CSV instead of a full scan) -- handy
for verifying a specific claim by hand without opening the web UI. Shares
the same config file / CLI flags as `pangenome_viewer.py`.

```
python3 pgv_lookup.py --build-index          # one-time, ~15-30s
python3 pgv_lookup.py group_4219 SOME_GENOME_ID
python3 pgv_lookup.py "some_cluster_name"    # row-level summary, no genome
```

## Design notes

`BUILD_NOTES_pangenome_viewer.md` has the fuller history of this tool's
design decisions (gutter markers, the genome-count toggle, the metadata
filter boxes, the refound-placeholder handling, etc.) from when it lived
inside the larger analysis project this was extracted from.
