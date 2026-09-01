# Pansynteny

A local, interactive gene-order viewer for a Panaroo pangenome. Search any
gene cluster by name or annotation, then build a gene synteny
chart on demand around this anchor gene, +/- a flanking window, sampled across carrier
genomes. For annotation, the script accepts a list of genes to highlight, plus a metadata file. 
A quick comparison bar shows whether each highlighted gene is visible in the window / present
on another contig or is truly absent in the strain.

Dataset-agnostic: a genome is viewable if it has a column in Panaroo's
`gene_presence_absence.csv` and a matching Prokka GFF -- nothing else is
required. Metadata is optional. 

## Requirements

Python 3, standard library only (no pip installs).

## Before you start: build the pangenome

Pansynteny reads the *output* of Prokka + Panaroo -- it doesn't run either
itself. If you don't already have these for your genome set:

1. **Annotate every assembly with Prokka**, one genome per output
   directory, each named after that genome (this is what becomes the
   `genome_id`):
   ```
   for fasta in /path/to/assemblies/*.fasta; do
     stem=$(basename "$fasta" .fasta)
     prokka --outdir prokka_out/"$stem" --prefix "$stem" "$fasta"
   done
   ```
2. **Run Panaroo across every genome's Prokka GFF** to build the
   pangenome:
   ```
   panaroo -i prokka_out/*/*.gff -o panaroo_out --clean-mode strict
   ```
   This produces `panaroo_out/gene_presence_absence.csv`, which is what
   `panaroo_csv` in the config points to.

## Setup

Every external data path is an absolute path, given either as a CLI flag
or in a config file -- nothing is inferred from where you put this folder.

`prokka_dir` must have the layout Prokka's own `--outdir`/`--prefix`
produces: one subdirectory per genome, named after that genome's stem,
containing at minimum `<stem>.gff` (required -- a genome without one
isn't viewable at all) and, optionally, `<stem>.faa`/`<stem>.ffn` (used
for the per-strain protein/nucleotide sequence panel; the tool works
without them, that panel just has nothing to show).
```
prokka_dir/
  GENOME_A/
    GENOME_A.gff
    GENOME_A.faa
    GENOME_A.ffn
    ... (Prokka's other output files, unused)
  GENOME_B/
    GENOME_B.gff
    ...
```
A genome is only viewable if its stem is both a genome-column header in
`panaroo_csv` *and* has a matching `<stem>.gff` under `prokka_dir` --
both are required, neither is inferred from the other.

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
