#!/usr/bin/env python3
"""
Local interactive Panaroo pangenome gene-order viewer.

Search for any Panaroo gene cluster by name or annotation, then build an
anchor-gene gene-order chart on demand: the anchor +/- a flanking window,
sampled across carrier genomes, RFE-important flanking genes colored,
gutter markers showing whether each top colored gene is visible/elsewhere/
absent, a genome-count toggle, and client-side metadata filter boxes.

Dataset-agnostic: every external input (Panaroo's gene_presence_absence.csv,
a Prokka GFF directory, and four optional extras -- genome metadata, an
RFE feature-importance list, an MGE-type table, a held-out test-score
table) is an absolute path given via CLI flag or a config file, never
inferred from this script's own location. A genome is viewable here if it
has both a column in the Panaroo CSV and a matching Prokka GFF -- nothing
else is required, so this works on any Panaroo+Prokka pangenome, not just
one that also happens to have a geNomad run behind it.

Config file: if --config isn't given, this script looks for
'pangenome_viewer.config' in its own directory (see
pangenome_viewer.config.example for the format -- simple 'key = value'
lines). CLI flags override the config file; the config file overrides
nothing else. Only --panaroo-csv/--prokka-dir (or their config-file
equivalents) are required; everything else is optional and the
corresponding feature just turns itself off if not supplied:
  - no --metadata-tsv: rows aren't grouped by serotype/host, just one
    section; row labels omit the source-type field; the client-side
    metadata filter boxes have nothing to offer.
  - no --rfe-features-txt: no gene gets colored, no RFE badges in search.
  - no --mge-genes-tsv: no MGE-type (plasmid/virus/provirus) labels, and
    multi-copy genomes fall back to picking a locus arbitrarily rather
    than preferring an MGE-classified copy.
  - no --test-predictions-csv: no "[test score X.XXX]" row annotation.

No pre-built index/database (deliberate -- see project plan): the search
index is built once at server startup with a single pass over the Panaroo
CSV, kept in memory; each chart build re-scans that same CSV plus one
Prokka GFF per displayed genome (~10-30s at the default 500-genome count;
scales roughly linearly with the genome-count toggle since GFF parsing
dominates -- "all" on a common gene can mean parsing thousands of GFFs).
Genome selection is a fixed shuffle of the carrier list (seed 42), sliced
to the requested count, so counts nest (1000 is 500's genomes plus 500
more) instead of each count drawing an unrelated random sample.

Usage:
    python3 pangenome_viewer.py --panaroo-csv /abs/path/gene_presence_absence.csv \\
        --prokka-dir /abs/path/prokka_out [--metadata-tsv ... --rfe-features-txt ...]
    # or with a config file (see pangenome_viewer.config.example):
    python3 pangenome_viewer.py

Then, since this machine has no GUI, reach it from a laptop via an SSH
tunnel: `ssh -L 8765:localhost:8765 <this-host>`, then open
http://localhost:8765/ there.
"""
import argparse
import csv
import json
import random
import re
import sys
import urllib.parse
from collections import Counter
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

BASE = Path(__file__).resolve().parent
TEMPLATE_HTML = BASE / "viewer_template.html"
CONFIG_FILENAME = "pangenome_viewer.config"

# Resolved once in main() via configure_globals(), before Context() is built
# or the server starts. Declared here (rather than threaded through every
# function signature) to keep the diff against the single-dataset version
# of this tool small -- everything below still just reads these names.
PANAROO_CSV = None
PROKKA_BASE = None
METADATA_TSV = None
METADATA_JOIN_COL = None
METADATA_O_ANTIGEN_COL = None
METADATA_H_ANTIGEN_COL = None
METADATA_SOURCE_COL = None
RFE_FEATURES_TXT = None
MGE_GENES_TSV = None
TEST_PREDICTIONS_CSV = None
TEST_PREDICTIONS_GENOME_COL = None
TEST_PREDICTIONS_SCORE_COL = None
STRIP_GENOME_SUFFIXES = ()

DEFAULT_COUNT = 500
COUNT_OPTIONS = [100, 200, 500, 1000, 2000, "all"]
TOP_N_COLORED = 16  # RFE features tracked + individually colored in the
                     # strain-comparison column
TOP_N_SEROTYPES = 10  # serotypes kept as their own row-group before "Other"
RANDOM_SEED = 42
WINDOW_OPTIONS_BP = [5000, 10000, 20000, 30000, 40000]
BUILD_WINDOW = WINDOW_OPTIONS_BP[-1]  # always build the widest superset;
                                       # client filters down for the toggle
DEFAULT_WINDOW = 20000
SEARCH_LIMIT = 25

csv.field_size_limit(10_000_000)


# ---------------------------------------------------------------------------
# CLI / config-file settings resolution.
# ---------------------------------------------------------------------------

CONFIG_KEYS = [
    "panaroo_csv", "prokka_dir",
    "metadata_tsv", "metadata_join_col", "metadata_o_antigen_col",
    "metadata_h_antigen_col", "metadata_source_col",
    "rfe_features_txt", "mge_genes_tsv",
    "test_predictions_csv", "test_predictions_genome_col", "test_predictions_score_col",
    "strip_genome_suffixes",
]


def build_arg_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--bind", default="127.0.0.1",
                     help="localhost-only by default; reach it via an SSH "
                          "tunnel, don't widen this to 0.0.0.0")
    ap.add_argument("--config", default=None,
                     help=f"path to a config file; if omitted, looks for "
                          f"'{CONFIG_FILENAME}' next to this script")
    ap.add_argument("--panaroo-csv", dest="panaroo_csv", default=None,
                     help="absolute path to Panaroo's gene_presence_absence.csv (required)")
    ap.add_argument("--prokka-dir", dest="prokka_dir", default=None,
                     help="absolute path to the Prokka output dir, one <stem>/<stem>.gff "
                          "per genome (required)")
    ap.add_argument("--metadata-tsv", dest="metadata_tsv", default=None,
                     help="absolute path to a tab-delimited genome metadata table (optional)")
    ap.add_argument("--metadata-join-col", dest="metadata_join_col", default=None,
                     help="column in --metadata-tsv matching Panaroo/Prokka genome stems "
                          "(required if --metadata-tsv is given)")
    ap.add_argument("--metadata-o-antigen-col", dest="metadata_o_antigen_col", default=None,
                     help="optional -- paired with --metadata-h-antigen-col to compute a "
                          "combined serotype column used for row grouping/coloring")
    ap.add_argument("--metadata-h-antigen-col", dest="metadata_h_antigen_col", default=None)
    ap.add_argument("--metadata-source-col", dest="metadata_source_col", default=None,
                     help="optional -- column shown in each row's label (e.g. host/source)")
    ap.add_argument("--rfe-features-txt", dest="rfe_features_txt", default=None,
                     help="absolute path to a 'feature'/'importance' TSV (optional -- "
                          "disables RFE-based coloring if omitted)")
    ap.add_argument("--mge-genes-tsv", dest="mge_genes_tsv", default=None,
                     help="absolute path to a pangenome_gene_cluster/genome_id/locus_tag/"
                          "mge_type TSV (optional -- disables MGE-type labels if omitted)")
    ap.add_argument("--test-predictions-csv", dest="test_predictions_csv", default=None,
                     help="absolute path to a held-out test-score CSV (optional)")
    ap.add_argument("--test-predictions-genome-col", dest="test_predictions_genome_col", default=None,
                     help="required if --test-predictions-csv is given")
    ap.add_argument("--test-predictions-score-col", dest="test_predictions_score_col", default=None,
                     help="required if --test-predictions-csv is given")
    ap.add_argument("--strip-genome-suffixes", dest="strip_genome_suffixes", default=None,
                     help="comma-separated suffixes to strip from Panaroo/Prokka stems for a "
                          "shorter display genome_id, e.g. '.result,.scaffold' (optional, cosmetic)")
    return ap


def load_config_file(path):
    """Simple 'key = value' lines, '#' comments, blank lines ignored. Never
    a hard requirement -- every key is also settable via CLI flag."""
    if not path.exists():
        return {}
    out = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, val = line.partition("=")
        out[key.strip()] = val.strip()
    return out


def resolve_settings(args):
    """CLI flag > config-file value > unset. No baked-in dataset-specific
    defaults for any path or column name -- see module docstring for which
    features silently disable themselves when their input is unset."""
    config_path = Path(args.config) if args.config else (BASE / CONFIG_FILENAME)
    config = load_config_file(config_path)

    resolved = {}
    for key in CONFIG_KEYS:
        cli_val = getattr(args, key)
        resolved[key] = cli_val if cli_val is not None else (config.get(key) or None)

    missing_required = [k for k in ("panaroo_csv", "prokka_dir") if not resolved[k]]
    if missing_required:
        sys.exit(f"error: missing required setting(s): {', '.join(missing_required)} "
                  f"-- pass --{missing_required[0].replace('_', '-')} or set it in "
                  f"{config_path} (see pangenome_viewer.config.example)")

    if resolved["metadata_tsv"] and not resolved["metadata_join_col"]:
        sys.exit("error: --metadata-tsv is set but --metadata-join-col is missing "
                  "(required together)")
    if resolved["test_predictions_csv"] and not (
            resolved["test_predictions_genome_col"] and resolved["test_predictions_score_col"]):
        sys.exit("error: --test-predictions-csv is set but --test-predictions-genome-col/"
                  "--test-predictions-score-col are missing (required together)")

    for key in ("panaroo_csv", "prokka_dir", "metadata_tsv", "rfe_features_txt",
                "mge_genes_tsv", "test_predictions_csv"):
        val = resolved[key]
        if val and not Path(val).exists():
            sys.exit(f"error: --{key.replace('_', '-')} path does not exist: {val}")

    resolved["strip_genome_suffixes"] = tuple(
        s.strip() for s in (resolved["strip_genome_suffixes"] or "").split(",") if s.strip()
    )
    return resolved


def configure_globals(resolved):
    global PANAROO_CSV, PROKKA_BASE, METADATA_TSV, METADATA_JOIN_COL
    global METADATA_O_ANTIGEN_COL, METADATA_H_ANTIGEN_COL, METADATA_SOURCE_COL
    global RFE_FEATURES_TXT, MGE_GENES_TSV
    global TEST_PREDICTIONS_CSV, TEST_PREDICTIONS_GENOME_COL, TEST_PREDICTIONS_SCORE_COL
    global STRIP_GENOME_SUFFIXES

    PANAROO_CSV = Path(resolved["panaroo_csv"])
    PROKKA_BASE = Path(resolved["prokka_dir"])
    METADATA_TSV = Path(resolved["metadata_tsv"]) if resolved["metadata_tsv"] else None
    METADATA_JOIN_COL = resolved["metadata_join_col"]
    METADATA_O_ANTIGEN_COL = resolved["metadata_o_antigen_col"]
    METADATA_H_ANTIGEN_COL = resolved["metadata_h_antigen_col"]
    METADATA_SOURCE_COL = resolved["metadata_source_col"]
    RFE_FEATURES_TXT = Path(resolved["rfe_features_txt"]) if resolved["rfe_features_txt"] else None
    MGE_GENES_TSV = Path(resolved["mge_genes_tsv"]) if resolved["mge_genes_tsv"] else None
    TEST_PREDICTIONS_CSV = Path(resolved["test_predictions_csv"]) if resolved["test_predictions_csv"] else None
    TEST_PREDICTIONS_GENOME_COL = resolved["test_predictions_genome_col"]
    TEST_PREDICTIONS_SCORE_COL = resolved["test_predictions_score_col"]
    STRIP_GENOME_SUFFIXES = resolved["strip_genome_suffixes"]


# ---------------------------------------------------------------------------
# Loaded once at startup, kept resident for the life of the process.
# ---------------------------------------------------------------------------

def build_genome_stem_map():
    """genome_id -> stem. stem = a genome-column header in the Panaroo CSV
    that also has a matching Prokka GFF (<stem>/<stem>.gff under
    PROKKA_BASE) -- a genome is viewable here if and only if both are true.
    No geNomad dependency at all: an earlier version of this tool borrowed
    genome discovery from a dereplication pipeline that gated on a
    completed geNomad run, which had nothing to do with what this viewer
    actually needs and made every new dataset require running geNomad
    first just to browse its pangenome. genome_id defaults to the raw
    stem; STRIP_GENOME_SUFFIXES optionally strips a cosmetic suffix for a
    shorter display ID (e.g. a dataset whose Panaroo/Prokka stems carry a
    trailing '.result'/'.scaffold' from the original FASTA filenames)."""
    with open(PANAROO_CSV, newline="") as f:
        header = next(csv.reader(f))
    out = {}
    for stem in header[3:]:
        gff = PROKKA_BASE / stem / f"{stem}.gff"
        if not gff.exists():
            continue
        genome_id = stem
        for suf in STRIP_GENOME_SUFFIXES:
            if genome_id.endswith(suf):
                genome_id = genome_id[: -len(suf)]
                break
        out[genome_id] = stem
    return out


def load_rfe_importance():
    """feature (Panaroo cluster name, possibly '~~~'-merged) -> importance
    str. Empty if RFE_FEATURES_TXT wasn't supplied -- RFE-based coloring
    just turns itself off, everything renders as a plain neutral gene."""
    if RFE_FEATURES_TXT is None:
        return {}
    out = {}
    with open(RFE_FEATURES_TXT, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            out[row["feature"]] = row["importance"]
    return out


def build_search_index(rfe_importance):
    """One pass over the Panaroo CSV: every cluster's name, merged
    gene-name variants, annotation, and genome-wide carrier count (free to
    compute here since the whole row is already being read)."""
    index = []
    with open(PANAROO_CSV, newline="") as f:
        r = csv.reader(f)
        header = next(r)
        n_genomes = len(header) - 3
        for row in r:
            cluster = row[0]
            carrier_count = sum(1 for cell in row[3:] if cell.strip())
            index.append({
                "cluster": cluster,
                "non_unique_name": row[1],
                "annotation": row[2],
                "carrier_count": carrier_count,
                "is_rfe": cluster in rfe_importance,
                "rfe_importance": rfe_importance.get(cluster),
                "_haystack": f"{cluster} {row[1]} {row[2]}".lower(),
            })
    print(f"[startup] indexed {len(index)} pangenome clusters across "
          f"{n_genomes} genomes", file=sys.stderr)
    return index


SYNTHETIC_SEROTYPE_COL = "Serotype (O:H combined, computed)"


def load_metadata():
    """Returns (meta, columns). meta: join-key -> {source_type, serotype,
    full}, where `full` is every raw column from the metadata TSV (for the
    generic column/value filter boxes) plus, if both antigen columns are
    configured and present, one synthetic column combining them. Empty if
    METADATA_TSV wasn't supplied. Each of the three semantic sub-features
    (join, computed serotype, source label) degrades independently and
    loudly (a stderr warning, not a silent no-op) if its configured column
    name isn't actually present in this file's header -- a typo here would
    otherwise look like "every genome has Unknown metadata" with no clue
    why."""
    if METADATA_TSV is None:
        return {}, []
    meta = {}
    with open(METADATA_TSV, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fieldnames = list(r.fieldnames)
        if METADATA_JOIN_COL not in fieldnames:
            print(f"[startup] WARNING: metadata_join_col {METADATA_JOIN_COL!r} not found in "
                  f"{METADATA_TSV} header -- metadata will not be loaded at all", file=sys.stderr)
            return {}, []
        has_sero = bool(METADATA_O_ANTIGEN_COL and METADATA_H_ANTIGEN_COL
                         and METADATA_O_ANTIGEN_COL in fieldnames
                         and METADATA_H_ANTIGEN_COL in fieldnames)
        if METADATA_O_ANTIGEN_COL and not has_sero:
            print("[startup] WARNING: metadata_o_antigen_col/metadata_h_antigen_col not both "
                  "found in the metadata header -- skipping the computed serotype column",
                  file=sys.stderr)
        has_source = bool(METADATA_SOURCE_COL and METADATA_SOURCE_COL in fieldnames)
        if METADATA_SOURCE_COL and not has_source:
            print(f"[startup] WARNING: metadata_source_col {METADATA_SOURCE_COL!r} not found in "
                  f"the metadata header -- row labels will show 'Unknown' for it", file=sys.stderr)
        columns = fieldnames + ([SYNTHETIC_SEROTYPE_COL] if has_sero else [])
        for row in r:
            key = row.get(METADATA_JOIN_COL)
            if not key:
                continue
            full = {k: (v or "").strip() for k, v in row.items()}
            if has_sero:
                o_ag = full.get(METADATA_O_ANTIGEN_COL, "")
                h_ag = full.get(METADATA_H_ANTIGEN_COL, "")
                sero = f"{o_ag}:{h_ag}".strip(":") or "Unknown"
                full[SYNTHETIC_SEROTYPE_COL] = sero
            else:
                sero = "Unknown"
            source_type = full.get(METADATA_SOURCE_COL, "Unknown") if has_source else "Unknown"
            meta[key] = {"source_type": source_type or "Unknown", "serotype": sero, "full": full}
    return meta, columns


def compute_metadata_value_totals(meta, columns):
    """column -> {value: genome count}, across ALL genomes in the metadata
    table (not restricted to any gene's carriers) -- the cohort-wide
    denominator for the filter boxes' "what % of <value> genomes carry this
    gene" stat. Anchor-independent, so computed once here at startup rather
    than per chart build. Naturally empty if no metadata was loaded."""
    totals = {c: Counter() for c in columns}
    for m in meta.values():
        full = m["full"]
        for c in columns:
            totals[c][full.get(c, "")] += 1
    return {c: dict(cnt) for c, cnt in totals.items()}


def load_test_scores():
    """Empty if TEST_PREDICTIONS_CSV wasn't supplied -- the "[test score
    X.XXX]" row annotation just doesn't appear."""
    if TEST_PREDICTIONS_CSV is None:
        return {}
    scores = {}
    with open(TEST_PREDICTIONS_CSV, newline="") as f:
        for row in csv.DictReader(f):
            scores[row[TEST_PREDICTIONS_GENOME_COL]] = float(row[TEST_PREDICTIONS_SCORE_COL])
    return scores


def parse_gff(stem):
    """Returns (genes, seqlens). genes: locus_tag -> (contig, start, end,
    strand, product). seqlens: contig -> length (from ##sequence-region)."""
    genes = {}
    seqlens = {}
    gff = PROKKA_BASE / stem / f"{stem}.gff"
    if not gff.exists():
        return genes, seqlens
    with open(gff) as f:
        for line in f:
            if line.startswith(">"):
                break
            if line.startswith("##sequence-region"):
                _, sid, s, e = line.split()
                seqlens[sid] = int(e)
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            m = re.search(r"locus_tag=([^;]+)", fields[8])
            if not m:
                continue
            pm = re.search(r"product=([^;]+)", fields[8])
            product = urllib.parse.unquote(pm.group(1)) if pm else "hypothetical protein"
            genes[m.group(1)] = (fields[0], int(fields[3]), int(fields[4]), fields[6], product)
    return genes, seqlens


# Runs of N shorter than this are single-base ambiguity calls, not scaffold
# gaps -- 33% of runs in this cohort are 1-9bp. 10bp filters those out while
# keeping every genuine assembly gap (observed gaps range ~10bp to ~1kb).
MIN_N_RUN = 10


def find_n_runs(stem, contig, min_len=MIN_N_RUN):
    """Scaffold-gap N-runs on ONE contig of a genome's assembly.

    The scaffold sequence lives in the ##FASTA block at the tail of the Prokka
    GFF. We read only the requested contig's sequence and stop -- because the
    anchor gene almost always sits on the first (largest) contig, this breaks
    out after a few KB rather than reading the whole ~5.5MB file. Returns a
    list of (start, end) 1-based inclusive coordinates for each run of >=
    min_len N's, in contig coordinates."""
    gff = PROKKA_BASE / stem / f"{stem}.gff"
    if not gff.exists():
        return []
    parts = []
    in_fasta = False
    in_target = False
    with open(gff) as f:
        for line in f:
            if not in_fasta:
                if line.startswith("##FASTA"):
                    in_fasta = True
                continue
            if line.startswith(">"):
                if in_target:
                    break  # captured the whole target contig; stop reading
                in_target = line[1:].split()[0] == contig
                continue
            if in_target:
                parts.append(line.strip())
    if not parts:
        return []
    seq = "".join(parts)
    return [(m.start() + 1, m.end())
            for m in re.finditer(r"[Nn]{%d,}" % min_len, seq)]


class Context:
    """Everything loaded once at server startup."""
    def __init__(self):
        print("[startup] loading RFE feature importances...", file=sys.stderr)
        self.rfe_importance = load_rfe_importance()
        self.rfe_features = set(self.rfe_importance)

        print("[startup] building genome stem map (Panaroo CSV header "
              "∩ Prokka GFFs)...", file=sys.stderr)
        self.genome_stem_map = build_genome_stem_map()
        self.stem_to_genome = {stem: g for g, stem in self.genome_stem_map.items()}
        print(f"[startup] {len(self.genome_stem_map)} genomes resolved", file=sys.stderr)

        print(f"[startup] building cluster search index (one pass over "
              f"{PANAROO_CSV.name}, this can take a while)...", file=sys.stderr)
        self.search_index = build_search_index(self.rfe_importance)
        self.cluster_lookup = {c["cluster"]: c for c in self.search_index}

        print("[startup] loading genome metadata + test-split scores...", file=sys.stderr)
        self.metadata, self.metadata_columns = load_metadata()
        self.metadata_value_totals = compute_metadata_value_totals(self.metadata, self.metadata_columns)
        self.test_scores = load_test_scores()
        print("[startup] ready.", file=sys.stderr)


# ---------------------------------------------------------------------------
# Chart builder.
# ---------------------------------------------------------------------------

def load_mge_types_for_cluster(anchor_cluster):
    """(genome_id, locus_tag) -> mge_type, for this anchor's gene instances
    only -- used for copy-selection when a genome carries >1 paralog, and
    for the row-label mge_type field. Empty if MGE_GENES_TSV wasn't
    supplied."""
    if MGE_GENES_TSV is None:
        return {}
    types = {}
    with open(MGE_GENES_TSV, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["pangenome_gene_cluster"] != anchor_cluster:
                continue
            types[(row["genome_id"], row["locus_tag"])] = row["mge_type"]
    return types


def build_chart_data(ctx, anchor_cluster, count=DEFAULT_COUNT, include_n_runs=False):
    if anchor_cluster not in ctx.cluster_lookup:
        return {"error": f"cluster {anchor_cluster!r} not found"}

    mge_types = load_mge_types_for_cluster(anchor_cluster)

    with open(PANAROO_CSV, newline="") as f:
        r = csv.reader(f)
        header = next(r)
        genome_cols = header[3:]
        anchor_row = None
        for row in r:
            if row[0] == anchor_cluster:
                anchor_row = row
                break
    if anchor_row is None:
        return {"error": f"cluster {anchor_cluster!r} not found in Panaroo CSV"}

    present = {}  # bare genome_id -> [locus_tag, ...]
    for stem, cell in zip(genome_cols, anchor_row[3:]):
        cell = cell.strip()
        if not cell:
            continue
        g = ctx.stem_to_genome.get(stem)
        if g is None:
            continue
        present[g] = [x.strip() for x in cell.split(";") if x.strip()]

    if not present:
        return {"error": f"cluster {anchor_cluster!r} has 0 carriers genome-wide"}

    # Genome-wide refound-placeholder stats -- free to compute here, no extra
    # file access: Panaroo's gene-refinding step (genes Prokka's own CDS pass
    # missed but Panaroo's re-search of the raw assembly found) writes
    # placeholder locus tags like "9_refound_2105" instead of a real Prokka
    # locus tag, and those never appear in the genome's own GFF -- the
    # locus-tag string already sitting in `present` (parsed from the same
    # CSV read above) is enough to detect this with a substring check, no
    # second scan of anything needed.
    carrier_loci_total = sum(len(loci) for loci in present.values())
    carrier_loci_refound = sum(1 for loci in present.values() for lt in loci if "refound" in lt)

    def lookup_meta(genome_id):
        m = ctx.metadata.get(genome_id)
        if not m:
            for suf in (".result", ".scaffold"):
                m = ctx.metadata.get(genome_id + suf)
                if m:
                    break
        return m

    # Cohort-wide-carrier crosstab: for every metadata column, how many of
    # THIS gene's full genome-wide carrier set (all of `present`, regardless
    # of the count toggle/sample) have each value. Paired client-side with
    # ctx.metadata_value_totals (total genomes with that value, gene
    # -independent, precomputed once at startup) to answer "what % of <value>
    # genomes carry this gene" -- a cohort-wide stat, deliberately
    # independent of both the genome-count toggle and the row-filter boxes
    # (which only narrow the displayed *sample*, not this).
    carrier_value_counts = {c: Counter() for c in ctx.metadata_columns}
    for g in present:
        m = lookup_meta(g)
        full = m["full"] if m else {}
        for c in ctx.metadata_columns:
            carrier_value_counts[c][full.get(c, "")] += 1
    carrier_value_counts = {c: dict(cnt) for c, cnt in carrier_value_counts.items()}

    # Fixed shuffle order (not random.sample per-count) so that different
    # counts are nested subsets of each other -- picking 1000 after 500 shows
    # the same 500 genomes plus 500 more, rather than a fresh unrelated draw.
    random.seed(RANDOM_SEED)
    shuffled = list(present.keys())
    random.shuffle(shuffled)
    sample_size = len(shuffled) if count == "all" else min(count, len(shuffled))
    sample_genome_ids = shuffled[:sample_size]

    stems = {g: ctx.genome_stem_map[g] for g in sample_genome_ids}
    gff_parsed = {g: parse_gff(stem) for g, stem in stems.items()}
    gff_genes = {g: genes for g, (genes, _) in gff_parsed.items()}
    gff_seqlens = {g: seqlens for g, (_, seqlens) in gff_parsed.items()}

    def pick_locus(g):
        loci = present[g]
        if len(loci) == 1:
            return loci[0]
        # Prefer a real Prokka-called locus over a Panaroo "refound"
        # placeholder -- a refound locus has no GFF entry to resolve
        # coordinates from at all, so picking one when a real locus is also
        # available for this genome would drop a row we could otherwise show.
        real = [lt for lt in loci if "refound" not in lt]
        candidates = real if real else loci
        classified = [lt for lt in candidates if (g, lt) in mge_types]
        return classified[0] if classified else candidates[0]

    anchor_locus = {g: pick_locus(g) for g in sample_genome_ids}

    col_for_genome = {g: stem for g, stem in stems.items()}
    reverse_map = {g: {} for g in stems}

    with open(PANAROO_CSV, newline="") as f:
        r = csv.reader(f)
        header = next(r)
        gene_idx = header.index("Gene")
        ann_idx = header.index("Annotation")
        col_idx = {}
        for g, stem in col_for_genome.items():
            if stem in header:
                col_idx[g] = header.index(stem)

        for row in r:
            cluster = row[gene_idx]
            annotation = row[ann_idx]
            for g, ci in col_idx.items():
                cell = row[ci].strip()
                if not cell:
                    continue
                loci = [x.strip() for x in cell.split(";") if x.strip()]
                for lt in loci:
                    reverse_map[g][lt] = (cluster, annotation)

    def sero_for(genome_id):
        m = lookup_meta(genome_id)
        return (m["serotype"] if m else "Unknown"), (m["source_type"] if m else "Unknown")

    sero_counts = Counter(sero_for(g)[0] for g in stems)
    top_seros = [s for s, _ in sero_counts.most_common(TOP_N_SEROTYPES)]

    def sero_bucket(genome_id):
        s, _ = sero_for(genome_id)
        return s if s in top_seros else "Other"

    rows_out = []
    mge_type_counts = Counter()
    rfe_cluster_counts = Counter()
    n_dropped_refound = 0
    n_dropped_other = 0

    for g in sample_genome_ids:
        if g not in stems:
            continue
        genes = gff_genes.get(g, {})
        lt = anchor_locus.get(g)
        if not lt or lt not in genes:
            if lt and "refound" in lt:
                n_dropped_refound += 1
            else:
                n_dropped_other += 1
            continue
        contig, e_start, e_end, e_strand, e_product = genes[lt]

        mge_type = mge_types.get((g, lt), "chromosome/unclassified")
        mge_type_counts[mge_type] += 1

        flip = e_strand == "-"
        anchor_ref = e_start if not flip else e_end
        win_lo, win_hi = e_start - BUILD_WINDOW, e_end + BUILD_WINDOW

        genes_out = []
        for glt, (g_contig, g_start, g_end, g_strand, g_product) in genes.items():
            if g_contig != contig:
                continue
            if g_end < win_lo or g_start > win_hi:
                continue
            if flip:
                rel_start = anchor_ref - g_end
                rel_end = anchor_ref - g_start
                disp_strand = "-" if g_strand == "+" else "+"
            else:
                rel_start = g_start - anchor_ref
                rel_end = g_end - anchor_ref
                disp_strand = g_strand

            cluster, annotation = reverse_map[g].get(glt, ("NA", None))
            product = urllib.parse.unquote(annotation.split(";")[0]) if annotation else g_product
            is_rfe = cluster in ctx.rfe_features
            if is_rfe and cluster != anchor_cluster:
                rfe_cluster_counts[cluster] += 1
            genes_out.append({
                "start": rel_start, "end": rel_end, "strand": disp_strand,
                "cluster": cluster, "product": product,
                "is_anchor": glt == lt,
                "is_rfe": is_rfe,
                "rfe_importance": None,
            })
        genes_out.sort(key=lambda x: x["start"])

        contiglen = gff_seqlens.get(g, {}).get(contig)
        if contiglen:
            if flip:
                contig_rel_start = anchor_ref - contiglen
                contig_rel_end = anchor_ref - 1
            else:
                contig_rel_start = 1 - anchor_ref
                contig_rel_end = contiglen - anchor_ref
        else:
            contig_rel_start = contig_rel_end = None

        # Scaffold-gap N-runs on the anchor's contig, mapped into the same
        # anchor-relative frame as the genes and clipped to the build window.
        # Off by default -- the FASTA reads add ~2s/build and gaps are sparse,
        # so this is opt-in via the client's "assembly gaps" checkbox.
        n_runs_out = []
        for ns, ne in (find_n_runs(stems[g], contig) if include_n_runs else ()):
            if ne < win_lo or ns > win_hi:
                continue
            if flip:
                r_start, r_end = anchor_ref - ne, anchor_ref - ns
            else:
                r_start, r_end = ns - anchor_ref, ne - anchor_ref
            n_runs_out.append({"start": r_start, "end": r_end})
        n_runs_out.sort(key=lambda x: x["start"])

        sero, source = sero_for(g)
        meta_row = lookup_meta(g)
        rows_out.append({
            "genome_id": g, "mge_type": mge_type,
            "serotype": sero, "sero_bucket": sero_bucket(g), "source_type": source,
            "test_score": ctx.test_scores.get(g),
            "contig_rel_start": contig_rel_start, "contig_rel_end": contig_rel_end,
            "genes": genes_out,
            "n_runs": n_runs_out,
            "_contig": contig,
            "metadata": meta_row["full"] if meta_row else {},
        })

    for row in rows_out:
        for gobj in row["genes"]:
            if gobj["is_rfe"]:
                gobj["rfe_importance"] = ctx.rfe_importance.get(gobj["cluster"])

    buckets_order = top_seros + (["Other"] if any(r["sero_bucket"] == "Other" for r in rows_out) else [])
    bucket_rank = {b: i for i, b in enumerate(buckets_order)}
    rows_out.sort(key=lambda r: (bucket_rank.get(r["sero_bucket"], 99), r["genome_id"]))

    top_rfe_clusters = [c for c, _ in rfe_cluster_counts.most_common(TOP_N_COLORED)]
    cluster_labels = {}
    for row in rows_out:
        for gobj in row["genes"]:
            if gobj["cluster"] in top_rfe_clusters and gobj["cluster"] not in cluster_labels:
                cluster_labels[gobj["cluster"]] = gobj["product"]

    # Per-row status marker for each of the tracked genes, whether or not
    # it's drawn in this row's window -- lets the client render one
    # strain-comparison column per gene so the whole chart is scannable for
    # presence even when synteny puts it out of frame.
    cluster_to_loci_cache = {}

    def cluster_to_loci(g):
        inv = cluster_to_loci_cache.get(g)
        if inv is None:
            inv = {}
            for lt, (cl, _ann) in reverse_map[g].items():
                inv.setdefault(cl, []).append(lt)
            cluster_to_loci_cache[g] = inv
        return inv

    for row in rows_out:
        g = row["genome_id"]
        row_contig = row.pop("_contig")
        visible = {gobj["cluster"] for gobj in row["genes"]}
        genes_g = gff_genes.get(g, {})
        markers = []
        for cl in top_rfe_clusters:
            if cl == anchor_cluster:
                continue
            if cl in visible:
                markers.append({"cluster": cl, "state": "visible"})
                continue
            loci = cluster_to_loci(g).get(cl, [])
            if not loci:
                continue  # truly absent -- no marker
            resolved = next((genes_g[lt] for lt in loci if lt in genes_g), None)
            if resolved is None:
                state = "refound_only"
            elif resolved[0] == row_contig:
                state = "same_contig_outside_window"
            else:
                state = "different_contig"
            markers.append({"cluster": cl, "state": state})
        row["gutter_markers"] = markers

    print(f"[build] {anchor_cluster}: {len(rows_out)}/{len(sample_genome_ids)} sampled "
          f"genomes shown; dropped {n_dropped_refound} (refound placeholder, no GFF "
          f"entry) + {n_dropped_other} (other lookup failure); genome-wide carrier "
          f"loci {carrier_loci_refound}/{carrier_loci_total} "
          f"({100 * carrier_loci_refound / carrier_loci_total:.1f}%) are refound",
          file=sys.stderr)

    meta = ctx.cluster_lookup[anchor_cluster]
    return {
        "group_id": f"{anchor_cluster}-anchored",
        "anchor_cluster": anchor_cluster,
        "anchor_annotation": meta["annotation"],
        "anchor_non_unique_name": meta["non_unique_name"],
        "total_group_members": len(present),
        "n_shown": len(rows_out),
        "n_dropped_refound": n_dropped_refound,
        "n_dropped_other": n_dropped_other,
        "carrier_loci_total": carrier_loci_total,
        "carrier_loci_refound": carrier_loci_refound,
        "sero_buckets": buckets_order,
        "top_rfe_clusters": top_rfe_clusters,
        "cluster_labels": cluster_labels,
        "min_n_run": MIN_N_RUN,
        "n_runs_included": include_n_runs,
        "max_window_bp": BUILD_WINDOW,
        "default_window_bp": DEFAULT_WINDOW,
        "window_options_bp": WINDOW_OPTIONS_BP,
        "count_options": COUNT_OPTIONS,
        "default_count": DEFAULT_COUNT,
        "requested_count": count,
        "sample_size": sample_size,
        "metadata_columns": ctx.metadata_columns,
        "metadata_source_col": METADATA_SOURCE_COL,
        "carrier_value_counts": carrier_value_counts,
        "rows": rows_out,
    }


# ---------------------------------------------------------------------------
# HTTP server.
# ---------------------------------------------------------------------------

def make_handler(ctx):
    template_bytes_cache = {}

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, fmt, *args):
            print(f"[http] {self.address_string()} {fmt % args}", file=sys.stderr)

        def _send_json(self, obj, status=200):
            body = json.dumps(obj).encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            parsed = urllib.parse.urlsplit(self.path)
            qs = urllib.parse.parse_qs(parsed.query)

            if parsed.path == "/":
                if "html" not in template_bytes_cache:
                    template_bytes_cache["html"] = TEMPLATE_HTML.read_bytes()
                body = template_bytes_cache["html"]
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return

            if parsed.path == "/api/search":
                q = (qs.get("q", [""])[0]).strip().lower()
                if not q:
                    self._send_json([])
                    return
                matches = [c for c in ctx.search_index if q in c["_haystack"]]
                matches.sort(key=lambda c: c["carrier_count"], reverse=True)
                out = [
                    {k: v for k, v in c.items() if not k.startswith("_")}
                    for c in matches[:SEARCH_LIMIT]
                ]
                self._send_json(out)
                return

            if parsed.path == "/api/metadata_totals":
                # Anchor-independent, fetched once by the client at page
                # load and reused for every subsequent chart build.
                self._send_json({
                    "columns": ctx.metadata_columns,
                    "totals": ctx.metadata_value_totals,
                })
                return

            if parsed.path == "/api/chart":
                anchor = (qs.get("anchor", [""])[0])
                if not anchor:
                    self._send_json({"error": "missing 'anchor' query param"}, status=400)
                    return
                raw_count = qs.get("count", [""])[0].strip()
                if not raw_count:
                    count = DEFAULT_COUNT
                elif raw_count == "all":
                    count = "all"
                else:
                    try:
                        count = int(raw_count)
                    except ValueError:
                        self._send_json({"error": f"invalid 'count' param {raw_count!r}"}, status=400)
                        return
                include_n_runs = qs.get("n_runs", [""])[0].strip() in ("1", "true", "yes")
                print(f"[build] anchor={anchor!r} count={count!r} n_runs={include_n_runs}", file=sys.stderr)
                data = build_chart_data(ctx, anchor, count, include_n_runs=include_n_runs)
                if "error" in data:
                    self._send_json(data, status=404)
                    return
                self._send_json(data)
                return

            self.send_response(404)
            self.end_headers()

    return Handler


def main():
    ap = build_arg_parser()
    args = ap.parse_args()
    resolved = resolve_settings(args)
    configure_globals(resolved)

    ctx = Context()

    server = ThreadingHTTPServer((args.bind, args.port), make_handler(ctx))
    print(f"[serve] listening on http://{args.bind}:{args.port}/ "
          f"(tunnel with: ssh -L {args.port}:localhost:{args.port} <this-host>)",
          file=sys.stderr)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[serve] shutting down", file=sys.stderr)


if __name__ == "__main__":
    main()
