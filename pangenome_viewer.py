#!/usr/bin/env python3
"""
Local interactive Panaroo pangenome gene-order viewer.

Search for any Panaroo gene cluster by name or annotation, then build an
anchor-gene gene-order chart on demand: the anchor +/- a flanking window,
sampled across carrier genomes, RFE-important flanking genes colored,
gutter markers showing whether each top colored gene is visible/elsewhere/
absent, a genome-count toggle, and client-side metadata filter boxes.

Dataset-agnostic: every external input (Panaroo's gene_presence_absence.csv,
a Prokka GFF directory, and five optional extras -- genome metadata, an
RFE feature-importance list, an MGE-type table, a held-out test-score
table, an stxtyper report) is an absolute path given via CLI flag or a
config file, never inferred from this script's own location. A genome is
viewable here if it has both a column in the Panaroo CSV and a matching
Prokka GFF -- nothing else is required, so this works on any Panaroo+Prokka pangenome, not just
one that also happens to have a geNomad run behind it.

Config file: if --config isn't given, this script looks for
'pangenome_viewer.config' in its own directory (see
pangenome_viewer.config.example for the format -- simple 'key = value'
lines). CLI flags override the config file; the config file overrides
nothing else. --setup builds that file interactively (tab-completing
paths, picking column names off the real header, and reporting how many
genomes/rows each join actually resolves) and is offered automatically
when a required setting is missing and stdin is a terminal -- under
nohup/systemd it stays the plain error message instead. Only --panaroo-csv/--prokka-dir (or their config-file
equivalents) are required; everything else is optional and the
corresponding feature just turns itself off if not supplied:
  - no --metadata-tsv: rows aren't grouped by serotype/host, just one
    section; row labels omit the source-type field; the client-side
    metadata filter boxes have nothing to offer.
  - no --rfe-features-txt: no gene gets colored, no RFE badges in search.
  - no --mge-genes-tsv: no MGE-type (plasmid/virus/provirus) labels, and
    multi-copy genomes fall back to picking a locus arbitrarily rather
    than preferring an MGE-classified copy.
  - no --stxtyper-tsv: no stx subtyping anywhere -- the stxA/stxB CDSs
    carry no stx_type/stx_operon/stx_identity, every row's stx_types is
    empty, and the "stx type (stxtyper)" column is absent from the filter
    boxes and row-label pickers.
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

The one exception: single-cluster/single-genome lookups that aren't part
of a chart build (the sequence panel, contig export) go through a small
on-disk byte-offset index (cluster_offsets.json, gitignored) instead of a
linear scan -- those are triggered by a UI click and need to feel instant,
unlike a chart build the user already expects to take a while. The index
records the mtime and size of the CSV it was built from and rebuilds
itself whenever either has changed: byte offsets are only valid for the
exact bytes they were computed over, and a stale one doesn't fail
cleanly -- it seeks into the middle of some other row and hands back a
different cluster's data as though it were the one asked for. Every seek
also asserts the row it landed on really is the requested cluster, so a
bad offset raises instead of lying. Same cache file dev/pgv_lookup.py's own
CLI uses, in the same format and with the same validation (see
build_cluster_offset_index()'s docstring for why the index-building code
itself isn't shared between the two). The chart-build scan itself is
untouched.

Usage:
    python3 pangenome_viewer.py --panaroo-csv /abs/path/gene_presence_absence.csv \\
        --prokka-dir /abs/path/prokka_out [--metadata-tsv ... --rfe-features-txt ...]
    # or with a config file (see pangenome_viewer.config.example):
    python3 pangenome_viewer.py
    # or build that config file interactively first:
    python3 pangenome_viewer.py --setup

Then, since this machine has no GUI, reach it from a laptop via an SSH
tunnel: `ssh -L 8765:localhost:8765 <this-host>`, then open
http://localhost:8765/ there.
"""
import argparse
import csv
import glob
import json
import os
import random
import re
import sys
import urllib.parse
from collections import Counter
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

BASE = Path(__file__).resolve().parent
TEMPLATE_HTML = BASE / "viewer_template.html"
# Same file dev/pgv_lookup.py's OFFSET_INDEX points at: that tool lives a
# directory down but resolves this path against the repo root, not its own
# -- whichever of the two builds it first, the other reuses it as-is.
# Deliberately not shared code between them (see
# build_cluster_offset_index()'s docstring for why).
CLUSTER_OFFSETS_PATH = BASE / "cluster_offsets.json"


def offset_index_stamp(path):
    """(mtime, size) of the file a byte-offset index was built from.

    Byte offsets are only meaningful for the exact bytes they were
    computed over, so an index has to be able to tell whether its source
    changed underneath it. Size is recorded alongside mtime because a file
    restored from a backup or written by a pipeline that preserves
    timestamps can change content without changing mtime."""
    st = Path(path).stat()
    return [st.st_mtime, st.st_size]
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
STXTYPER_TSV = None
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


def sniff_delimiter(path, default="\t"):
    """Tab or comma, decided from the header line.

    These tables are exported by hand from all sorts of places and arrive
    with either extension and either delimiter; reading a comma-delimited
    file as tab yields one giant column, so a lookup by column name fails
    with "no such column" about a column that is plainly there. Tab wins
    when both are present, since a TSV's free-text fields (annotations,
    descriptions) routinely contain commas while a CSV's rarely contain
    tabs."""
    try:
        with open(path, newline="") as f:
            first = f.readline()
    except OSError:
        return default
    if "\t" in first:
        return "\t"
    if "," in first:
        return ","
    return default


# ---------------------------------------------------------------------------
# CLI / config-file settings resolution.
# ---------------------------------------------------------------------------

CONFIG_KEYS = [
    "panaroo_csv", "prokka_dir",
    "metadata_tsv", "metadata_join_col", "metadata_o_antigen_col",
    "metadata_h_antigen_col", "metadata_source_col",
    "rfe_features_txt", "mge_genes_tsv", "stxtyper_tsv",
    "test_predictions_csv", "test_predictions_genome_col", "test_predictions_score_col",
    "strip_genome_suffixes",
]

# The only two settings without a "feature just turns itself off" fallback.
REQUIRED_KEYS = ("panaroo_csv", "prokka_dir")


def build_arg_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--bind", default="127.0.0.1",
                     help="localhost-only by default; reach it via an SSH "
                          "tunnel, don't widen this to 0.0.0.0")
    ap.add_argument("--config", default=None,
                     help=f"path to a config file; if omitted, looks for "
                          f"'{CONFIG_FILENAME}' next to this script")
    ap.add_argument("--setup", action="store_true",
                     help="interactively build a config file (prompting for each "
                          "input and validating it against the ones already given), "
                          "then exit; also offered automatically when a required "
                          "setting is missing and this is an interactive terminal")
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
    ap.add_argument("--stxtyper-tsv", dest="stxtyper_tsv", default=None,
                     help="absolute path to a combined stxtyper report -- one row per "
                          "operon call, with #name/stx_type/target_contig/target_start/"
                          "target_stop columns (optional -- disables stx subtyping if omitted)")
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


def missing_required_settings(args, config_path):
    """The required keys neither a CLI flag nor config_path supplies -- the
    same condition resolve_settings() exits on, checked separately so
    main() can offer --setup instead of just failing. Deliberately keyed on
    the settings being missing rather than on the config file being absent,
    so a legitimate flags-only run is never hijacked into a prompt."""
    config = load_config_file(config_path)
    return [k for k in REQUIRED_KEYS if not (getattr(args, k, None) or config.get(k))]


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

    missing_required = [k for k in REQUIRED_KEYS if not resolved[k]]
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
                "mge_genes_tsv", "stxtyper_tsv", "test_predictions_csv"):
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
    global RFE_FEATURES_TXT, MGE_GENES_TSV, STXTYPER_TSV
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
    STXTYPER_TSV = Path(resolved["stxtyper_tsv"]) if resolved["stxtyper_tsv"] else None
    TEST_PREDICTIONS_CSV = Path(resolved["test_predictions_csv"]) if resolved["test_predictions_csv"] else None
    TEST_PREDICTIONS_GENOME_COL = resolved["test_predictions_genome_col"]
    TEST_PREDICTIONS_SCORE_COL = resolved["test_predictions_score_col"]
    STRIP_GENOME_SUFFIXES = resolved["strip_genome_suffixes"]


# ---------------------------------------------------------------------------
# Interactive setup mode (--setup): writes a pangenome_viewer.config by
# filling values into pangenome_viewer.config.example, so the generated
# file keeps every explanatory comment the example carries rather than
# being a bare list of keys.
#
# Each prompt validates against the inputs already collected rather than
# just checking that a path exists: the prokka_dir answer is scored
# against the Panaroo CSV's genome columns, the metadata join column
# against the resulting genome IDs, and the RFE 'feature' column against
# the Panaroo cluster names. Those three joins are exactly the ones that
# otherwise fail silently at runtime (or, for a missing 'feature' column,
# crash mid-request), which is the whole reason this mode exists.
# ---------------------------------------------------------------------------

CONFIG_EXAMPLE_FILENAME = "pangenome_viewer.config.example"

# A join this weak almost always means the wrong column/file rather than a
# genuinely partial overlap, so it's worth an "are you sure" rather than a
# silent accept -- but it stays the user's call, since a deliberately
# partial metadata table is legitimate.
WEAK_JOIN_FRACTION = 0.10

# Config keys the walkthrough never prompts for; an existing value for one
# of these is preserved verbatim rather than blanked on a re-run.
UNPROMPTED_KEYS = ("test_predictions_csv", "test_predictions_genome_col",
                    "test_predictions_score_col")


def _use_color():
    """Colour only for a real terminal that hasn't opted out. Piping setup
    output to a file or a pager must not fill it with escape codes."""
    if not sys.stdout.isatty():
        return False
    if os.environ.get("NO_COLOR") is not None:
        return False
    return os.environ.get("TERM", "") not in ("", "dumb")


def _ok(text):
    """A check that passed -- green so the eye can find it while scrolling
    a long setup run. Falls back to a plain 'OK:' prefix, so the meaning
    survives when colour is off; the marker is never the only signal."""
    if _use_color():
        return f"\033[42;30m OK \033[0m {text}"
    return f"OK: {text}"


def _enable_path_completion():
    """Tab-completion over filesystem paths for the prompts below. The
    default completer delimiters include '/' and '-', which would hand the
    completer only the last path segment; clearing them down to whitespace
    is what makes completing a full absolute path work."""
    try:
        import readline
    except ImportError:  # pragma: no cover - readline is stdlib on Linux/macOS
        return

    def complete(text, state):
        stub = os.path.expanduser(text)
        matches = []
        for m in glob.glob(glob.escape(stub) + "*"):
            matches.append(m + "/" if os.path.isdir(m) else m)
        matches.sort()
        return matches[state] if state < len(matches) else None

    readline.set_completer_delims(" \t\n")
    readline.set_completer(complete)
    if "libedit" in (getattr(readline, "__doc__", "") or ""):
        readline.parse_and_bind("bind ^I rl_complete")
    else:
        readline.parse_and_bind("tab: complete")


def _prompt(text, default=None):
    """One line of input, or the default on an empty answer. Ctrl-D/Ctrl-C
    exits the whole setup rather than falling through with a half-built
    config."""
    suffix = f" [{default}]" if default else ""
    try:
        answer = input(f"{text}{suffix}: ").strip()
    except (EOFError, KeyboardInterrupt):
        sys.exit("\nsetup cancelled -- nothing was written")
    return answer or (default or "")


def _prompt_yes_no(text, default=True):
    hint = "Y/n" if default else "y/N"
    while True:
        answer = _prompt(f"{text} ({hint})").lower()
        if not answer:
            return default
        if answer in ("y", "yes"):
            return True
        if answer in ("n", "no"):
            return False
        print("  please answer y or n")


CLEAR_TOKEN = "none"


def _prompt_path(text, optional=False, must_be_dir=False, default=None):
    """An existing path, expanded and resolved to absolute (the rest of the
    tool requires absolute paths). Returns None when an optional prompt is
    skipped.

    With a `default` (a value carried over from an existing config), blank
    means "keep the default" -- so blank can no longer also mean "skip",
    and CLEAR_TOKEN is what unsets a path that was previously configured.
    A default whose path has since disappeared is shown but not silently
    accepted; the user has to retype or clear it."""
    hint = ""
    if optional:
        hint = f" ({CLEAR_TOKEN!r} to clear)" if default else " (blank to skip)"
    while True:
        answer = _prompt(f"{text}{hint}", default=default)
        if optional and answer.strip().lower() == CLEAR_TOKEN:
            return None
        if not answer:
            if optional:
                return None
            print("  a path is required here")
            continue
        path = Path(answer).expanduser()
        try:
            path = path.resolve()
        except OSError as exc:
            print(f"  cannot resolve that path: {exc}")
            continue
        if not path.exists():
            print(f"  no such path: {path}")
            continue
        if must_be_dir and not path.is_dir():
            print(f"  not a directory: {path}")
            continue
        return path


def _prompt_column(text, columns, optional=False, default=None):
    """Pick a column from a real header by number rather than typing it.
    Column names in these files routinely carry spaces and parentheses
    (e.g. 'Assembly barcode(Assembly stats)'), which is exactly the kind of
    value a hand-typed answer gets subtly wrong.

    `default` (a column name detected or carried over from an existing
    config) is offered as the blank answer, so the common case is one
    keypress."""
    print(f"\n{text}")
    default_num = None
    for i, col in enumerate(columns, 1):
        marker = ""
        if default is not None and col == default:
            default_num = i
            marker = "  <- default"
        print(f"  {i:>3}. {col}{marker}")
    if optional:
        print(f"    0. (skip -- leave this feature off)"
              f"{'  <- default' if default is None and optional else ''}")
    while True:
        answer = _prompt("  column number",
                          default=str(default_num) if default_num else None)
        if not answer:
            if optional:
                return None
            print("  pick a number from the list")
            continue
        try:
            choice = int(answer)
        except ValueError:
            print("  pick a number from the list")
            continue
        if optional and choice == 0:
            return None
        if 1 <= choice <= len(columns):
            return columns[choice - 1]
        print("  pick a number from the list")


def _read_header(path, delimiter=None):
    """First row of a delimited file, or None if it can't be read at all.
    Delimiter is sniffed unless one is given explicitly."""
    try:
        with open(path, newline="") as f:
            return next(csv.reader(f, delimiter=delimiter or sniff_delimiter(path)), None)
    except (OSError, csv.Error, UnicodeDecodeError) as exc:
        print(f"  cannot read {path}: {exc}")
        return None


def _setup_panaroo(default=None):
    """(csv path, genome-column stems). Accepts either the Panaroo output
    directory or the gene_presence_absence.csv inside it -- the directory
    is what people have in hand."""
    print("\n--- 1/6: Panaroo (required) ---")
    while True:
        path = _prompt_path("Panaroo output dir, or gene_presence_absence.csv directly",
                             default=default)
        if path.is_dir():
            candidate = path / "gene_presence_absence.csv"
            if not candidate.exists():
                print(f"  no gene_presence_absence.csv in {path}")
                continue
            path = candidate
        header = _read_header(path, delimiter=",")
        if header is None:
            continue
        # Panaroo's first three columns are Gene/Non-unique Gene name/
        # Annotation; every column after that is one genome.
        if len(header) <= 3:
            print(f"  {path} has no genome columns -- is this really Panaroo's "
                  "gene_presence_absence.csv?")
            continue
        stems = [s for s in header[3:] if s.strip()]
        print("  " + _ok(f"{len(stems)} genome columns "
                          f"(e.g. {', '.join(stems[:3])})"))
        return path, stems


def _setup_prokka(stems, default=None):
    """(prokka dir, stems that actually have a GFF). Mirrors
    build_genome_stem_map()'s rule exactly -- <stem>/<stem>.gff -- so the
    number reported here is the number of genomes the viewer will show."""
    print("\n--- 2/6: Prokka (required) ---")
    print("Expected layout: one subdirectory per genome, <stem>/<stem>.gff")
    while True:
        path = _prompt_path("Prokka output dir", must_be_dir=True, default=default)
        matched = [s for s in stems if (path / s / f"{s}.gff").exists()]
        if matched:
            print("  " + _ok(f"{len(matched)} of {len(stems)} Panaroo genome columns "
                              f"have a matching {path.name}/<stem>/<stem>.gff"))
        fraction = len(matched) / len(stems) if stems else 0
        if fraction >= WEAK_JOIN_FRACTION:
            return path, matched

        # A near-total miss is nearly always a naming mismatch between the
        # two tools' stems, and the only useful thing to show is both
        # sides of it side by side.
        print(f"  WARNING: only {len(matched)} of {len(stems)} Panaroo genome "
              "columns have a matching GFF.")
        subdirs = sorted(p.name for p in path.iterdir() if p.is_dir())[:3]
        print(f"    Panaroo genome columns look like: {', '.join(stems[:3]) or '(none)'}")
        print(f"    {path} subdirectories look like:  {', '.join(subdirs) or '(none)'}")
        print("    These stems have to match exactly (strip_genome_suffixes below is")
        print("    display-only and cannot rescue a mismatch here).")
        if matched and _prompt_yes_no("  Use this directory anyway?", default=False):
            return path, matched
        print("  Let's try another directory.")


# A suffix carried by at least this fraction of stems is a real assembler
# /pipeline suffix rather than part of a genome's name.
SUFFIX_MIN_FRACTION = 0.02


def detect_stem_suffixes(matched_stems):
    """Trailing '.xxx' suffixes shared by a meaningful share of the genome
    stems, longest first.

    Deliberately NOT "a suffix every stem shares": a real cohort is often
    assembled by more than one route, so the stems arrive as a mix (this
    dataset is 2995 '.result' + 254 '.scaffold'). Requiring one common
    suffix found nothing for such a mix, left the genome IDs with their
    suffixes attached, and took the metadata join to zero -- the failure
    looked like bad metadata rather than an unstripped suffix.

    Longest-first matters because build_genome_stem_map() strips at most
    one suffix per stem and breaks on the first match: with '.fasta'
    ordered before '.result.fasta', a stem 'X.result.fasta' would strip
    only '.fasta' and stop at 'X.result'."""
    counts = Counter()
    for stem in matched_stems:
        # Every trailing dotted segment, so '.result.fasta' contributes
        # both '.fasta' and '.result.fasta' as candidates.
        parts = stem.split(".")
        for i in range(1, len(parts)):
            counts["." + ".".join(parts[i:])] += 1
    if not matched_stems:
        return []
    floor = max(1, int(SUFFIX_MIN_FRACTION * len(matched_stems)))
    common = {suf for suf, n in counts.items() if n >= floor}

    # Keep only the longest form of any nested pair actually present, so
    # '.result.fasta' wins over the '.fasta' it contains -- emitting both
    # would let the shorter one match first and truncate the wrong amount.
    chosen = []
    for suf in sorted(common, key=len, reverse=True):
        if not any(longer.endswith(suf) for longer in chosen):
            chosen.append(suf)
    return sorted(chosen, key=len, reverse=True)


def _setup_strip_suffixes(matched_stems, default=None):
    """Cosmetic display-only suffix stripping -- but it decides what the
    genome IDs look like, so it also decides whether the metadata join in
    the next step can work at all."""
    suffixes = detect_stem_suffixes(matched_stems)
    if not suffixes:
        return default or ""
    counts = Counter()
    for stem in matched_stems:
        for suf in suffixes:
            if stem.endswith(suf):
                counts[suf] += 1
                break
    covered = sum(counts.values())
    print("\nDetected genome-stem suffix(es):")
    for suf in suffixes:
        print(f"  {suf:<20} on {counts[suf]} of {len(matched_stems)} stems")
    example = matched_stems[0]
    for suf in suffixes:
        if example.endswith(suf):
            example_stripped = example[: -len(suf)]
            break
    else:
        example_stripped = example
    print(f"  covers {covered}/{len(matched_stems)} stems; "
          f"{example!r} would display as {example_stripped!r}")
    print("Stripping is display-only -- it does not change which files are read --")
    print("but the stripped ID is what the metadata join below matches on.")
    if _prompt_yes_no("Strip these for display?", default=True):
        return ",".join(suffixes)
    return ""


# Column names worth trying for the antigen/source roles before falling
# back to asking -- Enterobase's export headers, which is what this tool
# is most often pointed at.
KNOWN_O_ANTIGEN_COLS = ("O Antigen(Serotype Prediction)",)
KNOWN_H_ANTIGEN_COLS = ("H Antigen(Serotype Prediction)",)
KNOWN_SOURCE_COLS = ("Source Type", "Source Niche")


def score_join_columns(path, header, genome_ids):
    """column -> how many genome_ids its values match, for every column at
    once. One pass over the file rather than one per column, since a
    metadata table can be 70+ columns wide.

    Returns (raw_scores, trimmed_scores, values). Raw is compared
    un-stripped, exactly as load_metadata() keys them, so a column's raw
    score is the number of genomes that will really join. The trimmed
    score exists only to tell "wrong column" apart from "right column,
    stray whitespace" when the raw score is zero."""
    wanted = set(genome_ids)
    seen = {c: set() for c in header}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter=sniff_delimiter(path)):
            for col in header:
                value = row.get(col)
                if value:
                    seen[col].add(value)
    raw = {c: len(v & wanted) for c, v in seen.items()}
    trimmed = {c: len({x.strip() for x in v if x.strip()} & wanted)
               for c, v in seen.items()}
    return raw, trimmed, seen


def _pick_known(header, candidates, configured=None):
    """The configured column if it's still in this header, else the first
    recognised name, else None."""
    if configured and configured in header:
        return configured
    for name in candidates:
        if name in header:
            return name
    return None


def _setup_metadata(genome_ids, defaults=None):
    """Values for the five metadata_* keys, or blanks if skipped.

    The join column is proposed by measuring every column's actual overlap
    with the genome IDs rather than by matching a name, so it works on any
    metadata table and the proposal is evidence rather than a guess. The
    user still confirms or overrides it."""
    defaults = defaults or {}
    print("\n--- 3/6: genome metadata TSV (optional) ---")
    print("Adds serotype/host row grouping and the client-side filter boxes.")
    blank = {k: "" for k in ("metadata_tsv", "metadata_join_col",
                             "metadata_o_antigen_col", "metadata_h_antigen_col",
                             "metadata_source_col")}
    while True:
        path = _prompt_path("Metadata table (TSV or CSV)", optional=True,
                             default=defaults.get("metadata_tsv"))
        if path is None:
            print("  skipped -- rows won't be grouped and the filter boxes stay empty")
            return blank
        header = _read_header(path)
        if not header:
            print("  that file has no readable header row")
            continue

        print("  (scoring each column against the genome IDs...)")
        scores, trimmed_scores, seen = score_join_columns(path, header, genome_ids)
        best = max(scores, key=lambda c: scores[c]) if scores else None
        if best and scores[best]:
            print("  " + _ok(f"{best!r} joins {scores[best]} of {len(genome_ids)} genomes"))
            proposed = best
        else:
            # Nothing joined on ANY column. Two very different causes, and
            # the user has to be told which: stray whitespace in an
            # otherwise-correct column, or (far more often) the stem-suffix
            # decision in step 2 leaving the genome IDs unstripped.
            proposed = _pick_known(header, (), defaults.get("metadata_join_col"))
            print("  WARNING: no column in this file joins to the genome IDs.")
            # Test the two hypotheses that actually explain a total miss,
            # in order of likelihood, so the message names a cause and a
            # fix instead of a column picked out of header order.
            probe = None
            suffixes_now = detect_stem_suffixes(sorted(genome_ids))
            if suffixes_now:
                stripped_ids = set()
                for gid in genome_ids:
                    for suf in suffixes_now:
                        if gid.endswith(suf):
                            gid = gid[: -len(suf)]
                            break
                    stripped_ids.add(gid)
                sfx = {c: len(v & stripped_ids) for c, v in seen.items()}
                best_sfx = max(sfx, key=lambda c: sfx[c])
                if sfx[best_sfx]:
                    print(f"    {best_sfx!r} matches {sfx[best_sfx]} of "
                          f"{len(genome_ids)} genomes once "
                          f"{'/'.join(suffixes_now)} is stripped from the stems.")
                    print("    That is the stem-suffix question in step 2 -- answer 'n'")
                    print("    here, re-run, and accept the suffix stripping.")
                    proposed = best_sfx
                    probe = best_sfx
            best_trimmed = (max(trimmed_scores, key=lambda c: trimmed_scores[c])
                            if trimmed_scores else None)
            if probe is not None:
                pass
            elif best_trimmed and trimmed_scores[best_trimmed]:
                # Naming the best-*trimmed* column matters: with every raw
                # score at zero, "best raw" is just whatever came first in
                # the header, which would point at an unrelated column.
                print(f"    {best_trimmed!r} matches {trimmed_scores[best_trimmed]} of "
                      f"{len(genome_ids)} genomes after trimming whitespace.")
                print("    The viewer does not trim this join key, so it will still not")
                print("    join -- fix the whitespace in the TSV.")
                proposed = best_trimmed
                probe = best_trimmed
            else:
                probe = header[0]
                print("    No column matches even after trimming whitespace or stripping")
                print("    stem suffixes -- this may simply be the wrong file.")
            sample_vals = sorted(seen.get(probe) or [])[:3]
            print(f"    genome IDs look like:            {sorted(genome_ids)[:3]}")
            print(f"    closest column {probe!r} looks like: {sample_vals}")
            if not _prompt_yes_no("  Continue picking a column anyway?", default=False):
                continue

        # Confirm the measured winner rather than opening with a 70-line
        # menu; the full list is one 'n' away for the case it got it wrong.
        if proposed and scores.get(proposed):
            if _prompt_yes_no(f"  Use {proposed!r} as the genome ID column?",
                               default=True):
                join_col = proposed
            else:
                join_col = _prompt_column("Which column holds the genome ID "
                                          "(matches the Panaroo/Prokka stems)?",
                                          header, default=proposed)
        else:
            join_col = _prompt_column("Which column holds the genome ID "
                                      "(matches the Panaroo/Prokka stems)?",
                                      header, default=proposed)
        # Keys are collected exactly as load_metadata() collects them --
        # un-stripped -- so this count is what the viewer will really see.
        # The stripped set is kept only to tell the two cases apart when
        # the raw join fails: wrong column vs. right column with stray
        # whitespace (Enterobase exports carry trailing spaces).
        keys = seen.get(join_col, set())
        trimmed = {v.strip() for v in keys if v.strip()}
        wanted = set(genome_ids)
        overlap = keys & wanted
        summary = (f"{len(overlap)} of {len(genome_ids)} viewable genomes matched "
                   f"({len(keys)} distinct values in {join_col!r})")
        print("  " + (_ok(summary) if overlap else summary))
        whitespace_only = (trimmed & wanted) - overlap
        if whitespace_only:
            print(f"  WARNING: another {len(whitespace_only)} value(s) match only after "
                  "trimming whitespace. The viewer does not trim this join key, so those "
                  "rows will NOT join -- fix the whitespace in the TSV.")
        if not overlap:
            example = next(iter(sorted(keys)), "(empty)")
            print(f"  WARNING: nothing joined. {join_col!r} looks like {example!r}, "
                  f"but genome IDs look like {sorted(genome_ids)[0]!r}.")
            if not _prompt_yes_no("  Keep this anyway?", default=False):
                continue
        elif len(overlap) < WEAK_JOIN_FRACTION * len(genome_ids):
            if not _prompt_yes_no("  That's a very partial join. Keep it?", default=False):
                continue

        print("\nThe O and H antigen columns are combined into one serotype used for")
        print("row grouping and coloring -- both are needed, or neither.")
        o_default = _pick_known(header, KNOWN_O_ANTIGEN_COLS,
                                 defaults.get("metadata_o_antigen_col"))
        h_default = _pick_known(header, KNOWN_H_ANTIGEN_COLS,
                                 defaults.get("metadata_h_antigen_col"))
        if o_default and h_default:
            print(f"  found {o_default!r} and {h_default!r}")
            if _prompt_yes_no("  Use these?", default=True):
                o_col, h_col = o_default, h_default
            else:
                o_col = _prompt_column("O antigen column?", header, optional=True,
                                        default=o_default)
                h_col = _prompt_column("H antigen column?", header, optional=True,
                                        default=h_default)
        else:
            # One recognised column is still a useful starting point, so
            # say which was found and offer it rather than reporting none.
            found = o_default or h_default
            if found:
                print(f"  found {found!r}, but both O and H are needed for the "
                      "computed serotype.")
            else:
                print("  no recognised O/H antigen columns in this file.")
            if _prompt_yes_no("  Pick them manually?", default=bool(found)):
                o_col = _prompt_column("O antigen column?", header, optional=True,
                                        default=o_default)
                h_col = _prompt_column("H antigen column?", header, optional=True,
                                        default=h_default)
            else:
                o_col = h_col = None
        if bool(o_col) != bool(h_col):
            print("  only one of the two given -- skipping the computed serotype column")
            o_col = h_col = None

        source_default = _pick_known(header, KNOWN_SOURCE_COLS,
                                      defaults.get("metadata_source_col"))
        source_prompt = "Host/source column shown in each row label?"
        if source_default:
            print(f"\n  found {source_default!r} for the host/source row label")
            if _prompt_yes_no("  Use it?", default=True):
                source_col = source_default
            else:
                source_col = _prompt_column(source_prompt, header, optional=True,
                                             default=source_default)
        else:
            print("\n  no recognised host/source column in this file.")
            if _prompt_yes_no("  Pick one manually?", default=False):
                source_col = _prompt_column(source_prompt, header, optional=True)
            else:
                source_col = None
        return {
            "metadata_tsv": str(path),
            "metadata_join_col": join_col,
            "metadata_o_antigen_col": o_col or "",
            "metadata_h_antigen_col": h_col or "",
            "metadata_source_col": source_col or "",
        }


def _panaroo_clusters(panaroo_csv):
    """Every cluster name in the Panaroo CSV -- the same first-column pass
    build_search_index() makes at startup, used here only to tell the user
    how much of their RFE file will actually resolve."""
    print("  (scanning Panaroo cluster names...)")
    clusters = set()
    with open(panaroo_csv, newline="") as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            if row:
                clusters.add(row[0])
    return clusters


def _setup_rfe(panaroo_csv, default=None):
    """rfe_features_txt, or blank if skipped. A missing 'feature' column is
    blocking rather than a warning: load_rfe_importance() subscripts
    row["feature"] directly, so accepting one here would just move the
    failure to a KeyError at startup."""
    print("\n--- 4/6: RFE / genes-of-interest TSV (optional) ---")
    print("Tab- or comma-delimited, with a 'feature' column of Panaroo cluster")
    print("names, plus optional 'importance' and 'Annotation' columns.")
    while True:
        path = _prompt_path("RFE features table", optional=True, default=default)
        if path is None:
            print("  skipped -- no gene coloring, no RFE badges in search")
            return ""
        header = _read_header(path)
        if not header:
            print("  that file has no readable header row")
            continue
        if "feature" not in header:
            print(f"  ERROR: no 'feature' column. Found: {', '.join(header)}")
            print("  This column is required -- the viewer reads it by name.")
            continue
        if "importance" in header:
            print("  'importance' column: yes")
        else:
            print("  'importance' column: no (features get colored but not ranked)")
        if "Annotation" in header:
            print("  'Annotation' column: yes")
        elif "annotation" in header:
            print("  WARNING: found lowercase 'annotation' -- the viewer only reads "
                  "capital-A 'Annotation', so this column will be ignored.")
        else:
            print("  'Annotation' column: no (Prokka annotations will be used)")

        clusters = _panaroo_clusters(panaroo_csv)
        # Same un-stripped keying as load_rfe_importance() -- see the note
        # in _setup_metadata() about why the trimmed set is tracked too.
        features, trimmed = set(), set()
        with open(path, newline="") as f:
            for row in csv.DictReader(f, delimiter=sniff_delimiter(path)):
                value = row.get("feature")
                if not value:
                    continue
                features.add(value)
                if value.strip():
                    trimmed.add(value.strip())
        overlap = features & clusters
        summary = f"{len(overlap)} of {len(features)} features matched a Panaroo cluster"
        print("  " + (_ok(summary) if overlap else summary))
        whitespace_only = (trimmed & clusters) - overlap
        if whitespace_only:
            print(f"  WARNING: another {len(whitespace_only)} feature(s) match only after "
                  "trimming whitespace. The viewer does not trim feature names, so those "
                  "genes will NOT be colored -- fix the whitespace in the file.")
        if not overlap:
            example = next(iter(sorted(features)), "(empty)")
            print(f"  WARNING: nothing matched -- {example!r} is not a cluster name "
                  "in this pangenome.")
            if not _prompt_yes_no("  Keep this file anyway?", default=False):
                continue
        return str(path)


def _setup_mge(default=None):
    """mge_genes_tsv, or blank if skipped."""
    print("\n--- 5/6: MGE-type TSV (optional) ---")
    print("Tab-delimited geNomad-style table with pangenome_gene_cluster/genome_id/")
    print("locus_tag/mge_type columns.")
    required = ["pangenome_gene_cluster", "genome_id", "locus_tag", "mge_type"]
    while True:
        path = _prompt_path("MGE genes TSV", optional=True, default=default)
        if path is None:
            print("  skipped -- no plasmid/virus/provirus labels")
            return ""
        header = _read_header(path)
        if not header:
            print("  that file has no readable header row")
            continue
        missing = [c for c in required if c not in header]
        if missing:
            print(f"  ERROR: missing required column(s): {', '.join(missing)}")
            print(f"  Found: {', '.join(header)}")
            continue
        print("  " + _ok("all four required columns present"))
        return str(path)


def _setup_stxtyper(genome_ids, default=None):
    """stxtyper_tsv, or blank if skipped. The five columns checked here are
    the ones load_stx_calls()/assign_stx_loci() read by name, so a missing
    one is blocking rather than a warning -- same reasoning as the RFE
    'feature' column: accepting the file here would only move the failure
    to startup.

    The overlap reported is against the genome IDs, because that join
    (stxtyper's --name value vs. the suffix-stripped Panaroo/Prokka stem)
    is the one thing that silently produces "no genome has stx" when the
    report itself is perfectly fine."""
    print("\n--- 6/6: stxtyper report (optional) ---")
    print("Tab- or comma-delimited stxtyper output, one row per operon call, with")
    print("'#name'/stx_type/target_contig/target_start/target_stop columns -- normally")
    print("the per-strain reports concatenated into one combined_stxtyper.tsv.")
    required = ("stx_type", "target_contig", "target_start", "target_stop")
    while True:
        path = _prompt_path("stxtyper report", optional=True, default=default)
        if path is None:
            print("  skipped -- no stx subtype labels, no stx filter/label column")
            return ""
        header = _read_header(path)
        if not header:
            print("  that file has no readable header row")
            continue
        name_col = _stx_name_col(header)
        missing = ([] if name_col else ["#name"]) + [c for c in required if c not in header]
        if missing:
            print(f"  ERROR: missing required column(s): {', '.join(missing)}")
            print(f"  Found: {', '.join(header)}")
            print("  These columns are required -- the viewer reads them by name.")
            continue
        names, n_calls = set(), 0
        with open(path, newline="") as f:
            for row in csv.DictReader(f, delimiter=sniff_delimiter(path)):
                if not (row.get("stx_type") or "").strip():
                    continue
                value = (row.get(name_col) or "").strip()
                if not value:
                    continue
                names.add(value)
                n_calls += 1
        overlap = names & set(genome_ids)
        summary = (f"{len(overlap)} of {len(genome_ids)} genomes matched a "
                    f"{name_col!r} value ({n_calls} operon calls from {len(names)} strains)")
        print("  " + (_ok(summary) if overlap else summary))
        if not overlap:
            example = next(iter(sorted(names)), "(empty)")
            print(f"  WARNING: nothing matched -- {example!r} is not one of this "
                  "pangenome's genome IDs. stxtyper's --name has to carry the same "
                  "strain name as the Panaroo/Prokka stem (after the suffix "
                  "stripping chosen above).")
            if not _prompt_yes_no("  Keep this file anyway?", default=False):
                continue
        return str(path)


def render_config(values):
    """The example file with each 'key = value' line rewritten in place, so
    the generated config keeps every comment block explaining what the keys
    mean. Keys the example doesn't mention are appended at the end."""
    example_path = BASE / CONFIG_EXAMPLE_FILENAME
    if not example_path.exists():
        return "".join(f"{k} = {values.get(k, '')}\n" for k in CONFIG_KEYS)

    lines = example_path.read_text().splitlines()
    seen = set()
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or "=" not in stripped:
            continue
        key = stripped.partition("=")[0].strip()
        if key in values:
            lines[i] = f"{key} = {values[key]}".rstrip()
            seen.add(key)

    extra = [k for k in CONFIG_KEYS if k not in seen and values.get(k)]
    if extra:
        lines += ["", "# --- added by --setup ---"]
        lines += [f"{k} = {values[k]}" for k in extra]
    return "\n".join(lines) + "\n"


def run_setup(config_path):
    """Walk the six required/optional inputs, validate each against the
    ones already given, and write config_path. Returns nothing -- the
    caller exits afterwards, since the settings the user just chose are
    read back from disk on the next run like any other config file."""
    print("=" * 70)
    print("pangenome_viewer setup")
    print("=" * 70)
    print(f"This will write {config_path}")
    print("Paths are tab-completable. Ctrl-C aborts without writing anything.")

    # An existing config becomes the source of defaults rather than just a
    # thing to refuse to clobber, so re-running to change one setting is a
    # walk through with Enter held down.
    defaults = {}
    if config_path.exists():
        defaults = load_config_file(config_path)
        print(f"\n{config_path} already exists.")
        print("Its current values are offered as defaults below -- press Enter to keep")
        print(f"one, or type {CLEAR_TOKEN!r} at an optional path to clear it.")
        if not _prompt_yes_no("Continue and overwrite it when done?", default=False):
            sys.exit("setup cancelled -- existing config left untouched")

    _enable_path_completion()

    panaroo_csv, stems = _setup_panaroo(defaults.get("panaroo_csv"))
    prokka_dir, matched = _setup_prokka(stems, defaults.get("prokka_dir"))
    strip_suffixes = _setup_strip_suffixes(
        matched, default=defaults.get("strip_genome_suffixes"))

    # The IDs the viewer will actually key metadata on: matched stems with
    # the cosmetic suffix removed, exactly as build_genome_stem_map() does.
    suffixes = tuple(s.strip() for s in strip_suffixes.split(",") if s.strip())
    genome_ids = []
    for stem in matched:
        for suffix in suffixes:
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        genome_ids.append(stem)

    values = {k: "" for k in CONFIG_KEYS}
    values["panaroo_csv"] = str(panaroo_csv)
    values["prokka_dir"] = str(prokka_dir)
    values["strip_genome_suffixes"] = strip_suffixes
    values.update(_setup_metadata(genome_ids, defaults))
    values["rfe_features_txt"] = _setup_rfe(panaroo_csv,
                                             defaults.get("rfe_features_txt"))
    values["mge_genes_tsv"] = _setup_mge(defaults.get("mge_genes_tsv"))
    values["stxtyper_tsv"] = _setup_stxtyper(genome_ids, defaults.get("stxtyper_tsv"))

    # Carry through only the settings this walkthrough never asks about,
    # so an existing test-score table isn't silently dropped. Scoped to
    # exactly those keys rather than "any blank value", because a blank
    # elsewhere is a deliberate clear and must not be undone here.
    for key in UNPROMPTED_KEYS:
        if defaults.get(key):
            values[key] = defaults[key]

    config_path.write_text(render_config(values))

    print("\n" + "=" * 70)
    print(f"Wrote {config_path}")
    for key in CONFIG_KEYS:
        if values[key]:
            print(f"  {key} = {values[key]}")
    print(f"\n{len(matched)} genomes are viewable. Start the viewer with:")
    print("  python3 pangenome_viewer.py")
    print("=" * 70)


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
    """feature (Panaroo cluster name, possibly '~~~'-merged) -> {"importance":
    str, "annotation": str}. Empty if RFE_FEATURES_TXT wasn't supplied --
    RFE-based coloring just turns itself off, everything renders as a plain
    neutral gene. "annotation" comes from an optional 'Annotation' column
    (a curated short name, e.g. "nleC") that -- when present for a given
    feature -- supersedes the Prokka annotation in the legend key, tooltip,
    and search results. Both "importance" and "annotation" degrade
    independently and silently to "" when their column is missing from the
    file or blank on a given row, so older two-column RFE files (just
    'feature'/'importance') keep working exactly as before, falling back to
    the Prokka annotation everywhere."""
    if RFE_FEATURES_TXT is None:
        return {}
    out = {}
    delim = sniff_delimiter(RFE_FEATURES_TXT)
    with open(RFE_FEATURES_TXT, newline="") as f:
        for row in csv.DictReader(f, delimiter=delim):
            out[row["feature"]] = {
                "importance": (row.get("importance") or "").strip(),
                "annotation": (row.get("Annotation") or "").strip(),
            }
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
            rfe_info = rfe_importance.get(cluster) or {}
            index.append({
                "cluster": cluster,
                "non_unique_name": row[1],
                "annotation": row[2],
                "carrier_count": carrier_count,
                "is_rfe": cluster in rfe_importance,
                "rfe_importance": rfe_info.get("importance") or None,
                "rfe_annotation": rfe_info.get("annotation") or None,
                "_haystack": f"{cluster} {row[1]} {row[2]}".lower(),
            })
    print(f"[startup] indexed {len(index)} pangenome clusters across "
          f"{n_genomes} genomes", file=sys.stderr)
    return index


SYNTHETIC_SEROTYPE_COL = "Serotype (O:H combined, computed)"

# A metadata table is often keyed by the raw assembly filename stem while
# genome_id has had STRIP_GENOME_SUFFIXES taken off it, so a genome's row
# has to be looked for under both spellings. Hardcoded rather than derived
# from STRIP_GENOME_SUFFIXES because the two settings answer different
# questions: what the display ID looks like vs. what the metadata file
# happens to be keyed by.
METADATA_KEY_SUFFIXES = (".result", ".scaffold")


def metadata_entry(meta, genome_id):
    """A genome's metadata row, or None -- the one place the genome_id ->
    metadata-key fallback lives, so anything that writes into a genome's
    metadata dict lands on the same row a chart build later reads out of."""
    m = meta.get(genome_id)
    if m:
        return m
    for suf in METADATA_KEY_SUFFIXES:
        m = meta.get(genome_id + suf)
        if m:
            return m
    return None


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
        r = csv.DictReader(f, delimiter=sniff_delimiter(METADATA_TSV))
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


# stxtyper reports coordinates against the ORIGINAL assembly contig name
# (SPAdes-style, "NODE_107_length_3367_cov_5.65253_ID_213"), which Prokka
# has since renamed to "gnl|Prokka|HJIKEAHO_33". The contig LENGTH is the
# only thing the two names still share, and the assembler wrote it into
# the name -- this is what digs it back out.
_STX_CONTIG_LEN_RE = re.compile(r"_length_(\d+)")

SYNTHETIC_STX_COL = "stx type (stxtyper)"


def _stx_name_col(fieldnames):
    """stxtyper's strain-name column, written '#name' -- its header line
    doubles as a comment line, so the '#' is part of the field name that
    csv.DictReader hands back. Matched by name (with any leading '#'
    ignored) rather than by position so a build that drops the '#', or
    reorders the columns, still reads."""
    for col in fieldnames or ():
        if col.lstrip("#").strip().lower() == "name":
            return col
    return None


def _stx_int(value):
    """A coordinate, or None if the file has 'NA'/blank/garbage there --
    stxtyper leaves fields unset on partial calls, and a call with no
    usable coordinates simply can't be placed on a locus."""
    try:
        return int((value or "").strip())
    except ValueError:
        return None


def load_stx_calls():
    """genome_id -> [call, ...], one entry per operon call in the stxtyper
    report. Empty if STXTYPER_TSV wasn't supplied -- every stx feature
    below then reduces to "no calls anywhere", the same way the other
    optional inputs turn themselves off.

    The join to genome_id is stxtyper's own `--name` value, which the
    run_stxtyper.sh convention sets to the suffix-stripped strain name --
    i.e. already exactly the Panaroo/Prokka stem this viewer keys on, so
    nothing is normalised here. A missing required column is a loud
    startup warning and an empty result rather than a KeyError mid-request:
    the report is optional, so a malformed one must not take the server
    down with it."""
    if STXTYPER_TSV is None:
        return {}
    calls = {}
    with open(STXTYPER_TSV, newline="") as f:
        r = csv.DictReader(f, delimiter=sniff_delimiter(STXTYPER_TSV))
        fieldnames = list(r.fieldnames or [])
        name_col = _stx_name_col(fieldnames)
        missing = [c for c in ("stx_type", "target_contig", "target_start", "target_stop")
                   if c not in fieldnames]
        if name_col is None:
            missing.insert(0, "#name")
        if missing:
            print(f"[startup] WARNING: {STXTYPER_TSV} is missing required column(s) "
                  f"{', '.join(missing)} -- stx subtyping will be off", file=sys.stderr)
            return {}
        for row in r:
            genome = (row.get(name_col) or "").strip()
            stx_type = (row.get("stx_type") or "").strip()
            # A blank stx_type is a row stxtyper emitted without actually
            # calling a subtype; it carries no information to display and
            # must not turn into a phantom "" subtype on a row label.
            if not genome or not stx_type:
                continue
            identity = (row.get("identity") or "").strip()
            m = _STX_CONTIG_LEN_RE.search(row.get("target_contig") or "")
            calls.setdefault(genome, []).append({
                "stx_type": stx_type,
                "operon": (row.get("operon") or "").strip() or None,
                "identity": identity if identity and identity != "NA" else None,
                "contig_len": int(m.group(1)) if m else None,
                "start": _stx_int(row.get("target_start")),
                "stop": _stx_int(row.get("target_stop")),
            })
    return calls


def attach_stx_metadata(meta, columns, genome_ids, stx_types):
    """Fold each genome's stx subtypes into its metadata row as one more
    column, so the existing client-side filter boxes and row-label pickers
    pick stx up for free instead of needing their own parallel plumbing --
    the same trick SYNTHETIC_SEROTYPE_COL plays with the antigen columns.

    The value is the subtypes joined ("stx1a + stx2c"), which is what makes
    an exact-match filter on a *combination* possible and what reads
    sensibly in a row label; the per-row stx_types array stays the source
    of truth for "carries stx2a at all" filtering. A genome with no call
    gets "", the same empty value every other metadata column uses for
    "nothing here", which the client already renders as (blank).

    Genomes absent from the metadata table (or with no metadata table at
    all -- this feature must not require one) get a row created for them,
    with the serotype/source fields sero_for() subscripts filled in as
    Unknown. Returns the updated column list."""
    if STXTYPER_TSV is None:
        return columns
    for genome_id in genome_ids:
        entry = metadata_entry(meta, genome_id)
        if entry is None:
            entry = {"source_type": "Unknown", "serotype": "Unknown", "full": {}}
            meta[genome_id] = entry
        entry["full"][SYNTHETIC_STX_COL] = " + ".join(stx_types.get(genome_id, ()))
    if SYNTHETIC_STX_COL in columns:
        return columns
    return list(columns) + [SYNTHETIC_STX_COL]


def assign_stx_loci(calls, seqlens, genes):
    """locus_tag -> the stx call it belongs to, for ONE genome, plus a
    Counter of what happened to each call.

    The chain: parse the contig length out of stxtyper's assembler contig
    name, look that length up among the GFF's ##sequence-region lines to
    recover Prokka's renamed contig, then take the CDSs on it that overlap
    the reported operon coordinates -- those are the stxA/stxB genes.

    Length is not a unique key (~6% of genomes here have two contigs the
    same length), so a tie is only broken by evidence: accept it when
    exactly one candidate contig actually has CDSs across the reported
    range, and otherwise drop the per-locus assignment entirely rather
    than guess. A wrong assignment would print "stx2a" on some unrelated
    gene, which is worse than printing nothing; the genome-level subtype
    set is read straight from the report, so it stays right regardless.

    Outcome counts: "resolved" (loci assigned), "no_cds" (contig found,
    nothing annotated over the range -- routine for the tiny
    PARTIAL_CONTIG_END contigs), "ambiguous" (several equally-supported
    candidate contigs), "no_contig" (no contig of that length in this GFF
    at all), "no_coords" (unusable coordinates in the report)."""
    mapping = {}
    stats = Counter()
    if not calls:
        return mapping, stats
    by_len = {}
    for contig, length in seqlens.items():
        by_len.setdefault(length, []).append(contig)
    # One pass over this genome's ~5k CDSs, keeping only the handful on
    # contigs some call could land on -- the alternative, re-scanning every
    # gene per call, is the same work multiplied by the operon count.
    wanted = set()
    for call in calls:
        wanted.update(by_len.get(call["contig_len"], ()))
    contig_genes = {}
    for locus_tag, (contig, start, end, _strand, _product) in genes.items():
        if contig in wanted:
            contig_genes.setdefault(contig, []).append((locus_tag, start, end))

    for call in calls:
        lo, hi = call["start"], call["stop"]
        if lo is None or hi is None:
            stats["no_coords"] += 1
            continue
        candidates = by_len.get(call["contig_len"], [])
        if not candidates:
            stats["no_contig"] += 1
            continue
        hits = {}
        for contig in candidates:
            overlapping = [lt for lt, s, e in contig_genes.get(contig, ()) if s <= hi and e >= lo]
            if overlapping:
                hits[contig] = overlapping
        if len(hits) == 1:
            for locus_tag in next(iter(hits.values())):
                mapping[locus_tag] = call
            stats["resolved"] += 1
        elif not hits:
            stats["no_cds"] += 1
        else:
            stats["ambiguous"] += 1
    return mapping, stats


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


def load_pan_genome_reference():
    """cluster -> representative CDS nucleotide sequence, from Panaroo's own
    pan_genome_reference.fa -- a real sequence Panaroo copied from one
    member genome during clustering (not a synthetic consensus), one per
    pangenome cluster, headers matching the same 'group_XXXX' (or gene
    name) convention used throughout this tool. Not a separate config
    setting: assumed to sit where Panaroo always puts it, next to
    PANAROO_CSV. Empty (feature silently off, same pattern as every other
    optional input) if that sibling file isn't there."""
    path = PANAROO_CSV.parent / "pan_genome_reference.fa"
    if not path.exists():
        return {}
    out = {}
    cluster = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cluster is not None:
                    out[cluster] = "".join(parts)
                cluster = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if cluster is not None:
        out[cluster] = "".join(parts)
    return out


# Standard genetic code (NCBI table 11, bacterial -- identical to table 1
# except for a couple of alternate start codons, handled separately below).
# Prokka/Prodigal-called CDSs are already in-frame from the first base, so
# translation here is just "read codons from position 0 until a stop".
_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
_ALT_START_CODONS = {"GTG", "TTG"}  # translate as Met when they're the first codon


def translate_dna(seq):
    """CDS nucleotide -> protein, stopping at (and dropping) the first stop
    codon. Trailing partial codon, if any, is ignored. An unrecognized
    codon (stray ambiguity code) becomes 'X' rather than aborting."""
    seq = seq.upper()
    aa = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i:i + 3]
        if i == 0 and codon in _ALT_START_CODONS:
            aa.append("M")
            continue
        res = _CODON_TABLE.get(codon, "X")
        if res == "*":
            break
        aa.append(res)
    return "".join(aa)


# Non-CDS feature types Prokka emits that are worth showing as a separate map
# layer -- structural/non-coding RNAs and repeats that otherwise render as
# blank track (the main gene layer draws CDS only). "gene" is skipped: it is
# the parent wrapper Prokka pairs with each CDS/tRNA/rRNA, not a feature.
NON_CDS_TYPES = {"tRNA", "rRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region", "CRISPR"}


def parse_gff(stem):
    """Returns (genes, seqlens, non_cds). genes: locus_tag -> (contig, start,
    end, strand, product). seqlens: contig -> length (from ##sequence-region).
    non_cds: list of (contig, start, end, strand, ftype, label) for the
    non-coding / repeat features in NON_CDS_TYPES -- read from the same pass
    that already scans every annotation line, so it costs no extra I/O."""
    genes = {}
    seqlens = {}
    non_cds = []
    gff = PROKKA_BASE / stem / f"{stem}.gff"
    if not gff.exists():
        return genes, seqlens, non_cds
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
            if len(fields) < 9:
                continue
            ftype = fields[2]
            if ftype == "CDS":
                m = re.search(r"locus_tag=([^;]+)", fields[8])
                if not m:
                    continue
                pm = re.search(r"product=([^;]+)", fields[8])
                product = urllib.parse.unquote(pm.group(1)) if pm else "hypothetical protein"
                genes[m.group(1)] = (fields[0], int(fields[3]), int(fields[4]), fields[6], product)
            elif ftype in NON_CDS_TYPES:
                pm = re.search(r"product=([^;]+)", fields[8])
                if pm:
                    label = urllib.parse.unquote(pm.group(1))
                else:
                    nm = re.search(r"(?:note|rpt_family|Name)=([^;]+)", fields[8])
                    label = urllib.parse.unquote(nm.group(1)) if nm else ftype
                non_cds.append((fields[0], int(fields[3]), int(fields[4]), fields[6], ftype, label))
    return genes, seqlens, non_cds


# Runs of N shorter than this are single-base ambiguity calls, not scaffold
# gaps -- 33% of runs in this cohort are 1-9bp. 10bp filters those out while
# keeping every genuine assembly gap (observed gaps range ~10bp to ~1kb).
MIN_N_RUN = 10


def extract_contig_fasta(stem, contig):
    """Full sequence of ONE contig, read from the ##FASTA block at the tail
    of a genome's Prokka GFF (Prokka's combined-output default). We read
    only the requested contig's records and stop -- because the anchor gene
    almost always sits on the first (largest) contig, this breaks out after
    a few KB rather than reading the whole ~5.5MB file. Returns "" if the
    GFF is missing, has no ##FASTA block, or doesn't contain this contig."""
    gff = PROKKA_BASE / stem / f"{stem}.gff"
    if not gff.exists():
        return ""
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
    return "".join(parts)


def _extract_locus_record(stem, locus_tag, ext):
    """One record's sequence, read straight out of Prokka's <stem>.<ext> by
    exact locus_tag match on the header's first whitespace-delimited token.
    Shared by extract_faa_record (ext='faa', protein) and
    extract_ffn_record (ext='ffn', per-gene nucleotide) -- neither file is
    required for the rest of this tool (only <stem>.gff is), so this is a
    soft dependency: returns "" if the file is missing, or the locus_tag
    isn't in it."""
    path = PROKKA_BASE / stem / f"{stem}.{ext}"
    if not path.exists():
        return ""
    parts = []
    in_target = False
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if in_target:
                    break  # captured the whole target record; stop reading
                in_target = line[1:].split()[0] == locus_tag
                continue
            if in_target:
                parts.append(line.strip())
    return "".join(parts)


def extract_faa_record(stem, locus_tag):
    """One protein sequence, read straight out of Prokka's own <stem>.faa --
    already translated by Prokka/Prodigal (correct start codon, correct
    genetic code), so this is preferred over re-translating from the GFF
    nucleotide whenever a per-strain (as opposed to Panaroo's pangenome
    reference) amino acid sequence is wanted."""
    return _extract_locus_record(stem, locus_tag, "faa")


def extract_ffn_record(stem, locus_tag):
    """One gene's own nucleotide (CDS) sequence, read straight out of
    Prokka's <stem>.ffn -- the per-strain nucleotide counterpart to
    extract_faa_record's protein."""
    return _extract_locus_record(stem, locus_tag, "ffn")


def find_n_runs(stem, contig, min_len=MIN_N_RUN):
    """Scaffold-gap N-runs on ONE contig of a genome's assembly. Returns a
    list of (start, end) 1-based inclusive coordinates for each run of >=
    min_len N's, in contig coordinates."""
    seq = extract_contig_fasta(stem, contig)
    if not seq:
        return []
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

        # stx subtypes are read once here, keyed by genome, and are the
        # authority for a row's stx_types -- a row label has to name the
        # right subtype whether or not any stx gene lands in the anchor's
        # window. The per-locus assignment (assign_stx_loci) is a separate,
        # per-build concern. Injected into the metadata dicts before the
        # value totals below are computed so the filter boxes' "% of X
        # genomes carry this gene" denominators cover the stx column too;
        # this deliberately does not depend on a metadata table existing.
        self.stx_calls = load_stx_calls()
        self.stx_types = {g: sorted({c["stx_type"] for c in calls})
                          for g, calls in self.stx_calls.items()}
        self.metadata_columns = attach_stx_metadata(
            self.metadata, self.metadata_columns, self.genome_stem_map, self.stx_types)
        if STXTYPER_TSV is not None:
            n_typed = sum(1 for g in self.genome_stem_map if self.stx_types.get(g))
            print(f"[startup] {sum(len(v) for v in self.stx_calls.values())} stx operon "
                  f"calls across {len(self.stx_calls)} strains; {n_typed}/"
                  f"{len(self.genome_stem_map)} viewable genomes typed", file=sys.stderr)

        self.metadata_value_totals = compute_metadata_value_totals(self.metadata, self.metadata_columns)
        self.test_scores = load_test_scores()

        print("[startup] loading pangenome reference sequences...", file=sys.stderr)
        self.pan_ref_nt = load_pan_genome_reference()
        self.pan_ref_aa_cache = {}  # cluster -> translated AA, filled lazily per request
        print(f"[startup] {len(self.pan_ref_nt)} reference sequences "
              f"{'loaded' if self.pan_ref_nt else '(pan_genome_reference.fa not found -- feature off)'}",
              file=sys.stderr)

        print("[startup] loading cluster byte-offset index (single-cluster "
              "lookups: sequence panel, contig export)...", file=sys.stderr)
        self.cluster_offsets = load_cluster_offsets()
        if len(self.cluster_offsets) != len(self.cluster_lookup):
            print(f"[startup] cluster_offsets.json looks stale "
                  f"({len(self.cluster_offsets)} vs {len(self.cluster_lookup)} "
                  f"clusters) -- rebuilding...", file=sys.stderr)
            self.cluster_offsets = build_cluster_offset_index()
        print(f"[startup] {len(self.cluster_offsets)} cluster offsets ready", file=sys.stderr)

        if MGE_GENES_TSV is not None:
            print("[startup] loading MGE-gene byte-offset index...", file=sys.stderr)
        self.mge_offsets = load_mge_offsets()
        if MGE_GENES_TSV is not None:
            print(f"[startup] {sum(len(v) for v in self.mge_offsets.values())} MGE rows "
                  f"across {len(self.mge_offsets)} clusters ready", file=sys.stderr)

        print("[startup] ready.", file=sys.stderr)


# ---------------------------------------------------------------------------
# Chart builder.
# ---------------------------------------------------------------------------

MGE_OFFSETS_PATH = BASE / "mge_offsets.json"


def build_mge_offset_index():
    """One pass over MGE_GENES_TSV, recording every row's start byte offset
    under its pangenome_gene_cluster. Unlike the Panaroo CSV (exactly one
    row per cluster), this file has one row per (genome, locus) instance
    and is ordered by genome, not by cluster -- so a cluster's rows are
    scattered throughout the file, and this maps to a LIST of offsets per
    cluster rather than a single one. Persisted to MGE_OFFSETS_PATH
    (gitignored) alongside the source file's mtime, so a later load can
    detect staleness (a repointed/regenerated mge_genes_tsv) without first
    re-scanning to find out."""
    if MGE_GENES_TSV is None:
        return {}
    offsets = {}
    with open(MGE_GENES_TSV, newline="") as f:
        header_line = f.readline()
        cluster_col = header_line.rstrip("\n").split("\t").index("pangenome_gene_cluster")
        pos = f.tell()
        line = f.readline()
        while line:
            parts = line.split("\t")
            if len(parts) > cluster_col:
                offsets.setdefault(parts[cluster_col], []).append(pos)
            pos = f.tell()
            line = f.readline()
    MGE_OFFSETS_PATH.write_text(json.dumps({
        "_source_mtime": MGE_GENES_TSV.stat().st_mtime,
        "offsets": offsets,
    }))
    n_rows = sum(len(v) for v in offsets.values())
    print(f"[startup] indexed {n_rows} MGE rows across {len(offsets)} clusters "
          f"-> {MGE_OFFSETS_PATH.name}", file=sys.stderr)
    return offsets


def load_mge_offsets():
    """MGE_OFFSETS_PATH's cached index if its recorded source mtime still
    matches MGE_GENES_TSV, else a fresh build. Empty (not an error) if
    MGE_GENES_TSV wasn't supplied at all -- load_mge_types_for_cluster
    already handles that same way."""
    if MGE_GENES_TSV is None:
        return {}
    if MGE_OFFSETS_PATH.exists():
        cached = json.loads(MGE_OFFSETS_PATH.read_text())
        if cached.get("_source_mtime") == MGE_GENES_TSV.stat().st_mtime:
            return cached["offsets"]
    return build_mge_offset_index()


_mge_header_cache = None


def load_mge_types_for_cluster(anchor_cluster, offsets=None):
    """(genome_id, locus_tag) -> mge_type, for this anchor's gene instances
    only -- used for copy-selection when a genome carries >1 paralog, and
    for the row-label mge_type field. Empty if MGE_GENES_TSV wasn't
    supplied.

    `offsets` (pass ctx.mge_offsets), when given, turns this into a
    handful of seeks against this cluster's known row offsets instead of a
    full DictReader scan of a several-hundred-MB file -- the difference
    between single-digit seconds and single-digit milliseconds. Falls back
    to the original full-scan behavior when omitted, so any caller that
    doesn't have a live ctx handy still works, just slower."""
    if MGE_GENES_TSV is None:
        return {}
    if offsets is None:
        types = {}
        with open(MGE_GENES_TSV, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row["pangenome_gene_cluster"] != anchor_cluster:
                    continue
                types[(row["genome_id"], row["locus_tag"])] = row["mge_type"]
        return types

    global _mge_header_cache
    positions = offsets.get(anchor_cluster, [])
    if not positions:
        return {}
    if _mge_header_cache is None:
        with open(MGE_GENES_TSV, newline="") as f:
            _mge_header_cache = next(csv.reader(f, delimiter="\t"))
    types = {}
    with open(MGE_GENES_TSV, newline="") as f:
        for pos in positions:
            f.seek(pos)
            row = dict(zip(_mge_header_cache, next(csv.reader(f, delimiter="\t"))))
            types[(row["genome_id"], row["locus_tag"])] = row["mge_type"]
    return types


def pick_locus(loci, mge_types, g, genes=None):
    """Pick a genome's locus tag for a gene cluster when Panaroo lists more
    than one paralog copy for it (semicolon-joined loci in the CSV cell).
    Prefer a real Prokka-called locus over a Panaroo "refound" placeholder
    -- a refound locus has no GFF entry to resolve coordinates from at all
    -- then, among real candidates, prefer one classified in `mge_types`
    (a (genome_id, locus_tag) -> mge_type dict for this anchor cluster).

    Any tie left after that is broken by position, given `genes` (the
    genome's parsed GFF): the copy furthest upstream in its own reading
    direction. Panaroo's order within the cell is arbitrary, so taking the
    first-listed copy anchored a tandem pair on copy 1 in some genomes and
    copy 2 in others, shifting those rows a gene-width apart and turning
    every link between them into a diagonal. Every caller passes `genes`
    so the chart, contig export and sequence panel agree on the copy."""
    if len(loci) == 1:
        return loci[0]
    real = [lt for lt in loci if "refound" not in lt]
    candidates = real if real else loci
    classified = [lt for lt in candidates if (g, lt) in mge_types]
    if classified:
        candidates = classified
    if genes and len(candidates) > 1 and all(lt in genes for lt in candidates):
        # genes[lt] = (contig, start, end, strand, product). On the minus
        # strand the upstream copy is the one with the larger end.
        candidates = sorted(candidates, key=lambda lt: (
            genes[lt][0],
            -genes[lt][2] if genes[lt][3] == "-" else genes[lt][1]))
    return candidates[0]


def build_chart_data(ctx, anchor_cluster, count=DEFAULT_COUNT, include_n_runs=False):
    if anchor_cluster not in ctx.cluster_lookup:
        return {"error": f"cluster {anchor_cluster!r} not found"}

    mge_types = load_mge_types_for_cluster(anchor_cluster, ctx.mge_offsets)

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
        return metadata_entry(ctx.metadata, genome_id)

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
    gff_genes = {g: genes for g, (genes, _s, _n) in gff_parsed.items()}
    gff_seqlens = {g: seqlens for g, (_g, seqlens, _n) in gff_parsed.items()}
    gff_noncds = {g: non_cds for g, (_g, _s, non_cds) in gff_parsed.items()}

    # locus_tag -> stx call, per sampled genome. Needs the parsed GFF above
    # (contig lengths to join on, CDS coordinates to overlap against), so it
    # can't sit up beside load_mge_types_for_cluster despite being the same
    # kind of per-(genome, locus) attribute.
    stx_loci = {}
    stx_stats = Counter()
    if STXTYPER_TSV is not None:
        for g in sample_genome_ids:
            mapping, call_stats = assign_stx_loci(
                ctx.stx_calls.get(g, ()), gff_seqlens.get(g, {}), gff_genes.get(g, {}))
            if mapping:
                stx_loci[g] = mapping
            stx_stats.update(call_stats)

    anchor_locus = {g: pick_locus(present[g], mge_types, g, gff_genes.get(g))
                    for g in sample_genome_ids}

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
            gene_out = {
                "start": rel_start, "end": rel_end, "strand": disp_strand,
                "cluster": cluster, "product": product,
                "is_anchor": glt == lt,
                "is_rfe": is_rfe,
                "rfe_importance": None,
                "rfe_annotation": None,
            }
            # Present only on the stxA/stxB CDSs of a called operon -- every
            # other gene simply doesn't carry these keys, so the client can
            # test for the key rather than for a sentinel value.
            stx_call = stx_loci.get(g, {}).get(glt)
            if stx_call:
                gene_out["stx_type"] = stx_call["stx_type"]
                gene_out["stx_operon"] = stx_call["operon"]
                gene_out["stx_identity"] = stx_call["identity"]
            genes_out.append(gene_out)
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

        # Non-CDS features (rRNA/tRNA/repeat) on the anchor contig, in-window,
        # in the same anchor-relative frame as genes. Always included (nearly
        # free -- already parsed); the client toggles their visibility.
        non_cds_out = []
        for nc_contig, nc_start, nc_end, nc_strand, nc_ftype, nc_label in gff_noncds.get(g, []):
            if nc_contig != contig:
                continue
            if nc_end < win_lo or nc_start > win_hi:
                continue
            if flip:
                nc_rel_start, nc_rel_end = anchor_ref - nc_end, anchor_ref - nc_start
                nc_disp_strand = "-" if nc_strand == "+" else "+"
            else:
                nc_rel_start, nc_rel_end = nc_start - anchor_ref, nc_end - anchor_ref
                nc_disp_strand = nc_strand
            non_cds_out.append({"start": nc_rel_start, "end": nc_rel_end,
                                "strand": nc_disp_strand, "ftype": nc_ftype,
                                "label": nc_label})
        non_cds_out.sort(key=lambda x: x["start"])

        sero, source = sero_for(g)
        meta_row = lookup_meta(g)
        rows_out.append({
            "genome_id": g, "mge_type": mge_type,
            "serotype": sero, "sero_bucket": sero_bucket(g), "source_type": source,
            "stx_types": ctx.stx_types.get(g, []),
            "test_score": ctx.test_scores.get(g),
            "contig_rel_start": contig_rel_start, "contig_rel_end": contig_rel_end,
            "genes": genes_out,
            "n_runs": n_runs_out,
            "non_cds": non_cds_out,
            "_contig": contig,
            "metadata": meta_row["full"] if meta_row else {},
        })

    for row in rows_out:
        for gobj in row["genes"]:
            if gobj["is_rfe"]:
                info = ctx.rfe_importance.get(gobj["cluster"]) or {}
                gobj["rfe_importance"] = info.get("importance") or None
                gobj["rfe_annotation"] = info.get("annotation") or None

    buckets_order = top_seros + (["Other"] if any(r["sero_bucket"] == "Other" for r in rows_out) else [])
    bucket_rank = {b: i for i, b in enumerate(buckets_order)}
    rows_out.sort(key=lambda r: (bucket_rank.get(r["sero_bucket"], 99), r["genome_id"]))

    top_rfe_clusters = [c for c, _ in rfe_cluster_counts.most_common(TOP_N_COLORED)]
    # Legend-key label per cluster: the curated RFE 'Annotation' (e.g. "nleC")
    # when the RFE features file supplies one for this cluster, else the
    # Prokka annotation as before -- the same fallback used in the tooltip
    # and search results below.
    cluster_labels = {}
    for row in rows_out:
        for gobj in row["genes"]:
            c = gobj["cluster"]
            if c in top_rfe_clusters and c not in cluster_labels:
                rfe_ann = (ctx.rfe_importance.get(c) or {}).get("annotation")
                cluster_labels[c] = rfe_ann if rfe_ann else gobj["product"]

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
    if STXTYPER_TSV is not None and stx_stats:
        print(f"[build] {anchor_cluster}: stx calls in these genomes -- "
              + ", ".join(f"{k} {v}" for k, v in sorted(stx_stats.items())), file=sys.stderr)

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
        "stx_enabled": STXTYPER_TSV is not None,
        "carrier_value_counts": carrier_value_counts,
        "rows": rows_out,
    }


def build_cluster_offset_index():
    """One pass over the Panaroo CSV, recording each cluster row's start
    byte offset -- written to CLUSTER_OFFSETS_PATH so future server starts
    (and dev/pgv_lookup.py's own CLI, which points at the same file) can just
    load it instead of rebuilding. Deliberately NOT shared code with
    dev/pgv_lookup.py's own near-identical build_index(): that tool imports
    this module (`import pangenome_viewer as pv`) and is meant to run as
    its own standalone entry point, where its `pv` reference is configured
    via its own _configure(); this module can't import it back the same
    way -- when pangenome_viewer.py itself is what's running (always, for
    the live server), `import pgv_lookup` would re-import this file as a
    second, freshly-reset module object rather than binding to the live
    __main__ one, so PANAROO_CSV etc. on that copy would be None. Small
    enough (~15 lines) that duplicating it here is simpler and more robust
    than fighting that.

    Uses explicit readline() (not `for line in f` / csv.reader on the file
    object) because file.tell() is only reliable after readline() in text
    mode -- iterating a file object internally buffers ahead and silently
    desyncs tell() from the logical line position. Splitting the first
    field on a raw comma (not csv.reader) is safe here because Panaroo
    cluster names never contain a comma or quote character."""
    offsets = {}
    with open(PANAROO_CSV, newline="") as f:
        f.readline()  # header
        pos = f.tell()
        line = f.readline()
        while line:
            cluster = line.split(",", 1)[0]
            offsets[cluster] = pos
            pos = f.tell()
            line = f.readline()
    CLUSTER_OFFSETS_PATH.write_text(json.dumps({
        "_source_stamp": offset_index_stamp(PANAROO_CSV),
        "offsets": offsets,
    }))
    print(f"[startup] indexed {len(offsets)} clusters -> {CLUSTER_OFFSETS_PATH.name}", file=sys.stderr)
    return offsets


def load_cluster_offsets():
    """CLUSTER_OFFSETS_PATH's cached index if it was built from the
    Panaroo CSV as it stands now, else a fresh build.

    An index whose source has changed doesn't just go missing -- its
    offsets point into the middle of other rows, so a lookup silently
    returns a different cluster's data. A cache written before this stamp
    existed has no '_source_stamp' key and is treated as stale, which
    costs one rebuild on first run after upgrading."""
    if CLUSTER_OFFSETS_PATH.exists():
        try:
            cached = json.loads(CLUSTER_OFFSETS_PATH.read_text())
        except (json.JSONDecodeError, OSError):
            cached = None
        if isinstance(cached, dict) and "offsets" in cached:
            if cached.get("_source_stamp") == offset_index_stamp(PANAROO_CSV):
                return cached["offsets"]
            print(f"[startup] {CLUSTER_OFFSETS_PATH.name} was built from a different "
                  f"{PANAROO_CSV.name} -- rebuilding", file=sys.stderr)
        else:
            print(f"[startup] {CLUSTER_OFFSETS_PATH.name} predates source-change "
                  "detection -- rebuilding once", file=sys.stderr)
    return build_cluster_offset_index()


_panaroo_header_cache = None


def get_cluster_row(cluster_name, offsets):
    """(header, row) for one cluster via an instant seek instead of a
    scan. Raises KeyError if cluster_name isn't in `offsets`."""
    global _panaroo_header_cache
    if cluster_name not in offsets:
        raise KeyError(f"{cluster_name!r} not in the cluster offset index "
                        f"(or the index is stale)")
    if _panaroo_header_cache is None:
        with open(PANAROO_CSV, newline="") as f:
            _panaroo_header_cache = next(csv.reader(f))
    with open(PANAROO_CSV, newline="") as f:
        f.seek(offsets[cluster_name])
        row = next(csv.reader(f))
    # Belt and braces behind the staleness check above: if an offset ever
    # does point at the wrong row, fail loudly here rather than return
    # another cluster's data as though it were this one.
    if not row or row[0] != cluster_name:
        raise KeyError(
            f"cluster offset index is stale: offset for {cluster_name!r} landed on "
            f"{(row[0] if row else '(empty line)')!r}. Delete "
            f"{CLUSTER_OFFSETS_PATH.name} (it rebuilds automatically).")
    return _panaroo_header_cache, row


def resolve_anchor_loci(ctx, anchor_cluster, genome_ids):
    """For each of the given genome_ids, this cluster's carrier locus tag(s)
    (semicolon-split Panaroo CSV cell). An instant seek via the byte-offset
    index (ctx.cluster_offsets), not a scan -- this is meant for a handful
    of specific genomes (export-selection checkboxes, a single
    sequence-viewer lookup) triggered directly by a UI click, so it needs
    to feel instant. A full linear scan of a several-hundred-MB Panaroo CSV
    per click was the previous approach here and took several seconds
    every time; build_chart_data's own scan is untouched (see the module
    docstring -- that one stays scan-based on purpose).

    Returns (present, error): present maps genome_id -> [locus_tag, ...]
    for whichever of genome_ids are carriers; error is set (present then
    empty) only if anchor_cluster itself doesn't resolve at all."""
    if anchor_cluster not in ctx.cluster_lookup:
        return {}, f"cluster {anchor_cluster!r} not found"

    try:
        header, row = get_cluster_row(anchor_cluster, ctx.cluster_offsets)
    except KeyError as e:
        return {}, str(e)

    wanted = set(genome_ids)
    present = {}
    for stem, cell in zip(header[3:], row[3:]):
        cell = cell.strip()
        if not cell:
            continue
        g = ctx.stem_to_genome.get(stem)
        if g is None or g not in wanted:
            continue
        present[g] = [x.strip() for x in cell.split(";") if x.strip()]
    return present, None


def build_contig_export(ctx, anchor_cluster, genome_ids):
    """Multi-FASTA text: for each requested genome_id, the WHOLE contig its
    anchor_cluster locus sits on, in raw assembly (GFF) orientation --
    genomes anchored on the '-' strand are NOT reverse-complemented to match
    the chart's display orientation, so this is the sequence exactly as it
    reads in that genome's own assembly.

    Returns (fasta_text, warnings): fasta_text is None if nothing could be
    exported at all; warnings lists one message per genome_id that couldn't
    be resolved (not a carrier, refound-only locus, or no ##FASTA block in
    that genome's GFF), in the caller's requested order."""
    present, error = resolve_anchor_loci(ctx, anchor_cluster, genome_ids)
    if error:
        return None, [error]

    mge_types = load_mge_types_for_cluster(anchor_cluster, ctx.mge_offsets)
    records, warnings = [], []
    for g in genome_ids:  # preserve the caller's (selection) order
        loci = present.get(g)
        if not loci:
            warnings.append(f"{g}: not a carrier of {anchor_cluster}")
            continue
        stem = ctx.genome_stem_map[g]
        genes, _seqlens, _non_cds = parse_gff(stem)
        locus = pick_locus(loci, mge_types, g, genes)
        if locus not in genes:
            warnings.append(f"{g}: locus {locus} has no GFF entry (refound placeholder)")
            continue
        contig = genes[locus][0]
        seq = extract_contig_fasta(stem, contig)
        if not seq:
            warnings.append(f"{g}: contig {contig!r} sequence not found (no ##FASTA block?)")
            continue
        header_line = f">{g} contig={contig} anchor={anchor_cluster} locus={locus} length={len(seq)}"
        wrapped = "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70))
        records.append(f"{header_line}\n{wrapped}")

    if not records:
        return None, warnings or ["no exportable genomes"]
    return "\n".join(records) + "\n", warnings


def build_sequence_response(ctx, anchor_cluster, genome_id=None):
    """Amino acid + nucleotide sequences for the gene-sequence panel:
    Panaroo's pangenome representative (a real member sequence -- the
    nucleotide is pan_genome_reference.fa as-is, the AA translated from it
    on first request per cluster and cached), plus -- if genome_id is given
    -- that specific genome's own copy, both read straight out of its
    Prokka .faa/.ffn (already correct, not re-derived). Returns a
    JSON-able dict with an "error" key on total failure (bad cluster);
    otherwise "representative" and "strain" are each either {aa, nt,
    length, ...} or None with a "note" explaining why (no
    pan_genome_reference.fa loaded, genome isn't a carrier, no .faa/.ffn on
    disk, etc.) -- never a hard error just because one of the two views
    isn't available."""
    if anchor_cluster not in ctx.cluster_lookup:
        return {"error": f"cluster {anchor_cluster!r} not found"}

    representative = None
    if anchor_cluster in ctx.pan_ref_nt:
        nt = ctx.pan_ref_nt[anchor_cluster]
        aa = ctx.pan_ref_aa_cache.get(anchor_cluster)
        if aa is None:
            aa = translate_dna(nt)
            ctx.pan_ref_aa_cache[anchor_cluster] = aa
        representative = {"aa": aa, "nt": nt, "length": len(aa)}
    else:
        representative = {"note": "pan_genome_reference.fa not available for this deployment"}

    strain = None
    if genome_id:
        present, error = resolve_anchor_loci(ctx, anchor_cluster, [genome_id])
        loci = present.get(genome_id) if not error else None
        if not loci:
            strain = {"note": f"{genome_id} is not a resolvable carrier of {anchor_cluster}"}
        else:
            mge_types = load_mge_types_for_cluster(anchor_cluster, ctx.mge_offsets)
            stem = ctx.genome_stem_map[genome_id]
            # Parsed only when there is a choice to make, so the common
            # single-copy case stays a pure .faa/.ffn lookup.
            genes = parse_gff(stem)[0] if len(loci) > 1 else None
            locus = pick_locus(loci, mge_types, genome_id, genes)
            aa = extract_faa_record(stem, locus)
            nt = extract_ffn_record(stem, locus)
            if not aa and not nt:
                strain = {"note": f"no .faa/.ffn entry for locus {locus} (missing files, or a refound placeholder)"}
            else:
                strain = {"genome_id": genome_id, "locus": locus, "aa": aa, "nt": nt, "length": len(aa)}

    return {"cluster": anchor_cluster, "representative": representative, "strain": strain}


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

        def _send_download(self, text, filename, warnings=(), status=200):
            body = text.encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "text/x-fasta; charset=utf-8")
            self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
            self.send_header("Content-Length", str(len(body)))
            if warnings:
                # Surfaced by the client even though the body is a raw FASTA
                # file, not JSON -- exposed via Access-Control-Expose-Headers
                # since the download is fetch()'d, not a plain <a href>, so
                # the browser hides custom response headers by default.
                joined = "; ".join(w.replace("\r", " ").replace("\n", " ") for w in warnings)
                self.send_header("X-Skipped-Genomes", joined)
                self.send_header("Access-Control-Expose-Headers", "X-Skipped-Genomes")
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

            if parsed.path == "/api/export_contigs":
                anchor = (qs.get("anchor", [""])[0])
                genome_ids = [g for g in (qs.get("genomes", [""])[0]).split(",") if g]
                if not anchor or not genome_ids:
                    self._send_json({"error": "requires 'anchor' and 'genomes' (comma-separated) query params"}, status=400)
                    return
                print(f"[export] anchor={anchor!r} genomes={len(genome_ids)}", file=sys.stderr)
                fasta_text, warnings = build_contig_export(ctx, anchor, genome_ids)
                if fasta_text is None:
                    self._send_json({"error": "; ".join(warnings)}, status=404)
                    return
                safe_anchor = re.sub(r"[^A-Za-z0-9_.-]", "_", anchor)
                self._send_download(fasta_text, f"{safe_anchor}_contigs.fasta", warnings=warnings)
                if warnings:
                    print(f"[export] {len(warnings)} genome(s) skipped: {'; '.join(warnings)}", file=sys.stderr)
                return

            if parsed.path == "/api/sequence":
                anchor = (qs.get("cluster", [""])[0])
                genome_id = (qs.get("genome", [""])[0]) or None
                if not anchor:
                    self._send_json({"error": "missing 'cluster' query param"}, status=400)
                    return
                data = build_sequence_response(ctx, anchor, genome_id)
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
    config_path = Path(args.config) if args.config else (BASE / CONFIG_FILENAME)

    if args.setup:
        run_setup(config_path)
        return

    # Offer setup instead of just failing -- but only with a real terminal
    # to prompt on. Under nohup/systemd/a pipe this has to stay the plain
    # sys.exit() message from resolve_settings(), or the server would hang
    # forever on input() instead of reporting what's missing.
    if missing_required_settings(args, config_path) and sys.stdin.isatty():
        print(f"No usable configuration found ({config_path}).", file=sys.stderr)
        if _prompt_yes_no("Run interactive setup now?", default=True):
            run_setup(config_path)
            return

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
