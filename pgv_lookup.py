#!/usr/bin/env python3
"""
Instant single-cluster/single-genome spot-checks,
without re-scanning a huge Panaroo CSV from scratch every time.

Builds a one-time byte-offset index (cluster_name -> file offset of its
row) so any future single-cluster lookup is an instant seek() instead of a
scan, and wraps the whole resolve-stem -> get-locus -> parse-GFF chain
into one function call.

Deliberately NOT wired into pangenome_viewer.py's HTTP server -- this is
scratch tooling for verifying things by hand, kept separate so it carries
none of the "no pre-indexing" tradeoffs that were a deliberate design
choice for the actual viewer (see BUILD_NOTES_pangenome_viewer.md). The
index is small on disk, rebuilds in well under a minute, and nothing
downstream depends on it existing.

Shares pangenome_viewer.py's CLI-flag/config-file settings resolution
(--panaroo-csv/--prokka-dir/etc., or pangenome_viewer.config next to this
script) rather than its own copy, so both tools always agree on where the
data lives.

Usage (CLI):
    python3 pgv_lookup.py --build-index          # one-time, ~15-30s
    python3 pgv_lookup.py group_4219 ESC_AA9696AA_AS
    python3 pgv_lookup.py "espK_2~~~espK_1~~~espK"   # row only, no genome
    # (all the usual --panaroo-csv/--prokka-dir/--config flags apply too;
    # see `python3 pgv_lookup.py --help`)

Usage (import, e.g. from a Bash heredoc for a quick check):
    from pgv_lookup import spot_check
    spot_check("group_4873", "ESC_AC5833AA_AS")
"""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pangenome_viewer as pv  # reuse parse_gff, build_genome_stem_map, settings resolution

BASE = Path(__file__).resolve().parent
OFFSET_INDEX = BASE / "cluster_offsets.json"

csv.field_size_limit(10_000_000)

_stem_map_cache = None
_header_cache = None


def build_index():
    """One pass over the CSV, recording each cluster's row start byte
    offset. Uses explicit readline() (not `for line in f` / csv.reader on
    the file object) because file.tell() is only reliable after readline()
    in text mode -- iterating a file object internally buffers ahead and
    silently desyncs tell() from the logical line position. Splitting the
    first field on a raw comma (not csv.reader) is safe here because
    Panaroo cluster names never contain a comma or quote character."""
    offsets = {}
    with open(pv.PANAROO_CSV, newline="") as f:
        f.readline()  # header
        pos = f.tell()
        line = f.readline()
        while line:
            cluster = line.split(",", 1)[0]
            offsets[cluster] = pos
            pos = f.tell()
            line = f.readline()
    OFFSET_INDEX.write_text(json.dumps(offsets))
    print(f"[pgv_lookup] indexed {len(offsets)} clusters -> {OFFSET_INDEX}", file=sys.stderr)
    return offsets


def _load_offsets():
    if not OFFSET_INDEX.exists():
        print("[pgv_lookup] no index yet, building one (~15-30s, one-time)...", file=sys.stderr)
        return build_index()
    return json.loads(OFFSET_INDEX.read_text())


def _header():
    global _header_cache
    if _header_cache is None:
        with open(pv.PANAROO_CSV, newline="") as f:
            _header_cache = next(csv.reader(f))
    return _header_cache


def get_cluster_row(cluster_name, offsets=None):
    """(header, row) for one cluster, via an instant seek instead of a
    full-file scan. Raises KeyError with a clear message if not found."""
    if offsets is None:
        offsets = _load_offsets()
    if cluster_name not in offsets:
        raise KeyError(f"{cluster_name!r} not in the Panaroo CSV (or index is stale -- "
                        f"rerun with --build-index if this cluster should exist)")
    header = _header()
    with open(pv.PANAROO_CSV, newline="") as f:
        f.seek(offsets[cluster_name])
        row = next(csv.reader(f))
    return header, row


def resolve_stem(genome_id):
    global _stem_map_cache
    if _stem_map_cache is None:
        _stem_map_cache = pv.build_genome_stem_map()
    return _stem_map_cache.get(genome_id)


def get_locus(cluster_name, genome_id, offsets=None):
    """List of locus tags for this cluster in this genome, or None if
    absent. Empty list never returned -- absence is always None."""
    stem = resolve_stem(genome_id)
    if stem is None:
        raise KeyError(f"genome_id {genome_id!r} not found in the genome stem map")
    header, row = get_cluster_row(cluster_name, offsets)
    if stem not in header:
        return None
    cell = row[header.index(stem)].strip()
    if not cell:
        return None
    return [x.strip() for x in cell.split(";") if x.strip()]


def get_gene_detail(genome_id, locus_tag):
    """(contig, start, end, strand, product) for one locus tag, or None if
    it's a refound placeholder / otherwise has no GFF entry."""
    stem = resolve_stem(genome_id)
    genes, _seqlens, _noncds = pv.parse_gff(stem)
    return genes.get(locus_tag)


def spot_check(cluster_name, genome_id=None):
    """Prints a quick summary of one cluster, optionally in one genome."""
    offsets = _load_offsets()
    if genome_id is None:
        header, row = get_cluster_row(cluster_name, offsets)
        n_carriers = sum(1 for cell in row[3:] if cell.strip())
        print(f"{cluster_name}: {n_carriers}/{len(header) - 3} genomes carry it "
              f"(non_unique_name={row[1]!r}, annotation={row[2]!r})")
        return
    loci = get_locus(cluster_name, genome_id, offsets)
    if loci is None:
        print(f"{cluster_name} is ABSENT from {genome_id}")
        return
    for lt in loci:
        detail = get_gene_detail(genome_id, lt)
        if detail is None:
            print(f"{cluster_name} in {genome_id}: locus {lt!r} -- no GFF entry "
                  f"(Panaroo 'refound' placeholder, Prokka never called it)")
        else:
            contig, start, end, strand, product = detail
            print(f"{cluster_name} in {genome_id}: locus {lt!r} on {contig} "
                  f"{start}-{end} ({strand}), {end - start + 1}bp -- {product}")


def _configure():
    """Reuses pangenome_viewer's argparse + config-file settings resolution
    (same --panaroo-csv/--prokka-dir/--config flags, same
    pangenome_viewer.config lookup) so this tool and the HTTP server always
    agree on where the data lives, plus this script's own positional args."""
    ap = pv.build_arg_parser()
    ap.add_argument("cluster", nargs="?", help="Panaroo cluster name")
    ap.add_argument("genome_id", nargs="?", help="optional -- omit to just see carrier count")
    ap.add_argument("--build-index", action="store_true", help="(re)build cluster_offsets.json")
    args = ap.parse_args()
    resolved = pv.resolve_settings(args)
    pv.configure_globals(resolved)
    return ap, args


if __name__ == "__main__":
    ap, args = _configure()
    if args.build_index:
        build_index()
    elif args.cluster:
        spot_check(args.cluster, args.genome_id)
    else:
        ap.print_help()
