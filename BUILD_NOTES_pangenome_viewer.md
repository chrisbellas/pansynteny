# How the local pangenome viewer was built

Not a published artifact — a local web tool, reached via SSH tunnel:
```
ssh -L 8765:localhost:8765 <this-host>
```
then `python3 pangenome_viewer.py --port 8765` here, open `http://localhost:8765/`
on your laptop.

Companion doc: the approved build plan, `~/.claude/plans/wiggly-snuggling
-heron.md` on this machine (two plans were written to that same file over
the course of the build — the original server+search+chart plan, then the
"present elsewhere" gutter-marker stretch goal; both fully executed).

## Why this exists, and how it differs from the 4 published charts

LEE/espK/sopA/derep_00002 (`BUILD_NOTES_lee_chart.md` etc.) are one
hand-written `build_<name>_chart_data.py` script per anchor gene, each
producing one static HTML artifact. This tool generalizes that into a
server: search any of the pangenome's 24,827 Panaroo clusters by name or
annotation, click a result, get the same chart style built on demand
(~10-30s, no pre-indexing). No SQLite/database — deliberately reuses the
existing scripts' full-CSV-scan approach as-is, per the user's explicit
choice (see the plan file).

## Files (all in `genomad_analysis_H2.2/`)

| File | Role |
|---|---|
| `pangenome_viewer.py` | stdlib-only (`http.server.ThreadingHTTPServer`, no Flask/pandas) server. Builds an in-memory search index once at startup, then `build_chart_data(ctx, anchor_cluster)` generates a chart's JSON on demand per request. |
| `viewer_template.html` | The page: search box, chart area, all JS. Served as-is (bytes cached in-process after first read); the chart's data comes from a live `fetch('/api/chart?...')`, not embedded JSON (this only works because it's a normal server, not a claude.ai artifact — the sandbox-blocks-`fetch()` constraint that forced inline-embedding on the other 4 charts doesn't apply here). |
| `pangenome_viewer.py.bak`, `.bak2`, `.bak3` / `viewer_template.html.bak`, `.bak2`, `.bak3` | Rollback checkpoints (no git repo in this directory — deliberate, see below). `.bak` = pre-any-marker-feature. `.bak2` = gutter markers for missing genes only. `.bak3` = all-8-genes-always-shown + legend click-to-highlight, pre-width-fix. Current working files are ahead of all three. |

Read-only reuse, no changes: `dereplicate_rfe_mges.build_genome_stem_map()`,
`gene_presence_absence.csv`, `important_feastures.txt`,
`mge_genes_all_genomes.tsv`, Prokka GFFs, `2026_combined_HC20`,
`human_all_gene_test_predictions.csv` — same sources as the 4 published
charts.

## Running it

```bash
cd genomad_analysis_H2.2
python3 pangenome_viewer.py --port 8765   # binds 127.0.0.1 only, by design
```
Restart after any edit — `template_bytes_cache` and the whole search index
are loaded once at process start, nothing is re-read from disk per request
except the two full CSV passes inside `build_chart_data()` itself.

## Key constants in `pangenome_viewer.py`

Same conventions as the 4 curated charts: `SAMPLE_SIZE = 500`,
`RANDOM_SEED = 42`, `TOP_N_COLORED = 8`. `WINDOW_OPTIONS_BP = [5000, 10000,
20000, 30000, 40000]`, always built at the widest (40kb) as a superset;
client-side toggle only filters visibility, never rebuilds. `flip = e_strand
== "-"` is the default (the LEE/espK/derep00002 convention) — sopA's
one-off inverted flip was never inherited here.

## Non-obvious design decisions worth knowing before modifying

1. **Robustness fixes the 4 curated scripts never needed, built in from the
   start**: `sample_size = min(SAMPLE_SIZE, len(present))` plus an explicit
   0-carrier early return, since arbitrary searched genes very often have
   <500 carriers (pangenome frequency is heavily rare-gene-skewed) —
   `random.sample(population, 500)` raises `ValueError` when `population <
   500`, which none of the 4 hand-picked anchors ever hit by luck.
2. **`reverse_map[g]` (built during the per-chart CSV pass, `pangenome_viewer.py:280-301`) is a *complete* `locus_tag -> cluster` map for each sampled genome's entire gene complement, genome-wide — not filtered to the drawn window.** This one fact is why the "present elsewhere" feature (point 4 below) cost zero extra file I/O: the data was already being loaded and then discarded before that feature existed.
3. **Panaroo "refound" placeholder loci** (format `"<n>_refound_<m>"`) —
   Panaroo's own gene-refinding re-search of a genome's raw assembly,
   assigned when enough *other* genomes' Prokka pass found a gene but this
   genome's own Prokka pass missed it. These never appear in that genome's
   own GFF (Prokka never called them), so any GFF-lookup-based tool
   silently drops them. Full discovery/prevalence numbers are in project
   memory (`group_22355` investigation, 2026-07-24) — full cohort: 3.9% of
   all gene instances, 15.0% of RFE-feature instances specifically, some
   clusters >90% refound. Two responses to this, both free (the refound
   -ness of a locus tag is derivable from the string itself, already in
   memory from the existing CSV read):
   - `pick_locus()` prefers a real (non-refound) locus over a refound one
     when a genome carries multiple copies of the anchor gene, before the
     existing geNomad-MGE-classification preference.
   - `n_dropped_refound`/`n_dropped_other` and `carrier_loci_total`/
     `carrier_loci_refound` are returned in the API response; the template
     surfaces a subtitle warning whenever >20% of a gene's genome-wide
     carriers are refound-only, so a low row count reads as a
     data-quality signal rather than implied rarity.
4. **"Present elsewhere" gutter markers (2026-07-24, user stretch-goal
   request), the main feature added after the base tool shipped.** For
   each of the 8 colored (top8_rfe_clusters) genes, every row gets a small
   square in a fixed-position gutter column to the right of the plotted
   gene track — one column per gene, so a column is scannable straight
   down all 500 rows for a consistency check. Went through two rounds:
   - **v1 (missing-only)**: a box appeared only for genes *not* drawn in
     that row's window, classified into 3 states via data already loaded
     (`reverse_map[g]` inverted to `cluster -> [locus_tags]`, cross
     -referenced against `gff_genes[g]`, which holds every CDS on every
     contig, not just the drawn one): `different_contig` (solid border),
     `same_contig_outside_window` (dashed), `refound_only` (dotted, no
     GFF entry to resolve at all). Truly absent genes got no box.
   - **v2 (all-8-always, current)**: per user request ("add them if
     visible or elsewhere... so I can scan across for consistency"),
     genes already drawn in the window now also get a box — a 4th state,
     `visible`, rendered as a plain filled square with no border (`stroke:
     none`). This makes "blank gutter slot" mean *only* "genuinely
     absent," never ambiguous with "present in the track." Field renamed
     `missing_top8` -> `top8_markers` to match (JSON key change, not
     backwards compatible, no other consumer existed).
   - Server-side classification block sits right after `top8_rfe_clusters`
     is computed (it has to — that list isn't known until after all rows'
     genes are drawn and counted) and before the function's final
     `print`/`return`, `pangenome_viewer.py` around line 407-440.
   - Client-side: `viewer_template.html`'s `draw()` widened `totalW` by a
     fixed `gutterW` (~110px, 8 squares x 12px pitch + left pad) and draws
     the squares in the same per-row loop as the gene arrows, using
     `clusterColor[cl]` (the same color the legend/track already use) for
     fill.
   - Click-to-highlight (`applyPin`/`pinnedCluster`) was extended to treat
     gutter squares exactly like track gene-arrows: `querySelectorAll(
     ".gene-arrow, .marker-sq")` toggles `.pinned` on whichever elements
     match `dataset.cluster`. Clicking a gutter box highlights the real
     gene copy on the track and vice versa.
   - **Legend swatches are also clickable (2026-07-24, same-day follow-up
     request)**: the 8 top-gene legend entries got `data-cluster`
     attributes and a click listener calling `applyPin`; `applyPin` also
     toggles a `.pinned-legend` class (box-shadow outline) on the matching
     legend item. So track, gutter, and legend are now one unified
     highlight group — click any one, all three light up together.
   - Verified against raw ground truth by hand, not just "it renders": one
     `refound_only` case (`group_4219` in `ESC_AA9696AA_AS`) confirmed the
     raw Panaroo cell was exactly `7_refound_1706_pseudo`, no real locus;
     one `different_contig` case (`group_4873` in `ESC_AC5833AA_AS`)
     confirmed via direct GFF lookup that it resolves to a different
     Prokka contig than espK's own locus in that same genome.
5. **Page width (2026-07-24, user follow-up: "why doesn't it stretch
   across my whole monitor").** `.wrap`'s CSS `max-width` was widened
   1520px -> 1800px — chosen specifically because the chart SVG's full
   pixel width (`labelW + plotW + gutterW + 20` = 1730px) now fits inside
   1800px without needing horizontal scroll on any monitor >=1778px wide.
   Deliberately did *not* also make `plotW` itself responsive/redraw-on
   -resize (offered as an option, user picked the simpler fix) — the gene
   track stays a fixed 1280px regardless of window size, only the outer
   page's blank margins shrink/grow.
6. **No git repo in this directory, deliberately** — rollback is via the
   numbered `.bak`/`.bak2`/`.bak3` file pairs (see the Files table above),
   not version control. If a future change should get finer-grained
   rollback than "restore one of the 3 checkpoints," consider `git init`
   at that point rather than continuing to hand-number backups.
