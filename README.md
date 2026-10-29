# processor-cleanet

Doublet removal for flow and mass cytometry with **Cleanet** (Ionita et al. 2025, Penn
Institute for Immunology and Immune Health), step 3 of the Immune Profiling & Differential
Abundance workflow: `fcs-prep → peacoqc → cleanet → flowsom → {cytometry-umap, diffcyt} → …`

A doublet is two cells that pass the laser together and are recorded as one event. Its
marker values are roughly the sum of the two cells', so a CD4 T cell stuck to a B cell shows
up as a "CD3+CD19+" event, and enough of them form small fake populations in clustering.
The usual remedy is a bivariate gate on FSC-A vs FSC-H (or DNA content on CyTOF), which
misses doublets of two small cells and can remove large single cells. Cleanet instead
**simulates** doublets: it adds the linear values of random pairs of events in the file
(as many simulated doublets as half the events), then looks at each real event's 15
nearest neighbours in marker space. An event whose neighbourhood is dominated by simulated
doublets (more than 5 of 15) is called a doublet. This finds heterotypic doublets (two
different cell types) well; homotypic ones (two cells of the same type) look like bright
singlets and are found only partly.

Each file is handled on its own. The app writes the singlets with the same file name,
channels and keywords, so later steps read them like `peacoqc`'s output.

## What Cleanet expects, and how the app feeds it

From the Cleanet source (commit `d2532096`) and the paper:

- **Linear values.** A simulated doublet is the sum of two events' *untransformed*
  values; Cleanet arcsinh-transforms internally (`cofactor`, one for all channels) before the
  neighbour search. The paper ran it on compensated, untransformed data. Our files arrive
  arcsinh-transformed by `fcs-prep` (x = asinh(raw / c), c from `panel.csv`), so the app
  undoes that first: raw = c · sinh(x). With `COFACTOR=panel` (default) it hands Cleanet
  raw / c with cofactor 1, so distances are computed on exactly the arcsinh scale the rest
  of the workflow uses; a numeric `COFACTOR` uses raw values with that one cofactor (Cleanet's
  defaults: 500 for flow, 5 for CyTOF). Without `panel.csv`, cofactors come from the
  `PENNSIEVE_COFACTORS` keyword `fcs-prep` writes into each file.
- **Which channels.** Protein markers only: `panel.csv` rows with `type == marker`. Scatter is
  not used for the neighbour search (the paper uses FSC/SSC only to deplete debris).
  On CyTOF include the DNA intercalator channels as markers (the Cleanet vignette uses
  DNA1/DNA2 together with the lineage markers).
- **Debris.** Debris plus a cell does not make a doublet, so Cleanet recommends keeping
  debris under ~10 % of events or excluding it from the simulation. `DEBRIS_SUM_MAX`
  applies Cleanet's flow rule (FSC-A + SSC-A below a cut-off, both positive); it is off by
  default because the workflow's inputs are gated (the demo files to CD45+ live cells).
- **Spectral flow.** The paper validates Cleanet on mass cytometry and on a 16-colour
  fluorescence dataset (it does not say whether conventional or spectral) and on imaging
  flow data; spectral data are not discussed explicitly. The method needs only that a
  doublet's signal is the sum of its two cells' signals, which holds for unmixed spectral
  data (unmixing is linear), so it applies; the demo results below are consistent with that.

## Ports

| Port | Kind | Files | Notes |
|---|---|---|---|
| `fcs` | input (subfolder) | `<sample_id>.fcs` | From `peacoqc` (or `fcs-prep`): FCS 3.1, float32, markers arcsinh-transformed. Found in `INPUT_DIR/fcs/`, else anywhere under `INPUT_DIR`. |
| `metadata` | input (subfolder) | `panel.csv` | Marker channels and their cofactors. Optional: without it the cofactors come from the `PENNSIEVE_COFACTORS` keyword and the markers are the channels with a `$PnS`. |
| `fcs` | output | `fcs/<sample_id>.fcs`, `fcs/provenance.json` | The singlets (`KEEP_DOUBLETS=true`: every event). Same channels, `$PnS` and keywords as the input, plus `PENNSIEVE_DOUBLETS=cleanet;doublets=<n>;pct_doublets=<p>;removed=<n>;thresh=5;cofactor=panel`. |
| `qc` | output | `qc/doublets.csv`, `qc/plots/<sample_id>.png`, `qc/provenance.json` | Per-sample table and diagnostic plot. |

### `qc/doublets.csv`

| Column | Meaning |
|---|---|
| `sample_id` | File name without `.fcs` |
| `events_in` | Events read |
| `doublets` | Events Cleanet called doublets |
| `pct_doublets` | 100 × doublets / events_in |
| `events_out` | Events written (events_in − doublets, or events_in with `KEEP_DOUBLETS=true`) |
| `sensitivity` | Fraction of the simulated doublets that Cleanet itself classifies as doublets: its internal check that doublets are separable in this panel (near 1 is good) |
| `debris_excluded` | Events left out of the simulation by `DEBRIS_SUM_MAX` (empty when off) |
| `n_channels` | Marker channels used |
| `status` | `ok`, `passed_through` (fewer than `MIN_EVENTS`; copied unchanged) or `failed` |
| `notes` | Warnings and remarks |

### Plots

Cleanet has no plotting function of its own for this (its debris helpers draw interactive
scatter plots), so the app draws one figure per sample: left, the distribution of the
number of simulated doublets among each event's 15 neighbours, for the real events and for
the simulated doublets, with the threshold; right, the calls on FSC-A vs FSC-H (not used by
Cleanet), where doublets should sit at higher FSC-A for their FSC-H.

## Parameters

| Parameter | Type | Default | What it does |
|---|---|---|---|
| `KEEP_DOUBLETS` | boolean | `false` | `true`: write every event (detection and report only). |
| `THRESH` | integer | `5` | Doublet when more than `THRESH` of the 15 neighbours are simulated doublets (Cleanet's `thresh`). Lower calls more. |
| `COFACTOR` | `panel` or number | `panel` | Distance scale; see above. |
| `DEBRIS_SUM_MAX` | number | `0` (off) | Exclude events with FSC-A + SSC-A below this from the simulation. |
| `MIN_EVENTS` | integer | `1000` | Smaller files are passed through unchanged. |
| `PLOTS` | boolean | `true` | Write the diagnostic plots. |
| `THREADS` | integer | `2` | Files processed at once (one R process each; Cleanet's neighbour search is single-threaded). |
| `SEED` | integer | `42` | Seed for the random pairs; results are reproducible for a given seed. |

## Results on the demo data

All 24 CytoBatchFlagR files (Koenig et al. 2025, Zenodo 15388817) after `fcs-prep` and
`peacoqc` (9.80 M events, 25 markers), default parameters.

**These files were already gated to single cells in FlowJo** (CD45+CD14− live singlets, by
scatter), so few doublets are expected: what Cleanet finds here are doublets the scatter
gate let through. Measured:

- 117,594 of 9,797,603 events called doublets: **1.20 %** overall; per file 0.41–2.28 %,
  median 1.25 %.
- Stimulated samples (anti-CD3/CD28 + IL-2, tcs3/tcs4) have more: 1.1–2.3 % against
  0.4–1.3 % for the unstimulated ones (tcu1/tcu2), in every batch; plausibly because
  activated T cells cluster together in culture.
- Sensitivity on simulated doublets 0.960–0.986 in every file: heterotypic doublets are well
  separated in this 25-marker panel.
- In the plots we inspected (e.g. batch06_tcs4) the called doublets sit mostly below the
  singlet diagonal on FSC-A vs FSC-H (higher area for their height), the classic doublet
  signature, at the edge of what the FlowJo singlet gate kept.

The full output is in `workflow-immune-cytometry/data/work/cleanet/` for the downstream apps.

## Runtime

`computeTypes: [standard]`, 2 vCPU / 8 GiB, 4 h timeout. The cost is the approximate
nearest-neighbour search (HNSW, single-threaded inside Cleanet) over 1.5 × the events of a
file; memory and time grow roughly linearly with events.

Measured on the 24 demo files (9.8 M events), amd64 image under emulation on an M-series Mac
(Docker Desktop), anonymous memory of the container's cgroup sampled every 0.5 s:

| | wall time | peak memory |
|---|---|---|
| largest file alone (600k events) | 4 min 25 s | 1.84 GB |
| all 24 files, `THREADS=2` (default) | 29 min 13 s – 33 min 30 s (two runs) | 3.8–3.9 GB |

Expect it to be several times faster on native amd64. 8 GiB holds two ~1M-event files at once.

Image: 1.88 GB unpacked, 474 MB compressed (`docker images`): `rocker/r-ver:4.5.3` plus a
412 MB layer of R packages (Cleanet, RcppHNSW, flowCore and their dependencies, built from
source, debug symbols stripped).

## Running it locally

```sh
pip install -r requirements-dev.txt   # flowio, numpy: only for the fixture script
make fixture                          # demo data -> data/fixture/{fcs,metadata}, 20k events/file
make run                              # build the image, run over data/fixture -> data/output
make test                             # unit tests inside the image
```

`make fixture` stands in for the upstream steps (see `sample-input/README.md`). The unit
tests build a tiny synthetic fixture in R — four populations plus planted doublets made by
adding two events' linear values, a file without doublets and one too small — and check
the calls, the output files and both input layouts. `make test` runs them in the image,
`make test-local` with a local R that has Cleanet, flowCore, jsonlite and testthat.

## Limitations

- Homotypic doublets (two cells of the same type) are close to bright singlets in marker
  space and are found only in part; a scatter or DNA gate remains complementary.
- Samples with much debris need `DEBRIS_SUM_MAX` (set for the cytometer's scatter scale).
- Cleanet's neighbour search is single-threaded; parallelism is across files only.
- flowCore's writer adds a non-standard `FCSversion=3` keyword and writes big-endian
  (`$BYTEORD` 4,3,2,1); the header says FCS3.1. Readers handle both.
- Doublet classification by cell type (Cleanet's `classify_doublets`) needs cluster labels
  and is not part of this app.

## Licences

Our code is Apache-2.0. The image runs Cleanet and RcppHNSW (GPL-3.0) and flowCore with
cytolib (AGPL-3.0-only) in the same R process, so the image is distributed under those terms.
See `THIRD_PARTY.md`. Cleanet is installed from GitHub at commit
`d2532096d00f69ab911080e008cc2e27f8cbcb45`; `provenance.json` records it.

## Citations

- Ionita M, McKeague ML, Painter MM, Mathew D, Pattekar A, Rezk A, Shwetank, Maseda D,
  Wherry EJ, Greenplate AR. Cleanet: Robust Doublet Detection in Cytometry Data Based on
  Protein Expression Patterns. *Cytometry A* 107(11):716–729 (2025). doi:10.1002/cyto.a.24961.
  PMID 41065042. Software: https://github.com/matei-ionita/Cleanet (GPL-3).
- Hahne F, LeMeur N, Brinkman RR, Ellis B, Haaland P, Sarkar D, Spidlen J, Strain E,
  Gentleman R. flowCore: a Bioconductor package for high throughput flow cytometry.
  *BMC Bioinformatics* 10:106 (2009). doi:10.1186/1471-2105-10-106. PMID 19358741.
- Demo data: Koenig Z, Andorf S, Tilburgs T, Eswar S. CytoBatchFlagR: Spectral flow
  cytometry dataset. Zenodo (2025). doi:10.5281/zenodo.15388817. CC BY 4.0.
