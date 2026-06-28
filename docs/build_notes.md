# Build Notes (instructor-facing — not for students)

## ⚠️ HANDOFF — READ THIS FIRST IF CONTINUING ON A NEW SERVER (2026-06-28)

**Why this section exists:** this cluster's CPU-only compute proved too slow/unpredictable
for full-scale validation in a reasonable time (see "Performance/hang issues" below). The
instructor is planning to move this project to a different server (possibly with GPU) and
asked me to leave clear notes so a fresh session can resume immediately without re-deriving
context. Read this whole section before doing anything else.

### Current state of the two solution notebooks
- `notebooks/level1/01_snrna_analysis_solution.ipynb`: **built and executed**, but against
  a **tiny 1,500-cell demo subsample** (`scratch_build/tiny_snrna_1500.h5ad`), NOT the real
  118,471-cell student dataset. This was a deliberate, explicit instructor decision to get a
  complete, real, working notebook produced quickly under a hard time limit — not a final
  deliverable. See "What MUST change for the real run" below.
- `notebooks/level2/02_spatial_cell2location_solution.ipynb`: build script
  (`scratch_build/build_solution_nb2.py`) is fully written (all 10 sections) but **not yet
  built/executed** — blocked on Level 1 finishing first (cell2location needs Level 1's
  saved annotated reference as input). Do this next once Level 1 is re-run at full scale.

### What MUST change before this is course-ready (do this on the new server)
1. **Data path** in `scratch_build/build_solution_nb.py`'s TASK 1.1 cell currently points at
   `scratch_build/tiny_snrna_1500.h5ad`. Change back to the real path:
   `/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/data/snRNA_seq/level1_prepared/gbm_l1_snrna_AT10_AT14_raw.h5ad`
   (118,471 cells, AT10+AT14, already QC-stripped of answer-key columns — see "Data prep"
   section below, that part is done and doesn't need to change).
2. **scVI epoch count**: TASK 4.3's code cell currently hardcodes `SCVI_MAX_EPOCHS = 5`
   (deliberately tiny, demo-only). On a faster machine (more cores, or GPU), re-benchmark
   with a *quick, bounded probe* (3-5 epochs, time it directly) and set a **fixed** number
   based on real timing — do NOT reintroduce `scvi.model._utils.get_max_epochs_heuristic()`,
   that is what caused a 2+ hour hang on this server's CPUs with no checkpoint and no
   warning (see "Performance/hang issues" below). A fixed, instructor-chosen number is
   strictly better for a teaching context regardless of hardware — the point is predictability.
3. **Re-execute** both notebooks end-to-end via:
   `conda run -n single_cell jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=<generous> notebooks/level1/01_snrna_analysis_solution.ipynb`
   (and same for level2 once built) — through Slurm (`srun --partition=fast --cpus-per-task=8 --mem=64G --time=<generous>`), not the login node, given the real data scale.
4. **cell2location FAST/FULL preset** (Level 2, TASK 5.2): currently `FAST` = 20 ref epochs /
   300 mapping epochs (~25 min on this server's CPUs, gene-filtered to ~16-21K genes). If
   the new server has GPU, switch the notebook's default to `FULL` (paper-exact: ref
   max_epochs=400, mapping max_epochs=6000) — should run in minutes on GPU vs. the ~8 CPU
   hours I measured here. Real per-epoch CPU numbers are in this file's cell2location section
   below if you need to re-derive a CPU fallback preset on the new hardware instead.
5. **Validate CellTypist + infercnvpy actually ran correctly** in the tiny-scale execution
   (check the executed notebook's own outputs for Sections 6/7) — these were KEPT in the
   notebook per explicit instructor instruction (do not remove them), but I was not able to
   independently re-confirm their behavior at full (117K-cell) scale before the time limit hit.
   Two placeholder cells (`cluster_to_celltype` in TASK 6.5, `REFERENCE_CELL_TYPES` in TASK
   7.2) were changed from hardcoded literals to **dynamic derivation** from the real
   CellTypist/marker output computed earlier in the same notebook run (`summary_df["celltypist"]`,
   and a `MALIGNANT_MIMIC_LABELS` exclusion set) specifically so they work correctly at ANY
   scale without needing pre-computed values — this should keep working as-is at full scale,
   but sanity-check the printed cluster counts make sense once you have real output.
6. **Derive the student notebook** (`01_snrna_analysis_student.ipynb`) from the solution
   once the solution is finalized at full scale — not done yet at all. Same for Level 2's
   student notebook. See the main project plan (`/shared/home/tp185005/.claude/plans/streamed-painting-thimble.md`)
   for the exact style-guide rules (TASK/HINT/QUESTION/CHECKPOINT conventions, coarse
   "Part"-level blanks, no answer-key leakage).

### Performance/hang issues hit on this server (context for why decisions were made)
- A full-scale (117,200-cell, post-QC) run using `scvi.model._utils.get_max_epochs_heuristic()`
  ran for **2+ hours without even finishing the integration stage** (no checkpoint saved) on
  an 8-core Slurm allocation, despite an adaptive capping scheme intended to target ~15 min —
  confirmed via `sstat` CPU-time deltas that it was genuinely computing the whole time, not
  hung/deadlocked, just far slower than the heuristic implied. Killed; root cause not fully
  diagnosed (suspect either the heuristic itself picked an extreme epoch count for this data
  size, or per-epoch cost on this server's CPUs was much higher than typical). **Do not
  reintroduce the heuristic** — use fixed, directly-timed epoch counts instead, on any server.
- `conda run -n single_cell python -u script.py` repeatedly buffered ALL stdout until process
  exit on this server, regardless of `-u`/`flush=True`, making progress monitoring unreliable
  through the harness's log capture. Workaround used: write progress directly to a dedicated
  file via explicit `open()/write()/flush()/close()` cycles, bypassing whatever pipe buffering
  was responsible. If the new server doesn't have this issue, the workaround is harmless but
  unnecessary — can simplify back to plain `print(..., flush=True)` if you confirm it streams.
- cell2location's `RegressionModel.export_posterior()` and `Cell2location.export_posterior()`
  use DIFFERENT real output-key conventions than commonly-quoted tutorial examples — confirmed
  by reading the cell2location 0.1.5 source directly (not docs): `RegressionModel` writes
  `adata.varm[f"{summary}_per_cluster_mu_fg"]`; `Cell2location` writes
  `adata.obsm[f"{summary}_cell_abundance_w_sf"]`. Both already fixed/used correctly in
  `scripts/03_benchmark_cell2location.py` and `scratch_build/build_solution_nb2.py`.
- cell2location's `RegressionModel`/`Cell2location` `setup_anndata()` need `layer="counts"`
  explicitly once `.X` has been normalized/log-transformed earlier in a notebook (easy to
  miss — the error is a cryptic GammaPoisson-support `ValueError`, not an obvious "wrong
  layer" message). Already fixed in both Level 1/2 scripts; if you add any NEW
  cell2location-adjacent code, check this every time.


Working notes from building this project, for future maintainers (or future-me).
Source paper: de Jong, Memi, Gracia, Lazareva et al., *"A spatiotemporal cancer cell
trajectory underlies glioblastoma heterogeneity"*, bioRxiv 2025.05.13.653495.
gbmspace.org. SpaceTree code: github.com/PMBio/spaceTree.

## Donor / sample selection

- 12 donors total in the cohort (AT3-AT15, no AT1/AT2/AT8 ever existed in this study).
  Every donor has *some* snRNA-seq + Visium coverage.
- **Used: AT10 + AT14 only**, all cells, no subsampling (118,471 cells: AT10=85,983,
  AT14=32,488). Originally considered AT10+AT14+AT15, but:
  - `data/snRNA_seq/donor_split/AT15.h5ad` is **HDF5-corrupted** ("bad object header
    version number", confirmed via h5py with multiple drivers — not a permissions issue).
    AT15's data is still recoverable by slicing the 62GB combined `GBM_space_snRNA.h5ad`,
    but wasn't worth the complexity for the gain. **TODO: flag this file for regeneration
    to whoever maintains the shared data folder** — independent of whether we ever use
    AT15, it's just broken and should be fixed at the source.
  - AT15's Visium sections use a different chemistry (CytAssist-style, ~half the gene
    panel, ~3x the spot count) than AT10/AT14's standard Visium — would've introduced a
    confusing technical confound for a teaching cell2location exercise.
  - AT15 is also heavily skewed toward NPC-neuronal-like (42% of its own malignant cells)
    — AT10+AT14 alone already cover all 18 `annotation_coarse` categories without that skew.
- Visium: `anndata_selected/` (9 sections total) is *already* curated down to exactly
  AT10 (5 sections) + AT14 (1 section) + AT15 (3 sections) — i.e. whoever prepared this
  folder for the course already encoded the same matched-donor logic. We use:
  - **Primary**: `AT10-BRA-5-FO-1_2` (3,999 spots, full cell2location/niche/histopath
    answer-key features present).
  - **Optional secondary**: `AT14-BRA-4-FO-2_1` (3,534 spots, no IvyGAP histopath overlap
    feature — flag this if it's used).

## Data quirks found during prep

- `.X` and `.raw.X` in the provided snRNA files are both **raw integer counts** (confirmed
  via h5py, not documented explicitly in the README) — standard normalize/log1p workflow
  applies as-is.
- snRNA `.obs` column names are **lowercase** (`cell_status`, `top_markers`) even though
  the README documents them title-cased (`Cell_status`, `Top_markers`) — caught this via
  a leaky first run of `scripts/01_prepare_snrna_subset.py` (those two columns slipped into
  the student file because the strip-list used the README's casing). Fixed in
  `src/gbmspace_utils/data.py::SNRNA_ANSWER_KEY_OBS_COLUMNS`. **Lesson: verify strip-lists
  against actual column names, never trust documentation casing.**
- **cell2location's `RegressionModel.export_posterior()` does NOT write `varm["q05_cell_abundance_w_sf"]`**
  (the key name some tutorials/old code quote) — confirmed by reading the cell2location 0.1.5
  source directly (`cell2location/models/reference/_reference_model.py`): the real keys are
  `f"{summary}_per_cluster_mu_fg"` for `summary` in `["means","stds","q05","q95"]`, i.e. use
  `adata_ref.varm["q05_per_cluster_mu_fg"]` as the reference signature. First benchmark run
  crashed on this (`IndexError`) — fixed in `scripts/03_benchmark_cell2location.py`.
- **cell2location reference-model training on the full 36,601-gene shared set costs ~170s/epoch**
  on CPU (8 cores, 118,471-cell reference) — far too slow. Applying the *standard* cell2location
  gene filter (`cell2location.utils.filtering.filter_genes`, default cutoffs) drops this to
  15,929 genes and **72s/epoch** (reference) / **3.9s/epoch** (spatial mapping, 3,999 spots).
  Real extrapolated CPU totals: paper-faithful epochs (ref 400 / mapping 6000) would take
  **~8 hours** — confirms GPU is genuinely required for full fidelity. Level 2 notebook uses
  a `FAST` (CPU, ~20 ref epochs / ~300 mapping epochs, ≈40 min total) vs `FULL` (GPU,
  paper-exact 400/6000) preset flag, same pattern as the `INTEGRATION_METHOD` flag in Level 1.
- **`harmonypy` 2.0.0's `Z_corr` is already shaped `(n_cells, n_pcs)`** — `scanpy.external.pp.harmony_integrate`'s
  wrapper still assumes the old `(n_pcs, n_cells)` convention and unconditionally transposes,
  which crashes (`ValueError`, shape mismatch) against this harmonypy version. Fix: call
  `harmonypy.run_harmony()` directly and only transpose if `Z_corr.shape[0] != n_obs` (defensive
  check, works regardless of which convention a given harmonypy version uses).
- **Live BioMart gene-position queries (`pybiomart`) failed** (malformed-XML server response) —
  too flaky to depend on for a live multi-student session anyway. Replaced with a one-time
  static download of the Ensembl GTF (`scripts/` — see `scratch_build/fetch_gene_positions_v2.py`
  for the approach to fold into a proper data-prep script) parsed locally for
  chromosome/start/end per gene symbol — no live API dependency at course time.
- Visium `anndata_selected/*.h5ad` files have ~10 duplicate gene symbols (multi-mapped
  Ensembl IDs, e.g. `TBCE`, `MATR3`) — handled with `var_names_make_unique()` in
  `scripts/02_prepare_visium_subset.py`.
- The paper's own Methods text and figure legends disagree in a few places — worth telling
  students about rather than silently picking one (good "real papers are imperfect"
  teaching moments):
  - Cluster-level malignant-call threshold: ≥3% (Methods text) vs. ≥5% (Ext. Data Fig. 3C/3E legend).
  - cell2fate QC thresholds: max latent time >20 / transition score >0.25 (Methods) vs. >25 / >0.2 (Ext. Data Fig. 17A legend).
  - Visium spot count: 338,481 (Results text) vs. 377,149 (Table S2 total) — likely post-QC vs. raw.

## Paper-faithful parameters reference (for whoever writes/audits notebook content)

- **snRNA QC** (nuclei): genes<500 removed, UMI<1000 removed, mito%>10 removed; Scrublet
  doublets + 2-step MAD filtering (FDR<0.05 cell-level, FDR<0.1 cluster-level); UMI>75,000
  removed post-doublet-calling.
- **Integration**: paper used scVI only (50 latent dims, 2 hidden layers, 1024 nodes/layer;
  batch key = 10x reaction; covariates = tumour ID, site, reaction date, cell-cycle phase).
  We additionally teach Harmony for comparison (paper didn't use it) — pedagogical choice,
  not a paper-fidelity one.
- **Malignant/TME split**: inferCNVpy, window=250 genes, reference = marker-clear TME
  clusters; CNA signal>0.02 AND CNA correlation>0.3 at cell level (see discrepancy above
  for the cluster-level threshold).
- **Malignant axis scoring**: `sc.tl.score_genes` per state (this is literally the paper's
  own method — see `src/gbmspace_utils/analysis.py::MALIGNANT_AXIS_MARKERS` for the exact
  gene sets used, transcribed from Methods/Table S5/S6), cross-checked in the paper via
  scPoli reference mapping (Braun et al. 2023 atlas) and decoupleR/MSigDB enrichment — we
  don't reproduce those two cross-checks, score_genes alone is the teaching-appropriate cut.
- **cell2location**: reference signature — max_epochs=400, batch_size=10000, lr=0.002,
  one reference per tumour. Spatial mapping — N_cells_per_location=30, detection_alpha=200,
  max_epochs=6000, batch_size≈25% of spots. See `scripts/03_benchmark_cell2location.py` for
  the CPU-timing-informed FAST preset used in the Level 2 notebook.
- **Niche analysis**: sklearn NMF, 16 factors/tumour in the paper (cross-tumour cohort);
  we scale down given far fewer spots in 1-2 sections — let students try a few factor counts.
- **Spatial proximity network**: pairwise minimum spot distance via k-d tree, 25th
  percentile summary (implemented as `gbmspace_utils.analysis.spatial_proximity_network`).
- **Spatial intermixing**: Shannon entropy of per-spot cell2location abundance vector.

## Explicitly out of scope (with rationale)

- **spaceTree**: needs a from-scratch clone-calling pipeline (infercnvpy + epiAneufinder on
  paired snATAC-seq, which we don't have students build) before the GNN is even
  applicable, plus a 1920-point hyperparameter grid search just to get the paper's own
  defaults, and no real ground truth to validate against in a classroom setting.
  Independent assessment (not just deferring to the instructor's hunch) confirmed this is
  reasonably out of scope. A cheap substitute (overlay infercnvpy-derived clone clusters
  on top of the already-built cell2location map, purely descriptive) is mentioned as a
  1-paragraph "further reading" pointer in Level 2, not built as an exercise.
- **cell2fate**: needs spliced/unspliced counts generated via STARsolo from raw FASTQs —
  not derivable from the processed h5ad files we have. Revisit if/when that data arrives.

## Environment

`single_cell` conda env (`/shared/projects/tp_2630_ubordeaux_neuromics_184418/envs/single_cell`)
is shared across several course groups (informally the deployed version of the planned
`neuromics-sc` env) and was mid-build when this project started — missing several packages
from the official `cluster_setup/neuromics-sc.yml` wishlist. Added for this project specifically:
`squidpy`, `celltypist`, `harmonypy`, `decoupler` (all via plain `pip install`, no version
pins needed — confirmed no breakage of the existing scanpy/anndata/cell2location/scvi-tools/
torch stack after each addition). `infercnvpy` added during Level 1 build if not already
present (check the actual install log / `pip list` if auditing later). Did NOT install the
rest of the official wishlist (jax, xgboost, snakemake, cooler, etc.) — those serve other
course groups' projects, not this one. Jupyter kernel registered as "Python (single_cell)".
