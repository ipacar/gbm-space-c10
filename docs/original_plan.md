# GBM-Space Summer School Project (C10) — Level 1 & Level 2 Build Plan

## Context

Project C10 for the CAJAL "Neuromics 2026" summer school needs two teaching modules built from real data in the GBM-Space paper (de Jong, Memi, Gracia, Lazareva et al., bioRxiv 2025, gbmspace.org): **Level 1** (snRNA-seq, 2 days) and **Level 2** (Visium spatial + cell2location mapping, 2 days). Each needs a student notebook (guided, blanks to fill) and an instructor solution notebook (fully executed), styled after the existing `cajal_comp_proj` example template in `C10/lederer/`. Level 3 (Xenium) and cell2fate (Fig. 6) are explicitly out of scope for now — Xenium data isn't available yet, and cell2fate needs spliced/unspliced counts we don't have.

Three research agents already completed deep-dives (template structure/pedagogy, paper methods, hands-on data profiling) — their findings drive every decision below. Already done: `lederer/` permissions locked to the instructor only; confirmed `single_cell` conda env (`/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/anaconda3/envs/single_cell`) already has scanpy 1.11.5, anndata 0.12.18, cell2location 0.1.5, scvi-tools 1.4.2, torch, leidenalg — but is missing squidpy, celltypist, harmonypy, decoupler (it's a pip-first, shared-but-instructor-owned env covering several course groups; ad hoc `pip install` additions are the documented, sanctioned way to extend it).

**Decisions confirmed with instructor this round:**
- **Donors: AT10 + AT14 only** (drop AT15 — its `donor_split/AT15.h5ad` is HDF5-corrupted, and its Visium sections use a different chemistry/gene panel than AT10/AT14). AT10 = 85,983 cells (49,653 malignant / 36,330 TME, 4 sites); AT14 = 32,488 cells (21,152 malignant / 11,336 TME, 2 sites); combined 118,471. **No subsampling** — use all cells from both donors as-is; real QC thresholds are the only filter applied (118K isn't far enough above the original ~100K target to justify the complexity/distortion risk of stratified subsampling).
- **cell2location compute:** benchmark CPU timing first on the actual chosen data to size epoch counts realistically; instructor will also pursue GPU partition access before the course — so the notebook should support both via a mode flag.
- **Integration (Level 1):** teach **both Harmony and scVI** side by side (UMAP + mixing-metric comparison), then students pick one via a `INTEGRATION_METHOD = "harmony" | "scvi"` flag for all downstream steps.

## Repository layout (new, under the private instructor folder)

Create `C10/lederer/gbm_space_proj/` (fresh git repo, mirroring `cajal_comp_proj`'s shape, adapted for conda instead of `uv`):
```
docs/
  project_background.md        # student-facing intro (paper withheld until Level 2 reveal point)
  build_notes.md                # instructor-facing: exact params/decisions used, for future maintainers
notebooks/
  level1/01_snrna_analysis_{student,solution}.ipynb
  level2/02_spatial_cell2location_{student,solution}.ipynb
src/gbmspace_utils/
  data.py        # load/subset helpers, answer-key stripping
  analysis.py    # shared malignant-axis scoring function (used in BOTH levels — paper-marker-gene-based sc.tl.score_genes), nhood composition helpers
  plotting.py    # consistent spatial plotting (tissue scatter, gene-on-tissue with 95th-pctile vmax cap)
scripts/
  01_prepare_snrna_subset.py      # builds student data + instructor answer-key from provided h5ads
  02_prepare_visium_subset.py     # strips cell2location/niche/histopath answer-key feature rows
  03_benchmark_cell2location.py   # CPU timing probe, decides fast-vs-full epoch defaults
README.md   # conda activate single_cell instead of uv; points at data already in C10/data/
```
Outputs split by audience:
- **Student-facing prepared data** → `C10/data/snRNA_seq/level1_prepared/` and `C10/data/visium/level2_prepared/` (shared location, same default course-group read access as rest of `C10/data/` — fine, since students must eventually read it).
- **Instructor-only answer keys** (stripped ground-truth columns/features) → `C10/lederer/answer_keys/` (private, same ACL treatment already applied to `lederer/`).

Final published notebooks for students move into the already-existing empty `C10/notebooks/` only when the instructor explicitly signs off — not part of this build.

## Phase A — Environment & data prep (prerequisite, do first)

1. **Complete `single_cell` env**: `pip install squidpy celltypist harmonypy decoupler` (not the full official `neuromics-sc.yml` list — no need for jax/xgboost/snakemake/cooler/etc., those serve other course groups). Add `liana` only if the optional cell-cell-communication stretch goal (Phase C) is confirmed. After each install, re-import scanpy/anndata/cell2location/scvi-tools/torch to confirm nothing broke. Register a Jupyter kernel for the env if missing.
2. **`01_prepare_snrna_subset.py`**: from `donor_split/AT10.h5ad` + `AT14.h5ad` (backed mode), concatenate **all** cells from both donors — 118,471 total, no subsampling. **Strip these `.obs` answer-key columns** before saving the student file: `cell_status, annotation_coarse, annotation_granular, neftel, celltypist, scPoli, ontology keywords, top_markers, CNV_signal_mean, cnv_corr, phase`. Also drop precomputed `.obsm['X_umap']` and `.obsp` neighbor graphs so students compute their own. Keep `donor_id, site_id, sample, n_genes_by_counts, total_counts, mt_frac, doublet_scores`. Save stripped columns + obs_names separately as the instructor answer key. **Verify first whether `.X` is raw or already-normalized counts** (check during this script, before Level 1 QC section is written — paper's own QC thresholds, §below, assume raw counts).
3. **`02_prepare_visium_subset.py`**: from `anndata_selected/AT10-BRA-5-FO-1_2.h5ad` (primary) and `AT14-BRA-4-FO-2_1.h5ad` (optional secondary), drop all `.var` rows where `feature_types != 'Gene Expression'` (removes the pre-computed cell2location "Cell state abundances", "Spatial niche abundances", "Histopath annotation overlap" rows) into the student file; save the removed rows separately as the instructor answer key/checkpoint reference.
4. **`03_benchmark_cell2location.py`**: time a small reference-signature fit + spatial-mapping run on CPU using the actual AT10+AT14 reference and AT10 Visium section, at reduced epoch counts. Use the result to set two epoch presets in the Level 2 notebook: `FAST` (CPU/demo, sized off real timing) and `FULL` (paper-faithful: ref `max_epochs=400, batch_size=10000, lr=0.002`; mapping `N_cells_per_location=30, detection_alpha=200, max_epochs=6000`), selectable via a flag — same pattern as the integration-method flag.

## Phase B — Level 1 notebook: Cellular Census (snRNA-seq, AT10+AT14)

Granular TASK-per-step structure for the solution (matching template precedent: ~8-9 numbered sections, TASK/HINT/QUESTION/CHECKPOINT markers, numeric-range checkpoints, no exact-value checkpoints); coarser "Part"-level blanks for the student version. Core pipeline:
1. **Load & explore** the prepared ~100K-cell subset.
2. **QC**: students propose thresholds; reference point from the paper (nuclei): genes<500, UMI<1000, mito%>10 removed; Scrublet doublets. Checkpoint as a range computed from our actual subset post-filtering (compute when building, don't guess).
3. **Normalize + HVG selection**.
4. **Integration — both methods**: Harmony (`harmonypy`, batch key `donor_id`, fast/simple, taught as primary) and scVI (already installed; batch key `donor_id`, ~paper-lite covariates) side by side — UMAP comparison + a batch-mixing metric. Student sets `INTEGRATION_METHOD` flag for everything downstream.
5. **Clustering** (Leiden) on the chosen embedding.
6. **Broad cell-type annotation**: marker-gene dotplots/scoring (microglia `P2RY12,CX3CR1`; macrophage/monocyte `CD163,STAB1,CD14,FCGR3A`; plus oligodendrocyte/astrocyte/neuron/endothelial/lymphocyte markers from the paper's Table S6) **and CellTypist** (`Developing_Human_Brain` model, matching what the paper itself used) as a cross-check/"reconstruction" tool per the instructor's request — reconcile any disagreements.
7. **Malignant vs. TME split**: lightweight `infercnvpy` pass (reference = clearly-marker-typed macrophage/oligodendrocyte/neuron/lymphocyte clusters), discuss the paper's actual threshold (CNA signal>0.02 AND correlation>0.3) as a reference point — explicitly tell students the paper's own Methods text (3%) and figure legend (5%) disagree on the cluster-level call threshold, as a "real papers are imperfect" moment.
8. **Malignant cell-state axis** (the core ask): score the 4 major classes / 9 `annotation_coarse` malignant subclasses using `sc.tl.score_genes` with the paper's exact marker sets (OPC-like `PDGFRA,OLIG1,SOX6`; NPC-neuronal-like `MYT1L,STMN2,SOX11,DCX`; AC-progenitor-like `SLC1A3,GFAP,EGFR,AQP4,ALDH1L1`; AC-gliosis-like `ITGAV,ITGB1,CDH2,ABCC3`; Gliosis-like `JAK2,STAT3,ANXA2,IL6R,SERPINE1,VEGFA,AKAP12`; Hypoxic `HILPDA,BNIP3L,VEGFA,JUN,FOS`; Proliferative cell-cycle set) — this is exactly the method the paper itself used (`sc.tl.score_genes`, not a black-box tool), so it's both pedagogically simple and paper-faithful. Bonus/optional: reproduce the paper's "myth-busting" finding that EMT genes (`SNAI1/2,TWIST1/2,ZEB1/2`) are flat across gliosis/hypoxia states.
9. **Composition & DE**: cell-type/state proportions per donor (AT10 vs AT14); DE between malignant states or donors.
10. **Publication-quality figure** + save processed h5ad (this becomes Level 2's cell2location reference input).

Flag clearly as **optional/stretch** (protect the 2-day budget): the 9-subclass granular scoring beyond the 4 major classes, and the EMT myth-busting digression.

## Phase C — Level 2 notebook: Spatial Context + cell2location (Visium, AT10 primary / AT14 optional)

1. **Load** prepared Visium (GEX-only, answer-key features already stripped).
2. **Spatial QC** + normalize/cluster directly on spot expression → "naive" spatial domains, fast first win.
3. **Spatial visualization basics** (squidpy `sq.pl.spatial_scatter`, gene-on-tissue).
4. **The cell-state axis in space**: reuse Level 1's marker sets (shared `gbmspace_utils.analysis` scoring function) directly on spots — shows the zonation gradient even pre-deconvolution using the paper's 4-gene minimal spatial panel (`AQP4`→`ABCC3`→`AKAP12`→`HILPDA`).
5. **cell2location**: build per-the-paper reference signatures from Level 1's integrated+annotated AT10+AT14 output, map onto the AT10 Visium section. Epoch counts from the `FAST`/`FULL` flag (Phase A.4 benchmark).
6. **Niche analysis**: NMF (`sklearn`, already installed) on the cell2location abundance matrix, students try a few factor counts — compare against the stripped-out original "Spatial niche abundances" as a checkpoint. Immediately followed by **spatial intermixing entropy** (Shannon entropy of each spot's cell2location abundance vector via `scipy.stats.entropy` — the paper's own "how mixed is this spot" metric, cheap to add here).
7. **Spatial neighborhood analysis, two ways**: squidpy (`nhood_enrichment`, `co_occurrence`) between cell types/niches, **and** the paper's own alternative — a pairwise minimum spot-distance network via `scipy.spatial.cKDTree`, summarized at the 25th percentile (their Fig. 2E/3C/6E/7E method) — explicit "same question, two implementations" comparison.
7b. **(Optional/stretch, needs `liana`)** Ligand-receptor / cell-cell communication: LIANA consensus scoring between TME states co-localized with dev-like niches vs. those co-localized with gliosis-hypoxia niches (using our own niche/cell2location output to define the groups). Skips the paper's cross-donor Tensor-cell2cell step — not meaningful with only 2 donors as "context." Clearly marked optional, scoped only if time/interest allow.
8. **Paper reveal** — placed consistently late (one fixed point, in both student and solution versions this time; the template had an inconsistency here — early in its student notebook, late in its solution — we pick one and apply it uniformly) — followed by guided comparison of the student's own cell2location/niche output against the stripped-out answer key and the paper's actual Fig. 2/3 panels.
9. **Write-up**: reproduce one specific paper panel.
10. (Optional discussion only, not an exercise) — 1-paragraph mention of **spaceTree** as further reading, explicitly not attempted (needs a clone-calling pipeline from paired snATAC-seq we don't have students build; confirmed reasonable to skip outright per the paper agent's independent assessment). **cell2fate** noted as future work pending spliced/unspliced counts.

If AT14's Visium section is used, scope it as an optional "does the same pattern hold in tumour 2?" comparison, not core (it has no IvyGAP histopathology overlay, unlike AT10).

## Style conventions (from template study, applied consistently this time)

🔬 TASK (one-liner, numbered `{section}.{n}`, directly above its code cell) · 💡 HINT (optional, more common early, can state an explicit default) · ❓ QUESTION (standalone reflective markdown, no code) · ⚠️ CHECKPOINT (numeric **range**, self-diagnostic, sparingly — Level 1 only, none in Level 2, matching the template's deliberate "progressive withdrawal of scaffolding"). Student code cells are bare one-line placeholders (`# Your ... here`), no `# TODO:` convention. Unlike the template, the `gbmspace_utils` helper library will actually be imported by the notebooks (at least the shared axis-scoring function) rather than left aspirational/unused.

## Verification plan

- After Phase A.1: re-import full package set in `single_cell`, confirm a kernel is selectable in Jupyter.
- After Phase A.2/A.3: assert no answer-key columns/features present in student files; spot-check counts against the plan.
- After Phase A.4: report actual wall-clock numbers before finalizing notebook epoch defaults.
- Each **solution** notebook: actually execute end-to-end in `single_cell` (`jupyter nbconvert --execute` or equivalent) so all reported numbers/checkpoints are real, not invented.
- Each **student** notebook: grep for answer-key column names / pre-reveal paper phrasing to confirm no leakage; confirm blanks match the style guide.

## Next steps after approval
Build Phase A first (env + data prep + benchmark), report back with real numbers (post-QC cell counts, actual epoch timings) before locking exact checkpoint values into the notebooks — then build Level 1 solution → derive Level 1 student → build Level 2 solution (depends on Level 1's real output) → derive Level 2 student.
