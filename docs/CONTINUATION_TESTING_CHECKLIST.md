# Continuation Checklist — Full-Scale Testing on a New Server

Written 2026-06-28 while building/testing on a CPU-limited cluster. Everything below is
true as of right now; check it against the actual notebook state before trusting it if
you're reading this much later. See also `docs/build_notes.md` (the "HANDOFF" section at
the top has the root-cause writeups for every bug found) and `docs/original_plan.md` (the
full approved plan this project was built from — read that first for overall intent).

## Current state in one line
**Level 1 solution notebook: built and fully executed (44/44 code cells, zero errors)** —
but against a **tiny 1,500-cell demo subsample**, not the real 118,471-cell dataset.
**Level 2 solution notebook: built, partially executed** (in progress as of this writing) —
against the real small Visium section (3,999 spots) but with cell2location capped at a tiny
5/20-epoch "DEMO" preset and consuming Level 1's tiny-scale demo output as its reference.
**No student notebooks built yet for either level.**

## Why everything ran at tiny/demo scale (context, not an excuse)
This server's CPUs made full-scale validation impractical in the available time:
a full 117,200-cell scVI integration step hung for 2+ hours using scvi-tools' adaptive
epoch heuristic (root-caused and removed — see build_notes.md); cell2location's reference
model costs ~72s/epoch even after gene filtering on the full reference. Under explicit time
pressure from the instructor, the decision was made to validate the **entire pipeline's
logic and code correctness** at tiny scale (fast iteration, several real bugs caught and
fixed this way — see "Bugs fixed" below) and defer **scientific/numerical re-validation at
real scale** to a faster server. This was a deliberate, explicit trade — not an oversight —
but it means the notebooks' actual printed numbers, thresholds, and even some qualitative
conclusions should NOT be trusted as scientifically final until re-run at full scale.

## Cell-by-cell: what needs re-testing at full scale

### Level 1 (`01_snrna_analysis_solution.ipynb`)
| Section | Status at tiny scale | What to check at full scale |
|---|---|---|
| 1. Load | OK | Just swap the data path (see below) |
| 2. QC | OK, but trivial (almost nothing filtered — see build_notes, the data was already pre-filtered) | Re-confirm the doublet-score cut is still the only meaningful lever at full scale |
| 3. Normalize/HVG | OK | None expected, just re-run |
| 4. Integration (Harmony + scVI) | Harmony OK. **scVI ran only 5 fixed epochs** ("demo", not meaningful convergence) | Re-time scVI with a quick probe (3-5 epochs, measure seconds/epoch directly), pick a **fixed** epoch count sized to your actual time budget — do NOT use `get_max_epochs_heuristic()` (root cause of the original hang). Re-check the neighbor-purity comparison numbers — they were computed but on a near-untrained scVI embedding, not meaningful yet |
| 5. Clustering | OK, but only **8 Leiden clusters** at 1,485 cells, resolution 0.5 — re-check resolution choice at full scale (template precedent suggests sweeping 0.3/0.5/1.0) |
| 6. Annotation (markers + CellTypist) | **Ran successfully** — CellTypist `Developing_Human_Brain` model produced real region-specific labels (e.g. "Hippocampus OPC", "Subcortex glioblast"). Mechanically fine. | Re-check whether `cluster_to_celltype = summary_df["celltypist"].to_dict()` (now a **dynamic** derivation, not hardcoded) gives sensible per-cluster calls at full scale — at 1,485 cells some clusters were very small and CellTypist's per-cell labels noisy |
| 7. Malignant/TME split (infercnvpy) | **Fallback triggered**: at tiny scale, *zero* clusters cleared the `SIGNAL_CUT>0.02 AND CORR_CUT>0.30` thresholds at the standard 20% cluster-fraction cut — the code fell back to "just take the single highest-CNA cluster," yielding only 31/1,485 "Malignant" cells. The paper's own finding is ~55-65% malignant in this cohort — **this 31-cell result is a tiny-scale artifact, not a real finding.** Re-run at full scale and confirm the *primary* threshold path is what actually fires (the fallback should NOT be needed with real statistical power) | This is the single most important "don't trust this number" flag in the whole project |
| 7. `REFERENCE_CELL_TYPES` (CNV reference) | Now a **dynamic keyword-based derivation** (`MALIGNANT_MIMIC_KEYWORDS = ("glioblast", "radial glia", "opc", "neural crest", "neuroblast")`) instead of a hardcoded list — worked at tiny scale (5 of 13 CellTypist labels classified as reference). Re-check this keyword list is comprehensive against the FULL CellTypist label vocabulary that emerges at 118K cells (more clusters → more distinct labels → keyword list may need additions) |
| 8. Malignant axis scoring | Ran on only 31 malignant cells (see above) — numbers are not meaningful. Re-run once the malignant/TME split is fixed at full scale |
| 8 (bonus). EMT myth-busting | Present in the notebook, ran without error, but on the same tiny malignant subset — re-check |
| 9. Composition/DE | Ran, but on tiny/skewed data — re-check |
| 10. Figure | Ran, cosmetically fine, content will look different at full scale |
| 11. Save | **Real bug-fixed code, works correctly** — writes `data/processed/gbm_l1_snrna_AT10_AT14_annotated.h5ad`. This is what Level 2 consumes — must be re-run at full scale before Level 2's cell2location section means anything for real |

### Level 2 (`02_spatial_cell2location_solution.ipynb`)
As of this writing, execution was in progress (14/24 code cells done, no errors) when this
checklist was written — **check the notebook's own cell `execution_count`/outputs directly
for current ground truth**, this table may already be stale.

| Section | Status | What to check |
|---|---|---|
| 1-4. Load/QC/naive clustering/axis-in-space | Should be fine — ran on the real 3,999-spot AT10 section, not synthetic data | Re-check Section 3's naive-cluster count and Section 4's axis-score gradients look sensible once Level 1's real (non-tiny) reference exists, in case anything downstream implicitly depends on cluster identity |
| 5. cell2location | Run at **"DEMO" preset (5 ref epochs / 20 mapping epochs)** — proves the code path works (the cell2location bugs below are real fixes, not guesses) but is nowhere near converged | Switch `C2L_MODE` to `"FAST"` (CPU, ~20/300 epochs, real timing already measured: ~72s/epoch ref @ ~16-21K genes, ~3.9-5.2s/epoch mapping — see build_notes.md) or `"FULL"` if GPU is available on the new server. **Also**: this section currently consumes Level 1's *tiny-scale* annotated output as its reference — must be re-run after Level 1 is re-run at full scale |
| 6. Niche analysis (NMF) + intermixing entropy | **Not yet confirmed executed** at time of writing — check directly | If it ran: sanity-check niche count choices (5/8/12 tried) make sense; the entropy computation is simple/robust and shouldn't have scale-dependent bugs |
| 7. Neighborhood (squidpy) + proximity network | **Not yet confirmed executed** | Check for errors; the custom `spatial_proximity_network` helper (`gbmspace_utils.analysis`) was unit-tested only implicitly through this one run |
| 8. Optional CCC (LIANA) | **Not implemented** — the cell is literally a comment placeholder (`# import liana as li ... implement if time allows`). This matches the *original plan's* "optional/stretch, only if time allows" framing, so it's not a regression, but it's genuinely empty right now if anyone goes looking for it |
| 9. Paper reveal | **Markdown placeholder not filled in** — currently reads `[INSTRUCTOR: insert full citation + 3-4 bullet summary of key findings here...]`. The real citation and findings ARE available (an earlier research pass extracted them in full) — see `docs/build_notes.md`'s header for the citation; **this needs to be written before the notebook is student-ready**, it's a TODO not a bug |
| 10. Write-up | Template/placeholder code cell (`# Your figure-reproduction code here.`) — this is *meant* to be open-ended (matches the plan), not a gap |

## Gaps vs. the original plan (`docs/original_plan.md`)

Comparing what's actually in the notebooks against Phase B/C of the original plan:

1. **Paper reveal content (Level 2, Section 9)** — planned to have a real citation + findings
   summary; currently a placeholder. Needs filling in before this notebook is student-ready.
2. **LIANA cell-cell-communication (Level 2, Section 7b)** — planned as explicitly optional/
   stretch; currently a comment stub with no real implementation. Consistent with the plan's
   own "only if time allows" framing, but flagging since "implemented" and "stubbed" look
   different from a quick glance at the notebook.
3. **Kernel/style consistency** — both notebooks use the `single_cell` kernel and the
   TASK/HINT/QUESTION/CHECKPOINT style correctly (verified). No gap here.
4. **Everything else in Phase B and Phase C (sections 1-7, 9 partial, 10, 11) is present
   and was at least mechanically exercised** — the gaps are specifically the two items above,
   plus the full-scale numerical re-validation already covered in the table above.
5. **Student notebooks (both levels) — not started at all.** This was always the last step
   in the plan's own sequencing ("...then build Level 1 solution → derive Level 1 student →
   build Level 2 solution → derive Level 2 student") and the plan was interrupted before
   reaching it. Style rules for deriving them are in `docs/original_plan.md`'s "Style
   conventions" section and in `docs/build_notes.md`.

## Bugs found and fixed this session (won't recur, but know they existed)
Full root-cause writeups are in `docs/build_notes.md`'s HANDOFF section. Quick index:
- scVI adaptive epoch heuristic caused a 2+ hour hang on full-scale data — removed, replaced
  with fixed epoch counts everywhere.
- cell2location needs `layer="counts"` explicitly once `.X` is normalized — fixed in both
  notebooks and `scripts/03_benchmark_cell2location.py`.
- cell2location's real output keys differ from commonly-quoted tutorial names — fixed
  (`varm[f"{summary}_per_cluster_mu_fg"]` for `RegressionModel`, `obsm[f"{summary}_cell_abundance_w_sf"]`
  for `Cell2location`).
- infercnvpy needs `cnv.tl.pca()` → `cnv.pp.neighbors()` → `cnv.tl.leiden()` before
  `cnv.tl.cnv_score()` — the inherited draft code was missing all three; fixed.
- Level 2's PCA-copy cell was missing `adata.uns["pca"] = adata_hvg.uns["pca"]` (needed by
  `sc.pl.pca_variance_ratio`) — fixed (bug in this session's own draft, not inherited).
- `jupyter nbconvert --execute` hung on this server's Jupyter-kernel startup for unknown
  reasons (ZMQ/IPC related, not diagnosed further given time pressure) — worked around with
  `scratch_build/direct_execute_nb.py`, a from-scratch in-process cell executor (exec() each
  code cell sequentially in a persistent namespace, captures stdout + matplotlib figures,
  writes outputs back into the notebook JSON). **Try plain `nbconvert --execute` first on the
  new server — this workaround may not be needed there.**
- This server's login node kills long-running heavy compute processes that aren't submitted
  through Slurm (manifested as background jobs being silently "stopped" partway through
  cell2location training, with no error, after roughly 8-10 minutes) — fixed by always
  wrapping execution in `srun --partition=fast --cpus-per-task=8 --mem=32G --time=<generous>`.
  **Confirm whether the new server has the same restriction before assuming you can run
  directly on its login node.**

## Recommended order of operations on the new server
1. Re-confirm `single_cell` env is available/equivalent (or rebuild — see `docs/build_notes.md`'s "Environment" section for exactly what was added: squidpy, celltypist, harmonypy, decoupler, infercnvpy).
2. Swap Level 1's data path back to the real file (see build_notes.md item 1) and re-time scVI with a quick probe; pick fixed epochs.
3. Re-run Level 1 solution notebook end-to-end at full scale via Slurm. Confirm the malignant/TME split's *primary* threshold path fires (not the fallback) and sanity-check the malignant axis distribution against the paper's own (~55-65% malignant; OPC-like/NPC-neuronal-like/AC-gliosis-hypoxia/Proliferative all represented).
4. Re-run Level 2 solution notebook (cell2location FAST or FULL preset depending on GPU availability) against Level 1's real output.
5. Fill in the Level 2 paper-reveal placeholder with the real citation/findings.
6. Derive both student notebooks from the now-finalized solutions (style rules in `docs/original_plan.md`).
7. Run the final verification pass from the original plan: confirm no answer-key leakage in student notebooks, confirm checkpoint ranges in Level 1 reflect real full-scale numbers (not the tiny-scale ones currently in the demo run).
