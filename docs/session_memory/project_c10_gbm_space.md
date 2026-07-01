---
name: project-c10-gbm-space
description: "GBM-Space C10 teaching project (CAJAL Neuromics 2026) — location, cross-cluster migration, and current full-scale-rerun status"
metadata: 
  node_type: memory
  type: project
  originSessionId: e2ce4505-eebd-4dd3-8fbe-be5b78067320
---

Project lives at `/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/C10/lederer/gbm_space_proj/` (docs/, notebooks/, src/gbmspace_utils/, scripts/, scratch_build/). Data lives one level up at `C10/data/` (snRNA_seq/level1_prepared, visium/level2_prepared, xenium/) and answer keys at `C10/lederer/answer_keys/`.

**Why:** Built in an earlier session on a different cluster (path signature `/shared/projects/tp_2630_ubordeaux_neuromics_184418/...`, likely a Bordeaux-hosted course allocation), then rsync'd to this CHUV/UNIL Curnagl cluster (account `rgottar1_spatial`) because that cluster's CPUs couldn't run full-scale workloads (118,471-cell scVI, paper-exact cell2location epochs). Both solution notebooks ran clean only at toy scale there; full-scale validation was deferred to this cluster, which has GPU (`gpu` partition, 7 nodes x 2xA100).

**How to apply:** The old-cluster path has been repointed in all load-bearing files. The full execution plan including all bugs found, fixes applied, and current run status is at `docs/full_scale_run_plan.md` — always read the "Session 3 update" section at the bottom before re-planning. See [[hpc-no-container-native-conda]], [[hpc-job-scheduling-walltime]], and [[feedback-pip-no-user]].

**Current state as of 2026-06-30:**
- Level 1 solution: fully executed at full scale ✅. Level 1 student: re-derived, leakage check passed ✅.
- Level 2 solution: **fully finalized ✅** (job 61863377 completed clean 2026-06-30; AT14 cell2location ~2.1h; both `vis_AT10_FULL.h5ad`/`vis_AT14_FULL.h5ad` checkpointed). AT14 maps verified spatially coherent (median CV 0.93, vs the original broken run's flat ~0.05) via `scratch_build/check_at14_coherence.py`. Write-up cell filled with real TASK 9.1 correlations (Gliosis r=0.70, Vasculature r=0.73) and TASK 10.3 fractions (AC-gliosis-hypoxia 99.6%/97.9% in AT10/AT14).
- Level 2 student: **re-derived ✅** (66 cells), leakage check passed. NOTE: instructor answer text that must not leak is now wrapped in `<!-- INSTRUCTOR-ONLY -->…<!-- /INSTRUCTOR-ONLY -->` sentinels in the solution; `derive_student_nb.py` strips those regions from markdown (the script otherwise keeps all markdown, blanks all code). The paper citation legitimately appears in the Section 9 "Revealing the paper" reveal cell — that is by design, not leakage.
- Level 3 notebooks: `notebooks/level3/03_xenium_organoid_analysis.ipynb` (Xenium organoid analysis, fully executed) and `04_xenium_segmentation_basis_solution.ipynb` (BASIS segmentation, partially executed — skipped for now). No Level 3 student notebook yet. Xenium data: 3 organoid samples (AT410, patient GBM organoids, day 7 and day 14) in `C10/data/xenium/`. BASIS v0.3.2 installed in `single_cell` env; local copy at `gbm_space_proj/BASIS/`.
