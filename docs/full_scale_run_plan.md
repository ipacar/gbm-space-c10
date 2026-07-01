# Full-scale re-run of GBM-Space C10 notebooks on new cluster (native conda + Slurm GPU)

## Context

This project (CAJAL "Neuromics 2026" course, project C10) was built in a previous session on a different cluster (path signature `/shared/projects/tp_2630_ubordeaux_neuromics_184418/...`) and has now been rsync'd to this cluster at
`/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/C10/`. Both solution notebooks execute end-to-end with zero logic errors, but only at toy scale (Level 1: 1,500-cell demo subsample; Level 2: 3,999-spot AT10 section with cell2location at 5/20 "DEMO" epochs) because the old cluster's CPUs made full-scale runs impractical. The goal now is to re-run both at full scale using this cluster's GPU, with the same discipline that fixed every bug last time: probe before committing to a long run, monitor with real evidence, never let anything run on a hope.

Confirmed facts that shape this plan:
- Real full-scale data is already in place: `C10/data/snRNA_seq/level1_prepared/gbm_l1_snrna_AT10_AT14_raw.h5ad` (4.0 GB, 118,471 cells) and `C10/data/visium/level2_prepared/AT10-BRA-5-FO-1_2_student.h5ad` (3,999 spots) + `AT14-BRA-4-FO-2_1_student.h5ad`. Answer keys are present at `C10/lederer/answer_keys/`.
- **No Singularity container needed.** The project's own builder scripts already invoke `conda run -n single_cell ...` directly; the user already has a working native mamba/anaconda3 install on this cluster (`/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/anaconda3`, mamba 2.4.0) with several other bioinformatics envs built the same way. `run.sh`/`run_gpu.sh`'s Singularity usage is specific to launching the R-focused `rserver_4.4.2.sif` image for interactive RStudio/Jupyter — not something this project's Python stack needs. User confirmed: build `single_cell` natively and skip Singularity entirely.
- `docs/single_cell_environment.yml` is fully self-contained: only `python=3.11.15` + basic build tools come from conda; everything else (torch 2.12.1, scvi-tools, cell2location, all `nvidia-*` CUDA 13 wheels) is pip-installed. No system CUDA/container dependency — just needs a compatible NVIDIA driver on the GPU node.
- This cluster's Slurm: account `rgottar1_spatial`; `gpu` partition = 7 nodes × 2×A100 (`gpu:A100:2`), `MaxTime=3-00:00:00`, `DefaultTime=00:15:00` (must always pass `--time` explicitly); `cpu` partition similarly 3-day max. No login-node compute-killer was found in this exploration (unlimited ulimits) — but jobs will be submitted via `sbatch`/`srun` regardless, never run unsupervised on the login node.
- Every project script/notebook/doc still has the **old cluster's literal path** baked in: `grep -rl tp_2630_ubordeaux .` hits `README.md`, `docs/build_notes.md`, `docs/original_plan.md`, `scripts/01_prepare_snrna_subset.py`, `scripts/02_prepare_visium_subset.py`, `scripts/03_benchmark_cell2location.py`, and both solution `.ipynb` files (12 and 15 occurrences respectively) — including inside `scratch_build/build_solution_nb.py` and `build_solution_nb2.py`, the canonical generators for those notebooks. These must be repointed to `/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/C10/...` before regenerating.
- `scratch_build/direct_execute_nb.py` (in-process cell executor, used because `nbconvert --execute`'s kernel hung on the old server) is path-agnostic (`NB_PATH = sys.argv[1]`) — reusable as-is. Try plain `jupyter nbconvert --execute` first on this cluster per the build notes' own suggestion; fall back to this script if it hangs.
- The malignant/TME split fallback (`build_solution_nb.py` TASK 7.4, ~line 509-519) triggers only if **no** Leiden cluster's malignant-cell fraction exceeds 20%. At full scale with paper-expected ~55-65% malignant nuclei this should clear easily via the primary path — but must be verified from real printed output, not assumed.
- cell2location epoch presets already exist as a dict in `build_solution_nb2.py` (~line 199-202): `{"DEMO": ..., "FAST": ..., "FULL": {"ref": 400, "map": 6000}}`. Switching to `"FULL"` is a one-line change once GPU is confirmed working.

## Scheduling note

Request shorter walltimes wherever the real timed probes support it (aim for a few hours, no more than 8-10h) rather than padding `--time` generously "just in case" — shorter requests queue faster on this cluster's shared `gpu`/`cpu` partitions. Size `--time` from the actual per-epoch/per-cell timing measured in each probe, with modest headroom, not from a pessimistic ceiling.

## Step 1 — Build the `single_cell` conda env natively (no container)

1. `mamba create -n single_cell python=3.11.15 pip=26.1.2 -c conda-forge -y` (or directly via `mamba env create -f docs/single_cell_environment.yml` if mamba handles the embedded pip block — confirm pip section installs correctly; if not, create the env then `pip install -r <(extract pip list from yml)`).
2. Sanity-check **before any long job**: `salloc -p gpu --gres=gpu:1 --time=00:15:00 --account=rgottar1_spatial` → inside, run `nvidia-smi` (confirm driver/A100 visible) and `source activate single_cell && python -c "import torch, scvi, cell2location; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"`. This single check resolves the only real risk in skipping a container: whether the pip-bundled CUDA 13 wheels are compatible with this node's driver. If incompatible, the fallback is a driver-compatible torch/cuda-toolkit pip pin — diagnose from the actual error before changing anything.
3. Confirm `scanpy`, `harmonypy`, `celltypist` (with the `Developing_Human_Brain` model downloaded/cached), `infercnvpy`, `squidpy` import cleanly too.

## Step 2 — Repoint every stale path

Across the 8 files found by `grep -rl tp_2630_ubordeaux .` inside `gbm_space_proj/`, replace the root
`/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10` → `/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/C10`. Prioritize the load-bearing ones first since they're what actually gets executed:
- `scratch_build/build_solution_nb.py`: `OUT`, `sys.path.insert(...)`, the `DATA = ".../tiny_snrna_1500.h5ad"` line (TASK 1.1, ~line 71) → point at the real `level1_prepared/gbm_l1_snrna_AT10_AT14_raw.h5ad`, and the gene-position parquet path (TASK 7.1, ~line 451).
- `scratch_build/build_solution_nb2.py`: `OUT`, `sys.path.insert(...)`, `VISIUM` path (~line 70), `ref_path` (~line 189), `answer_key` path (~line 371).
- Docs/README (`README.md`, `docs/build_notes.md`, `docs/original_plan.md`) — cosmetic, lower priority, but fix for consistency since future sessions will read them.
- Leave the already-built (stale, demo-scale) `.ipynb` files alone — they get regenerated fresh from the corrected builder scripts, not hand-edited.

## Step 3 — Level 1 full-scale run

1. In `build_solution_nb.py`, set `SCVI_MAX_EPOCHS` based on a **real timed probe on the actual 118,471-cell data**, not a guess: run a short script that does `scvi.model.SCVI.setup_anndata(...)`, `model.train(max_epochs=3, early_stopping=False)` on GPU, time it, compute per-epoch cost, then pick a fixed epoch count that fits a teaching-reasonable GPU runtime (tens of minutes at most, since GPU obviates the CPU 5-epoch compromise). Do **not** call `get_max_epochs_heuristic()` for the actual training decision (it can stay as a printed discussion value, as it already is in the notebook, since it's never used to set `max_epochs`).
2. Submit via a new sbatch script (e.g. `scratch_build/sbatch_level1.sh`) modeled on `run_gpu.sh`'s SBATCH header (`--account rgottar1_spatial`, `--partition gpu --gres gpu:1 --gres-flags enforce-binding`, generous `--cpus-per-task`/`--mem`, and a `--time` budget sized off the probe from step 3.1) but stripped of all Jupyter/port/SSH-tunnel logic — just:
   ```bash
   source ~/anaconda3/etc/profile.d/conda.sh   # or mamba.sh
   conda activate single_cell
   python -u scratch_build/build_solution_nb.py     # regenerate notebook with corrected paths/epochs
   python -u scratch_build/direct_execute_nb.py notebooks/level1/01_snrna_analysis_solution.ipynb 1800
   ```
   Write a dedicated progress log via explicit `open()/write()/flush()` (the project's own lesson — don't trust buffered stdout through Slurm capture); `direct_execute_nb.py` already does per-cell timing+flush prints and saves the notebook after every cell, which doubles as a checkpoint file — `cat` it mid-run, or watch the `.ipynb`'s mtime, instead of just checking "still running."
3. Check in against evidence: per-cell timestamps advancing in the saved notebook / log, and `sstat -j <jobid> --format=AveCPU,MaxRSS,Elapsed` showing CPU time still accruing. A stall longer than the slowest individual step so far (scVI training, expected to be the longest single cell) means something is hung — investigate via real evidence (read scvi-tools source, check `nvidia-smi` GPU utilization) before re-launching.
4. Once complete, verify directly from the printed output (don't assume): TASK 7.4's `frac_mal`/`malignant_clusters` print confirms the **primary** 20%-cluster-fraction path fired (not the "no cluster exceeded" fallback warning), and the malignant fraction lands roughly in the paper's 55–65% range (CHECKPOINT in the notebook already states this). If the fallback fires anyway at full scale, diagnose root cause (clustering resolution too coarse/fine, reference cluster selection) before accepting the result.
5. Confirm `data/processed/gbm_l1_snrna_AT10_AT14_annotated.h5ad` is freshly written (Level 2's input).

## Step 3b — Scope upgrade: LIANA CCC and AT14 secondary comparison now in-scope

User decided both of `original_plan.md`'s previously-optional Level 2 extensions should actually be built into the notebooks, not left as stubs:

- **LIANA cell-cell communication** (Section 8, currently `# import liana as li` commented out). Confirmed via [liana-py](https://github.com/saezlab/liana-py) / readthedocs basic-usage docs: pip package name is `liana` (`import liana as li`), and the standard consensus call is
  ```python
  li.mt.rank_aggregate(adata, groupby="<group_key>", resource_name="consensus", expr_prop=0.1, use_raw=True)
  # results land in adata.uns["liana_res"]; li.pl.dotplot(adata, colour="magnitude_rank", size="specificity_rank", inverse_size=True, inverse_colour=True) to plot
  ```
  `use_raw=True` expects log1p-normalized counts in `.raw` (Visium `vis` object already has this from Section 2/QC). LIANA's spot-based/spatial resolution is natively supported, so `groupby` can be a per-spot category — exactly the plan's original framing: define two spot groups from the niche output (TASK 8.1 already sketches this) — `"dev-like-dominant"` vs `"gliosis/hypoxia-dominant"` niches, using the niche→cell-type loading heatmap from Section 6 — then run `rank_aggregate` comparing them and plot/report the top-ranked interactions. Run this on the AT10 primary section (matches the original plan's framing and keeps the LIANA addition scoped to where the rest of the narrative already lives); no separate epoch/heuristic concerns since LIANA is a single deterministic scoring pass, not iterative training — just wrap in the same per-cell logging discipline as everything else.
- **AT14 secondary Visium comparison** ("does the same pattern hold in tumour 2?", `AT14-BRA-4-FO-2_1_student.h5ad`, 3,534 spots, no IvyGAP histopath overlay). Add as a new section before the existing Write-up (renumbering Write-up to follow it): repeat the condensed core pipeline — spatial QC + normalize + naive clustering, axis-in-space scoring (reuse `score_axis`/`ZONATION_PANEL`), and cell2location mapping (reusing the **same** Level 1 reference signatures already fit for AT10 — only the spatial-mapping model needs separate training per section, not the reference model) and niche NMF — then a comparison subsection: do the same niche structure and `AQP4→ABCC3→AKAP12→HILPDA` axis gradient appear in AT14, and how does its niche composition compare to AT10's. Explicitly scoped *condensed*, not a full duplicate of AT10's teaching depth: skip re-doing the squidpy-vs-k-d-tree neighborhood method comparison and LIANA for AT14 (those stay AT10-only) to avoid bloating Level 2 — flag this scoping choice back if full parity is wanted instead. Since AT14 has no IvyGAP/cell2location answer key, this comparison is AT10-vs-AT14 cross-tumor only, not vs. ground truth (TASK 9.1's answer-key comparison stays AT10-only, unaffected).
- **Env**: add `liana` to the Step 1 pip install list (not in the current `single_cell_environment.yml` — confirmed by grep — needs adding).
- **Scheduling**: per the "shorter `--time` requests queue faster" preference, run AT10's and AT14's cell2location mapping as two separate, shorter sbatch submissions rather than one long combined job — AT14 has fewer spots (3,534 vs AT10's 3,999) so a comparably-sized or shorter probe-timed budget applies.
- **Student derivation impact**: confirm `scratch_build/derive_student_nb.py` (Level 2's generic blanking script) handles the new LIANA and AT14 sections the same generic way as everything else (blank code cells, keep TASK/HINT/QUESTION markdown) — it should, since it isn't section-specific, but verify after running it rather than assuming.

## Step 4 — Level 2 full-scale run (cell2location FULL)

1. In `build_solution_nb2.py`, set `C2L_MODE = "FULL"` (paper-exact: ref `max_epochs=400`, mapping `max_epochs=6000` — already wired via the existing preset dict, no heuristic involved here since these counts are fixed by the user's explicit instruction, not adaptive).
2. Before committing to the full submission, time a short probe on GPU (3-5 epochs each for the reference and mapping models, same spirit as the scVI probe) to compute real per-epoch GPU cost and set a realistic `--time` budget — the old CPU numbers (~72s/epoch ref, ~3.9s/epoch mapping) are CPU-only and will not transfer; GPU should be far faster but must be measured, not assumed.
3. Submit via a similar sbatch script (`scratch_build/sbatch_level2.sh`), `--gres gpu:1`, regenerate via `build_solution_nb2.py` then execute via `direct_execute_nb.py` with the same explicit-flush logging/checkpoint-watching discipline as Level 1.
4. Fill in the Level 2 "paper reveal" markdown placeholder (`[INSTRUCTOR: insert citation + 3-4 bullet summary...]`) with the real citation (de Jong, Memi, Gracia, Lazareva et al., bioRxiv 2025.05.13.653495) and findings, per `build_notes.md`'s HANDOFF note that this still needs filling.
4b. **Cross-checked against `docs/original_plan.md` and confirmed two real gaps unique to the Level 2 solution** (Level 1's solution has every task fully worked through with real output; Level 2's does not, for these two): TASK 9.1 (~line 371 of `build_solution_nb2.py`) currently only loads the answer key and leaves the actual niche/axis comparison as a commented-out suggestion — replace with real comparison code (e.g. correlate `niche_loadings`/`vis.obs['niche']` against `answer_key`'s "Spatial niche abundances" rows) and real printed/plotted output. TASK 10.1 (~line 383, "Write-up") is an empty `# Your figure-reproduction code here.` placeholder — fill with an actual reproduction of one specific published panel (e.g. the zonation gradient or niche map) plus a short comparison paragraph, so the answer key is a complete exemplar, not a student-style blank. Do this before Step 5 derives the student notebook (the derivation script blanks code cells anyway, so filling the solution here doesn't leak anything).
5. Verify cell2location's known footguns are still respected after the path edits: `layer="counts"` passed explicitly to both `RegressionModel.setup_anndata` and `Cell2location.setup_anndata` (already present, ~lines 216/227 — just confirm the edits didn't disturb them), and that output keys (`varm[f"{summary}_per_cluster_mu_fg"]`, `obsm[f"{summary}_cell_abundance_w_sf"]`) are read correctly downstream.

## Step 5 — Re-derive student notebooks

1. Update the same stale-path root in `scratch_build/derive_student_nb1.py` (Level 1) and its Level 2 counterpart (find via `grep -rl tp_2630_ubordeaux scratch_build/derive_student*`), including the `HINT_SOFTENING` dict's literal `DATA = "..."` student-facing path replacement (~line 30 of `derive_student_nb1.py`) — point it at the real `level1_prepared` path, not the old cluster's.
2. Run the derivation scripts against the now-finalized full-scale solution notebooks.
3. Re-verify zero answer-key leakage: grep both student notebooks for the answer-key column names (`src/gbmspace_utils/data.py::SNRNA_ANSWER_KEY_OBS_COLUMNS`), for the paper's name/citation/authors, and for any of the hardcoded computed values (epoch counts, malignant fractions, cluster counts) that should be blanks rather than revealed answers — the derivation script's existing `HINT_SOFTENING` mechanism handles known cases, but a fresh full-scale run may print new specific numbers in markdown HINTs that weren't anticipated by the demo-scale derivation logic.

## Status as of this update

- **Step 1 (env)**: done. `single_cell` built natively, GPU/CUDA confirmed (Quadro RTX 8000 on `gpu-rtx` partition — far less contended than `gpu`/A100 nodes, jobs start in seconds; switched all sbatch scripts to `--partition gpu-rtx`).
- **Step 2 (paths)**: done, all load-bearing files repointed.
- **Step 3 (Level 1)**: done and verified. Along the way, found and fixed two real bugs: `infercnvpy`'s `window_size` had drifted from the documented value (100→250, matching `build_notes.md`'s own "paper-faithful parameters" — and matching the paper's actual Methods, confirmed: "window size of 250 genes"), and — bigger — the malignant/TME split's absolute signal threshold (`cnv_score > 0.02`) didn't transfer to this pipeline: `cnv.tl.cnv_score()`'s cluster-broadcast design compressed signal into a narrow 0.005–0.02 band for *all* cells. Fixed by computing genuine per-cell signal directly (`mean(|X_cnv|)`, bypassing the `cnv_leiden` clustering step entirely — also ~20min faster per run) and calibrating `SIGNAL_CUT` as the 75th percentile of the known-diploid reference cells' *own* signal distribution, rather than trusting an external absolute number. Result: **59.4% malignant**, primary threshold path, matching the paper's expected 55-65% range. Documented in `build_notes.md`.
- **Step 3b (LIANA + AT14)**: code written, verified against actual installed library signatures (`li.mt.rank_aggregate`, `li.pl.dotplot`) via `inspect.signature`, not guessed.
- **Step 4 (Level 2)**: ran successfully end-to-end after fixing one bug (a missing import caught after a costly ~3h21min retrain — added real model checkpointing to `build_solution_nb2.py` afterward so a future failure past the training stages doesn't repeat that cost: `inf_aver` and each `vis`/`vis14` get saved to `scratch_build/checkpoints/`, keyed by `C2L_MODE`, and reused if present). **However**, inspecting the actual output surfaced a deeper problem (see below) — the results from this run should not be treated as final.

## Step 4b — Cell2location produced spatially flat, uninformative output: root cause found, fix planned

**The problem, found by actually looking at the output, not just checking for errors:** the per-spot cell-type abundance maps and the NMF niche map are pure speckle — no spatial coherence at all (visually confirmed from the embedded figures). Cross-checked: training ELBO loss dropped substantially (1.22e8 → 5.7e7) over the 6000 epochs, so the model *did* converge to something — it just isn't spatially meaningful. The within-spot column-wise std of the abundance matrix is only ~3-8% of the mean, far too tight for genuine biological heterogeneity.

**Root cause, confirmed with real evidence (not assumed):**
1. Computed the actual pairwise correlation between the 15 `cell_type` reference categories' mean expression profiles (cheap, via sparse one-hot matrix multiply on the already-saved Level 1 h5ad, no retraining needed). The 10 "malignant-mimic" categories (CellTypist's region-specific labels like "Hypothalamus glioblast", "Pons OPC", "Dorsal midbrain glioblast") are **mean 0.80 correlated with each other, up to 0.96-0.98 for some pairs** — they are not biologically distinct cell types, just the same underlying malignant population relabeled by which brain region CellTypist's `Developing_Human_Brain` model happened to match it to. cell2location has essentially no orthogonal signal to resolve between them.
2. Read the paper's actual Methods section directly (`data/paper/GBM-Space-Paper.pdf`) to check what reference categories and parameters they really used. Confirmed: "We mapped **malignant and TME clusters** from each patient's snRNA-seq profile to their matched Visium data using cell2location" — they used genuine, distinct malignant cell-*state* clusters (their 4-class/9-subclass hierarchy) and TME clusters, never the kind of region-confounded developmental-brain-atlas mimic labels CellTypist produces. This is exactly what our own `malignant_state` (9 categories, matching the paper's 9 subclasses almost exactly: OPC-like, OPC-NPC-like, OPC-neuronal-like, NPC-neuronal-like, AC-progenitor-like, AC-gliosis-like, Gliosis-like, Hypoxic, Proliferative) already represents — it was sitting unused in Level 1's saved output the whole time.
3. Same Methods section also confirms the mapping model's `batch_size` should be **~25% of total spots per tumour** — we used full-batch (100%) instead, reasoning at build time that full-batch was "defensible" for a small (~4,000-spot) section. Combined with `detection_alpha=200`'s strong regularization, full-batch (deterministic gradient every step) plausibly converges to a smoother/flatter solution than the paper's stochastic minibatching. Every other parameter (`max_epochs=400/6000`, `N_cells_per_location=30`, `detection_alpha=200`, `lr=0.002`) already matched the paper exactly — confirmed by reading the Methods text, not just trusting `build_notes.md`'s transcription of it.

**Fix — implemented and validated (approved by user):**
1. `build_solution_nb2.py` TASK 5.1b: combined reference label on `adata_ref` — `malignant_state` (9 marker-defined axis states) for malignant cells, `cell_type` for TME cells, used as `labels_key` instead of plain `cell_type`.
2. **Refinement found while validating**: naively dropping non-malignant nuclei whose `cell_type` looks malignant-mimic (CNV signal below threshold) would have discarded ~26k real cells — confirmed via `MBP/PLP1/MOG/MOBP/ST18` marker scoring that ~23k of those are genuine oligodendrocytes (CellTypist's `Developing_Human_Brain` model has no clean adult-oligo category, so it mismatched them). Relabelled those `"Oligodendrocyte"` instead of dropping; only ~3,300 truly ambiguous cells get dropped.
3. Mapping model `batch_size` changed from `vis.n_obs` (100%, full-batch) to `int(0.25 * vis.n_obs)`, matching the paper's Methods exactly.
4. **Validated cheaply in FAST mode** (20 ref / 300 mapping epochs, ~12 min on GPU, via `scratch_build/validate_c2l_fix.py`) before committing to the full ~5-6h FULL rerun: abundance-map coefficient of variation went from ~3-8% (flat, broken) to ~20-36%, and the spatial maps show real, coherent regional structure (visually confirmed), not speckle. Ran this validation twice — once before the oligodendrocyte fix (already showed strong improvement) and once after (confirmed the fix doesn't regress anything, Oligodendrocyte's own CV sits mid-pack at 0.255, neither flat nor dominant).
5. **Now running for real**: `C2L_MODE = "FULL"`. No stale checkpoints to worry about (the broken first run predates the checkpointing code, so `scratch_build/checkpoints/` was empty going in).
6. Once complete: re-verify TASK 9.1's answer-key niche correlation and the AT10-vs-AT14 comparison with real numbers — the previous run's captured numbers were provisional/wrong and must not be reused in the write-up.

## Cluster GPU congestion mid-run: had to switch launch method

Right when the FULL fix was ready to submit, the whole cluster's GPU pool hit a `QOSGrpGRES` cap (a cluster-wide concurrent-GPU limit, not specific to any one partition — confirmed by seeing *other users'* `sbatch` jobs stuck pending with the identical reason on `gpu`, `gpu-rtx`, and `gpu-h100` simultaneously). Resolution, found by testing directly rather than guessing:
- `sbatch` jobs on `gpu`/`gpu-rtx`/`gpu-h100` were all stuck behind this cap.
- The `interactive` partition (node `dnagpu001`, MIG-sliced A100s: `gpu:3g.20gb:4`) is **exempt** from that cap and grants instantly.
- **Gotcha**: `salloc --gres=gpu:3g.20gb:1 <command>` alone does **not** bind the GPU device to the command on this cluster (`CUDA_VISIBLE_DEVICES` empty, `/dev/nvidia*` missing) — you must nest an explicit `srun --gres=gpu:3g.20gb:1 <command>` *inside* the `salloc` allocation for the device to actually bind. Confirmed working pattern:
  ```bash
  salloc --partition=interactive --gres=gpu:3g.20gb:1 --account=rgottar1_spatial --time=08:00:00 --cpus-per-task=4 --mem=30gb \
    srun --gres=gpu:3g.20gb:1 bash -c "source .../conda.sh && conda activate single_cell && cd .../gbm_space_proj && python -u scratch_build/build_solution_nb2.py && python -u scratch_build/direct_execute_nb.py notebooks/level2/02_spatial_cell2location_solution.ipynb 14400" \
    > scratch_build/level2_run_slurm.out 2> scratch_build/level2_run_slurm.err
  ```
  Run via Bash with `run_in_background: true` (it's a long-lived foreground process under `salloc`, not a detached `sbatch` job — there's no separate job-monitoring needed beyond watching the log file and `squeue -u alederer`).
- **`interactive` QOS resource caps** (from `sacctmgr show qos interactive`): `MaxTRESMins`: `cpu=1920`, `gres/gpu=480`, `mem=15T` (per-job CPU-minutes / GPU-minutes / memory-minutes products). At `--time=08:00:00` (480 min, the partition's `MaxWall`), this caps `--cpus-per-task` at 4 (1920/480) and `--mem` at ~30gb (15,000,000 MB-min / 480 min ≈ 31,250 MB). Hit both caps by trial before landing on `--cpus-per-task=4 --mem=30gb`.
- **Risk accepted**: the MIG slice (`3g.20gb`, ~3/7 of an A100, 20GB) is slower per-epoch than the full RTX 8000/A100 GPUs benchmarked on, and `interactive`'s hard `MaxWall=08:00:00` is tighter than the ~5h41m the previous (broken-label) FULL run actually took on a full GPU. If this run gets cut off by the 8h wall-time before finishing, the checkpointing in `build_solution_nb2.py` (saves `inf_aver` and each `vis`/`vis14` to `scratch_build/checkpoints/`, keyed by `C2L_MODE`) means a second `salloc`+`srun` leg with the identical command will skip already-completed training stages rather than redoing them — check `ls scratch_build/checkpoints/` before relaunching to confirm what's already saved.
- Job 61851814 launched and confirmed RUNNING on `dnagpu001` as of this writing; monitor `scratch_build/level2_run_slurm.out` for per-cell progress (same `direct_execute_nb.py` logging as before).

## Verification

- After Step 1: `nvidia-smi` + `torch.cuda.is_available()` inside an interactive GPU salloc session — confirms the no-container approach actually works on this cluster before committing to long batch jobs.
- After Step 3: read the actual printed TASK 7.4 output from the executed notebook (cluster malignant fractions, which path fired) and the CHECKPOINT's expected ~55-65% range — both must be confirmed from real output, not assumed from the plan.
- After Step 4: confirm cell2location converged sensibly (loss curve / ELBO not still dropping sharply at the last epoch) and that niche analysis / spatial sections downstream still run without the cryptic GammaPoisson error.
- After Step 3b's additions: confirm LIANA's `adata.uns["liana_res"]` is populated with real ranked interactions (not empty/all-NaN from a mismatched `groupby` key or wrong `.raw` state), and that the AT14 secondary section's cell2location mapping actually converges (same ELBO sanity check as AT10) before writing the cross-tumor comparison.
- After Step 5: grep-based leakage check on both student notebooks, plus a quick scan that all the original plan's section structure (TASK/HINT/QUESTION/CHECKPOINT) survived intact.
- Throughout: every long-running step gets its own log file with explicit `open()/write()/flush()`, a wall-clock `--time` budget set from a real timed probe (never a guess or an heuristic), and an explicit check-in plan (CPU-time delta via `sstat`, or the saved notebook's per-cell progress) — if either check shows no progress for materially longer than the slowest prior step took, kill and diagnose before relaunching.

---

## Session 3 update (2026-06-30)

### What happened in the previous (Session 2) Level 2 run

The FULL fix run (job 61851814, `salloc interactive`, `dnagpu001`) got through:
- Cell 30: AT10 reference model training — OK, ~42.6 min ✅
- Cell 32: AT10 mapping model training — OK, ~2.53h ✅ → saved to `scratch_build/checkpoints/vis_AT10_FULL.h5ad`
- Cells 34–58: NMF niches, LIANA, answer-key loading (TASK 9.1), AT14 QC/normalize — all OK ✅
- Cell 60: AT14 cell2location mapping — **FAILED** with `IndexError: Dimension specified as 0 but tensor has no dimensions`

Also saved: `scratch_build/checkpoints/inf_aver_FULL.parquet`

### Root cause of cell 60 failure and fix

AT14's spot count after QC produced `n_obs % int(0.25 * n_obs) == 1` — the last mini-batch has exactly 1 spot. PyTorch squeezes a 1-element batch dimension to a 0-d scalar; cell2location then tries `tensor[0]` on it and crashes. AT10 avoids this because 3,999 spots gives a remainder ≥ 2.

**Fix applied** in `scratch_build/build_solution_nb2.py` ~line 575:
```python
_bs14 = max(int(0.25 * vis14.n_obs), 2)
if vis14.n_obs % _bs14 == 1:  # last mini-batch of 1 sample → 0-d tensor crash in PyTorch
    _bs14 += 1
sp_model14.train(max_epochs=MAP_EPOCHS, batch_size=_bs14)
```

### Missing conda env dependencies discovered and fixed

Running `direct_execute_nb.py` from the login node (not from the GPU node) surfaced missing packages in the `single_cell` env that somehow worked before (likely because the GPU node had them system-wide):
- `zarr` (required by both `anndata` and `squidpy`)
- `pooch` (required by `squidpy.datasets`)
- `dask` (required by `squidpy`)
- `spatialdata` (required by `squidpy.datasets._downloader`)

Installed all four to the correct location (`anaconda3/envs/single_cell/lib/python3.11/site-packages`) using `pip install --no-user` (see below for the `--no-user` rule). All verified with `pip show <pkg> | grep Location`.

**Important pip rule**: always use `--no-user` when installing into a conda env on this cluster. Installing without that flag can write to `~/.local/lib/python3.11/site-packages/`, which is shared across ALL environments and caused breakage in other envs in a prior session.

### Level 3 notebooks discovered and their status

Two Level 3 notebooks were found in `notebooks/level3/`:
- `03_xenium_organoid_analysis.ipynb` (2.8 MB, fully executed 17/17 code cells): GBM organoid Xenium analysis — load 3 sections (AT410, patient-derived organoids, day 7 and day 14), QC, clustering, malignant-state-axis scoring, spatial viz, cross-modality comparison (snRNA/Visium/Xenium), nuclear vs cytoplasmic spillover (TASKS 7.1-7.5). No separate student version.
- `04_xenium_segmentation_basis_solution.ipynb` (solution-only, partially executed ec=1-7): BASIS (Bayesian Assignment for Spatial Inference and Segmentation) run on a single Xenium section. Preprocessing complete; `basis.tl.segment` call and all evaluation cells (6.1-6.4) not yet run. Skipped for now.
- BASIS package: installed in `single_cell` env (v0.3.2), local copy at `gbm_space_proj/BASIS/`.
- Xenium data: `C10/data/xenium/` — 3 h5ad files (pre-annotated) and 3 raw Xenium Ranger output folders for organoid samples AT410-BRA-5-ORG-{E30/e25}-MOI0_25-{D7/D14}.

### Student notebook status

- **Level 1 student** (`notebooks/level1/01_snrna_analysis_student.ipynb`): **re-derived this session** from the full-scale Level 1 solution. 88 cells (44 md, 44 code). Leakage check passed (only `117,200` appears in a HINT explaining calibration methodology — acceptable context, not an analytical answer). `scratch_build/derive_student_nb1.py` was updated: ROOT and DATA path both fixed from old cluster path → new cluster path.
- **Level 2 student** (`notebooks/level2/02_spatial_cell2location_student.ipynb`): **still stale** — derived Jun 28 from an earlier solution version (53 cells, missing LIANA/AT14/write-up sections). Must be re-derived after Level 2 solution is finalized.

### Current run

Job 61863377 (`salloc interactive`, `dnagpu001`, MIG slice `3g.20gb`) launched 2026-06-30 ~11:54 and is **actively running**. Monitor with:
```bash
tail -f scratch_build/level2_run_slurm.out   # per-cell timing
cat scratch_build/level2_run_slurm.err        # errors / salloc messages
squeue -u alederer                             # job still alive?
ls -lh scratch_build/checkpoints/             # checkpoint files saved so far
```

AT10 loads from `scratch_build/checkpoints/vis_AT10_FULL.h5ad` (skip ~2.5h retraining). AT14 training is the remaining expensive step (~2h on MIG slice, extrapolating from AT10's 2.53h × 3534/3999 spots). The run will self-checkpoint AT14 to `scratch_build/checkpoints/vis_AT14_FULL.h5ad` on completion.

**If the run gets cut off by the 8h wall limit:** check `ls scratch_build/checkpoints/`. If `vis_AT14_FULL.h5ad` is NOT there, AT14 didn't finish — re-launch with the identical salloc command and it will load both AT10 and the AT10 `inf_aver` from checkpoints, skip straight to AT14 training.

The salloc command to re-launch:
```bash
PROJECT=/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/C10/lederer/gbm_space_proj
salloc --partition=interactive --gres=gpu:3g.20gb:1 --account=rgottar1_spatial --time=08:00:00 --cpus-per-task=4 --mem=30gb \
  srun --gres=gpu:3g.20gb:1 bash -c "
    source /work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/anaconda3/etc/profile.d/conda.sh && \
    conda activate single_cell && \
    cd ${PROJECT} && \
    python -u scratch_build/build_solution_nb2.py && \
    python -u scratch_build/direct_execute_nb.py notebooks/level2/02_spatial_cell2location_solution.ipynb 14400
  " \
  > ${PROJECT}/scratch_build/level2_run_slurm.out 2> ${PROJECT}/scratch_build/level2_run_slurm.err
```

### After Level 2 run completes: remaining steps

1. **Fill Level 2 write-up placeholder** (last markdown cell in `02_spatial_cell2location_solution.ipynb`, currently `[INSTRUCTOR: write-up paragraph filled in after the full-scale execution...]`): write 2-3 sentences comparing the spatial malignant-class map to the paper's Fig. 2/3 zonation claim, using the real TASK 9.1 correlation values and TASK 10.3 class fractions printed by the notebook. These numbers are in the notebook's cell outputs — read them directly, don't invent.

2. **Re-derive Level 2 student notebook**: run `scratch_build/derive_student_nb.py` (the Level 2 script — note it's `derive_student_nb.py` not `derive_student_nb2.py`; verify the exact filename first). Before running, check for any old-cluster path references inside it (`grep -n tp_2630_ubordeaux scratch_build/derive_student_nb.py`). Expected output: ~53 → ~66 cells to match the solution (LIANA, AT14, write-up sections now included).

3. **Leakage check Level 2 student**: grep for answer-key column names (`Spatial niche abundances`, `Cell state abundances`), paper citation/authors (`de Jong`, `Lazareva`, `bioRxiv`, `gbmspace.org`), and any specific computed numbers (TASK 9.1 correlation values, AT14 class fractions). Same pattern as Level 1 check above.

4. **Level 3 student notebook (optional)**: `03_xenium_organoid_analysis.ipynb` has TASK/HINT/QUESTION/CHECKPOINT structure with full solution code — a student version would blank code cells. Not started; confirm with user before creating.

### Session 3 — run completed and post-run checklist DONE (2026-06-30 PM)

Job 61863377 finished cleanly: cell 60 (AT14 cell2location) ran 7551s (~2.1h), then all cells passed — `=== ALL CELLS EXECUTED SUCCESSFULLY ===`. `vis_AT14_FULL.h5ad` saved 14:09. The batch_size-of-1 fix held.

Post-run checklist (Steps 1–3 of "remaining steps" above) now complete:

1. **AT14 spatial-coherence sanity check** (`scratch_build/check_at14_coherence.py`): per-class coefficient of variation across spots, both sections. AT14 median CV **0.931** (range 0.55–1.57), AT10 median **1.060** — far above the broken first run's flat 0.03–0.08 speckle, confirming the label/batch_size fix produced genuinely heterogeneous, spatially-coherent maps. Most-variable classes (Hypoxic, Brain vascular cells) are biologically sensible in both tumours.

2. **Write-up filled** (solution cell 65 / `build_solution_nb2.py` ~line 649), using real printed numbers — *not* invented:
   - TASK 9.1: AT10 NMF niches vs answer-key niches over 3,928 shared spots — Gliosis r=0.70, Vasculature r=0.73, Immune/TAM r≈0.45–0.50, plus dev-like(OPC) / grey-/white-matter matches.
   - TASK 10.3: dominant per-spot malignant class is AC-gliosis-hypoxia in 99.6% of AT10 and 97.9% of AT14 spots; both resolve into 8 NMF niches.
   - The write-up is wrapped in `<!-- INSTRUCTOR-ONLY -->…<!-- /INSTRUCTOR-ONLY -->` sentinels, preceded by a student-facing **TASK 11.1** prompt. `derive_student_nb.py` was extended to strip INSTRUCTOR-ONLY regions from markdown (generic, mechanical) — and `build_solution_nb2.py` updated to emit the prompt + sentinel-wrapped placeholder so a future regen keeps the structure (the answer itself stays a hand-filled-after-run placeholder in the builder, since it depends on post-execution numbers).

3. **Level 2 student notebook re-derived** (`02_spatial_cell2location_student.ipynb`): 66 cells (36 md, 30 code), matching the solution. `derive_student_nb.py` had no stale paths.

4. **Leakage check passed**: no answer-key column names (`Spatial niche abundances`, `Cell state abundances`), no computed numbers (correlations / class fractions), instructor write-up block stripped, all 30 code cells blanked to the placeholder. Paper citation/authors appear ONLY in the intended Section 9 "Revealing the paper" reveal cell (and a TASK 7.2 method reference) — by design, identical to the solution, not a leak.

**Level 2 is now fully finalized.** Only remaining optional item: Level 3 student notebook (needs user go-ahead).
