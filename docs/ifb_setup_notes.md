# IFB Core setup notes — Session 4 (2026-07-01)

Context for whoever picks this up next. The project was developed on Curnagl (UNIL/CHUV),
then transferred to the **IFB Core cluster** (Bordeaux course allocation) where the course
actually runs. This file records what is true about the IFB environment and what was changed
this session. Read this together with `CLAUDE.md` and `full_scale_run_plan.md`.

## Where things live on IFB
- Project root: `/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10`
- Repo: `.../C10/lederer/gbm_space_proj`
- Data: `.../C10/data/{snRNA_seq,visium,xenium,paper}`; answer keys `.../C10/lederer/answer_keys`
- Student package: **`.../C10/starting_materials/`** (new this session — see below)
- Login host: `core.cluster.france-bioinformatique.fr` (fqdn seen: `core-login2…`)
- IFB user docs: https://ifb-elixirfr.gitlab.io/cluster/doc/

## Conda on IFB (important gotcha)
- conda + mamba live at `/shared/software/miniconda` (conda 25.11).
- `module load conda` **only adds conda to PATH** — it does NOT enable `conda activate`
  (fresh users get *"Run 'conda init' before 'conda activate'"*). The working recipe is:
  - `source /shared/software/miniconda/etc/profile.d/conda.sh` (per shell), **or**
  - `/shared/software/miniconda/bin/conda init bash` once (persistent).
- The `single_cell` env is at `/shared/projects/tp_2630_ubordeaux_neuromics_184418/envs/single_cell`.
  It is **not** on the default `envs_dirs`; either activate by full path or
  `conda config --append envs_dirs /shared/projects/tp_2630_ubordeaux_neuromics_184418/envs`.
- **Permissions:** the env dir carries an ACL granting `group:tp_2630_ubordeaux_neuromics_184418:rwx`
  (verified propagated to deep files, incl. `bin/python`, `site-packages`), so all project-group
  students can activate and import it. No per-student rebuild needed.

## Environment reconciliation done this session
- `single_cell` already had the full stack (scanpy 1.11.5, scvi-tools 1.4.2, cell2location 0.1.5,
  celltypist 1.7.1, infercnvpy 0.6.1, squidpy 1.8.2, harmonypy 2.0.0, torch **2.12.1+cu130**,
  decoupler, scikit-misc, leidenalg) but was **missing `liana`**, which the Level 2 solution
  imports. Installed `liana==1.7.3` (`pip install --no-user`) to match `docs/single_cell_environment.yml`.
  All L1+L2 libraries now import cleanly.

## Paths repointed (Curnagl → IFB)
- The old Curnagl root `/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/C10` was
  replaced with the IFB root in **all code**: the 4 notebooks, `scripts/*.py`, the level3
  notebooks, and 15 `scratch_build/*.py|*.sh` builder/derive/probe scripts. Docs/markdown were
  left as historical record. Backups of edited files are under the session scratch dir.
- **Do not** re-run `build_solution_nb*.py` unless you intend to regenerate; the executed
  notebooks are already repointed and finalized. If you regenerate, the builders are repointed too.

## GPU situation on IFB (as of 2026-07-01)
- The course account `tp_2630_ubordeaux_neuromics_184418` currently has associations for the
  **`fast`** and **`long`** CPU partitions only. The **`gpu`** partition rejects jobs
  ("Invalid account or account/partition combination"). GPU access is to be sorted out later.
- GPU-accelerated steps are exactly two: **scVI** integration (Level 1) and **cell2location**
  (Level 2). Both have escape hatches:
  - Level 1 also computes **Harmony** (CPU) and downstream already uses `INTEGRATION_METHOD="harmony"`.
  - Level 2 cell2location results are **checkpointed** on disk:
    `scratch_build/checkpoints/{vis_AT10_FULL.h5ad,vis_AT14_FULL.h5ad,inf_aver_FULL.parquet}`.

## Intended student workflow (confirmed with instructor)
- **CPU steps run on the FULL dataset.** GPU steps: students do a **proof-of-concept on a
  subsample**, then **load the full-scale output** (instructor is generating those on GPUs
  separately). The student notebooks already point at the full input files.
- TODO (pending instructor input): a stable path for the instructor's GPU-generated full outputs
  (L1 annotated `.h5ad` / scVI latent; L2 c2l result), and whether to add a small
  "PoC-on-subsample → else load full output from `<path>`" scaffold into the notebooks.

## Running jobs on IFB
- **Always via Slurm**, never the login node. Partition `fast` (≤24h), `long` (≤30d).
- **Slurm artifacts (script, --output, data) must be on the shared filesystem**, not the
  session `/tmp/...` scratch — `/tmp` is node-local, so a compute node can't see login-node
  paths and the job fails immediately (learned the hard way: job 410780 failed for exactly this;
  re-run from `scratch_build/ifb_smoke/` on shared FS worked).
- Working scratch for validation lives at `scratch_build/ifb_smoke/` (patch scripts + sbatch +
  logs for the CPU validation runs).

## Validation status (this session)
- **Level 1 subsample smoke** (1,500 cells, scVI skipped, Harmony): job 409959 — **PASSED**,
  44/44 code cells, 0 errors, inferCNV primary threshold fired.
- **Level 1 full-data CPU** (118k cells, scVI skipped, Harmony), job 410287: heavy cells
  (integration, inferCNV) passed with 0 errors; confirms students can run the full-data CPU
  workload on `fast` (used 16 CPU / 200 GB). Outputs redirected to `ifb_smoke/processed` so the
  real `data/processed/gbm_l1_snrna_AT10_AT14_annotated.h5ad` is untouched.
- **Level 2 CPU** (full data, load c2l checkpoints): not yet run.

## `starting_materials/` (new — the student hand-out)
`/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/starting_materials/`
- `setup_c10.sh` — one-shot: copies materials to `~/C10`, `conda init`, verifies env, registers kernel.
- `README.md`/`.html`, `INSTALL.md`/`.html` (conda + VS Code Remote-SSH + Jupyter-on-compute-node
  + Claude + Slurm). **Markdown is always also rendered to `.html`** (instructor preference) via
  `scratch_build/md2html.py`.
- `notebooks/level{1,2}/*_student.ipynb` (leak-checked), `gbmspace_utils/`, `reference/grch38_gene_positions.parquet`,
  `single_cell_environment.yml`.
- Deliberately **no `data/` copy or symlink** (the big files stay shared read-only; `data/paper/`
  holds the paper PDF which must not leak before the Level 2 reveal). Input paths are documented in README.
