---
name: hpc-execution-preferences
description: "How this user wants Slurm/HPC jobs run on the rgottar1_spatial Curnagl cluster — no container, short walltimes"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: e2ce4505-eebd-4dd3-8fbe-be5b78067320
---

Use native conda/mamba environments + plain `sbatch`/`srun`, not Singularity containers, for Python batch/compute work on this cluster (account `rgottar1_spatial`).

**Why:** User confirmed this directly when asked (see [[project-c10-gbm-space]]) after I found their own project's builder scripts already call `conda run -n single_cell ...` and they already maintain several other native conda envs on this cluster (anaconda3 with mamba, multiple `envs/`). The `run.sh`/`run_gpu.sh` Singularity templates on this cluster are specifically for launching an R-focused RStudio image interactively — not a general requirement for batch Python jobs.

**How to apply:** Default to native conda env creation + direct `sbatch` scripts for any new batch/compute work on this cluster unless the workload specifically needs container isolation (e.g. R version pinning via the rserver image). Don't assume Singularity is mandatory just because interactive-session templates use it.

---

When requesting Slurm `--time`, prefer shorter, tightly-measured walltimes (aim for a few hours, no more than 8-10h) over generous padding.

**Why:** User stated explicitly: shorter time requests queue faster on this cluster's shared `gpu`/`cpu` partitions than long ones, even if the job would actually finish well within a longer window.

**How to apply:** Always size `--time` from an actual timed probe (e.g. a few epochs/iterations measured directly) with modest headroom, not a pessimistic ceiling "just in case." If the probe doesn't support a confident estimate under ~8-10h, break the work into checkpointed stages rather than requesting one long allocation.
