# GBM-Space Computational Mini-Project (C10)

Instructor working repo for CAJAL "Neuromics 2026" project **C10**, built from real
snRNA-seq + Visium spatial transcriptomics data from the GBM-Space study
(de Jong, Memi, Gracia, Lazareva et al., bioRxiv 2025; [gbmspace.org](https://www.gbmspace.org/)).

Structured the same way as the `cajal_comp_proj` example template (see
`../cajal_comp_proj/`) — see that repo's README for the general pedagogical philosophy
(progressive withdrawal of scaffolding, paper withheld in Level 1, students write all the
code). This README only covers what differs for this specific project.

## Levels

- **Level 1** (`notebooks/level1/`, ~2 days): standard snRNA-seq analysis (scanpy) on
  donors AT10+AT14 (118,471 cells, no subsampling) — QC, normalization, integration
  (Harmony *and* scVI, compared), clustering, cell type annotation (markers + CellTypist),
  malignant/TME calling, and the paper's 4-state malignant cell axis (OPC-NPC-like /
  NPC-neuronal-like / AC-gliosis-hypoxia / Proliferative) via marker-gene scoring.
- **Level 2** (`notebooks/level2/`, ~2.5 days): Visium spatial analysis on the matched AT10
  (+ optional AT14) section(s) — spatial QC/clustering, the cell-state axis in space,
  cell2location mapping of the Level 1 reference onto the tissue, NMF niche analysis,
  spatial neighborhood/proximity analysis, paper reveal, figure reproduction.
- **Level 3** (Xenium) and **cell2fate** (Fig. 6): out of scope for now — Xenium data and
  spliced/unspliced counts aren't available yet.

## Environment

Uses the shared `single_cell` conda environment (not `uv`, unlike the template):
```bash
module load conda   # if needed
conda activate single_cell
jupyter lab
```
`single_cell` lives at `/work/PRTNR/CHUV/DIR/rgottar1/single_cell_all/users/alederer/anaconda3/envs/single_cell`
— it's an instructor-owned environment shared across several course groups (the official
planning name is `neuromics-sc`). It already ships scanpy, anndata, cell2location,
scvi-tools, torch; squidpy, celltypist, harmonypy, and decoupler were added for this
project specifically (see `docs/build_notes.md` for exact versions/dates).

## Data

All source data already lives in `../../data/` (snRNA-seq, Visium, the paper PDFs) — see
`docs/build_notes.md` for exact donor/section selection rationale. Prepared (answer-key
stripped) student-facing datasets are written to:
- `../../data/snRNA_seq/level1_prepared/`
- `../../data/visium/level2_prepared/`

Instructor-only answer keys (stripped ground-truth columns/features, used for solution
notebooks and grading) live in `../answer_keys/` — **private**, do not copy into anything
students can read.

## Repo layout

```
docs/project_background.md   # student-facing intro (paper withheld until Level 2 reveal)
docs/build_notes.md          # instructor-facing: exact params/decisions, for future maintainers
notebooks/level1/            # 01_snrna_analysis_{student,solution}.ipynb
notebooks/level2/            # 02_spatial_cell2location_{student,solution}.ipynb
src/gbmspace_utils/          # shared helpers, actually imported by the notebooks
scripts/                     # data prep (01, 02) and cell2location CPU/GPU benchmark (03)
```
