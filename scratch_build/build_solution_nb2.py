"""Construct the Level 2 SOLUTION notebook (02_spatial_cell2location_solution.ipynb) via
nbformat. Sections 1-5 use REAL numbers from scratch_build/explore_level2_visium.py.
Sections 6+ depend on Level 1's annotated reference (data/processed/gbm_l1_snrna_AT10_AT14_annotated.h5ad)
-- code is drafted and will be validated/filled with real numbers once that exists, then
the whole notebook gets executed once via nbconvert as the final verification pass.
"""
from pathlib import Path
import nbformat as nbf

OUT = Path("/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/lederer/gbm_space_proj/notebooks/level2/02_spatial_cell2location_solution.ipynb")

cells = []
def md(s): cells.append(("md", s))
def code(s): cells.append(("code", s))

# ============================================================ TITLE
md(r"""# Level 2 — Spatial Context: Mapping the Glioma Microenvironment

## CAJAL "Neuromics 2026" — Computational Mini-Project C10 (Level 2)

**Estimated time:** ~2.5 days

**Learning objectives**
- Load, QC, and explore Visium spatial transcriptomics data
- Get a first ("naive") spatial domain map directly from spot expression, before any deconvolution
- See the malignant cell-state axis directly in space, even before deconvolving spots
- Map your Level 1 reference onto tissue with **cell2location**
- Identify spatial niches (tissue domains) from the deconvolved cell-state map
- Quantify spatial organization: neighborhood enrichment, a proximity network, and spatial intermixing
- Compare your own results to the published figures (the paper is revealed in this notebook)

**Dataset:** Visium spatial transcriptomics from the **same two donors** as Level 1 — `AT10`
(primary, full feature set) and `AT14` (optional secondary section). Each spot covers
multiple cells (~1-10), so spot expression is a *mixture*, unlike Level 1's single nuclei.

> The paper is still not named yet. You'll recognize the cell-type and cell-state language
> from Level 1 — that continuity is the point. The reveal happens partway through this notebook.

---""")

# ============================================================ 0. SETUP
md(r"""## 0. Setup""")
code(r"""import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, "/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/lederer/gbm_space_proj/src")
from gbmspace_utils.analysis import (
    MALIGNANT_AXIS_MARKERS, MAJOR_CLASS_OF, ZONATION_PANEL, score_axis,
    assign_dominant_state, spatial_proximity_network,
)
from gbmspace_utils.plotting import plot_gene_on_tissue, plot_spatial_categories

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(5, 4))
%matplotlib inline

print("scanpy", sc.__version__, "| squidpy", sq.__version__)""")

# ============================================================ 1. LOAD
md(r"""## 1. Load and explore the spatial data

🔬 **TASK 1.1:** Load the AT10 Visium section and inspect the object — note how it differs from Level 1's AnnData (spatial coordinates, a tissue image, no per-nucleus QC metrics yet).""")
code(r"""VISIUM = "/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/data/visium/level2_prepared/AT10-BRA-5-FO-1_2_student.h5ad"
adata = sc.read_h5ad(VISIUM)
print(adata)
print(f"\n{adata.n_obs} spots x {adata.n_vars} genes")
print(f".obsm: {list(adata.obsm.keys())}  |  .uns: {list(adata.uns.keys())}")
lib_id = list(adata.uns["spatial"].keys())[0]
print(f"Library: {lib_id}, images: {list(adata.uns['spatial'][lib_id]['images'].keys())}")""")

code(r"""fig, ax = plt.subplots(figsize=(6, 6))
sq.pl.spatial_scatter(adata, color=None, ax=ax, size=1.3)
ax.set_title(f"{lib_id} — H&E + spot grid ({adata.n_obs} spots)")""")

md(r"""❓ **QUESTION:** Each Visium spot is ~55 µm in diameter. Given typical nucleus/cell sizes, roughly how many cells might a single spot cover? What does that imply about interpreting any single spot's gene expression?""")

# ============================================================ 2. SPATIAL QC
md(r"""## 2. Spatial quality control

Same idea as Level 1 — total counts, genes detected, %mito — but spots, not nuclei, and no doublet score (a "doublet" concept doesn't apply the same way to multi-cell spots).

🔬 **TASK 2.1:** Compute QC metrics and look at their distributions.""")
code(r"""adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
print(adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe().round(1))""")

code(r"""fig, axes = plt.subplots(1, 3, figsize=(13, 4))
axes[0].hist(adata.obs["total_counts"], bins=60); axes[0].set_title("Total counts")
axes[1].hist(adata.obs["n_genes_by_counts"], bins=60); axes[1].set_title("Genes detected")
axes[2].hist(adata.obs["pct_counts_mt"], bins=60); axes[2].set_title("% mitochondrial")
plt.tight_layout(); plt.show()""")

md(r"""💡 **HINT:** Visium spots are much "deeper" than single nuclei (each spot pools several cells), so don't reuse Level 1's per-nucleus thresholds verbatim — look at *these* distributions. A light touch is usually right: drop near-empty spots (very low counts/genes), keep almost everything else.

🔬 **TASK 2.2:** Apply QC and a minimum-cells gene filter. Report spots/genes remaining.""")
code(r"""n0 = adata.n_obs
adata = adata[(adata.obs["total_counts"] >= 500) & (adata.obs["n_genes_by_counts"] >= 250)].copy()
sc.pp.filter_genes(adata, min_cells=3)
print(f"Spots: {n0} -> {adata.n_obs}")
print(f"Genes (min_cells=3): {adata.n_vars}")""")

md(r"""⚠️ **CHECKPOINT:** This section should remove very few spots — on the order of **1-2% of spots** (this dataset is high quality; expect roughly **3,900-3,950 spots** remaining out of ~4,000, and somewhere around **24,000-25,000 genes** after the gene filter). If you lost a large fraction of spots, your thresholds are too strict for Visium-scale counts.""")

# ============================================================ 3. NORMALIZE + CLUSTER (NAIVE)
md(r"""## 3. Normalization and a *naive* spatial domain map

Before any deconvolution, let's see what plain clustering of spot expression gives us — a "naive" map of spatial domains, mixing whatever cell types happen to co-occur in each spot.

🔬 **TASK 3.1:** Normalize, log-transform, select HVGs, and run PCA.""")
code(r"""adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", layer="counts")
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=30)
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata.uns["pca"] = adata_hvg.uns["pca"]
sc.pl.pca_variance_ratio(adata, n_pcs=30, log=True)""")

md(r"""🔬 **TASK 3.2:** Build the neighbor graph, UMAP, and cluster at a couple of resolutions. Pick one.""")
code(r"""sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
for res in [0.5, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_r{res}", flavor="igraph", n_iterations=2)
    print(f"resolution {res}: {adata.obs[f'leiden_r{res}'].nunique()} naive domains")
adata.obs["leiden"] = adata.obs["leiden_r1.0"]""")

md(r"""💡 At resolution 0.5 you should see roughly **10-12** domains; at 1.0, roughly **16-20**. We carry resolution 1.0 forward as `leiden` — finer domains are easier to relate to the cell-state axis later, but feel free to use 0.5 instead.

🔬 **TASK 3.3:** Plot the naive domains both on the tissue and on UMAP, side by side.""")
code(r"""fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
sq.pl.spatial_scatter(adata, color="leiden", ax=axes[0], size=1.3, legend_fontsize=6)
axes[0].set_title("Naive spatial domains (tissue)")
sc.pl.umap(adata, color="leiden", ax=axes[1], show=False, title="Naive spatial domains (UMAP)")
plt.tight_layout(); plt.show()""")

md(r"""❓ **QUESTION:** Do the naive domains form spatially coherent regions (contiguous patches of one color), or are they speckled/scattered across the tissue? What would each pattern imply about whether expression-based clustering alone is picking up real tissue architecture?""")

# ============================================================ 4. AXIS IN SPACE
md(r"""## 4. The malignant cell-state axis, in space (before deconvolution)

Level 1's marker-gene scoring works the same way here — `score_genes` doesn't care whether
an observation is a nucleus or a multi-cell spot. The catch: a spot's score is a *blend* of
whatever cell states happen to be in that spot, not a clean call.

🔬 **TASK 4.1:** Score every spot against the Level 1 malignant-state marker sets, using the shared `score_axis()` helper.""")
code(r"""state_scores = score_axis(adata, MALIGNANT_AXIS_MARKERS, use_raw=True)
for col in state_scores.columns:
    adata.obs[f"score_{col}"] = state_scores[col].values
print(adata.obs.groupby("leiden")[[f"score_{c}" for c in state_scores.columns]].mean().round(3))""")

md(r"""🔬 **TASK 4.2:** Plot the paper's minimal 4-gene spatial zonation panel (`ZONATION_PANEL`: dev-like → gliosis → hypoxia) directly on tissue.""")
code(r"""present = [g for g in ZONATION_PANEL if g in adata.raw.var_names]
fig2, axes2 = plt.subplots(1, len(present), figsize=(4.2 * len(present), 4.2))
for ax, gene in zip(axes2, present):
    expr = np.asarray(adata[:, gene].X.todense()).flatten()
    sca = ax.scatter(adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1], c=expr, cmap="Reds", s=8,
                      vmax=np.percentile(expr[expr > 0], 95) if (expr > 0).any() else None)
    ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_title(gene); ax.axis("off")
    fig2.colorbar(sca, ax=ax, shrink=0.7)
plt.tight_layout(); plt.show()""")

md(r"""❓ **QUESTION:** Do you see any spatial gradient across the four zonation genes — e.g. do `AQP4`-high and `HILPDA`-high regions occupy *different, non-overlapping* areas of the tissue? Compare this to the per-cluster axis scores above. This blended, spot-level picture is exactly the limitation **cell2location** is designed to address — keep this figure in mind for comparison once you've deconvolved.

---

## 5. Mapping single cells onto space with cell2location

So far every spot has been treated as one observation, even though it's really a mixture.
**cell2location** uses your Level 1 reference (cell types learned from single nuclei) to
estimate *how many cells of each type* are in every spot — turning "this spot's expression
looks bulk hypoxic" into "this spot is ~60% Hypoxic-state malignant cells + ~20% Macrophage + ...".

🔬 **TASK 5.1:** Load your saved, annotated Level 1 reference, and build the label cell2location
will actually deconvolve.""")
code(r"""ref_path = "/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/lederer/gbm_space_proj/data/processed/gbm_l1_snrna_AT10_AT14_annotated.h5ad"
adata_ref = sc.read_h5ad(ref_path)
print(f"Reference: {adata_ref.n_obs} nuclei, cell types: {adata_ref.obs['cell_type'].nunique()}")
print(adata_ref.obs['cell_type'].value_counts())""")

md(r"""💡 **HINT — which label should cell2location deconvolve?** `cell_type` (Section 6's CellTypist
majority call per cluster) is a reasonable broad label for TME cells, but for malignant cells it's
a region-mimic label (e.g. "Hypothalamus glioblast" vs "Striatum glioblast") — the same malignant
population matched to whichever normal-brain-region profile CellTypist happened to land on, not a
real biological distinction (their average expression profiles are correlated >0.9). Per the
paper's own Methods, cell2location was run on genuine **malignant cell-state clusters**, not
region-mimic labels. Build a combined label: `malignant_state` (Section 8's 9 marker-defined axis
states) for malignant cells, `cell_type` for TME cells.

⚠️ Some non-malignant clusters carry a malignant-mimic CellTypist label too (their CNV signal
didn't clear the malignant threshold, so they're correctly TME — but CellTypist still gave them a
"glioblast"/"OPC"-sounding name). Don't just drop these: check their actual marker-gene expression
first. A quick oligodendrocyte marker score (`MBP, PLP1, MOG, MOBP, ST18`) reveals most of them
*are* real oligodendrocytes (mean score ~120-155 vs. ~1-6 in true TME/malignant cells) — CellTypist's
`Developing_Human_Brain` model has no clean adult-oligodendrocyte category, so it matched them to
the closest thing it had (developing OPC/neural-crest programmes). Relabel those properly instead
of throwing away a real, abundant TME population; only drop the few small clusters that show
*neither* a malignant-state signature *nor* oligodendrocyte markers (truly ambiguous, ~2% of cells).

🔬 **TASK 5.1b:** Build `cell_state_for_c2l` and use it as the reference label.""")
code(r"""MALIGNANT_MIMIC_KEYWORDS = ("glioblast", "radial glia", "opc", "neural crest", "neuroblast")
OLIGODENDROCYTE_MARKERS = ["MBP", "PLP1", "MOG", "MOBP", "ST18"]

is_malignant = (adata_ref.obs["cell_status_derived"] == "Malignant").to_numpy()
is_mimic_but_tme = (~is_malignant) & adata_ref.obs["cell_type"].str.lower().str.contains(
    "|".join(MALIGNANT_MIMIC_KEYWORDS)).to_numpy()

present_markers = [g for g in OLIGODENDROCYTE_MARKERS if g in adata_ref.var_names]
counts = adata_ref.layers["counts"]
libsize = np.asarray(counts.sum(axis=1)).flatten()
marker_sum = np.asarray(counts[:, [adata_ref.var_names.get_loc(g) for g in present_markers]].sum(axis=1)).flatten()
oligo_score = marker_sum / np.maximum(libsize, 1) * 1e4

is_real_oligo = is_mimic_but_tme & (oligo_score > 20)  # >>10x the ~1-6 baseline in true TME/malignant cells
is_truly_ambiguous = is_mimic_but_tme & ~is_real_oligo
print(f"Of {int(is_mimic_but_tme.sum())} mimic-labelled-but-TME nuclei: "
      f"{int(is_real_oligo.sum())} are real oligodendrocytes (relabelled), "
      f"{int(is_truly_ambiguous.sum())} stay ambiguous (dropped)")

adata_ref = adata_ref[~is_truly_ambiguous].copy()
is_malignant = (adata_ref.obs["cell_status_derived"] == "Malignant").to_numpy()
is_real_oligo = is_real_oligo[~is_truly_ambiguous]

cell_type_relabelled = adata_ref.obs["cell_type"].astype(str).to_numpy()
cell_type_relabelled[is_real_oligo] = "Oligodendrocyte"
adata_ref.obs["cell_state_for_c2l"] = np.where(
    is_malignant, adata_ref.obs["malignant_state"].astype(str), cell_type_relabelled)
print("\nReference labels for cell2location:")
print(adata_ref.obs["cell_state_for_c2l"].value_counts())""")

md(r"""💡 **HINT — runtime.** cell2location is two models: a reference **signature** model (NB regression on your Level 1 single-cell counts) and a **spatial mapping** model (maps that signature onto spots). Both are slow on CPU. We benchmarked this on the actual data: training on all shared genes costs **~170s/epoch**; the standard cell2location gene filter (`cell2location.utils.filtering.filter_genes`, defaults) cuts that to **~15,900 genes** and **~72s/epoch** (reference) / **~3.9s/epoch** (spatial mapping, this Visium section's spot count). Paper-faithful epoch counts (ref 400 / mapping 6000) would take **~8 hours total on CPU** — only realistic on GPU. Set the mode below accordingly.

🔬 **TASK 5.2:** Set the compute mode, filter genes, and train the reference signature model.

💡 **HINT:** cell2location's reference model assumes a negative-binomial (GammaPoisson) likelihood over **raw integer counts** — but by this point your Level 1 reference's `.X` is log-normalized (from Level 1 Section 3). Point `setup_anndata` at the raw-counts layer explicitly (`layer="counts"`), or training will crash with a cryptic "value... not within the support of GammaPoisson" error the first time it tries to evaluate a likelihood on fractional log-values.""")
code(r"""import os
C2L_MODE = "FULL"   # "DEMO" (tiny, seconds-minutes; for fast iteration) / "FAST" (CPU, ~25 min) / "FULL" (GPU, paper-exact)
REF_EPOCHS = {"DEMO": 5, "FAST": 20, "FULL": 400}[C2L_MODE]
MAP_EPOCHS = {"DEMO": 20, "FAST": 300, "FULL": 6000}[C2L_MODE]
print(f"Mode={C2L_MODE}: reference {REF_EPOCHS} epochs, mapping {MAP_EPOCHS} epochs")

# Checkpoint dir for the expensive trained-model stages below (reference signature: ~1h;
# each spatial mapping: ~2-4h at FULL). If a later, unrelated cell fails (e.g. a typo three
# cells later), re-running this notebook from scratch would otherwise repeat hours of already-
# correct training. Checkpoint files are keyed by C2L_MODE so a DEMO/FAST run's cache is never
# mistaken for a FULL run's.
CKPT_DIR = "/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/lederer/gbm_space_proj/scratch_build/checkpoints"
os.makedirs(CKPT_DIR, exist_ok=True)""")

code(r"""from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel, Cell2location

shared = sorted(set(adata_ref.var_names) & set(adata.var_names))
ref = adata_ref[:, shared].copy()
vis = adata.copy()[:, shared].copy()

selected = filter_genes(ref, cell_count_cutoff=15, cell_percentage_cutoff2=0.05, nonz_mean_cutoff=1.12)
ref = ref[:, selected].copy()
vis = vis[:, [g for g in selected if g in vis.var_names]].copy()
print(f"Genes after filtering: {ref.n_vars}")

REF_CKPT = f"{CKPT_DIR}/inf_aver_{C2L_MODE}.parquet"
if os.path.exists(REF_CKPT):
    inf_aver = pd.read_parquet(REF_CKPT)
    print(f"Loaded cached reference signature from {REF_CKPT} (skipped {REF_EPOCHS}-epoch training)")
else:
    RegressionModel.setup_anndata(ref, layer="counts", batch_key="donor_id", labels_key="cell_state_for_c2l")
    ref_model = RegressionModel(ref)
    ref_model.train(max_epochs=REF_EPOCHS, batch_size=10000)
    ref = ref_model.export_posterior(ref, sample_kwargs={"num_samples": 100, "batch_size": 10000})
    inf_aver = ref.varm["q05_per_cluster_mu_fg"]
    inf_aver.to_parquet(REF_CKPT)
    print(f"Saved reference signature checkpoint -> {REF_CKPT}")
print(f"Reference signature: {inf_aver.shape} (genes x cell types)")""")

md(r"""🔬 **TASK 5.3:** Train the spatial mapping model and export cell-type abundance per spot.""")
code(r"""vis = vis[:, [g for g in inf_aver.index if g in vis.var_names]].copy()
inf_aver_aligned = inf_aver.loc[vis.var_names]

MAP_CKPT = f"{CKPT_DIR}/vis_AT10_{C2L_MODE}.h5ad"
if os.path.exists(MAP_CKPT):
    vis = sc.read_h5ad(MAP_CKPT)
    print(f"Loaded cached AT10 mapping result from {MAP_CKPT} (skipped {MAP_EPOCHS}-epoch training)")
else:
    Cell2location.setup_anndata(vis, layer="counts", batch_key="sample_name" if "sample_name" in vis.obs else None)
    sp_model = Cell2location(vis, cell_state_df=inf_aver_aligned, N_cells_per_location=30, detection_alpha=200)
    # Paper used batch_size ~= 25% of spots per tumour (Methods), not full-batch.
    sp_model.train(max_epochs=MAP_EPOCHS, batch_size=max(int(0.25 * vis.n_obs), 1))
    vis = sp_model.export_posterior(vis, sample_kwargs={"num_samples": 100, "batch_size": vis.n_obs})
    vis.write_h5ad(MAP_CKPT)
    print(f"Saved AT10 mapping checkpoint -> {MAP_CKPT}")

abundance = vis.obsm["q05_cell_abundance_w_sf"] if "q05_cell_abundance_w_sf" in vis.obsm else \
            vis.obs[[c for c in vis.obs.columns if c.startswith("q05")]]
print(f"Cell-type abundance per spot: {abundance.shape}")
print(abundance.describe().T[["mean", "std", "max"]].round(2))""")

md(r"""🔬 **TASK 5.4:** Plot a few cell-type abundance maps on tissue.""")
code(r"""top_types = abundance.mean().nlargest(4).index.tolist()
fig, axes = plt.subplots(1, len(top_types), figsize=(4.2 * len(top_types), 4.2))
for ax, ct in zip(axes, top_types):
    vals = abundance[ct].to_numpy()
    coords = vis.obsm["spatial"]
    sca = ax.scatter(coords[:, 0], coords[:, 1], c=vals, cmap="viridis", s=8)
    ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_title(ct); ax.axis("off")
    fig.colorbar(sca, ax=ax, shrink=0.7)
plt.tight_layout(); plt.show()""")

md(r"""❓ **QUESTION:** Compare these deconvolved abundance maps to the blended Section 4 scores. Are the malignant-state spatial patterns sharper now? Does any region's *dominant* cell type surprise you given what the naive expression clustering (Section 3) suggested was there?

---""")

# ============================================================ 6. NICHE ANALYSIS
md(r"""## 6. Spatial niches and intermixing

Individual cell-type abundances are noisy spot-by-spot. **NMF** on the abundance matrix finds
recurring *co-occurrence patterns* — niches — the way the original study does.

🔬 **TASK 6.1:** Run NMF on the cell2location abundance matrix with a few different factor counts.""")
code(r"""from sklearn.decomposition import NMF

for n_factors in [5, 8, 12]:
    nmf = NMF(n_components=n_factors, init="nndsvda", random_state=0, max_iter=500)
    W = nmf.fit_transform(abundance.clip(lower=0))
    print(f"n_factors={n_factors}: reconstruction error = {nmf.reconstruction_err_:.1f}")
# Carry forward one choice for the rest of the notebook:
N_NICHES = 8
nmf = NMF(n_components=N_NICHES, init="nndsvda", random_state=0, max_iter=500)
niche_loadings = nmf.fit_transform(abundance.clip(lower=0))
vis.obs["niche"] = pd.Categorical(niche_loadings.argmax(axis=1).astype(str))
print(vis.obs["niche"].value_counts())""")

md(r"""🔬 **TASK 6.2:** Plot niches on tissue, and look at which cell types load most strongly onto each niche factor.""")
code(r"""fig, ax = plt.subplots(figsize=(6, 6))
plot_spatial_categories(vis, "niche", spatial_key="spatial", ax=ax)
ax.set_title(f"{N_NICHES} NMF-derived niches"); plt.show()

components = pd.DataFrame(nmf.components_, columns=abundance.columns,
                           index=[f"niche_{i}" for i in range(N_NICHES)])
fig, ax = plt.subplots(figsize=(10, 6))
sns.heatmap(components.div(components.max(axis=1), axis=0), cmap="viridis", ax=ax)
ax.set_title("Cell-type loading per niche (row-normalized)")
plt.tight_layout(); plt.show()""")

md(r"""🔬 **TASK 6.3 — spatial intermixing.** For each spot, compute the Shannon entropy of its cell-type abundance distribution — a high-entropy spot has many cell types evenly mixed; a low-entropy spot is dominated by one.""")
code(r"""from scipy.stats import entropy

props = abundance.clip(lower=0)
props = props.div(props.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
vis.obs["intermixing_entropy"] = props.apply(lambda row: entropy(row + 1e-12), axis=1)

fig, ax = plt.subplots(figsize=(6, 6))
coords = vis.obsm["spatial"]
sca = ax.scatter(coords[:, 0], coords[:, 1], c=vis.obs["intermixing_entropy"], cmap="magma", s=8)
ax.invert_yaxis(); ax.set_aspect("equal"); ax.set_title("Spatial intermixing (entropy)"); ax.axis("off")
fig.colorbar(sca, ax=ax, shrink=0.7); plt.show()
print(vis.obs.groupby("niche")["intermixing_entropy"].mean().sort_values())""")

md(r"""❓ **QUESTION:** Which niche has the lowest intermixing entropy (most "pure")? Which has the highest? Does low entropy correspond to a niche you'd expect to be compositionally homogeneous (e.g. a dense malignant core) versus a mixed immune/stromal region?

---""")

# ============================================================ 7. NEIGHBORHOOD / PROXIMITY
md(r"""## 7. Spatial neighborhood analysis, two ways

🔬 **TASK 7.1 — squidpy.** Build the spatial neighbor graph and compute neighborhood enrichment between niches.""")
code(r"""sq.gr.spatial_neighbors(vis, coord_type="generic", n_neighs=6)
sq.gr.nhood_enrichment(vis, cluster_key="niche")
sq.pl.nhood_enrichment(vis, cluster_key="niche", figsize=(6, 5), annotate=False)""")

md(r"""🔬 **TASK 7.2 — the paper's own method.** Build the same kind of "which niches are near which" picture a different way: pairwise minimum spot-distance via a k-d tree, summarized at the 25th percentile (`gbmspace_utils.spatial_proximity_network` — this mirrors the paper's actual Fig. 2E/3C/6E/7E method, an alternative to squidpy's enrichment z-scores).""")
code(r"""prox = spatial_proximity_network(vis, cluster_key="niche", spatial_key="spatial", percentile=25)
fig, ax = plt.subplots(figsize=(7, 6))
sns.heatmap(prox, cmap="viridis_r", ax=ax)  # reversed: closer (smaller distance) = brighter
ax.set_title("Niche-niche proximity (25th-pctile nearest distance, smaller=closer)")
plt.tight_layout(); plt.show()""")

md(r"""❓ **QUESTION:** Do squidpy's neighborhood-enrichment z-scores and the proximity-distance heatmap agree on which niches are spatial neighbors? Where (if anywhere) do they disagree, and why might two reasonable methods for "spatial closeness" give different answers?

---""")

# ============================================================ 8. CCC
md(r"""## 8. Cell-cell communication with LIANA

The published study runs **LIANA** (a consensus of several ligand-receptor methods) comparing
TME states co-localized with dev-like niches vs. those co-localized with gliosis/hypoxia
niches. We skip their cross-donor Tensor-cell2cell step (not meaningful with only 2 donors) but
reproduce the core comparison directly on our own spots.

🔬 **TASK 8.1:** Define two spot groups directly from the per-spot malignant-axis scores
already computed in Section 4 (`state_scores`, on the full normalized `adata` — deliberately
*not* the heavily gene-filtered `vis` from cell2location, so LIANA sees the full gene panel
rather than just the ~15,900 genes that survived cell2location's deconvolution-specific
filter) — "dev-like-dominant" (OPC-like / OPC-NPC-like / OPC-neuronal-like / NPC-neuronal-like
/ AC-progenitor-like) vs. "gliosis/hypoxia-dominant" (AC-gliosis-like / Gliosis-like /
Hypoxic) — leaving Proliferative-dominant spots out of the comparison (neither end of the
trajectory).""")
code(r"""DEV_LIKE_STATES = ["OPC-like", "OPC-NPC-like", "OPC-neuronal-like", "NPC-neuronal-like", "AC-progenitor-like"]
GLIOSIS_HYPOXIA_STATES = ["AC-gliosis-like", "Gliosis-like", "Hypoxic"]

dominant_axis = assign_dominant_state(state_scores)  # from Section 4, computed on `adata`

adata.obs["niche_group"] = np.select(
    [dominant_axis.isin(DEV_LIKE_STATES).to_numpy(), dominant_axis.isin(GLIOSIS_HYPOXIA_STATES).to_numpy()],
    ["dev-like-dominant", "gliosis-hypoxia-dominant"],
    default="other",
)
print(adata.obs["niche_group"].value_counts())""")

md(r"""🔬 **TASK 8.2:** Run LIANA's consensus ligand-receptor scoring with `niche_group` as the
grouping variable — this scores every source-group -> target-group pair, so with two groups you
get all four directions, including the two cross-group ("between niche") directions we actually
care about.""")
code(r"""import liana as li

vis_ccc = adata[adata.obs["niche_group"] != "other"].copy()
vis_ccc.obs["niche_group"] = vis_ccc.obs["niche_group"].astype("category")
print(f"Spots used for LIANA: {vis_ccc.n_obs} ({vis_ccc.obs['niche_group'].value_counts().to_dict()})")

li.mt.rank_aggregate(vis_ccc, groupby="niche_group", resource_name="consensus", expr_prop=0.1, use_raw=True, verbose=True)
liana_res = vis_ccc.uns["liana_res"]
print(f"\nLIANA interactions scored: {len(liana_res)}")
print(liana_res["source"].astype(str).str.cat(liana_res["target"].astype(str), sep=" -> ").value_counts())""")

md(r"""🔬 **TASK 8.3:** Look at the top cross-group interactions in each direction (dev-like spots
signalling to gliosis/hypoxia spots, and vice versa), ranked by LIANA's aggregate
`magnitude_rank` (lower = stronger consensus across methods).""")
code(r"""for src, tgt in [("dev-like-dominant", "gliosis-hypoxia-dominant"), ("gliosis-hypoxia-dominant", "dev-like-dominant")]:
    cross = liana_res[(liana_res["source"] == src) & (liana_res["target"] == tgt)].sort_values("magnitude_rank")
    print(f"\nTop 10 {src} -> {tgt}:")
    print(cross[["ligand_complex", "receptor_complex", "magnitude_rank", "specificity_rank"]].head(10).to_string(index=False))

li.pl.dotplot(adata=vis_ccc, colour="magnitude_rank", size="specificity_rank",
              inverse_size=True, inverse_colour=True,
              source_labels=["dev-like-dominant", "gliosis-hypoxia-dominant"],
              target_labels=["dev-like-dominant", "gliosis-hypoxia-dominant"],
              top_n=15, orderby="magnitude_rank", orderby_ascending=True, figure_size=(9, 7))
plt.show()""")

md(r"""❓ **QUESTION:** Which ligand-receptor pairs are specific to one direction (e.g. dev-like
spots signalling to gliosis/hypoxia spots) rather than the reverse? Do any involve genes you
already recognise from the hypoxia/gliosis marker sets in `MALIGNANT_AXIS_MARKERS` (e.g.
`VEGFA`)? What would that suggest about how these two ends of the trajectory interact spatially,
rather than just co-occurring?

---""")

# ============================================================ 9. PAPER REVEAL
md(r"""## 9. Revealing the paper

📄 **de Jong, Memi, Gracia, Lazareva et al. "A spatiotemporal cancer cell trajectory
underlies glioblastoma heterogeneity." bioRxiv 2025.05.13.653495.** Companion website:
[gbmspace.org](https://www.gbmspace.org/). The data you have been working with (AT10 and
AT14, snRNA-seq + Visium) are two of the 12 IDH-wildtype glioblastoma tumours profiled in
this study (1,025,329 nuclei total across the cohort).

**Key findings, for comparison against your own results:**
- Malignant cells occupy a **continuous trajectory**, not discrete subtypes: from
  developmental-like states (OPC-like, NPC-neuronal-like, AC-progenitor-like) through a
  **gliosis** and **hypoxia** axis — what was historically called "mesenchymal-like" (MES1/2)
  in the Neftel et al. 2019 framework, but the authors show classical EMT regulators
  (`SNAI1/2`, `TWIST1/2`, `ZEB1/2`) are *not* specifically enriched there, arguing against an
  EMT interpretation.
- This trajectory maps onto **spatial zonation**: AC-progenitor-like cells dominate near the
  tumour core / infiltrating edge; gliosis and hypoxic states concentrate deep in the tumour,
  around and within necrotic regions — exactly the `AQP4` → `ABCC3` → `AKAP12` → `HILPDA`
  gradient you looked for in Section 4.
- Spatial **niches** were derived the same way you just did it: NMF on cell2location
  cell-state abundances (16 factors per tumour in the paper, clustered into ~14-16 recurrent
  niches across the cohort), cross-validated against pathologist-annotated IvyGAP regions
  (leading edge, infiltrating tumour, cellular tumour, necrosis, pseudopalisading cells,
  perinecrotic zone, microvascular proliferation).
- A major caveat the authors flag explicitly: **single-biopsy sampling can be misleading** —
  one of their tumours (AT10, the same donor in your data!) showed a different dominant
  malignant state in each of its 4 sampled sites, challenging the idea of one fixed
  "subtype" per tumour.

🔬 **TASK 9.1:** Now that you know the source, compare your own niche map and axis scores
against the paper's actual cell2location/niche outputs for this exact section (these were
withheld from your input data — load them now from the answer-key file for comparison only).""")
code(r"""answer_key = sc.read_h5ad("/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/lederer/answer_keys/AT10-BRA-5-FO-1_2_answer_key.h5ad")
print(answer_key.var["feature_types"].value_counts())

def extract_answer_block(ak, feature_type):
    mask = (ak.var["feature_types"] == feature_type).to_numpy()
    block = ak.X[:, mask]
    block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
    return pd.DataFrame(block, index=ak.obs_names, columns=ak.var_names[mask])

niche_answer = extract_answer_block(answer_key, "Spatial niche abundances")
state_answer = extract_answer_block(answer_key, "Cell state abundances")
print(f"\nAnswer-key niche columns ({niche_answer.shape[1]}): {list(niche_answer.columns)}")
print(f"Answer-key cell-state columns ({state_answer.shape[1]}): {list(state_answer.columns)}")

common_spots = vis.obs_names.intersection(niche_answer.index)
print(f"\nSpots in common with answer key: {len(common_spots)} / {vis.n_obs}")

niche_onehot = pd.get_dummies(vis.obs.loc[common_spots, "niche"])
joint = np.hstack([niche_onehot.to_numpy(), niche_answer.loc[common_spots].to_numpy()])
full_corr = np.corrcoef(joint, rowvar=False)
corr = pd.DataFrame(full_corr[:niche_onehot.shape[1], niche_onehot.shape[1]:],
                     index=niche_onehot.columns, columns=niche_answer.columns)
print("\nCorrelation: your NMF niche (rows) vs the paper's niche abundance (cols):")
print(corr.round(2))
print("\nBest-matching paper niche per your niche (by |correlation|):")
print(corr.abs().idxmax(axis=1))

fig, ax = plt.subplots(figsize=(8, 6))
sns.heatmap(corr, cmap="RdBu_r", center=0, ax=ax)
ax.set_title("Your niches vs the paper's niche abundances (correlation)")
plt.tight_layout(); plt.show()""")

md(r"""❓ **QUESTION:** Where does your independent analysis agree with the published result, and where does it diverge? For any divergence, what's your best hypothesis — different gene filtering, different epoch budget (FAST vs FULL), different niche factor count, or something about how the reference was built in Level 1?

---

## 10. Secondary check: does the same pattern hold in tumour 2? (AT14)

Everything above used AT10 only. The paper's own caveat from Section 9 — that a single biopsy
can be misleading, and that AT10 itself showed different dominant states across its 4 sampled
sites — is exactly why a second, independent tumour is worth checking. `AT14`'s Visium section
has **no IvyGAP histopathology overlay**, so this is a cross-tumour sanity check against AT10,
not a check against ground truth. It's also intentionally **condensed**: the reference
signature from Level 1 is reused as-is (only a new spatial-mapping model is trained), and we
don't repeat the squidpy-vs-k-d-tree neighborhood comparison or LIANA here — those stay AT10-only.

🔬 **TASK 10.1:** Load, QC, normalize, and score the malignant-state axis on AT14's spots, the
same way as Sections 2-4.""")
code(r"""VISIUM_AT14 = "/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/data/visium/level2_prepared/AT14-BRA-4-FO-2_1_student.h5ad"
adata14 = sc.read_h5ad(VISIUM_AT14)
adata14.var["mt"] = adata14.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata14, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

n0 = adata14.n_obs
adata14 = adata14[(adata14.obs["total_counts"] >= 500) & (adata14.obs["n_genes_by_counts"] >= 250)].copy()
sc.pp.filter_genes(adata14, min_cells=3)
print(f"AT14 spots: {n0} -> {adata14.n_obs}  |  genes: {adata14.n_vars}")

adata14.layers["counts"] = adata14.X.copy()
sc.pp.normalize_total(adata14, target_sum=1e4)
sc.pp.log1p(adata14)
adata14.raw = adata14

state_scores14 = score_axis(adata14, MALIGNANT_AXIS_MARKERS, use_raw=True)
dominant14 = assign_dominant_state(state_scores14)
adata14.obs["malignant_class"] = dominant14.map(MAJOR_CLASS_OF)
print("\nAT14 dominant malignant-class fractions:")
print(adata14.obs["malignant_class"].value_counts(normalize=True).round(3))""")

md(r"""🔬 **TASK 10.2:** Map cell2location onto AT14 — reusing the **same** reference signature
(`inf_aver`) fit once in Section 5, training only a new spatial-mapping model — then derive
niches the same way as Section 6.""")
code(r"""shared14 = sorted(set(inf_aver.index) & set(adata14.var_names))
vis14 = adata14.copy()[:, shared14].copy()
inf_aver_aligned14 = inf_aver.loc[shared14]
print(f"AT14 genes shared with the Level 1 reference signature: {len(shared14)}")

MAP_CKPT14 = f"{CKPT_DIR}/vis_AT14_{C2L_MODE}.h5ad"
if os.path.exists(MAP_CKPT14):
    vis14 = sc.read_h5ad(MAP_CKPT14)
    print(f"Loaded cached AT14 mapping result from {MAP_CKPT14} (skipped {MAP_EPOCHS}-epoch training)")
else:
    Cell2location.setup_anndata(vis14, layer="counts", batch_key="sample_name" if "sample_name" in vis14.obs else None)
    sp_model14 = Cell2location(vis14, cell_state_df=inf_aver_aligned14, N_cells_per_location=30, detection_alpha=200)
    _bs14 = max(int(0.25 * vis14.n_obs), 2)
    if vis14.n_obs % _bs14 == 1:  # last mini-batch of 1 sample causes 0-d tensor error in PyTorch
        _bs14 += 1
    sp_model14.train(max_epochs=MAP_EPOCHS, batch_size=_bs14)
    vis14 = sp_model14.export_posterior(vis14, sample_kwargs={"num_samples": 100, "batch_size": vis14.n_obs})
    vis14.write_h5ad(MAP_CKPT14)
    print(f"Saved AT14 mapping checkpoint -> {MAP_CKPT14}")

abundance14 = vis14.obsm["q05_cell_abundance_w_sf"] if "q05_cell_abundance_w_sf" in vis14.obsm else \
              vis14.obs[[c for c in vis14.obs.columns if c.startswith("q05")]]
print(f"AT14 cell-type abundance per spot: {abundance14.shape}")

nmf14 = NMF(n_components=N_NICHES, init="nndsvda", random_state=0, max_iter=500)
niche_loadings14 = nmf14.fit_transform(abundance14.clip(lower=0))
vis14.obs["niche"] = pd.Categorical(niche_loadings14.argmax(axis=1).astype(str))

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
plot_spatial_categories(vis14, "niche", spatial_key="spatial", ax=axes[0])
axes[0].set_title(f"AT14: {N_NICHES} NMF-derived niches")
plt.tight_layout(); plt.show()
print(vis14.obs["niche"].value_counts())""")

md(r"""🔬 **TASK 10.3:** Compare AT10 vs AT14 directly — does the dev-like -> gliosis -> hypoxia
axis, and a comparable niche structure, show up in both tumours?""")
code(r"""adata.obs["malignant_class"] = dominant_axis.map(MAJOR_CLASS_OF)  # dominant_axis from Section 8

cross_tumor = pd.DataFrame({
    "AT10": adata.obs["malignant_class"].value_counts(normalize=True),
    "AT14": adata14.obs["malignant_class"].value_counts(normalize=True),
}).fillna(0).round(3)
print("Dominant malignant-class fraction per spot, AT10 vs AT14:")
print(cross_tumor)

fig, ax = plt.subplots(figsize=(7, 4))
cross_tumor.plot(kind="bar", ax=ax)
ax.set(ylabel="fraction of spots", title="Malignant-class composition: AT10 vs AT14")
plt.tight_layout(); plt.show()

print(f"\nAT10 niches: {N_NICHES} (from {vis.n_obs} spots) | AT14 niches: {N_NICHES} (from {vis14.n_obs} spots)")
print("Note: AT14 has no IvyGAP/cell2location answer key, so this is a cross-tumour check, "
      "not a check against ground truth (TASK 9.1's answer-key comparison stays AT10-only).")""")

md(r"""❓ **QUESTION:** Is the malignant-class composition similar between AT10 and AT14, or does
one tumour skew much more dev-like / gliosis-hypoxia than the other? Given the paper's own
caveat about single-biopsy sampling (Section 9), how confident should you be that either
section alone represents "the" malignant-state distribution of its tumour?

---

## 11. Write-up

🔬 **TASK 11.1:** Reproduce one specific published figure panel using your own pipeline, and
write a short paragraph comparing your result to the original. We reproduce the spatial
malignant-class map (the panel underlying the paper's Fig. 2/3 zonation claim) directly from
our own AT10 pipeline.""")
code(r"""fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
classes = sorted(adata.obs["malignant_class"].dropna().unique())
palette = dict(zip(classes, sns.color_palette("Set2", len(classes))))
coords = adata.obsm["spatial"]
colors = adata.obs["malignant_class"].map(palette)
axes[0].scatter(coords[:, 0], coords[:, 1], c=colors, s=8)
axes[0].invert_yaxis(); axes[0].set_aspect("equal"); axes[0].axis("off")
axes[0].set_title("Our spatial malignant-class map (AT10)")
handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[c], label=c, markersize=8)
           for c in classes]
axes[0].legend(handles=handles, fontsize=7, loc="lower left")

plot_spatial_categories(vis, "niche", spatial_key="spatial", ax=axes[1])
axes[1].set_title(f"Our {N_NICHES}-niche map (AT10), for comparison")
plt.tight_layout()
fig.savefig("/shared/projects/tp_2630_ubordeaux_neuromics_184418/projects/C10/lederer/gbm_space_proj/notebooks/level2/level2_summary_figure.png",
            dpi=200, bbox_inches="tight")
plt.show()""")

md(r"""🔬 **TASK 11.1 — Write-up.** In 2-3 sentences, compare your spatial malignant-class map and the niches you recovered to the paper's zonation claim. Back it with your own numbers: the TASK 9.1 niche-correlation values (your NMF niches vs the answer key) and the TASK 10.3 AT10-vs-AT14 dominant-class fractions.

<!-- INSTRUCTOR-ONLY -->
**[INSTRUCTOR: write-up paragraph filled in by hand in the executed notebook after the full-scale
run, using the real correlation values and class fractions printed in TASK 9.1/10.3 above — not
invented before the notebook was actually run. This whole INSTRUCTOR-ONLY block is stripped from
the derived student notebook by scratch_build/derive_student_nb.py.]**
<!-- /INSTRUCTOR-ONLY -->

---

## Summary

You have:
1. ✅ QC'd and explored real Visium spatial data
2. ✅ Built a naive (pre-deconvolution) spatial domain map
3. ✅ Seen the malignant cell-state axis directly in space
4. ✅ Mapped your Level 1 reference onto tissue with **cell2location**
5. ✅ Identified spatial niches via NMF, and quantified spatial intermixing
6. ✅ Compared two methods for spatial neighborhood/proximity analysis
7. ✅ Run LIANA cell-cell communication scoring between dev-like and gliosis/hypoxia spots
8. ✅ Compared your independent results against the published findings
9. ✅ Checked whether the same axis/niche pattern holds in a second tumour (AT14)

**Further reading, not built here:** the paper also describes **spaceTree** (joint cell-type
+ genetic-clone mapping) and **cell2fate** (RNA-velocity-based temporal ordering of malignant
states) — both require data/pipelines beyond this course's current scope (paired snATAC-seq
clone-calling, and spliced/unspliced counts, respectively).""")

# ============================================================ BUILD
nb = nbf.v4.new_notebook()
nb.cells = [nbf.v4.new_markdown_cell(s) if k == "md" else nbf.v4.new_code_cell(s) for k, s in cells]
nb.metadata["kernelspec"] = {"display_name": "Python (single_cell)", "language": "python", "name": "single_cell"}
nb.metadata["language_info"] = {"name": "python"}
OUT.parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, OUT)
print(f"Wrote {OUT} with {len(nb.cells)} cells "
      f"({sum(1 for k,_ in cells if k=='md')} md, {sum(1 for k,_ in cells if k=='code')} code)")
