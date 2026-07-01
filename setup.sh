#!/bin/bash
# =============================================================================
# C10 — student setup (run AFTER you have cloned this repo on the IFB cluster)
#
#   git clone <REPO_URL> ~/gbm-space-c10
#   cd ~/gbm-space-c10
#   bash setup.sh
#
# What it does (it does NOT copy anything — you already have the clone):
#   1. enables conda for all future shells (conda init)
#   2. makes the shared `single_cell` environment discoverable and verifies it
#   3. registers the `single_cell` Jupyter kernel for VS Code / Jupyter
#
# NOTE: this only sets things up. Run heavy compute via Slurm (see INSTALL.md §5),
#       never on the login node.
# =============================================================================
set -uo pipefail

ENVROOT=/shared/projects/tp_2630_ubordeaux_neuromics_184418/envs
CONDA_SH=/shared/software/miniconda/etc/profile.d/conda.sh
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "==============================================================="
echo " C10 setup   (repo: $REPO_DIR)"
echo "==============================================================="

# --- 1. enable conda for all future shells ------------------------------------
echo
echo "[1/3] Enabling conda"
if [ ! -f "$CONDA_SH" ]; then
  echo "      ERROR: cannot find conda at $CONDA_SH — is this the IFB cluster?" >&2
  exit 1
fi
source "$CONDA_SH"
/shared/software/miniconda/bin/conda init bash >/dev/null 2>&1 || true
conda config --set auto_activate_base false 2>/dev/null || true
echo "      conda ready (added to ~/.bashrc for future shells)."

# --- 2. make single_cell discoverable + verify --------------------------------
echo
echo "[2/3] Locating and verifying the 'single_cell' environment"
if ! conda config --show envs_dirs 2>/dev/null | grep -q "$ENVROOT"; then
  conda config --append envs_dirs "$ENVROOT"
fi
conda activate single_cell 2>/dev/null || conda activate "$ENVROOT/single_cell"
if [ "${CONDA_PREFIX:-}" != "$ENVROOT/single_cell" ]; then
  echo "      ERROR: could not activate single_cell (CONDA_PREFIX=${CONDA_PREFIX:-unset})" >&2
  exit 1
fi
python - <<'PY'
import importlib.util, sys
mods = ["scanpy","anndata","scvi","cell2location","celltypist",
        "infercnvpy","squidpy","harmonypy","liana","decoupler","torch"]
missing = [m for m in mods if importlib.util.find_spec(m) is None]
if missing:
    print("      MISSING packages:", missing); sys.exit(1)
print(f"      env OK — python {sys.version.split()[0]}, all {len(mods)} key packages import.")
PY
[ $? -ne 0 ] && exit 1

# --- 3. register the Jupyter kernel -------------------------------------------
echo
echo "[3/3] Registering the 'single_cell' Jupyter kernel"
python -m ipykernel install --user --name single_cell \
       --display-name "Python (single_cell)" >/dev/null 2>&1 \
  && echo "      kernel registered." \
  || echo "      (kernel registration skipped — not fatal; pick the env in VS Code instead)"

echo
echo "==============================================================="
echo " Setup complete."
echo "   • Activate anytime:   conda activate single_cell"
echo "   • Add helpers to path in a notebook:"
echo "       import sys; sys.path.insert(0, '$REPO_DIR/src')"
echo "   • Open notebooks/level1/01_snrna_analysis_student.ipynb in VS Code (INSTALL.md)."
echo "   • Pull updates later with:  git -C $REPO_DIR pull"
echo "   • Reminder: heavy compute via Slurm, NOT the login node (INSTALL.md §5)."
echo "==============================================================="
