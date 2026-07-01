# IFB cluster quickstart — VS Code, passwordless login, conda test

A short, project-agnostic setup for working on the **IFB Core cluster** from **VS Code**. Do this
once; it applies to any course project. ~15 minutes.

Cluster login host: `core.cluster.france-bioinformatique.fr` (you need an IFB account).

---

## 1. Install VS Code + extensions

1. Install **VS Code** on your laptop: <https://code.visualstudio.com/>
2. Open the Extensions panel (`Ctrl/Cmd+Shift+X`) and install **Remote - SSH** (Microsoft).
3. The **Python**, **Jupyter**, and (optional) **Claude** extensions are installed *later*, **on
   the cluster**, once you're connected — VS Code shows an "Install in SSH: …" button for each.

---

## 2. Passwordless login with an SSH key

Instead of typing your password every time, register your laptop with an SSH key.

**a. Make a key** (skip if `~/.ssh/id_ed25519` already exists). In a laptop terminal:
```bash
ssh-keygen -t ed25519          # press Enter through the prompts
```

**b. Copy the public key to the cluster** (asks for your password one last time):
```bash
# macOS / Linux:
ssh-copy-id <your-username>@core.cluster.france-bioinformatique.fr

# Windows PowerShell (if ssh-copy-id is missing):
type $env:USERPROFILE\.ssh\id_ed25519.pub | ssh <your-username>@core.cluster.france-bioinformatique.fr "mkdir -p ~/.ssh && cat >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys"
```

**c. Tell your SSH config to use the key.** Open `~/.ssh/config` (in VS Code: `F1` →
"Remote-SSH: Open SSH Configuration File") and add:
```
Host ifb
    HostName core.cluster.france-bioinformatique.fr
    User <your-username>
    IdentityFile ~/.ssh/id_ed25519
```

**d. Test — should log in with no password:**
```bash
ssh ifb
```
> If it still asks for a password: on the cluster run `chmod 700 ~/.ssh && chmod 600 ~/.ssh/authorized_keys`.
> If you gave the key a passphrase, you'll be asked for *that* (not your account password); run
> `ssh-add ~/.ssh/id_ed25519` once per session to cache it.

---

## 3. Connect VS Code to the cluster

1. `F1` → **Remote-SSH: Connect to Host…** → pick **`ifb`**. No password prompt.
2. **File → Open Folder…** to open your working directory on the cluster.
3. Install the **Python** and **Jupyter** extensions when VS Code offers "Install in SSH".

---

## 4. Minimal test — can you build a conda environment?

This confirms conda works for your account (independent of any project).

**a. Enable conda once** (in a cluster terminal — the VS Code terminal is fine):
```bash
/shared/software/miniconda/bin/conda init bash
```
Then **open a new terminal** (or `source ~/.bashrc`) so `conda` is available.

**b. Create a tiny test environment and verify it:**
```bash
conda create -y -n envtest python=3.11 numpy
conda activate envtest
python -c "import numpy; print('conda works — numpy', numpy.__version__)"
```
You should see `conda works — numpy X.Y.Z`. 🎉

**c. Clean up** (optional):
```bash
conda deactivate
conda env remove -y -n envtest
```

> **Rule of thumb:** the login node is for editing and light setup only. Anything heavy
> (training, large data, long jobs) must go through **Slurm** — e.g. grab a compute node with
> `srun --account=<your-project-account> --partition=fast --cpus-per-task=4 --mem=16G --time=02:00:00 --pty bash`.
> Your project's setup guide lists its Slurm account and any shared environment to use instead of
> building your own.

---

That's it — VS Code connects passwordlessly, extensions are in place, and conda works. You're
ready to open a specific project and follow its own setup notes.
