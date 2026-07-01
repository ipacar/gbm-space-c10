---
name: feedback-pip-no-user
description: Always use --no-user when running pip install in conda envs; installing to ~/.local broke other environments
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d29c3aa3-e8f4-4b3a-a213-0683aceb9983
---

Always add `--no-user` to every `pip install` call inside a conda environment on this cluster.

**Why:** A previous session installed packages without `--no-user`, which caused some packages to land in `~/.local/lib/python3.11/site-packages` instead of the conda env's site-packages. The `~/.local` path is shared across ALL conda environments (Python shares the same `~/.local` regardless of which env is active), so packages installed there can silently override or conflict with packages in other envs. This broke other unrelated environments on the cluster — a severe enough incident that the user flagged it explicitly.

**How to apply:** Whenever calling pip in this project (or any other conda env on this cluster), always use:
```bash
/path/to/env/bin/pip install --no-user <package>
```
or equivalently via conda run:
```bash
conda run -n env_name pip install --no-user <package>
```
Never use bare `pip install <package>` which respects PYTHONUSERSITE and can write to `~/.local`.
