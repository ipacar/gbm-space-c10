# ⛔ do_not_use — instructor solutions

**Students: do not open these files.** This folder holds the completed **solution**
notebooks for the C10 project. Using them defeats the entire point of the project — you learn
by working the analysis out yourself (and by discussing with instructors and classmates).

The whole design of Level 1 is that you reach the biology through your **own** clustering,
markers, and reasoning — including discovering the dataset's structure without being told the
source paper up front. Reading the solutions short-circuits that.

Contents:
- `01_snrna_analysis_solution.ipynb` — Level 1 completed solution
- `02_spatial_cell2location_solution.ipynb` — Level 2 completed solution
- `03_xenium_organoid_analysis.ipynb` — Level 3 (Xenium) completed solution
- `level1_summary_figure.png`, `level2_summary_figure.png` — the completed-analysis result
  figures (i.e. what your own final figure should roughly look like — don't peek)

**Exception —** `grch38_gene_positions.parquet` **is NOT a solution.** It's a reference table of
gene genomic coordinates that the inferCNV step (Level 1) needs as *input*. It's parked in this
folder only so a `git clone` of the repo has it. You may load it:
```python
gene_pos = pd.read_parquet(".../do_not_use/grch38_gene_positions.parquet")
```

These are kept in the repository only for **instructors** (teaching, grading, reference). If you
are a student and opened this by accident, close it and go back to your `*_student.ipynb`
notebooks. 🙂
