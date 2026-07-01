"""Derive a student notebook from a completed, executed solution notebook: keep all
markdown (TASK/HINT/QUESTION/CHECKPOINT guidance), blank out every code cell to a single
placeholder comment, clear all outputs/execution_counts. Mechanical and safe by construction.

One exception to "keep all markdown": instructor-only answer text (e.g. the Section 11
write-up, which necessarily embeds real computed numbers) is wrapped in the solution between
`<!-- INSTRUCTOR-ONLY -->` and `<!-- /INSTRUCTOR-ONLY -->` HTML-comment sentinels. Those
regions are stripped here so the student keeps the surrounding TASK prompt but not the answer.
"""
import re
import sys

import nbformat as nbf

SRC, DST = sys.argv[1], sys.argv[2]

INSTRUCTOR_BLOCK = re.compile(
    r"\n?<!--\s*INSTRUCTOR-ONLY\s*-->.*?<!--\s*/INSTRUCTOR-ONLY\s*-->\n?",
    re.DOTALL,
)


def strip_instructor(md):
    return INSTRUCTOR_BLOCK.sub("\n", md)

# One short, generic placeholder comment per section -- matches the template's bare
# "# Your ... here" convention (no `# TODO:`), inferred from the nearest preceding TASK line.
nb = nbf.read(SRC, as_version=4)
out_cells = []
last_task_text = "this step"

for cell in nb.cells:
    if cell.cell_type == "markdown":
        m = re.search(r"TASK\s+[\d.]+:\*\*\s*(.+?)(?:\n|$)", cell.source)
        if m:
            last_task_text = m.group(1).strip().rstrip(".")
        out_cells.append(nbf.v4.new_markdown_cell(strip_instructor(cell.source)))
    else:
        placeholder = f"# Your code for: {last_task_text}\n"
        out_cells.append(nbf.v4.new_code_cell(placeholder))

nb.cells = out_cells
nbf.write(nb, DST)
print(f"Wrote student notebook: {DST} ({len(out_cells)} cells)")
