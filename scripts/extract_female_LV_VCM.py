#!/usr/bin/env python3
"""
Extract selected columns from the full Matrix Market count matrix.

Input
-----
data_raw/multi_study_snRNA_Seq_counts_mm_file.txt
results/tables/female_LV_VCM_matrix_column_indices.csv

Output
------
data/female_LV_VCM_counts.mtx

The R-generated column indices are 1-based, matching Matrix Market
coordinate indices. The script streams the ~15 GB source file and
never materializes the full matrix in memory.
"""

from pathlib import Path
import csv
import shutil
import tempfile

SOURCE = Path("data_raw/multi_study_snRNA_Seq_counts_mm_file.txt")
INDEX_FILE = Path("results/tables/female_LV_VCM_matrix_column_indices.csv")
OUTPUT = Path("data/female_LV_VCM_counts.mtx")

if not SOURCE.exists():
    raise FileNotFoundError(f"Missing source count matrix: {SOURCE}")
if not INDEX_FILE.exists():
    raise FileNotFoundError(f"Missing target column index file: {INDEX_FILE}")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with INDEX_FILE.open(newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh)
    target_cols = [int(row["column_index"]) for row in reader]

if not target_cols:
    raise RuntimeError("No target columns were supplied.")

# Map original 1-based column index -> reduced 1-based column index.
col_map = {old: new for new, old in enumerate(target_cols, start=1)}

with SOURCE.open("r", encoding="utf-8") as src:
    header = src.readline()
    if not header.startswith("%%MatrixMarket"):
        raise RuntimeError("Input does not look like a Matrix Market file.")

    comments = []
    line = src.readline()
    while line.startswith("%"):
        comments.append(line)
        line = src.readline()

    n_rows, n_cols, _ = map(int, line.split())

    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, encoding="utf-8", newline=""
    ) as body:
        body_path = Path(body.name)
        kept_nnz = 0

        for raw in src:
            if not raw.strip():
                continue
            row, col, value = raw.split()
            col = int(col)

            new_col = col_map.get(col)
            if new_col is None:
                continue

            body.write(f"{row} {new_col} {value}\n")
            kept_nnz += 1

with OUTPUT.open("w", encoding="utf-8", newline="") as out:
    out.write(header)
    for comment in comments:
        out.write(comment)
    out.write(f"{n_rows} {len(target_cols)} {kept_nnz}\n")

    with body_path.open("r", encoding="utf-8") as body:
        shutil.copyfileobj(body, out, length=1024 * 1024)

body_path.unlink(missing_ok=True)

print(f"Wrote: {OUTPUT}")
print(f"Rows: {n_rows}")
print(f"Columns retained: {len(target_cols)}")
print(f"Non-zero entries retained: {kept_nnz}")
