#!/usr/bin/env python3

import csv
import statistics as st
from collections import Counter
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from docx import Document
from docx.shared import Inches, Pt

BASE = Path("... base path ...")
DATA = BASE / " ... "
OUT = BASE / "03_results/01_Neoantigen_Burden"
OUT.mkdir(parents=True, exist_ok=True)

SOURCE = DATA / "03_TCGA_BRCA_Neoantigen_Patient_Summary.tsv"

with SOURCE.open(newline="") as f:
    patients = list(csv.DictReader(f, delimiter="\t"))

assert len(patients) == 80
assert len({r["patient_id"] for r in patients}) == 80

def classify(value):
    v = value.strip().lower().replace("-", "_")
    if "missense" in v:
        return "Missense"
    if "frameshift" in v or v == "fs":
        return "Frameshift"
    if "inframe" in v or "in_frame" in v:
        return "Inframe"
    return "Other"

def write_tsv(path, rows, fields):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

records = []
validation = []

for patient in patients:
    pid = patient["patient_id"]
    files = list((DATA / pid / "MHC_Class_I").glob("*.MHC_I.filtered.tsv"))

    if len(files) != 1:
        raise RuntimeError(f"{pid}: expected one filtered file, found {len(files)}")

    with files[0].open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = reader.fieldnames or []

        peptide_col = next((x for x in ["MT Epitope Seq", "MT Epitope"] if x in columns), None)
        hla_col = next((x for x in ["HLA Allele", "HLA"] if x in columns), None)
        type_col = next((x for x in ["Variant Type", "Mutation Type", "variant_type"] if x in columns), None)

        if not all([peptide_col, hla_col, type_col]):
            raise RuntimeError(f"{pid}: missing peptide, HLA, or variant-type column: {columns}")

        pairs = {}
        raw_rows = 0

        for row in reader:
            raw_rows += 1
            variant = tuple(row[k] for k in [
                "Chromosome", "Start", "Stop", "Reference", "Variant"
            ])
            key = (variant, row[peptide_col].strip(), row[hla_col].strip())
            category = classify(row[type_col])

            if not all(key[1:]):
                raise RuntimeError(f"{pid}: missing peptide or HLA value")

            if key in pairs and pairs[key] != category:
                raise RuntimeError(f"{pid}: conflicting variant classes for {key}")

            pairs[key] = category

    expected_rows = int(patient["neo_filtered_rows"])
    expected_pairs = int(patient["neo_peptide_hla_count"])

    if raw_rows != expected_rows or len(pairs) != expected_pairs:
        raise RuntimeError(
            f"{pid}: row/pair mismatch; observed {raw_rows}/{len(pairs)}, "
            f"expected {expected_rows}/{expected_pairs}"
        )

    counts = Counter(pairs.values())

    record = {
        "patient_id": pid,
        "peptide_hla_pairs": expected_pairs,
        "unique_peptides": int(patient["neo_unique_peptide_count"]),
        "neoantigen_variants": int(patient["neo_variant_count"]),
        "missense_pairs": counts["Missense"],
        "frameshift_pairs": counts["Frameshift"],
        "inframe_pairs": counts["Inframe"],
        "other_pairs": counts["Other"],
    }

    assert sum(counts.values()) == expected_pairs

    records.append(record)
    validation.append({
        "patient_id": pid,
        "filtered_rows": raw_rows,
        "unique_pairs": len(pairs),
        "status": "PASS"
    })

records.sort(key=lambda r: (-r["peptide_hla_pairs"], r["patient_id"]))

for rank, row in enumerate(records, 1):
    row["rank"] = rank

patient_fields = [
    "rank", "patient_id", "peptide_hla_pairs",
    "unique_peptides", "neoantigen_variants",
    "missense_pairs", "frameshift_pairs",
    "inframe_pairs", "other_pairs"
]

write_tsv(
    OUT / "02_Neoantigen_Burden_Patient_Ranking.tsv",
    records,
    patient_fields
)

write_tsv(
    OUT / "03_Neoantigen_Burden_Validation.tsv",
    validation,
    ["patient_id", "filtered_rows", "unique_pairs", "status"]
)

metrics = [
    ("peptide_hla_pairs", "Peptide-HLA pairs"),
    ("unique_peptides", "Unique mutant peptides"),
    ("neoantigen_variants", "Neoantigen-generating variants"),
    ("missense_pairs", "Missense-derived pairs"),
    ("frameshift_pairs", "Frameshift-derived pairs"),
    ("inframe_pairs", "In-frame indel-derived pairs"),
    ("other_pairs", "Other-derived pairs")
]

summary = []

for field, label in metrics:
    values = sorted(r[field] for r in records)
    q = st.quantiles(values, n=4, method="inclusive")

    summary.append({
        "metric": label,
        "n": len(values),
        "mean": round(st.mean(values), 2),
        "median": round(st.median(values), 2),
        "q1": round(q[0], 2),
        "q3": round(q[2], 2),
        "min": min(values),
        "max": max(values),
        "total": sum(values)
    })

summary_fields = [
    "metric", "n", "mean", "median", "q1", "q3", "min", "max", "total"
]

write_tsv(
    OUT / "01_Neoantigen_Burden_Cohort_Summary.tsv",
    summary,
    summary_fields
)

plt.rcParams.update({"font.size": 10})

values = np.array([r["peptide_hla_pairs"] for r in records])
median = np.median(values)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.hist(values, bins=20, color=plt.cm.viridis(0.65), edgecolor="black", linewidth=0.5)
ax.set_xlabel("Predicted peptide-HLA pairs per patient")
ax.set_ylabel("Number of patients")
ax.set_title("Neoantigen burden")
fig.tight_layout()
fig.savefig(OUT / "01_Neoantigen_Burden_Distribution.png", dpi=300)
plt.close(fig)

categories = [
    ("missense_pairs", "Missense", "#4DBBD5"),
    ("frameshift_pairs", "Frameshift", "#E64B35"),
    ("inframe_pairs", "In-frame indel", "#00A087"),
    ("other_pairs", "Other", "#3C5488")
]

fig, ax = plt.subplots(figsize=(11, 4.5))
x = np.arange(len(records))
bottom = np.zeros(len(records))

for field, label, color in categories:
    y = np.array([r[field] for r in records])
    if y.sum() == 0:
        continue
    ax.bar(x, y, bottom=bottom, width=0.9, color=color, label=label, edgecolor="black", linewidth=0.5)
    bottom += y

ax.set_xlabel("Patients ranked by neoantigen burden")
ax.set_ylabel("Predicted peptide-HLA pairs")
ax.set_title("Neoantigen burden by variant class")
ax.set_xlim(-0.7, len(records) - 0.3)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "02_Neoantigen_Burden_Patient_Ranking.png", dpi=300)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.5))
bins = np.geomspace(values.min(), values.max(), 20)
ax.hist(values, bins=bins, color=plt.cm.viridis(0.65), edgecolor="black", linewidth=0.5)
ax.set_xscale("log")
ax.set_xlabel("Predicted peptide-HLA pairs per patient (log scale)")
ax.set_ylabel("Number of patients")
ax.set_title("Neoantigen burden distribution")
fig.tight_layout()
fig.savefig(OUT / "03_Neoantigen_Burden_Log_Distribution.png", dpi=300)
plt.close(fig)

doc = Document()
section = doc.sections[0]
section.page_width = Inches(11.7)
section.page_height = Inches(8.3)
section.top_margin = Inches(0.65)
section.bottom_margin = Inches(0.65)
section.left_margin = Inches(0.65)
section.right_margin = Inches(0.65)

normal = doc.styles["Normal"]
normal.font.name = "Arial"
normal.font.size = Pt(9)

def add_table(headers, rows):
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"

    for cell, value in zip(table.rows[0].cells, headers):
        cell.text = str(value)

    for record in rows:
        cells = table.add_row().cells
        for cell, value in zip(cells, record):
            cell.text = str(value)

    for row in table.rows:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.name = "Arial"
                    run.font.size = Pt(8)

doc.add_heading("Neoantigen Burden - TCGA-BRCA (n=80)", 0)

doc.add_paragraph(
    "Burden is defined using pVACseq-filtered mutant peptide-HLA predictions. "
    "The original pVACseq run used a binding threshold of 50,000 nM "
    "and an exploratory percentile-threshold strategy. "
    "These counts should not be interpreted as strong-binding neoantigens."
)

doc.add_heading("Table 1. Cohort Summary", level=2)

add_table(
    ["Metric", "N", "Mean", "Median", "Q1", "Q3", "Min", "Max", "Total"],
    [[r[k] for k in summary_fields] for r in summary]
)

doc.add_heading("Table 2. Patient-Level Neoantigen Burden", level=2)

display_fields = [
    "rank", "patient_id", "peptide_hla_pairs",
    "unique_peptides", "neoantigen_variants",
    "missense_pairs", "frameshift_pairs",
    "inframe_pairs", "other_pairs"
]

add_table(
    ["Rank", "Patient", "pMHC", "Peptides", "Variants",
     "Missense", "Frameshift", "In-frame", "Other"],
    [[r[k] for k in display_fields] for r in records]
)

doc.save(OUT / "04_Neoantigen_Burden_Tables.docx")

print("Patients:", len(records))
print("Total peptide-HLA pairs:", sum(r["peptide_hla_pairs"] for r in records))
print("Median burden:", median)
print("Missense pairs:", sum(r["missense_pairs"] for r in records))
print("Frameshift pairs:", sum(r["frameshift_pairs"] for r in records))
print("In-frame indel pairs:", sum(r["inframe_pairs"] for r in records))
print("Other pairs:", sum(r["other_pairs"] for r in records))
print("Validation: PASS")
print("Outputs:", OUT)

