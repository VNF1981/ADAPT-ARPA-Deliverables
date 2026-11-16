#!/usr/bin/env python3
import csv
import math
import os
import statistics as st
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from docx import Document
from docx.shared import Inches, Pt, RGBColor

BASE = Path(os.environ.get("BRCA80_BASE", "...path..."))
DATA = BASE / "... additional path ..."
OUT = BASE / "03_results/02_Peptide_HLA_Binding"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE = DATA / "03_TCGA_BRCA_Neoantigen_Patient_Summary.tsv"

with SOURCE.open(newline="") as f:
    patients = list(csv.DictReader(f, delimiter="\t"))
if len(patients) != 80 or len({p["patient_id"] for p in patients}) != 80:
    raise RuntimeError("Expected 80 distinct patients in the cohort summary")

AFFINITY_CHOICES = ["Best MT IC50 Score", "Median MT IC50 Score", "MT IC50 Score"]
PERCENTILE_CHOICES = ["Best MT Percentile", "Median MT Percentile", "MT Percentile"]
CLASSES = ["Strong (≤50 nM)", "Intermediate (50–500 nM)", "Weak (>500 nM)", "Unavailable"]
LOCI = ["HLA-A", "HLA-B", "HLA-C"]
COLORS = [plt.cm.viridis(v) for v in (0.18, 0.48, 0.78)]

def choose(cols, choices, required=False):
    name = next((name for name in choices if name in cols), None)
    if required and name is None:
        raise RuntimeError(f"Required column missing; expected one of {choices}; found {cols}")
    return name

def number(value):
    try:
        x = float(str(value).strip())
        return x if math.isfinite(x) and x >= 0 else None
    except (ValueError, TypeError):
        return None

def bind_class(value):
    if value is None:
        return "Unavailable"
    if value <= 50:
        return CLASSES[0]
    if value <= 500:
        return CLASSES[1]
    return CLASSES[2]

def locus_of(allele):
    v = allele.upper().strip()
    for locus in LOCI:
        if v.startswith(locus + "*") or v.startswith(locus + ":"):
            return locus
    raise RuntimeError(f"Unrecognized HLA-A/B/C allele: {allele}")

def fmt(x, places=2):
    return "NA" if x is None else f"{x:.{places}f}"

def pct(a, b):
    return (100 * a / b) if b else None

def median(xs):
    return st.median(xs) if xs else None

def write_tsv(path, rows, fields):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

patient_rows = []
all_candidates = []
col_usage = Counter()
raw_total = 0

for p in patients:
    pid = p["patient_id"]
    files = list((DATA / pid / "MHC_Class_I").glob("*.MHC_I.filtered.tsv"))
    if len(files) != 1:
        raise RuntimeError(f"{pid}: expected exactly one filtered pVACseq file; found {len(files)}")
    with files[0].open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        cols = reader.fieldnames or []
        peptide_col = choose(cols, ["MT Epitope Seq", "MT Epitope"], True)
        hla_col = choose(cols, ["HLA Allele", "HLA"], True)
        affinity_col = choose(cols, AFFINITY_CHOICES, True)
        percentile_col = choose(cols, PERCENTILE_CHOICES)
        col_usage[(affinity_col, percentile_col or "Unavailable")] += 1
        required = ["Chromosome", "Start", "Stop", "Reference", "Variant"]
        if any(k not in cols for k in required):
            raise RuntimeError(f"{pid}: missing mutation identity fields")
        candidates = {}
        row_count = 0
        for row in reader:
            row_count += 1
            mutation = tuple(row[k].strip() for k in required)
            peptide = row[peptide_col].strip()
            hla = row[hla_col].strip()
            if not peptide or not hla:
                raise RuntimeError(f"{pid}: empty peptide or HLA allele")
            locus = locus_of(hla)
            key = (mutation, peptide, hla)
            affinity = number(row[affinity_col])
            percentile = number(row[percentile_col]) if percentile_col else None
            if key not in candidates:
                candidates[key] = {"patient_id": pid, "hla_locus": locus, "ic50": affinity, "percentile": percentile}
            else:
                rec = candidates[key]
                if affinity is not None and (rec["ic50"] is None or affinity < rec["ic50"]):
                    rec["ic50"] = affinity
                if percentile is not None and (rec["percentile"] is None or percentile < rec["percentile"]):
                    rec["percentile"] = percentile
    expected_rows = int(p["neo_filtered_rows"])
    expected_pairs = int(p["neo_peptide_hla_count"])
    if row_count != expected_rows or len(candidates) != expected_pairs:
        raise RuntimeError(f"{pid}: observed {row_count} rows and {len(candidates)} candidates; expected {expected_rows} and {expected_pairs}")
    raw_total += row_count
    values = list(candidates.values())
    all_candidates.extend(values)
    classes = Counter(bind_class(r["ic50"]) for r in values)
    assessed = len(values) - classes["Unavailable"]
    affinities = [r["ic50"] for r in values if r["ic50"] is not None]
    percentiles = [r["percentile"] for r in values if r["percentile"] is not None]
    patient_rows.append({
        "patient_id": pid, "candidates_total": len(values), "affinity_assessed": assessed,
        "strong_n": classes[CLASSES[0]], "intermediate_n": classes[CLASSES[1]],
        "weak_n": classes[CLASSES[2]], "affinity_unavailable_n": classes["Unavailable"],
        "strong_pct_assessed": fmt(pct(classes[CLASSES[0]], assessed)),
        "intermediate_pct_assessed": fmt(pct(classes[CLASSES[1]], assessed)),
        "weak_pct_assessed": fmt(pct(classes[CLASSES[2]], assessed)),
        "median_ic50_nM": fmt(median(affinities)),
        "median_percentile_rank": fmt(median(percentiles)),
        "percentile_available_n": len(percentiles)
    })

if sum(r["candidates_total"] for r in patient_rows) != 5533:
    raise RuntimeError("Candidate total differs from validated Analysis 1 total (5,533)")

cohort_counts = Counter(bind_class(r["ic50"]) for r in all_candidates)
assessed_total = len(all_candidates) - cohort_counts["Unavailable"]
affinity_values = [r["ic50"] for r in all_candidates if r["ic50"] is not None]
percentile_values = [r["percentile"] for r in all_candidates if r["percentile"] is not None]
cohort_rows = [
    {"metric": "Patients", "value": 80, "denominator": "", "percentage": ""},
    {"metric": "Filtered pVACseq rows", "value": raw_total, "denominator": "", "percentage": ""},
    {"metric": "Distinct mutation–peptide–HLA candidates", "value": len(all_candidates), "denominator": "", "percentage": ""},
    {"metric": "Candidates with available IC50", "value": assessed_total, "denominator": len(all_candidates), "percentage": fmt(pct(assessed_total, len(all_candidates)))},
]
for label in CLASSES:
    value = cohort_counts[label]
    denom = assessed_total if label != "Unavailable" else len(all_candidates)
    cohort_rows.append({"metric": label, "value": value, "denominator": denom, "percentage": fmt(pct(value, denom))})
cohort_rows.extend([
    {"metric": "Median IC50 (nM; available)", "value": fmt(median(affinity_values)), "denominator": len(affinity_values), "percentage": ""},
    {"metric": "Candidates with percentile rank", "value": len(percentile_values), "denominator": len(all_candidates), "percentage": fmt(pct(len(percentile_values), len(all_candidates)))},
    {"metric": "Median percentile rank (available)", "value": fmt(median(percentile_values)), "denominator": len(percentile_values), "percentage": ""}
])

locus_rows = []
for locus in LOCI:
    subset = [r for r in all_candidates if r["hla_locus"] == locus]
    counts = Counter(bind_class(r["ic50"]) for r in subset)
    assessed = len(subset) - counts["Unavailable"]
    locus_rows.append({
        "hla_locus": locus, "candidates_total": len(subset), "affinity_assessed": assessed,
        "strong_n": counts[CLASSES[0]], "intermediate_n": counts[CLASSES[1]],
        "weak_n": counts[CLASSES[2]], "affinity_unavailable_n": counts["Unavailable"],
        "strong_pct_assessed": fmt(pct(counts[CLASSES[0]], assessed)),
        "intermediate_pct_assessed": fmt(pct(counts[CLASSES[1]], assessed)),
        "weak_pct_assessed": fmt(pct(counts[CLASSES[2]], assessed)),
        "median_ic50_nM": fmt(median([r["ic50"] for r in subset if r["ic50"] is not None])),
        "median_percentile_rank": fmt(median([r["percentile"] for r in subset if r["percentile"] is not None]))
    })

patient_rows.sort(key=lambda r: (-r["candidates_total"], r["patient_id"]))
cohort_fields = ["metric", "value", "denominator", "percentage"]
patient_fields = list(patient_rows[0])
locus_fields = list(locus_rows[0])
write_tsv(OUT / "01_Peptide_HLA_Binding_Cohort_Summary.tsv", cohort_rows, cohort_fields)
write_tsv(OUT / "02_Peptide_HLA_Binding_Patient_Summary.tsv", patient_rows, patient_fields)
write_tsv(OUT / "03_Peptide_HLA_Binding_HLA_Locus_Summary.tsv", locus_rows, locus_fields)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

positive = np.array([x for x in affinity_values if x > 0])
zero_count = len(affinity_values) - len(positive)
if not len(positive):
    raise RuntimeError("No positive IC50 values available for affinity histogram")
bins = np.geomspace(max(positive.min() * 0.99, 1e-6), positive.max() * 1.01, 31)
fig, ax = plt.subplots(figsize=(7.5, 4.5))
ax.hist(positive, bins=bins, color=plt.cm.viridis(0.65), edgecolor="black", linewidth=0.5)
ax.set_xscale("log")
ax.axvline(50, color="black", linestyle="--", linewidth=0.8)
ax.axvline(500, color="black", linestyle=":", linewidth=0.8)
ax.set_xlabel("Predicted mutant peptide–HLA IC50 (nM, log scale)")
ax.set_ylabel("Candidate combinations")
ax.set_title("Predicted peptide–HLA binding affinity")
fig.tight_layout()
fig.savefig(OUT / "01_Peptide_HLA_Binding_Affinity_Distribution.png", dpi=300)
plt.close(fig)

fig, ax = plt.subplots(figsize=(12, 4.5))
x = np.arange(len(patient_rows))
bottom = np.zeros(len(patient_rows))
for col, label, color in zip(["strong_pct_assessed", "intermediate_pct_assessed", "weak_pct_assessed"], CLASSES[:3], COLORS):
    heights = np.array([float(r[col]) if r[col] != "NA" else 0 for r in patient_rows])
    ax.bar(x, heights, bottom=bottom, width=0.9, color=color, edgecolor="black", linewidth=0.35, label=label)
    bottom += heights
ax.set_xlim(-0.7, len(patient_rows) - 0.3)
ax.set_ylim(0, 103)
ax.set_xlabel("Patients ranked by candidate neoantigen burden")
ax.set_ylabel("Percentage of candidates with available IC50")
ax.set_title("Predicted binding categories by patient")
ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.16))
fig.tight_layout()
fig.savefig(OUT / "02_Peptide_HLA_Binding_Patient_Categories.png", dpi=300, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 4.5))
x = np.arange(3)
bottom = np.zeros(3)
for col, label, color in zip(["strong_pct_assessed", "intermediate_pct_assessed", "weak_pct_assessed"], CLASSES[:3], COLORS):
    heights = np.array([float(r[col]) if r[col] != "NA" else 0 for r in locus_rows])
    ax.bar(x, heights, bottom=bottom, width=0.68, color=color, edgecolor="black", linewidth=0.5, label=label)
    bottom += heights
ax.set_xticks(x, LOCI)
ax.set_ylim(0, 103)
ax.set_ylabel("Percentage of candidates with available IC50")
ax.set_title("Predicted binding categories by HLA locus")
ax.legend(frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.16))
fig.tight_layout()
fig.savefig(OUT / "03_Peptide_HLA_Binding_HLA_Locus.png", dpi=300, bbox_inches="tight")
plt.close(fig)

doc = Document()
section = doc.sections[0]
section.page_width = Inches(11.7)
section.page_height = Inches(8.3)
section.top_margin = Inches(0.65)
section.bottom_margin = Inches(0.65)
section.left_margin = Inches(0.65)
section.right_margin = Inches(0.65)
for name in ("Normal", "Title", "Heading 1", "Heading 2"):
    style = doc.styles[name]
    style.font.name = "Times New Roman"
    style.font.size = Pt(9)
    style.font.color.rgb = RGBColor(0, 0, 0)
for name in ("Title", "Heading 1", "Heading 2"):
    doc.styles[name].paragraph_format.space_before = Pt(5)
    doc.styles[name].paragraph_format.space_after = Pt(3)

def add_table(title, fields, rows, headers=None):
    p = doc.add_paragraph()
    r = p.add_run(title)
    r.bold = True
    t = doc.add_table(rows=1, cols=len(fields))
    t.style = "Table Grid"
    for cell, value in zip(t.rows[0].cells, headers or fields):
        cell.text = str(value)
    for record in rows:
        cells = t.add_row().cells
        for cell, field in zip(cells, fields):
            cell.text = str(record[field])
    for row in t.rows:
        for cell in row.cells:
            for para in cell.paragraphs:
                para.paragraph_format.space_after = Pt(0)
                for run in para.runs:
                    run.font.name = "Times New Roman"
                    run.font.size = Pt(9)
                    run.font.color.rgb = RGBColor(0, 0, 0)
    for cell in t.rows[0].cells:
        for para in cell.paragraphs:
            for run in para.runs:
                run.bold = True

doc.add_paragraph("Peptide–HLA Binding Affinity — TCGA-BRCA (n=80)")
doc.add_paragraph("Each candidate is a distinct somatic mutation–mutant peptide–HLA combination, consistent with Analysis 1. Binding categories use predicted mutant-peptide IC50: strong ≤50 nM; intermediate >50–500 nM; weak >500 nM. Percentages use candidates with available IC50 as denominator. These are computational predictions, not measured presentation or immunogenicity.")
add_table("Table 1. Cohort binding summary", cohort_fields, cohort_rows, ["Metric", "Value", "Denominator", "%"])
add_table("Table 2. HLA locus summary", locus_fields, locus_rows)
add_table("Table 3. Patient-level binding summary", patient_fields, patient_rows)
doc.save(OUT / "04_Peptide_HLA_Binding_Tables.docx")

print("Patients:", len(patient_rows))
print("Filtered pVACseq rows:", raw_total)
print("Distinct candidates:", len(all_candidates))
print("IC50 assessed:", assessed_total)
print("Strong (≤50 nM):", cohort_counts[CLASSES[0]])
print("Intermediate (50–500 nM):", cohort_counts[CLASSES[1]])
print("Weak (>500 nM):", cohort_counts[CLASSES[2]])
print("Unavailable IC50:", cohort_counts["Unavailable"])
print("Zero IC50 excluded from log histogram:", zero_count)
print("Affinity/percentile columns:", dict(col_usage))
print("Validation: PASS")
print("Outputs:", OUT)
