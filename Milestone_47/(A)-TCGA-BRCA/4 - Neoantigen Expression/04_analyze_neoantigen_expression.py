import csv
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from docx import Document
from docx.shared import Inches, Pt, RGBColor

BASE = Path(os.environ.get("BRCA80_BASE", "Path to parent dir"))
DATA = BASE / "02_data/02_pVACseq"
OUT = BASE / "03_results/04_Neoantigen_Expression"
OUT.mkdir(parents=True, exist_ok=True)
SUMMARY = DATA / "03_TCGA_BRCA_Neoantigen_Patient_Summary.tsv"
MANIFEST = BASE / "02_data/01_TCGA_BRCA_80_Input_Manifest.tsv"
TPM_THRESHOLD = 1.0
RNA_ALT_THRESHOLD = 3
BINDING_THRESHOLD = 500.0

EXP_COLS = {
    "transcript_tpm": ["Transcript Expression", "Transcript Expression (TPM)", "Transcript TPM", "transcript_expression", "transcript_tpm"],
    "gene_tpm": ["Gene Expression", "Gene Expression (TPM)", "Gene TPM", "gene_expression", "gene_tpm"],
    "rna_vaf": ["Tumor RNA VAF", "Tumor RNA Variant Allele Frequency", "tumor_rna_vaf"],
    "rna_depth": ["Tumor RNA Depth", "Tumor RNA Coverage", "tumor_rna_depth"],
    "rna_alt": ["Tumor RNA Alt Read Support", "Tumor RNA Alt Reads", "Tumor RNA ALT Reads", "tumor_rna_alt_reads"],
    "ic50": ["Best MT IC50 Score", "Median MT IC50 Score", "MT IC50 Score"]
}

def choose(columns, aliases):
    normalized = {"".join(c.lower().split()): c for c in columns}
    for name in aliases:
        if name in columns:
            return name
        result = normalized.get("".join(name.lower().split()))
        if result:
            return result
    return None

def value_number(value):
    if value is None or str(value).strip().lower() in {"", ".", "na", "n/a", "nan", "none", "null"}:
        return None
    try:
        number = float(str(value).strip())
        return number if math.isfinite(number) and number >= 0 else None
    except ValueError:
        return None

def read_table(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_table(name, rows, columns):
    with (OUT / name).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def percent(n, d):
    return round(100.0 * n / d, 2) if d else None

def display(value, places=2):
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.{places}f}"
    return str(value)

def median(values):
    return float(np.median(values)) if values else None

def tpm_label(record):
    if record["transcript_tpm"] is not None:
        return record["transcript_tpm"], "Transcript"
    if record["gene_tpm"] is not None:
        return record["gene_tpm"], "Gene proxy"
    return None, "Unavailable"

def make_doc(cohort, patient_rows, source_rows):
    doc = Document()
    sec = doc.sections[0]
    sec.top_margin = Inches(.65)
    sec.bottom_margin = Inches(.65)
    sec.left_margin = Inches(.7)
    sec.right_margin = Inches(.7)
    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(9)
    normal.font.color.rgb = RGBColor(0, 0, 0)
    for name in ["Title", "Heading 1", "Heading 2"]:
        style = doc.styles[name]
        style.font.name = "Times New Roman"
        style.font.size = Pt(9)
        style.font.color.rgb = RGBColor(0, 0, 0)
        style.font.bold = True
    doc.add_paragraph("Analysis 4 — Neoantigen Expression", style="Heading 1")
    doc.add_paragraph("RNA transcript expression is used when available; gene-level expression is labeled as a proxy. Expression-positive denotes measured TPM ≥1. RNA alternate-read evidence (≥3 reads), if available, is reported separately. Missing measurements are not treated as absence of expression. Predicted HLA binding and transcript expression do not establish presentation or T-cell recognition.")
    def add_table(title, headers, data):
        doc.add_paragraph(title, style="Heading 2")
        table = doc.add_table(rows=1, cols=len(headers))
        table.style = "Table Grid"
        for cell, label in zip(table.rows[0].cells, headers):
            cell.text = label
        for values in data:
            row = table.add_row()
            for cell, item in zip(row.cells, values):
                cell.text = display(item)
        for row in table.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    for run in paragraph.runs:
                        run.font.name = "Times New Roman"
                        run.font.size = Pt(8)
                        run.font.color.rgb = RGBColor(0, 0, 0)
    add_table("Table 1. Cohort summary", ["Metric", "Value", "Denominator", "Percentage"], [[r["metric"], r["value"], r["denominator"], r["percentage"]] for r in cohort])
    cols = ["patient_id", "subtype", "candidates_total", "tpm_assessed", "expressed_n", "expressed_pct_assessed", "unavailable_n", "median_available_tpm", "rna_alt_assessed", "rna_alt_supported_n", "expressed_and_binding_500_n"]
    add_table("Table 2. Patient-level summary", ["Patient", "Subtype", "Candidates", "TPM assessed", "Expressed", "Expressed %", "TPM missing", "Median TPM", "RNA ALT assessed", "RNA ALT ≥3", "Expressed & IC50 ≤500"], [[r[c] for c in cols] for r in patient_rows])
    add_table("Table 3. Expression data availability", ["Source field", "Patient files containing column", "Rows with numeric values"], [[r["field"], r["files_with_column"], r["rows_with_value"]] for r in source_rows])
    doc.save(OUT / "04_Neoantigen_Expression_Tables.docx")

summary = read_table(SUMMARY)
manifest = read_table(MANIFEST)
if len(summary) != 80 or len({r["patient_id"] for r in summary}) != 80:
    raise RuntimeError("Expected exactly 80 unique patients in the pVACseq summary")
if len(manifest) != 80 or {r["patient_id"] for r in manifest} != {r["patient_id"] for r in summary}:
    raise RuntimeError("The 80-patient manifest does not match the pVACseq summary")
subtypes = {r["patient_id"]: r.get("subtype", "") for r in manifest}
records = []
source_counts = Counter()
source_files = Counter()
raw_total = 0
for patient in summary:
    pid = patient["patient_id"]
    files = list((DATA / pid / "MHC_Class_I").glob("*.MHC_I.filtered.tsv"))
    if len(files) != 1:
        raise RuntimeError(f"{pid}: expected exactly one filtered pVACseq file, found {len(files)}")
    candidates = {}
    raw_rows = 0
    with files[0].open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = reader.fieldnames or []
        peptide_col = choose(columns, ["MT Epitope Seq", "MT Epitope"])
        hla_col = choose(columns, ["HLA Allele", "HLA"])
        required = ["Chromosome", "Start", "Stop", "Reference", "Variant"]
        if not peptide_col or not hla_col or any(c not in columns for c in required):
            raise RuntimeError(f"{pid}: candidate identifier columns missing")
        selected = {k: choose(columns, choices) for k, choices in EXP_COLS.items()}
        for key, column in selected.items():
            if column:
                source_files[key] += 1
        for row in reader:
            raw_rows += 1
            mutation = tuple(str(row[c]).strip() for c in required)
            peptide = str(row[peptide_col]).strip()
            hla = str(row[hla_col]).strip()
            if not peptide or not hla:
                raise RuntimeError(f"{pid}: missing peptide/HLA in candidate")
            key = (mutation, peptide, hla)
            measurements = {k: value_number(row[column]) if column else None for k, column in selected.items()}
            for k in EXP_COLS:
                if measurements[k] is not None:
                    source_counts[k] += 1
            if key not in candidates:
                candidates[key] = {"patient_id": pid, "subtype": subtypes[pid], "chromosome": mutation[0], "start": mutation[1], "stop": mutation[2], "reference": mutation[3], "variant": mutation[4], "mt_epitope": peptide, "hla_allele": hla, **measurements}
            else:
                existing = candidates[key]
                for k in EXP_COLS:
                    incoming = measurements[k]
                    previous = existing[k]
                    if incoming is not None and (previous is None or (k == "ic50" and incoming < previous)):
                        existing[k] = incoming
    if raw_rows != int(patient["neo_filtered_rows"]) or len(candidates) != int(patient["neo_peptide_hla_count"]):
        raise RuntimeError(f"{pid}: candidate counts do not match validated pVACseq summary")
    raw_total += raw_rows
    records.extend(candidates.values())
if len(records) != 5533 or raw_total != 5596:
    raise RuntimeError(f"Candidate count mismatch: {len(records)} unique candidates, {raw_total} filtered rows")

for record in records:
    tpm, source = tpm_label(record)
    record["expression_tpm"] = tpm
    record["expression_source"] = source
    record["expression_status"] = "Not assessed" if tpm is None else "Expressed" if tpm >= TPM_THRESHOLD else "Below threshold"
    alt = record["rna_alt"]
    record["rna_alt_status"] = "Not assessed" if alt is None else "Supported" if alt >= RNA_ALT_THRESHOLD else "Below threshold"
    record["binding_500"] = "NA" if record["ic50"] is None else "Yes" if record["ic50"] <= BINDING_THRESHOLD else "No"

fields = ["patient_id", "subtype", "chromosome", "start", "stop", "reference", "variant", "mt_epitope", "hla_allele", "transcript_tpm", "gene_tpm", "expression_tpm", "expression_source", "expression_status", "rna_vaf", "rna_depth", "rna_alt", "rna_alt_status", "ic50", "binding_500"]
write_table("03_Neoantigen_Expression_Candidates.tsv", records, fields)
source_rows = [{"field": key, "files_with_column": source_files[key], "rows_with_value": source_counts[key]} for key in EXP_COLS]
write_table("05_Expression_Field_Audit.tsv", source_rows, ["field", "files_with_column", "rows_with_value"])

if not any(r["expression_tpm"] is not None for r in records):
    print("Expression inputs unavailable: none of the 5,533 candidates have numeric transcript/gene TPM in the pVACseq filtered files.")
    print("Field audit:", OUT / "05_Expression_Field_Audit.tsv")
    raise SystemExit("Cannot perform RNA expression analysis without actual RNA expression data; no expression plots or Word summary were generated.")

patient_rows = []
for pid in sorted({r["patient_id"] for r in records}):
    subset = [r for r in records if r["patient_id"] == pid]
    with_tpm = [r for r in subset if r["expression_tpm"] is not None]
    positive = [r for r in with_tpm if r["expression_status"] == "Expressed"]
    with_alt = [r for r in subset if r["rna_alt"] is not None]
    patient_rows.append({
        "patient_id": pid, "subtype": subtypes[pid], "candidates_total": len(subset),
        "tpm_assessed": len(with_tpm), "expressed_n": len(positive),
        "expressed_pct_assessed": percent(len(positive), len(with_tpm)),
        "unavailable_n": len(subset) - len(with_tpm),
        "transcript_tpm_assessed": sum(r["transcript_tpm"] is not None for r in subset),
        "gene_proxy_assessed": sum(r["expression_source"] == "Gene proxy" for r in subset),
        "median_available_tpm": median([r["expression_tpm"] for r in with_tpm]),
        "rna_alt_assessed": len(with_alt),
        "rna_alt_supported_n": sum(r["rna_alt_status"] == "Supported" for r in with_alt),
        "expressed_and_binding_500_n": sum(r["expression_status"] == "Expressed" and r["binding_500"] == "Yes" for r in subset)
    })
patient_rows.sort(key=lambda r: (-r["expressed_n"], r["patient_id"]))
write_table("02_Neoantigen_Expression_Patient_Summary.tsv", patient_rows, list(patient_rows[0]))

assessed = [r for r in records if r["expression_tpm"] is not None]
expressed = [r for r in assessed if r["expression_status"] == "Expressed"]
transcript_assessed = sum(r["expression_source"] == "Transcript" for r in records)
gene_proxy = sum(r["expression_source"] == "Gene proxy" for r in records)
rna_alt_assessed = sum(r["rna_alt"] is not None for r in records)
cohort = [
    ("Patients", 80, "", ""),
    ("Distinct mutation–peptide–HLA candidates", len(records), "", ""),
    ("Candidates with expression TPM", len(assessed), len(records), percent(len(assessed), len(records))),
    ("Candidates with transcript TPM", transcript_assessed, len(records), percent(transcript_assessed, len(records))),
    ("Candidates with gene-level TPM proxy", gene_proxy, len(records), percent(gene_proxy, len(records))),
    ("Expression-positive candidates (TPM ≥1)", len(expressed), len(assessed), percent(len(expressed), len(assessed))),
    ("Below expression threshold (TPM <1)", len(assessed) - len(expressed), len(assessed), percent(len(assessed) - len(expressed), len(assessed))),
    ("Expression unavailable", len(records) - len(assessed), len(records), percent(len(records) - len(assessed), len(records))),
    ("Expression-positive and IC50 ≤500 nM", sum(r["binding_500"] == "Yes" for r in expressed), len(expressed), percent(sum(r["binding_500"] == "Yes" for r in expressed), len(expressed))),
    ("Candidates with RNA alternate-read count", rna_alt_assessed, len(records), percent(rna_alt_assessed, len(records))),
    ("Candidates with RNA alternate reads ≥3", sum(r["rna_alt_status"] == "Supported" for r in records), rna_alt_assessed, percent(sum(r["rna_alt_status"] == "Supported" for r in records), rna_alt_assessed)),
    ("Patients with any expression measurements", sum(r["tpm_assessed"] > 0 for r in patient_rows), 80, percent(sum(r["tpm_assessed"] > 0 for r in patient_rows), 80)),
    ("Median TPM among assessed candidates", median([r["expression_tpm"] for r in assessed]), len(assessed), "")
]
cohort_rows = [{"metric": metric, "value": display(value), "denominator": display(denominator) if denominator != "" else "", "percentage": display(percentage) if percentage != "" else ""} for metric, value, denominator, percentage in cohort]
write_table("01_Neoantigen_Expression_Cohort_Summary.tsv", cohort_rows, ["metric", "value", "denominator", "percentage"])

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(7.5, 4.5))
values = np.log10(1 + np.array([r["expression_tpm"] for r in assessed]))
ax.hist(values, bins=30, color=plt.cm.viridis(.65), edgecolor="black", linewidth=.5)
ax.axvline(np.log10(1 + TPM_THRESHOLD), color="black", linestyle="--", linewidth=1, label="TPM = 1")
ax.set(xlabel="log10(1 + expression TPM)", ylabel="Candidate combinations", title="Available neoantigen-associated RNA expression")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "01_Neoantigen_Expression_Distribution.png", dpi=300)
plt.close(fig)

fig, ax = plt.subplots(figsize=(11.5, 4.6))
ordered = sorted(patient_rows, key=lambda r: (-r["expressed_pct_assessed"] if r["expressed_pct_assessed"] is not None else float("inf"), r["patient_id"]))
x = np.arange(len(ordered))
positive_pc = np.array([100 * r["expressed_n"] / r["candidates_total"] for r in ordered])
below_pc = np.array([100 * (r["tpm_assessed"] - r["expressed_n"]) / r["candidates_total"] for r in ordered])
missing_pc = np.array([100 * r["unavailable_n"] / r["candidates_total"] for r in ordered])
for heights, bottom, color, label in [
    (positive_pc, np.zeros(len(ordered)), plt.cm.viridis(.78), "Expressed (TPM ≥1)"),
    (below_pc, positive_pc, plt.cm.viridis(.38), "Measured TPM <1"),
    (missing_pc, positive_pc + below_pc, plt.cm.viridis(.12), "TPM unavailable")
]:
    ax.bar(x, heights, bottom=bottom, width=.9, color=color, edgecolor="black", linewidth=.3, label=label)
ax.set(xlabel="Patients ranked by expression-positive percentage", ylabel="Percentage of all candidates (%)", title="RNA expression classification by patient")
ax.set_ylim(0, 103)
ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(.5, -.15))
fig.tight_layout()
fig.savefig(OUT / "02_Neoantigen_Expression_Patient_Percentages.png", dpi=300)
plt.close(fig)

fig, ax = plt.subplots(figsize=(11.5, 4.5))
ordered = sorted(patient_rows, key=lambda r: (-r["expressed_n"], r["patient_id"]))
ax.bar(np.arange(len(ordered)), [r["expressed_n"] for r in ordered], width=.85, color=plt.cm.viridis(.65), edgecolor="black", linewidth=.4)
ax.set(xlabel="Patients ranked by expression-positive candidate burden", ylabel="Candidates with TPM ≥1", title="Expression-supported candidate burden by patient")
fig.tight_layout()
fig.savefig(OUT / "03_Neoantigen_Expression_Patient_Burden.png", dpi=300)
plt.close(fig)

make_doc(cohort_rows, patient_rows, source_rows)
with (OUT / "06_Analysis_Metadata.json").open("w") as f:
    json.dump({"patients": 80, "candidate_rows": raw_total, "distinct_candidates": len(records), "expression_tpm_threshold": TPM_THRESHOLD, "rna_alt_read_threshold": RNA_ALT_THRESHOLD, "binding_ic50_threshold_nm": BINDING_THRESHOLD, "expression_assessed": len(assessed), "expression_positive": len(expressed), "transcript_tpm_records": transcript_assessed, "gene_proxy_tpm_records": gene_proxy, "rna_alt_count_assessed": rna_alt_assessed, "note": "TPM reflects RNA expression of a transcript or gene, not necessarily expression of the mutant allele. RNA mutant allele support is reported separately when available. Missing expression values are not zeros."}, f, indent=2)
print("Patients:", len(patient_rows))
print("Distinct candidates:", len(records))
print("Expression TPM assessed:", len(assessed))
print("Transcript TPM:", transcript_assessed)
print("Gene-level TPM proxy:", gene_proxy)
print("Expressed (TPM ≥1):", len(expressed))
print("RNA ALT read counts assessed:", rna_alt_assessed)
print("Validation: PASS")
print("Outputs:", OUT)
