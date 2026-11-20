import csv
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from docx import Document
from docx.shared import Inches, Pt, RGBColor

BASE = Path(os.environ.get("BRCA80_BASE", Path(__file__).resolve().parents[2]))
DATA = BASE / "02_data/02_pVACseq"
OUT = BASE / "03_results/04_Neoantigen_Expression"
ST = OUT / "00_StringTie_Quantification"
SUMMARY = DATA / "03_TCGA_BRCA_Neoantigen_Patient_Summary.tsv"
MANIFEST = BASE / "02_data/01_TCGA_BRCA_80_Input_Manifest.tsv"
TPM_THRESHOLD = 1.0
BINDING_THRESHOLD = 500.0
RNA_ALT_THRESHOLD = 3
OUT.mkdir(parents=True, exist_ok=True)

EXP_COLS = {
    "rna_vaf": ["Tumor RNA VAF", "Tumor RNA Variant Allele Frequency", "tumor_rna_vaf"],
    "rna_depth": ["Tumor RNA Depth", "Tumor RNA Coverage", "tumor_rna_depth"],
    "rna_alt": ["Tumor RNA Alt Read Support", "Tumor RNA Alt Reads", "Tumor RNA ALT Reads", "tumor_rna_alt_reads"],
    "ic50": ["Best MT IC50 Score", "Median MT IC50 Score", "MT IC50 Score"]
}

def choose(columns, aliases):
    normalized = {"".join(x.lower().split()): x for x in columns}
    for alias in aliases:
        if alias in columns:
            return alias
        found = normalized.get("".join(alias.lower().split()))
        if found:
            return found
    return None

def numeric(value):
    try:
        result = float(value)
        return result if math.isfinite(result) and result >= 0 else None
    except (TypeError, ValueError):
        return None

def read_table(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_table(name, rows, fields):
    with (OUT / name).open("w", newline="") as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

def pct(a, b):
    return round(100 * a / b, 2) if b else None

def clean_id(value):
    return str(value).split(".")[0].strip()

def ids(value, prefix):
    return set(re.findall(r"\b" + prefix + r"\d+(?:\.\d+)?\b", str(value)))

def gtf_attributes(text):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))

def load_stringtie(pid, wanted_transcripts):
    folder = ST / pid
    gene_file = folder / f"{pid}.gene_abundance.tsv"
    transcript_file = folder / f"{pid}.transcripts.gtf"
    done_file = folder / f"{pid}.complete"
    if not all(f.is_file() and f.stat().st_size > 0 for f in (gene_file, transcript_file, done_file)):
        raise RuntimeError(f"{pid}: incomplete StringTie outputs")
    gene_by_id = {}
    gene_names = defaultdict(list)
    with gene_file.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or not {"Gene ID", "Gene Name", "TPM"}.issubset(reader.fieldnames):
            raise RuntimeError(f"{pid}: unexpected gene TPM header: {reader.fieldnames}")
        for row in reader:
            tpm = numeric(row["TPM"])
            if tpm is None:
                raise RuntimeError(f"{pid}: invalid gene TPM")
            gid = clean_id(row["Gene ID"])
            gene_by_id[gid] = tpm
            gene_names[row["Gene Name"]].append((gid, tpm))
    wanted = {clean_id(x) for x in wanted_transcripts}
    transcripts = {}
    with transcript_file.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t", 8)
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            attrs = gtf_attributes(parts[8])
            tid = clean_id(attrs.get("transcript_id", ""))
            if tid not in wanted:
                continue
            tpm = numeric(attrs.get("TPM"))
            if tpm is None:
                raise RuntimeError(f"{pid}: invalid transcript TPM for {tid}")
            transcripts[tid] = max(tpm, transcripts.get(tid, 0))
    return transcripts, gene_by_id, gene_names

def get_expression(record, transcript_tpm, gene_tpm, gene_names):
    matches = [(t, transcript_tpm[clean_id(t)]) for t in record["transcript_ids"] if clean_id(t) in transcript_tpm]
    if matches:
        tid, tpm = max(matches, key=lambda item: item[1])
        return tpm, "Transcript", tid
    matches = [(g, gene_tpm[clean_id(g)]) for g in record["gene_ids"] if clean_id(g) in gene_tpm]
    if matches:
        gid, tpm = max(matches, key=lambda item: item[1])
        return tpm, "Gene proxy", gid
    matches = []
    for symbol in record["gene_names"]:
        linked = gene_names.get(symbol, [])
        if len({g for g, _ in linked}) == 1:
            matches.append((symbol, linked[0][1]))
    if matches:
        symbol, tpm = max(matches, key=lambda item: item[1])
        return tpm, "Gene proxy", symbol
    return None, "Unavailable", ""

def fmt(x):
    return "NA" if x is None else f"{x:.2f}" if isinstance(x, float) else str(x)

def make_doc(cohort, patients, audit):
    doc = Document()
    section = doc.sections[0]
    section.top_margin = Inches(.65)
    section.bottom_margin = Inches(.65)
    section.left_margin = Inches(.7)
    section.right_margin = Inches(.7)
    for name in ("Normal", "Title", "Heading 1", "Heading 2"):
        style = doc.styles[name]
        style.font.name = "Times New Roman"
        style.font.size = Pt(9)
        style.font.color.rgb = RGBColor(0, 0, 0)
    doc.add_paragraph("Analysis 4 — Neoantigen Expression")
    doc.add_paragraph("StringTie 3.0.3 quantification of tumor genomic RNA-seq alignments using GENCODE v47. Transcript TPM is used when a pVACseq transcript matches; otherwise gene TPM is reported as a proxy. When multiple candidate-associated transcripts match, the maximum TPM is used. Expression-positive means TPM ≥1. Transcript or gene expression does not establish mutant-allele expression. RNA alternate-read evidence, if available, is presented separately.")
    def table(title, headers, rows):
        doc.add_paragraph(title)
        tab = doc.add_table(rows=1, cols=len(headers))
        tab.style = "Table Grid"
        for cell, value in zip(tab.rows[0].cells, headers):
            cell.text = value
        for values in rows:
            for cell, value in zip(tab.add_row().cells, values):
                cell.text = fmt(value)
        for row in tab.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    for run in paragraph.runs:
                        run.font.name = "Times New Roman"
                        run.font.size = Pt(8)
                        run.font.color.rgb = RGBColor(0, 0, 0)
    table("Table 1. Cohort summary", ["Metric", "Value", "Denominator", "Percentage"], [[r[x] for x in ("metric", "value", "denominator", "percentage")] for r in cohort])
    columns = ["patient_id", "subtype", "candidates_total", "tpm_assessed", "expressed_n", "expressed_pct_assessed", "unavailable_n", "median_available_tpm", "rna_alt_assessed", "rna_alt_supported_n", "expressed_and_binding_500_n"]
    table("Table 2. Patient-level summary", ["Patient", "Subtype", "Candidates", "TPM assessed", "Expressed", "Expressed %", "Missing", "Median TPM", "RNA ALT assessed", "RNA ALT ≥3", "Expressed & IC50 ≤500"], [[r[x] for x in columns] for r in patients])
    table("Table 3. Expression matching audit", ["Source", "Count"], [[k, v] for k, v in audit.items()])
    doc.save(OUT / "04_Neoantigen_Expression_Tables.docx")

summary = read_table(SUMMARY)
manifest = read_table(MANIFEST)
patient_ids = {r["patient_id"] for r in summary}
if len(summary) != 80 or len(patient_ids) != 80:
    raise RuntimeError("Expected 80 unique patients in pVACseq summary")
if len(manifest) != 80 or {r["patient_id"] for r in manifest} != patient_ids:
    raise RuntimeError("80-patient input manifest mismatch")
subtypes = {r["patient_id"]: r.get("subtype", "") for r in manifest}
records = []
raw_total = 0
for patient in summary:
    pid = patient["patient_id"]
    files = list((DATA / pid / "MHC_Class_I").glob("*.MHC_I.filtered.tsv"))
    if len(files) != 1:
        raise RuntimeError(f"{pid}: expected one pVACseq file, got {len(files)}")
    candidates = {}
    with files[0].open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        columns = reader.fieldnames or []
        peptide_col = choose(columns, ["MT Epitope Seq", "MT Epitope"])
        hla_col = choose(columns, ["HLA Allele", "HLA"])
        transcript_col = choose(columns, ["Transcript", "Transcript ID", "Transcript Id"])
        gene_id_col = choose(columns, ["Gene ID", "Ensembl Gene ID", "Ensembl Gene Id"])
        gene_name_col = choose(columns, ["Gene Name", "Gene Symbol", "Hugo Symbol"])
        fields = ["Chromosome", "Start", "Stop", "Reference", "Variant"]
        if not peptide_col or not hla_col or any(c not in columns for c in fields):
            raise RuntimeError(f"{pid}: candidate identifying columns missing")
        chosen = {key: choose(columns, vals) for key, vals in EXP_COLS.items()}
        count = 0
        for row in reader:
            count += 1
            mutation = tuple(str(row[x]).strip() for x in fields)
            peptide = row[peptide_col].strip()
            hla = row[hla_col].strip()
            if not peptide or not hla:
                raise RuntimeError(f"{pid}: empty peptide/HLA")
            key = (mutation, peptide, hla)
            measure = {k: numeric(row[alias]) if alias else None for k, alias in chosen.items()}
            tids = ids(row.get(transcript_col, ""), "ENST")
            gids = ids(row.get(gene_id_col, ""), "ENSG")
            gene_label = str(row.get(gene_name_col, "")).strip() if gene_name_col else ""
            gids.update(ids(gene_label, "ENSG"))
            symbols = {gene_label} if gene_label and not gene_label.startswith("ENSG") and gene_label not in {"NA", "."} else set()
            if key not in candidates:
                candidates[key] = dict(patient_id=pid, subtype=subtypes[pid], chromosome=mutation[0], start=mutation[1], stop=mutation[2], reference=mutation[3], variant=mutation[4], mt_epitope=peptide, hla_allele=hla, transcript_ids=set(), gene_ids=set(), gene_names=set(), **measure)
            rec = candidates[key]
            rec["transcript_ids"].update(tids)
            rec["gene_ids"].update(gids)
            rec["gene_names"].update(symbols)
            for k, value in measure.items():
                if value is not None and (rec[k] is None or (k == "ic50" and value < rec[k])):
                    rec[k] = value
    if count != int(patient["neo_filtered_rows"]) or len(candidates) != int(patient["neo_peptide_hla_count"]):
        raise RuntimeError(f"{pid}: pVACseq candidate counts differ from validated summary")
    raw_total += count
    records.extend(candidates.values())
if len(records) != 5533 or raw_total != 5596:
    raise RuntimeError(f"Candidate total mismatch: {len(records)} distinct, {raw_total} raw")

by_patient = defaultdict(list)
for rec in records:
    by_patient[rec["patient_id"]].append(rec)
for pid, candidates in by_patient.items():
    wanted = set().union(*(r["transcript_ids"] for r in candidates))
    transcripts, genes, symbols = load_stringtie(pid, wanted)
    for rec in candidates:
        tpm, source, match = get_expression(rec, transcripts, genes, symbols)
        rec["expression_tpm"] = tpm
        rec["expression_source"] = source
        rec["matched_identifier"] = match
        rec["transcript_tpm"] = tpm if source == "Transcript" else None
        rec["gene_tpm"] = tpm if source == "Gene proxy" else None
        rec["expression_status"] = "Not assessed" if tpm is None else "Expressed" if tpm >= TPM_THRESHOLD else "Below threshold"
        rec["rna_alt_status"] = "Not assessed" if rec["rna_alt"] is None else "Supported" if rec["rna_alt"] >= RNA_ALT_THRESHOLD else "Below threshold"
        rec["binding_500"] = "NA" if rec["ic50"] is None else "Yes" if rec["ic50"] <= BINDING_THRESHOLD else "No"
        rec["transcript_ids"] = ";".join(sorted(rec["transcript_ids"]))
        rec["gene_ids"] = ";".join(sorted(rec["gene_ids"]))
        rec["gene_names"] = ";".join(sorted(rec["gene_names"]))

cols = ["patient_id", "subtype", "chromosome", "start", "stop", "reference", "variant", "mt_epitope", "hla_allele", "transcript_ids", "gene_ids", "gene_names", "matched_identifier", "transcript_tpm", "gene_tpm", "expression_tpm", "expression_source", "expression_status", "rna_vaf", "rna_depth", "rna_alt", "rna_alt_status", "ic50", "binding_500"]
write_table("03_Neoantigen_Expression_Candidates.tsv", records, cols)
audit = Counter(r["expression_source"] for r in records)
audit["Patients with TPM files"] = len(by_patient)
audit["Candidates with transcript annotations"] = sum(bool(r["transcript_ids"]) for r in records)
audit["Candidates with gene annotations"] = sum(bool(r["gene_ids"] or r["gene_names"]) for r in records)
write_table("05_Expression_Field_Audit.tsv", [{"field": k, "value": v} for k, v in audit.items()], ["field", "value"])
if not any(r["expression_tpm"] is not None for r in records):
    raise RuntimeError("No expression matches: inspect pVACseq transcript/gene identifiers")

patient_rows = []
for pid in sorted(by_patient):
    sub = by_patient[pid]
    assessed = [r for r in sub if r["expression_tpm"] is not None]
    expressed = [r for r in assessed if r["expression_tpm"] >= TPM_THRESHOLD]
    alt = [r for r in sub if r["rna_alt"] is not None]
    patient_rows.append(dict(patient_id=pid, subtype=subtypes[pid], candidates_total=len(sub), tpm_assessed=len(assessed), expressed_n=len(expressed), expressed_pct_assessed=pct(len(expressed), len(assessed)), unavailable_n=len(sub)-len(assessed), transcript_tpm_assessed=sum(r["expression_source"] == "Transcript" for r in sub), gene_proxy_assessed=sum(r["expression_source"] == "Gene proxy" for r in sub), median_available_tpm=median([r["expression_tpm"] for r in assessed]) if assessed else None, rna_alt_assessed=len(alt), rna_alt_supported_n=sum(r["rna_alt"] >= RNA_ALT_THRESHOLD for r in alt), expressed_and_binding_500_n=sum(r["binding_500"] == "Yes" for r in expressed)))
patient_rows.sort(key=lambda r: (-r["expressed_n"], r["patient_id"]))
write_table("02_Neoantigen_Expression_Patient_Summary.tsv", patient_rows, list(patient_rows[0]))
assessed = [r for r in records if r["expression_tpm"] is not None]
expressed = [r for r in assessed if r["expression_tpm"] >= TPM_THRESHOLD]
values = [
    ("Patients", 80, "", ""),
    ("Distinct mutation–peptide–HLA candidates", len(records), "", ""),
    ("Candidates with expression TPM", len(assessed), len(records), pct(len(assessed), len(records))),
    ("Candidates with transcript TPM", audit["Transcript"], len(records), pct(audit["Transcript"], len(records))),
    ("Candidates with gene-level TPM proxy", audit["Gene proxy"], len(records), pct(audit["Gene proxy"], len(records))),
    ("Expression-positive candidates (TPM ≥1)", len(expressed), len(assessed), pct(len(expressed), len(assessed))),
    ("Below expression threshold (TPM <1)", len(assessed)-len(expressed), len(assessed), pct(len(assessed)-len(expressed), len(assessed))),
    ("Expression unavailable", len(records)-len(assessed), len(records), pct(len(records)-len(assessed), len(records))),
    ("Expression-positive and IC50 ≤500 nM", sum(r["binding_500"] == "Yes" for r in expressed), len(expressed), pct(sum(r["binding_500"] == "Yes" for r in expressed), len(expressed))),
    ("Candidates with RNA alternate-read count", sum(r["rna_alt"] is not None for r in records), len(records), pct(sum(r["rna_alt"] is not None for r in records), len(records))),
    ("Patients with expression measurements", sum(r["tpm_assessed"] > 0 for r in patient_rows), 80, pct(sum(r["tpm_assessed"] > 0 for r in patient_rows), 80)),
    ("Median TPM among assessed candidates", median([r["expression_tpm"] for r in assessed]), len(assessed), "")
]
cohort = [dict(metric=a, value=fmt(b), denominator=fmt(c) if c != "" else "", percentage=fmt(d) if d != "" else "") for a,b,c,d in values]
write_table("01_Neoantigen_Expression_Cohort_Summary.tsv", cohort, ["metric", "value", "denominator", "percentage"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(7.5,4.5))
ax.hist(np.log10(1+np.array([r["expression_tpm"] for r in assessed])), bins=30, color=plt.cm.viridis(.65), edgecolor="black", linewidth=.5)
ax.axvline(np.log10(2), color="black", linestyle="--", linewidth=1, label="TPM = 1")
ax.set(xlabel="log10(1 + expression TPM)", ylabel="Candidate combinations", title="Neoantigen-associated RNA expression")
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(OUT / "01_Neoantigen_Expression_Distribution.png", dpi=300)
plt.close(fig)
ordered = sorted(patient_rows, key=lambda r: (-r["expressed_n"] / r["tpm_assessed"] if r["tpm_assessed"] else float("inf"),r["patient_id"]))
x = np.arange(len(ordered))
total = np.array([r["candidates_total"] for r in ordered])
pos = 100*np.array([r["expressed_n"] for r in ordered])/total
below = 100*np.array([r["tpm_assessed"]-r["expressed_n"] for r in ordered])/total
missing = 100*np.array([r["unavailable_n"] for r in ordered])/total
fig, ax = plt.subplots(figsize=(max(16,len(ordered)*.23),6))
for heights, bottom, color, label in [(pos,np.zeros(len(x)),plt.cm.viridis(.78),"TPM ≥1"),(below,pos,plt.cm.viridis(.38),"TPM <1"),(missing,pos+below,plt.cm.viridis(.12),"Unavailable")]:
    ax.bar(x,heights,bottom=bottom,width=.88,color=color,edgecolor="black",linewidth=.3,label=label)
ax.set(ylabel="Percentage of all candidates (%)",title="RNA expression classification by patient",ylim=(0,105))
ax.set_xticks(x)
ax.set_xticklabels([r["patient_id"] for r in ordered],rotation=90,fontsize=7)
ax.legend(frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.5,1.10))
fig.subplots_adjust(bottom=.32,top=.87,left=.055,right=.99)
fig.savefig(OUT / "02_Neoantigen_Expression_Patient_Percentages.png",dpi=300)
plt.close(fig)
ordered = sorted(patient_rows,key=lambda r:(-r["expressed_n"],r["patient_id"]))
x = np.arange(len(ordered))
fig, ax = plt.subplots(figsize=(max(16,len(ordered)*.23),6))
ax.bar(x,[r["expressed_n"] for r in ordered],width=.85,color=plt.cm.viridis(.65),edgecolor="black",linewidth=.4)
ax.set(ylabel="Candidates with TPM ≥1",title="Expression-supported candidate burden by patient")
ax.set_xticks(x)
ax.set_xticklabels([r["patient_id"] for r in ordered],rotation=90,fontsize=7)
fig.subplots_adjust(bottom=.32,top=.91,left=.055,right=.99)
fig.savefig(OUT / "03_Neoantigen_Expression_Patient_Burden.png",dpi=300)
plt.close(fig)
make_doc(cohort,patient_rows,audit)
with (OUT / "06_Analysis_Metadata.json").open("w") as f:
    json.dump(dict(patients=80,raw_rows=raw_total,distinct_candidates=len(records),expression_tpm_threshold=TPM_THRESHOLD,binding_threshold_nm=BINDING_THRESHOLD,expression_assessed=len(assessed),expression_positive=len(expressed),transcript_tpm_records=audit["Transcript"],gene_proxy_tpm_records=audit["Gene proxy"],unavailable_tpm_records=audit["Unavailable"],expression_source="StringTie 3.0.3; GENCODE v47; tumor genomic RNA-seq BAM",caution="Transcript or gene TPM does not establish expression of the mutant allele. Highest TPM of matched candidate transcripts used."),f,indent=2)
print("Patients:",len(patient_rows))
print("Unique mutation-peptide-HLA candidates:",len(records))
print("TPM assessed:",len(assessed))
print("Transcript TPM matched:",audit["Transcript"])
print("Gene TPM proxy:",audit["Gene proxy"])
print("TPM unavailable:",audit["Unavailable"])
print("TPM >=1:",len(expressed))
print("Validation: PASS")
print("Outputs:",OUT)
