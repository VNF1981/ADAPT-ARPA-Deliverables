import csv
from collections import Counter
from pathlib import Path

base = Path("/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/05_Downstream_Analysis/02_80_patients")
root = base / "03_results/07_HLA_LOH_Production"
manifest = root / "production_manifest.tsv"

with manifest.open() as handle:
    patients = [row["patient_id"] for row in csv.DictReader(handle, delimiter="\t")]

records = []
missing = []
malformed = []
expected = {"A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"}

for patient in patients:
    path = root / patient / "03_LOH/merge.hla.copy.txt"

    if not path.is_file() or path.stat().st_size == 0:
        missing.append(patient)
        continue

    with path.open() as handle:
        lines = [line.split() for line in handle if line.strip()]

    if not lines:
        malformed.append((patient, "empty"))
        continue

    header = lines[0]
    required = {"Sample", "HLA", "Allele1", "Allele2", "copyratio",
                "KeptHLA", "LossHLA", "Freq1", "Freq2",
                "Purity", "Het_num", "LOH"}

    if not required.issubset(header):
        malformed.append((patient, "invalid header"))
        continue

    rows = []
    valid = True

    for line in lines[1:]:
        if len(line) != len(header):
            valid = False
            break
        row = dict(zip(header, line))
        if row["LOH"] not in {"Y", "N"}:
            valid = False
            break
        row["patient_id"] = patient
        row.pop("Sample", None)
        rows.append(row)

    if not valid or set(r["HLA"] for r in rows) != expected or len(rows) != len(expected):
        malformed.append((patient, "incomplete or malformed loci"))
        continue

    records.extend(rows)

fields = ["patient_id", "HLA", "Allele1", "Allele2", "copyratio",
          "KeptHLA", "LossHLA", "Freq1", "Freq2", "Purity", "Het_num", "LOH"]

summary = root / "13_HLA_LOH_All_Loci.tsv"
with summary.open("w") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(records)

loh = [r for r in records if r["LOH"] == "Y"]
class_i = [r for r in records if r["HLA"] in {"A", "B", "C"}]
class_i_loh = [r for r in class_i if r["LOH"] == "Y"]
homozygous = [r for r in class_i if r["Allele1"] == r["Allele2"]]
low_het = [r for r in class_i if int(r["Het_num"]) < 5]
patients_loh = set(r["patient_id"] for r in loh)
patients_class_i_loh = set(r["patient_id"] for r in class_i_loh)

print("===== COHORT COMPLETENESS =====")
print("Expected patients:", len(patients))
print("Patients with complete results:", len(set(r["patient_id"] for r in records)))
print("Missing results:", len(missing))
print("Malformed/incomplete results:", len(malformed))
print("Total valid locus records:", len(records))

print("\n===== HLA-LOH SUMMARY =====")
print("Patients with LOH at any locus:", len(patients_loh))
print("Patients with HLA-I LOH (A/B/C):", len(patients_class_i_loh))
print("Total LOH-positive loci:", len(loh))
print("HLA-I LOH-positive loci:", len(class_i_loh))

print("\n===== LOH BY LOCUS =====")
for locus in ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]:
    subset = [r for r in records if r["HLA"] == locus]
    positive = sum(r["LOH"] == "Y" for r in subset)
    print(f"{locus:5s}: {positive:3d} / {len(subset):3d}")

print("\n===== HLA-I QUALITY INDICATORS =====")
print("Homozygous locus calls:", len(homozygous))
print("Locus calls with fewer than 5 heterozygous variants:", len(low_het))
print("These categories can overlap.")

print("\n===== HLA-I LOH-POSITIVE RECORDS =====")
for r in class_i_loh:
    print(f'{r["patient_id"]}\tHLA-{r["HLA"]}\t'
          f'Allele1={r["Allele1"]}\tAllele2={r["Allele2"]}\t'
          f'LossHLA={r["LossHLA"]}\t'
          f'Freq1={r["Freq1"]}\tFreq2={r["Freq2"]}\t'
          f'Het_num={r["Het_num"]}')

if missing:
    print("\nMISSING:", ", ".join(missing))
if malformed:
    print("\nMALFORMED:", malformed)

with (root / "13_HLA_LOH_Positive_Loci.tsv").open("w") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()
    writer.writerows(loh)

print("\nOutput:", summary)
print("Output:", root / "13_HLA_LOH_Positive_Loci.tsv")

