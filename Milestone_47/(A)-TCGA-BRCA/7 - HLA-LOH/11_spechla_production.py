import argparse
import csv
import math
import subprocess
import shutil
from pathlib import Path

BASE = Path("/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/05_Downstream_Analysis/02_80_patients")
MASTER = Path("/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/01_data_audit/02_results/47_BRCA_FINAL_MASTER_MANIFEST_PAM50_TNBC_ABSOLUTE_EVOFLUx.tsv")
SUBSET = BASE / "02_data/01_TCGA_BRCA_80_Input_Manifest.tsv"
OUT = BASE / "03_results/07_HLA_LOH_Production"
MANIFEST = OUT / "production_manifest.tsv"

def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def valid_number(value):
    try:
        number = float(value)
        return math.isfinite(number) and number > 0
    except (ValueError, TypeError):
        return False

def prepare():
    patients = [row["patient_id"].strip() for row in read_table(SUBSET)]
    if len(patients) != 80 or len(set(patients)) != 80:
        raise ValueError("Subset must contain 80 unique patients")

    master = {}
    for row in read_table(MASTER):
        patient = row["patient_id"].strip()
        if patient in patients:
            if patient in master:
                raise ValueError(f"Duplicate master record: {patient}")
            master[patient] = row

    if set(patients) != set(master):
        raise ValueError(f"Missing master patients: {set(patients) - set(master)}")

    eligible = []
    excluded = []

    for patient in patients:
        row = master[patient]
        purity = row["ABSOLUTE_purity"]
        ploidy = row["ABSOLUTE_ploidy"]

        if not valid_number(purity) or not valid_number(ploidy):
            excluded.append((patient, purity, ploidy))
            continue

        tumor = Path(row["tumor_wxs_bam"])
        normal = Path(row["normal_wxs_bam"])

        for path in (tumor, normal):
            if not path.is_file():
                raise FileNotFoundError(path)

        eligible.append({
            "patient_id": patient,
            "purity": purity,
            "ploidy": ploidy,
            "tumor_bam": str(tumor),
            "normal_bam": str(normal)
        })

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "logs").mkdir(exist_ok=True)

    if MANIFEST.exists():
        existing = read_table(MANIFEST)
        if existing != eligible:
            raise RuntimeError("Production manifest already exists with different contents")
    else:
        with MANIFEST.open("w") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(eligible[0]))
            writer.writeheader()
            writer.writerows(eligible)

    print(f"Selected patients: {len(patients)}")
    print(f"Eligible patients: {len(eligible)}")
    print(f"Excluded patients: {len(excluded)}")
    for patient, purity, ploidy in excluded:
        print(f"EXCLUDED: {patient}, purity={purity}, ploidy={ploidy}")
    print(f"Manifest: {MANIFEST}")

def run_command(command):
    print("COMMAND:", " ".join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), check=True)

def check_typing(directory):
    required = [
        directory / "hla.result.txt",
        directory / "HLA_A_freq.txt",
        directory / "HLA_B_freq.txt",
        directory / "HLA_C_freq.txt"
    ]
    return all(path.is_file() and path.stat().st_size > 0 for path in required)

def execute(index):
    rows = read_table(MANIFEST)
    if index < 1 or index > len(rows):
        raise ValueError(f"Invalid array index: {index}")

    row = rows[index - 1]
    patient = row["patient_id"]
    root = OUT / patient
    extraction = root / "01_extracted_reads"
    typing = root / "02_typing"
    loh = root / "03_LOH"

    print(f"PATIENT: {patient}", flush=True)

    for label in ("normal", "tumor"):
        sample = f"{patient}_{label}"
        bam = row[f"{label}_bam"]
        destination = extraction / label
        destination.mkdir(parents=True, exist_ok=True)

        prefix = destination / f"{sample}_extract"
        r1 = Path(str(prefix) + "_1.fq.gz")
        r2 = Path(str(prefix) + "_2.fq.gz")

        if not (r1.is_file() and r1.stat().st_size > 0 and r2.is_file() and r2.stat().st_size > 0):
            run_command([
                "spechla-extract-hla-reads",
                "-s", sample,
                "-b", bam,
                "-r", "hg38",
                "-o", destination
            ])

        if not (r1.is_file() and r1.stat().st_size > 0 and r2.is_file() and r2.stat().st_size > 0):
            raise RuntimeError(f"Missing extracted paired FASTQs: {sample}")

        result_directory = typing / label / sample
        if not check_typing(result_directory):
            if result_directory.exists():
                shutil.rmtree(result_directory)

            (typing / label).mkdir(parents=True, exist_ok=True)

            run_command([
                "spechla",
                "-n", sample,
                "-1", r1,
                "-2", r2,
                "-o", typing / label,
                "-u", "1",
                "-j", "4"
            ])

        if not check_typing(result_directory):
            raise RuntimeError(f"SpecHLA typing outputs missing: {sample}")

    tumor = typing / "tumor" / f"{patient}_tumor"
    normal = typing / "normal" / f"{patient}_normal"

    loh.mkdir(parents=True, exist_ok=True)
    freq_list = loh / "tumor_freq.list"
    frequencies = sorted(tumor.glob("HLA_*_freq.txt"))

    if not {"HLA_A", "HLA_B", "HLA_C"}.issubset({
        path.name.replace("_freq.txt", "") for path in frequencies
    }):
        raise RuntimeError(f"Required tumor HLA frequencies missing: {patient}")

    freq_list.write_text("\n".join(map(str, frequencies)) + "\n")
    result = loh / "merge.hla.copy.txt"

    if not result.is_file() or result.stat().st_size == 0:
        run_command([
            "spechla-loh",
            "-S", patient,
            "-C", "5",
            "-purity", row["purity"],
            "-ploidy", row["ploidy"],
            "-F", freq_list,
            "-T", normal / "hla.result.txt",
            "-O", loh
        ])

    if not result.is_file() or result.stat().st_size == 0:
        raise RuntimeError(f"Missing LOH result: {patient}")

    with result.open() as handle:
        lines = [line.split() for line in handle if line.strip()]

    required_columns = {"Sample", "HLA", "copyratio", "Het_num", "LOH"}
    header = lines[0] if lines else []

    if not required_columns.issubset(set(header)):
        raise RuntimeError(f"Invalid LOH result header: {patient}")

    locus_index = header.index("HLA")
    loh_index = header.index("LOH")
    loci = {line[locus_index] for line in lines[1:] if len(line) == len(header)}

    if not {"A", "B", "C"}.issubset(loci):
        raise RuntimeError(f"Incomplete HLA-A/B/C results: {patient}")

    for line in lines[1:]:
        if len(line) != len(header) or line[loh_index] not in {"Y", "N"}:
            raise RuntimeError(f"Malformed LOH result: {patient}")

    print(f"COMPLETED: {patient}", flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--index", type=int)
    args = parser.parse_args()

    if args.prepare:
        prepare()
    elif args.index is not None:
        execute(args.index)
    else:
        parser.error("Provide --prepare or --index")

