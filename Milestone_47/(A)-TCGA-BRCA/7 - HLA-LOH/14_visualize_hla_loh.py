import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from docx import Document
from docx.shared import Inches, Pt, RGBColor

LOCI = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
CLASS_I = {"A", "B", "C"}
NPG_BLUE = "#3C5488"
NPG_TEAL = "#00A087"

def add_table(doc, heading, columns, rows):
    doc.add_paragraph(heading)
    table = doc.add_table(rows=1, cols=len(columns))
    table.style = "Table Grid"
    for cell, value in zip(table.rows[0].cells, columns):
        cell.text = str(value)
    for values in rows:
        cells = table.add_row().cells
        for cell, value in zip(cells, values):
            cell.text = str(value)
    for row in table.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                for run in p.runs:
                    run.font.name = "Times New Roman"
                    run.font.size = Pt(10)
                    run.font.color.rgb = RGBColor(0, 0, 0)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cohort-size", type=int, default=80)
    args = parser.parse_args()

    output = Path(args.output).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input, sep="\t", dtype={"patient_id": str, "HLA": str, "LOH": str})
    required = {"patient_id", "HLA", "LOH"}
    if not required.issubset(df.columns):
        raise ValueError(f"Missing required columns: {sorted(required - set(df.columns))}")
    if df.empty or df[["patient_id", "HLA", "LOH"]].isna().any().any():
        raise ValueError("Input has empty or missing patient/locus/call values")
    if not set(df.HLA).issubset(LOCI) or not set(df.LOH).issubset({"Y", "N"}):
        raise ValueError("Unexpected HLA locus or LOH value")
    if df.duplicated(["patient_id", "HLA"]).any():
        raise ValueError("Duplicated patient/locus combinations")
    patients = df.patient_id.nunique()
    if len(df) != patients * len(LOCI) or not df.groupby("patient_id").HLA.nunique().eq(len(LOCI)).all():
        raise ValueError("Incomplete patient-level HLA locus coverage")

    df["Class"] = np.where(df.HLA.isin(CLASS_I), "Class I", "Class II")
    df["LOH_positive"] = df.LOH.eq("Y")
    patient = df.groupby("patient_id", sort=True).agg(
        positive_loci=("LOH_positive", "sum"),
        positive_class_I=("LOH_positive", lambda s: int(s[df.loc[s.index, "Class"].eq("Class I")].sum())),
        positive_class_II=("LOH_positive", lambda s: int(s[df.loc[s.index, "Class"].eq("Class II")].sum()))
    ).reset_index()
    patient["any_LOH"] = patient.positive_loci.gt(0)
    patient["class_I_LOH"] = patient.positive_class_I.gt(0)
    patient["class_II_LOH"] = patient.positive_class_II.gt(0)
    patient.to_csv(output / "14_HLA_LOH_Patient_Summary.tsv", sep="\t", index=False)

    locus = df.groupby(["HLA", "Class"], sort=False).agg(
        patients=("patient_id", "nunique"),
        positive=("LOH_positive", "sum")
    ).reset_index()
    locus["HLA"] = pd.Categorical(locus.HLA, LOCI, ordered=True)
    locus = locus.sort_values("HLA").reset_index(drop=True)
    locus["percentage"] = 100 * locus.positive / locus.patients
    locus.to_csv(output / "14_HLA_LOH_Locus_Summary.tsv", sep="\t", index=False)

    n_any = int(patient.any_LOH.sum())
    n_i = int(patient.class_I_LOH.sum())
    n_ii = int(patient.class_II_LOH.sum())
    n_both = int((patient.class_I_LOH & patient.class_II_LOH).sum())
    n_i_only = int((patient.class_I_LOH & ~patient.class_II_LOH).sum())
    n_ii_only = int((~patient.class_I_LOH & patient.class_II_LOH).sum())
    n_none = int((~patient.any_LOH).sum())
    class_i_rows = df[df.Class.eq("Class I")]
    class_ii_rows = df[df.Class.eq("Class II")]
    n_i_records = int(class_i_rows.LOH_positive.sum())
    n_ii_records = int(class_ii_rows.LOH_positive.sum())

    overview = [
        ("Patients in selected cohort", args.cohort_size, ""),
        ("Patients with complete SpecHLA results", patients, f"{100 * patients / args.cohort_size:.1f}%"),
        ("Evaluated HLA locus records", len(df), ""),
        ("Patients LOH-positive at any locus", n_any, f"{100 * n_any / patients:.1f}%"),
        ("Patients LOH-positive at class I loci", n_i, f"{100 * n_i / patients:.1f}%"),
        ("Patients LOH-positive at class II loci", n_ii, f"{100 * n_ii / patients:.1f}%"),
        ("Patients LOH-positive in both classes", n_both, f"{100 * n_both / patients:.1f}%"),
        ("Patients LOH-positive in class I only", n_i_only, f"{100 * n_i_only / patients:.1f}%"),
        ("Patients LOH-positive in class II only", n_ii_only, f"{100 * n_ii_only / patients:.1f}%"),
        ("Patients without LOH-positive calls", n_none, f"{100 * n_none / patients:.1f}%"),
        ("Class I LOH-positive locus records", n_i_records, f"{100 * n_i_records / len(class_i_rows):.1f}%"),
        ("Class II LOH-positive locus records", n_ii_records, f"{100 * n_ii_records / len(class_ii_rows):.1f}%")
    ]
    pd.DataFrame(overview, columns=["Measure", "Count", "Percentage"]).to_csv(
        output / "14_HLA_LOH_Cohort_Summary.tsv", sep="\t", index=False
    )

    doc = Document()
    section = doc.sections[0]
    section.top_margin = Inches(0.7)
    section.bottom_margin = Inches(0.7)
    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(10)
    normal.font.color.rgb = RGBColor(0, 0, 0)
    doc.add_paragraph("Table 7. SpecHLA-reported HLA loss-of-heterozygosity classifications across TCGA-BRCA patients.")
    add_table(doc, "Cohort summary", ["Measure", "Count", "Percentage"], overview)
    add_table(
        doc,
        "Locus-specific summary",
        ["HLA locus", "Class", "Evaluated patients", "LOH-positive calls", "Percentage"],
        [(f"HLA-{r.HLA}", r.Class, r.patients, r.positive, f"{r.percentage:.1f}%") for r in locus.itertuples(index=False)]
    )
    doc.add_paragraph(
        "Note. Percentages for patient-level outcomes use the 79 analyzed patients; locus-level percentages use the "
        "number of patients evaluated at each locus. A patient can have LOH-positive calls at multiple loci. "
        "LOH-positive denotes the software-reported classification."
    )
    doc.save(output / "Table_7_HLA_LOH.docx")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "savefig.dpi": 300,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(8.2, 4.9))
    xpos = np.arange(len(locus))
    colors = [NPG_BLUE if cls == "Class I" else NPG_TEAL for cls in locus["Class"]]
    bars = ax.bar(xpos, locus.percentage, color=colors, width=0.72)
    for rect, r in zip(bars, locus.itertuples(index=False)):
        ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + 1.1,
                f"{int(r.positive)}/{int(r.patients)}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(xpos)
    ax.set_xticklabels([f"HLA-{x}" for x in locus.HLA])
    ax.set_ylabel("Patients with SpecHLA LOH-positive call (%)")
    ax.set_ylim(0, max(65, float(locus.percentage.max()) + 11))
    ax.legend(handles=[Patch(facecolor=NPG_BLUE, label="Class I"), Patch(facecolor=NPG_TEAL, label="Class II")],
              frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(output / "Figure_9_HLA_LOH.png", dpi=300)
    fig.savefig(output / "Figure_9_HLA_LOH.pdf")
    plt.close(fig)
    print(f"Complete patients: {patients}; records: {len(df)}")
    print(f"Any-locus positive: {n_any}; class I: {n_i}; class II: {n_ii}")
    print(f"Saved: {output}")

if __name__ == "__main__":
    main()
