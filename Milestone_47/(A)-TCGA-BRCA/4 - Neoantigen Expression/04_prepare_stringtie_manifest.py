#!/usr/bin/env python3
import csv
from pathlib import Path

BASE = Path('/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/05_Downstream_Analysis/02_80_patients')
INPUT = Path('/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/07_BRCA_NextNEOpi/02_subset_TCGA_Patients/08_BRCA_80_NextNEOpi_BAM_Batch.csv')
GTF = Path('/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/06_NextNEOpi/nextNEOpi/resources/references/hg38/annotation/gencode.v47.primary_assembly.annotation.gtf')
OUT = BASE / '03_results/04_Neoantigen_Expression'
OUT.mkdir(parents=True, exist_ok=True)

if not INPUT.is_file() or not GTF.is_file():
    raise FileNotFoundError(f'Missing manifest or annotation: {INPUT}, {GTF}')
with INPUT.open(newline='') as handle:
    entries = list(csv.DictReader(handle))
selected = [r for r in entries if r['sampleType'] == 'tumor_RNA']
if len(selected) != 80:
    raise RuntimeError(f'Expected exactly 80 tumor_RNA records; found {len(selected)}')
seen = set()
rows = []
for r in selected:
    pid = r['sampleName'].strip()
    bam = Path(r['bam'].strip())
    if pid in seen:
        raise RuntimeError(f'Duplicate tumor RNA patient: {pid}')
    seen.add(pid)
    if '.genomic.' not in bam.name:
        raise RuntimeError(f'Not a genomic RNA BAM for {pid}: {bam}')
    if not bam.is_file() or bam.stat().st_size == 0:
        raise FileNotFoundError(f'Missing/empty BAM for {pid}: {bam}')
    rows.append((pid, str(bam)))
rows.sort()
manifest = OUT / '00_StringTie_80_Genomic_RNA_BAMs.tsv'
with manifest.open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['patient_id', 'bam'])
    writer.writerows(rows)
print('Validated patients:', len(rows))
print('Annotation:', GTF)
print('Manifest:', manifest)
print('Pilot patient (array task 1):', rows[0][0])
