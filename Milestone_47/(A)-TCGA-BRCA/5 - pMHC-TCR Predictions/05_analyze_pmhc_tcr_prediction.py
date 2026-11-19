import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document
from docx.shared import Inches, Pt, RGBColor

ROOT = Path('/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics')
BASE = ROOT / '05_Downstream_Analysis/02_80_patients'
PMT = ROOT / '08_pMHC_TCR_screening/01_old_pipeline/03_pMTnet/05_downstream_analysis/02_results'
SOURCE = PMT / '01_TCGA_BRCA_80_pMTnet_Combined_Predictions.tsv'
MANIFEST = BASE / '02_data/01_TCGA_BRCA_80_Input_Manifest.tsv'
OUT = BASE / '03_results/05_pMHC_TCR_Prediction'
OUT.mkdir(parents=True, exist_ok=True)

if not SOURCE.is_file():
    raise FileNotFoundError(f'Missing completed pMTnet combined predictions: {SOURCE}')
if not MANIFEST.is_file():
    raise FileNotFoundError(f'Missing patient manifest: {MANIFEST}')

manifest = pd.read_csv(MANIFEST, sep='\t', dtype=str, keep_default_na=False)
if 'patient_id' not in manifest.columns or len(manifest) != 80 or manifest['patient_id'].nunique() != 80:
    raise ValueError('Expected 80 unique patient IDs in input manifest')
patient_ids = set(manifest['patient_id'])

original = pd.read_csv(SOURCE, sep='\t', dtype=str, keep_default_na=False, low_memory=False)
required = ['patient_id', 'CDR3', 'Antigen', 'HLA', 'Rank_raw', 'Rank_for_summary']
missing = [x for x in required if x not in original.columns]
if missing:
    raise ValueError(f'Missing columns in combined pMTnet output: {missing}; observed: {list(original.columns)}')

for col in ['patient_id', 'CDR3', 'Antigen', 'HLA']:
    original[col] = original[col].astype(str).str.strip()

if original[required[:4]].eq('').any().any():
    raise ValueError('Combined predictions include blank patient, TCR, peptide or HLA identifiers')
if not set(original['patient_id']).issubset(patient_ids):
    raise ValueError(f'Unknown patients in predictions: {sorted(set(original["patient_id"]) - patient_ids)}')

original['Rank_raw'] = pd.to_numeric(original['Rank_raw'], errors='coerce')
original['Rank_for_summary'] = pd.to_numeric(original['Rank_for_summary'], errors='coerce')
valid = original['Rank_for_summary'].notna() & original['Rank_for_summary'].between(0, 1, inclusive='both')
invalid = int((~valid).sum())
if invalid:
    raise ValueError(f'{invalid:,} invalid Rank_for_summary values: must be numeric from 0 to 1')

raw_negative = int((original['Rank_raw'] < 0).sum())
raw_missing = int(original['Rank_raw'].isna().sum())
changed_rank = int((~np.isclose(original['Rank_raw'], original['Rank_for_summary'], rtol=0, atol=1e-12, equal_nan=True)).sum())

data = original[required].copy()
data['Rank'] = data['Rank_for_summary']
data['pMHC'] = data['Antigen'] + ' | ' + data['HLA']

def count_le(g, t):
    return int((g['Rank'] <= t).sum())

def fraction(n, d):
    return float(n / d) if d else 0.0

counts = {t: count_le(data, t) for t in [0.01, 0.05, 0.10]}
if len(data) != 345374:
    raise ValueError(f'Expected 345,374 completed original pMTnet predictions, observed {len(data):,}')
if data['patient_id'].nunique() != 78:
    raise ValueError(f'Expected 78 patients with predictions, observed {data["patient_id"].nunique()}')

pairs = data[['patient_id', 'Antigen', 'HLA']].drop_duplicates()
tested_by_pair = data.groupby(['patient_id', 'Antigen', 'HLA'], sort=False).agg(
    tested_tcrs=('Rank', 'size'), minimum_rank=('Rank', 'min'), median_rank=('Rank', 'median')
).reset_index()

for t, name in [(0.01, '001'), (0.05, '005'), (0.10, '010')]:
    tmp = data.loc[data['Rank'] <= t].groupby(['patient_id', 'Antigen', 'HLA']).size().rename(f'interactions_le_{name}')
    tested_by_pair = tested_by_pair.join(tmp, on=['patient_id', 'Antigen', 'HLA'])
    tested_by_pair[f'interactions_le_{name}'] = tested_by_pair[f'interactions_le_{name}'].fillna(0).astype(int)
    tested_by_pair[f'fraction_le_{name}'] = tested_by_pair[f'interactions_le_{name}'] / tested_by_pair['tested_tcrs']

tested_by_pair = tested_by_pair.sort_values(['patient_id', 'minimum_rank', 'Antigen', 'HLA']).reset_index(drop=True)
tested_by_pair.to_csv(OUT / '03_Peptide_HLA_Recognition_Summary.tsv', sep='\t', index=False)

rows = []
for patient in manifest['patient_id']:
    g = data[data['patient_id'] == patient]
    gp = tested_by_pair[tested_by_pair['patient_id'] == patient]
    n = len(g)
    row = {
        'patient_id': patient,
        'subtype': manifest.loc[manifest['patient_id'] == patient, 'subtype'].iloc[0] if 'subtype' in manifest else '',
        'status': 'SCREENED' if n else 'NO_USABLE_TRB_OR_PREDICTIONS',
        'tested_interactions': n,
        'unique_TCR_CDR3beta': int(g['CDR3'].nunique()),
        'unique_peptide_HLA': len(gp),
        'minimum_rank': round(float(g['Rank'].min()), 6) if n else '',
        'median_rank': round(float(g['Rank'].median()), 6) if n else ''
    }
    for t, name in [(0.01, '001'), (0.05, '005'), (0.10, '010')]:
        k = count_le(g, t)
        m = int((gp[f'interactions_le_{name}'] > 0).sum())
        row[f'interactions_rank_le_{name}'] = k
        row[f'fraction_interactions_rank_le_{name}'] = round(fraction(k, n), 6) if n else ''
        row[f'peptide_HLA_with_rank_le_{name}'] = m
        row[f'fraction_peptide_HLA_with_rank_le_{name}'] = round(fraction(m, len(gp)), 6) if n else ''
    rows.append(row)

patient = pd.DataFrame(rows)
patient.to_csv(OUT / '02_pMHC_TCR_Patient_Summary.tsv', sep='\t', index=False)

best_tcr = data.sort_values(['Rank', 'patient_id', 'Antigen', 'HLA']).drop_duplicates(['patient_id', 'Antigen', 'HLA'])
top = tested_by_pair.merge(best_tcr[['patient_id', 'Antigen', 'HLA', 'CDR3', 'Rank']].rename(columns={'CDR3':'best_CDR3beta', 'Rank':'best_TCR_rank'}), on=['patient_id', 'Antigen', 'HLA'], how='left', validate='one_to_one')
top = top.sort_values(['minimum_rank', 'fraction_le_001', 'interactions_le_001', 'patient_id'], ascending=[True, False, False, True]).head(100)
top.to_csv(OUT / '04_Top_100_Peptide_HLA_Candidates.tsv', sep='\t', index=False)

screened = patient[patient['status'] == 'SCREENED']
summary = [
    ('Patients in cohort', len(patient)),
    ('Patients with pMTnet predictions', len(screened)),
    ('Patients without usable predictions', len(patient)-len(screened)),
    ('Evaluated TCR-peptide-HLA interactions', len(data)),
    ('Unique TCR CDR3beta sequences across cohort', data['CDR3'].nunique()),
    ('Distinct patient-peptide-HLA combinations', len(tested_by_pair)),
    ('Interactions with pMTnet rank <= 0.01', counts[0.01]),
    ('Interactions with pMTnet rank <= 0.05', counts[0.05]),
    ('Interactions with pMTnet rank <= 0.10', counts[0.10]),
    ('Fraction of interactions with rank <= 0.01 (%)', round(100*fraction(counts[0.01], len(data)), 2)),
    ('Fraction of interactions with rank <= 0.05 (%)', round(100*fraction(counts[0.05], len(data)), 2)),
    ('Fraction of interactions with rank <= 0.10 (%)', round(100*fraction(counts[0.10], len(data)), 2)),
    ('Patient-peptide-HLA combinations with at least one rank <= 0.01', int((tested_by_pair['interactions_le_001']>0).sum())),
    ('Patient-peptide-HLA combinations with at least one rank <= 0.05', int((tested_by_pair['interactions_le_005']>0).sum())),
    ('Patient-peptide-HLA combinations with at least one rank <= 0.10', int((tested_by_pair['interactions_le_010']>0).sum())),
    ('Raw pMTnet ranks below zero', raw_negative),
    ('Missing raw pMTnet ranks', raw_missing),
    ('Rows with different raw versus summary rank', changed_rank),
    ('Patients with receptor subtype available', int((screened['subtype'] != '').sum()))
]
summary_df = pd.DataFrame(summary, columns=['metric', 'value'])
summary_df.to_csv(OUT / '01_pMHC_TCR_Cohort_Summary.tsv', sep='\t', index=False)

audit = pd.DataFrame([
    ('Source combined predictions', str(SOURCE)),
    ('Source patient manifest', str(MANIFEST)),
    ('Expected interactions', '345374'),
    ('Observed interactions', str(len(data))),
    ('Expected screened patients', '78'),
    ('Observed screened patients', str(len(screened))),
    ('Missing patients', ', '.join(patient.loc[patient['status'] != 'SCREENED', 'patient_id'])),
    ('Rank interpretation', 'Lower percentile rank indicates higher predicted pMHC-TCR recognition'),
    ('Rank thresholds', '0.01, 0.05, 0.10: exploratory prioritization, not calibrated binding probability'),
    ('Rank transformation', 'Uses Rank_for_summary from original completed pMTnet downstream analysis'),
    ('Raw negative ranks', str(raw_negative)),
    ('No negative-control model runs', 'Do not interpret unscreened or high-rank pairs as experimentally confirmed negatives'),
    ('TCR source', 'Bulk RNA-seq CDR3beta; no alpha-beta pairing')
], columns=['item', 'value'])
audit.to_csv(OUT / '05_pMHC_TCR_Input_Validation.tsv', sep='\t', index=False)

plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False, 'savefig.dpi':300})
viridis = plt.get_cmap('viridis')
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.hist(data['Rank'], bins=np.linspace(0, 1, 51), color=viridis(.62), edgecolor='black', linewidth=.35)
for x, ls in [(0.01, '--'), (0.05, ':'), (0.10, '-.')]:
    ax.axvline(x, color='black', linewidth=.8, linestyle=ls, label=f'Rank ≤ {x:g}')
ax.set(xlabel='pMTnet percentile rank (lower = stronger prediction)', ylabel='Predicted interactions', title='pMHC–TCR recognition rank distribution')
ax.set_ylim(0, 10000)
ax.legend(loc='upper right', bbox_to_anchor=(0.99, 0.99), frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / '01_pMHC_TCR_Rank_Distribution.png')
plt.close(fig)

ordered = screened.sort_values('fraction_interactions_rank_le_001', ascending=False).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(20, 6.5))
ax.bar(np.arange(len(ordered)), 100*ordered['fraction_interactions_rank_le_001'], color=viridis(.65), edgecolor='black', linewidth=.35, width=.88)
ax.set(xlabel='Patients ranked by fraction of high-priority interactions', ylabel='Interactions with rank ≤ 0.01 (%)', title='Patient-level pMHC–TCR screening results')
ax.title.set_fontsize(24)
ax.xaxis.label.set_fontsize(20)
ax.yaxis.label.set_fontsize(20)
ax.tick_params(axis='y', labelsize=17)
ax.set_xlim(-.7, len(ordered)-.3)
ax.set_xticks(np.arange(len(ordered)), ordered['patient_id'], rotation=90, fontsize=11)
fig.tight_layout()
fig.savefig(OUT / '02_pMHC_TCR_Patient_Fractions.png')
plt.close(fig)

ordered = screened.sort_values('fraction_peptide_HLA_with_rank_le_001', ascending=False).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(20, 6.5))
ax.bar(np.arange(len(ordered)), 100*ordered['fraction_peptide_HLA_with_rank_le_001'], color=viridis(.50), edgecolor='black', linewidth=.35, width=.88)
ax.set(xlabel='Patients ranked by fraction of candidate pMHCs with at least one prioritized TCR', ylabel='Candidate peptide–HLA combinations (%)', title='Peptide–HLA complexes with predicted TCR recognition')
ax.title.set_fontsize(24)
ax.xaxis.label.set_fontsize(20)
ax.yaxis.label.set_fontsize(20)
ax.tick_params(axis='y', labelsize=17)
ax.set_xlim(-.7, len(ordered)-.3)
ax.set_xticks(np.arange(len(ordered)), ordered['patient_id'], rotation=90, fontsize=11)
fig.tight_layout()
fig.savefig(OUT / '03_pMHC_TCR_Patient_pMHC_Coverage.png')
plt.close(fig)

doc = Document()
sec = doc.sections[0]
sec.page_width = Inches(11.7)
sec.page_height = Inches(8.3)
sec.left_margin = sec.right_margin = Inches(.65)
sec.top_margin = sec.bottom_margin = Inches(.65)
for sty in ['Normal', 'Title', 'Heading 1', 'Heading 2']:
    obj = doc.styles[sty]
    obj.font.name = 'Times New Roman'
    obj.font.size = Pt(9)
    obj.font.color.rgb = RGBColor(0,0,0)


def add_table(title, table, cols, display=None):
    p = doc.add_paragraph(title)
    p.runs[0].bold = True
    display = display or cols
    t = doc.add_table(rows=1, cols=len(cols))
    t.style = 'Table Grid'
    for c, val in zip(t.rows[0].cells, display):
        c.text = str(val)
    for row in table[cols].itertuples(index=False, name=None):
        cells = t.add_row().cells
        for cell, val in zip(cells, row):
            cell.text = '' if pd.isna(val) else str(val)
    for row in t.rows:
        for cell in row.cells:
            for p in cell.paragraphs:
                for run in p.runs:
                    run.font.name = 'Times New Roman'
                    run.font.size = Pt(8)
                    run.font.color.rgb = RGBColor(0,0,0)


doc.add_paragraph('Analysis 5 — Predicted pMHC–TCR recognition in TCGA-BRCA (80 patients)')
doc.add_paragraph('Original pMTnet predictions using productive TCR CDR3β, mutant peptides and HLA alleles. Lower percentile ranks indicate stronger predicted recognition. Thresholds (0.01, 0.05, 0.10) are exploratory prioritization cutoffs, not experimentally validated binding classifications. Patients without usable productive TRB are retained in the patient table as unscreened.')
add_table('Table 1. Cohort summary', summary_df, ['metric','value'], ['Metric','Value'])
doc.add_page_break()
cols = ['patient_id','subtype','status','tested_interactions','unique_peptide_HLA','interactions_rank_le_001','fraction_interactions_rank_le_001','peptide_HLA_with_rank_le_001']
add_table('Table 2. Patient-level summary', patient, cols, ['Patient','Subtype','Status','Interactions','pMHCs','Rank ≤0.01','Fraction ≤0.01','pMHCs with hit'])
doc.add_page_break()
add_table('Table 3. Top 20 patient-peptide-HLA candidates', top.head(20), ['patient_id','Antigen','HLA','tested_tcrs','minimum_rank','interactions_le_001','best_CDR3beta'], ['Patient','Peptide','HLA','TCRs','Min rank','≤0.01 hits','Best CDR3β'])
doc.add_paragraph('Technical limitations: These predictions do not demonstrate measured binding or T-cell activation. The top candidates are chosen by computational scores; the minimum rank and presence of any hit are influenced by the number of TCRs tested. TRA–TRB pairs were not reconstructed from bulk RNA-seq. Results are based on existing original pMTnet outputs; pMTnet-Omni results are not included.')
doc.save(OUT / '06_pMHC_TCR_Summary_Tables.docx')

print('===== ANALYSIS 5 — pMHC–TCR PREDICTION =====')
print('Patients:', len(patient))
print('Patients screened:', len(screened))
print('Patients without predictions:', ', '.join(patient.loc[patient['status']!='SCREENED', 'patient_id']))
print('Predicted interactions:', len(data))
print('Distinct patient-peptide-HLA combinations:', len(tested_by_pair))
print('Ranks <= 0.01:', counts[0.01])
print('Ranks <= 0.05:', counts[0.05])
print('Ranks <= 0.10:', counts[0.10])
print('Raw negative ranks:', raw_negative)
print('Validation: PASS')
print('Outputs:', OUT)
