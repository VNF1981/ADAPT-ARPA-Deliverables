import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.stats import wilcoxon
from docx import Document
from docx.shared import Pt, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH

PALETTE = plt.get_cmap('viridis')
NA_VALUES = ['NA', 'N/A', 'None', 'null', '.', '']
SUBTYPES = ['TNBC', 'HR+/HER2-', 'HR+/HER2+', 'HR-/HER2+']


def clean_columns(df):
    df.columns = [str(x).strip() for x in df.columns]
    return df


def find_column(df, options, required=False):
    lookup = {re.sub(r'[^a-z0-9]', '', c.lower()): c for c in df.columns}
    for option in options:
        k = re.sub(r'[^a-z0-9]', '', option.lower())
        if k in lookup:
            return lookup[k]
    if required:
        raise ValueError('Missing required column. Expected one of: ' + ', '.join(options))
    return None


def numeric(df, options):
    col = find_column(df, options)
    return pd.to_numeric(df[col], errors='coerce') if col else pd.Series(np.nan, index=df.index, dtype=float)


def strings(df, options):
    col = find_column(df, options)
    return df[col].astype('string').str.strip() if col else pd.Series(pd.NA, index=df.index, dtype='string')


def categorize(values):
    return np.select([values <= 50, values <= 500, values > 500], ['Strong (<=50 nM)', 'Intermediate (50-500 nM)', 'Weak (>500 nM)'], default='Unavailable')


def subtype_of(row):
    for key, value in row.items():
        k = re.sub(r'[^a-z0-9]', '', str(key).lower())
        if k in ('subtype', 'breastcancersubtype', 'ihcsubtype', 'receptorsubtype', 'molecularsubtype') and pd.notna(value):
            s = str(value).upper().replace(' ', '').replace('−', '-').replace('–', '-')
            if 'TNBC' in s or 'TRIPLENEGATIVE' in s:
                return 'TNBC'
            if 'HER2' in s:
                hr = 'HR+' if 'HR+' in s else 'HR-' if 'HR-' in s else None
                her2 = 'HER2+' if 'HER2+' in s else 'HER2-' if 'HER2-' in s else None
                if hr and her2:
                    return hr + '/' + her2
    return 'Unknown'


def patient_of(row):
    for k, v in row.items():
        name = re.sub(r'[^a-z0-9]', '', str(k).lower())
        if name in ('patientid', 'patient', 'caseid', 'tcgapatient', 'submitterid') and pd.notna(v):
            found = re.search(r'TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}', str(v))
            return found.group(0) if found else str(v).strip()
    return None


def resolve_path(raw, roots):
    if pd.isna(raw) or str(raw).strip() in ('', '.', 'NA'):
        return None
    p = Path(str(raw).strip()).expanduser()
    if p.is_file():
        return p
    for root in roots:
        q = root / p
        if q.is_file():
            return q
    return None


def candidate_files(row, roots):
    found = []
    for key, raw in row.items():
        if pd.isna(raw):
            continue
        k = str(key).lower()
        v = str(raw)
        if ('pvac' in k or 'epitope' in k or '.tsv' in v.lower()) and ('filter' in v.lower() or 'epitope' in v.lower()):
            p = resolve_path(v, roots)
            if p:
                found.append(p)
    return found


def choose_files(files):
    files = sorted(set(files))
    filtered = [f for f in files if '.filtered.tsv' in f.name.lower() and 'aggregated' not in f.name.lower()]
    all_epitopes = [f for f in files if '.all_epitopes.tsv' in f.name.lower() and 'aggregated' not in f.name.lower()]
    use = filtered if filtered else all_epitopes
    if any('combined' in f.name.lower() for f in use):
        use = [f for f in use if 'combined' in f.name.lower()]
    else:
        class_i = [f for f in use if 'mhc_i.' in f.name.lower() or 'class_i' in str(f).lower()]
        if class_i:
            use = class_i
    return use, 'filtered' if filtered else 'all_epitopes' if all_epitopes else 'unavailable'


def extract_file(path, patient, subtype):
    df = clean_columns(pd.read_csv(path, sep='\t', low_memory=False, na_values=NA_VALUES))
    mt = numeric(df, ['Best MT IC50 Score', 'Best MT Score'])
    wt = numeric(df, ['Corresponding WT IC50 Score', 'Corresponding WT Score'])
    mt_median = numeric(df, ['Median MT IC50 Score', 'Median MT Score'])
    wt_median = numeric(df, ['Median WT IC50 Score', 'Median WT Score'])
    needed = find_column(df, ['Best MT IC50 Score', 'Best MT Score'])
    paired = find_column(df, ['Corresponding WT IC50 Score', 'Corresponding WT Score'])
    if not needed or not paired:
        raise ValueError(f'{path}: Required paired IC50 columns unavailable. Columns: {list(df.columns)}')
    out = pd.DataFrame({
        'Patient': patient,
        'Subtype': subtype,
        'HLA Allele': strings(df, ['HLA Allele']),
        'MT Epitope Seq': strings(df, ['MT Epitope Seq']),
        'WT Epitope Seq': strings(df, ['WT Epitope Seq']),
        'Gene Name': strings(df, ['Gene Name', 'Gene']),
        'Variant Type': strings(df, ['Variant Type']),
        'Mutation': strings(df, ['Mutation', 'HGVSp']),
        'Chromosome': strings(df, ['Chromosome']),
        'Start': strings(df, ['Start']),
        'Stop': strings(df, ['Stop']),
        'Transcript': strings(df, ['Transcript']),
        'Best MT IC50 (nM)': mt,
        'Corresponding WT IC50 (nM)': wt,
        'Median MT IC50 (nM)': mt_median,
        'Median WT IC50 (nM)': wt_median,
        'Best Method': strings(df, ['Best MT IC50 Score Method', 'Best MT Score Method']),
        'Source File': str(path)
    })
    return out


def pct(num, den):
    return round(100 * num / den, 2) if den else np.nan


def summary_counts(g):
    valid = g[g['Paired Valid']]
    return {
        'Total candidate rows': len(g),
        'Paired valid rows': len(valid),
        'Paired valid (%)': pct(len(valid), len(g)),
        'Missing or invalid WT/MT': len(g) - len(valid),
        'MT stronger (ratio >1)': int((valid['WT/MT ratio'] > 1).sum()),
        'MT >=2-fold stronger': int((valid['WT/MT ratio'] >= 2).sum()),
        'MT >=10-fold stronger': int((valid['WT/MT ratio'] >= 10).sum()),
        'WT stronger (ratio <1)': int((valid['WT/MT ratio'] < 1).sum()),
        'Equal binding (ratio =1)': int((valid['WT/MT ratio'] == 1).sum()),
        'MT strong <=50 nM': int((valid['Best MT IC50 (nM)'] <= 50).sum()),
        'MT <=500 nM': int((valid['Best MT IC50 (nM)'] <= 500).sum()),
        'WT <=500 nM': int((valid['Corresponding WT IC50 (nM)'] <= 500).sum()),
        'MT-only binder <=500 nM': int(((valid['Best MT IC50 (nM)'] <= 500) & (valid['Corresponding WT IC50 (nM)'] > 500)).sum()),
        'WT-only binder <=500 nM': int(((valid['Best MT IC50 (nM)'] > 500) & (valid['Corresponding WT IC50 (nM)'] <= 500)).sum()),
        'Both bind <=500 nM': int(((valid['Best MT IC50 (nM)'] <= 500) & (valid['Corresponding WT IC50 (nM)'] <= 500)).sum()),
        'Neither binds <=500 nM': int(((valid['Best MT IC50 (nM)'] > 500) & (valid['Corresponding WT IC50 (nM)'] > 500)).sum()),
        'Median WT/MT ratio': round(valid['WT/MT ratio'].median(), 4) if len(valid) else np.nan,
        'Median MT IC50 (nM)': round(valid['Best MT IC50 (nM)'].median(), 2) if len(valid) else np.nan,
        'Median WT IC50 (nM)': round(valid['Corresponding WT IC50 (nM)'].median(), 2) if len(valid) else np.nan,
        'Median log2(WT/MT)': round(valid['log2(WT/MT)'].median(), 4) if len(valid) else np.nan
    }


def render_doc(path, cohort, by_patient, by_subtype, stats, sources):
    doc = Document()
    sec = doc.sections[0]
    sec.top_margin = Inches(.7)
    sec.bottom_margin = Inches(.7)
    normal = doc.styles['Normal']
    normal.font.name = 'Times New Roman'
    normal.font.size = Pt(10)
    normal.font.color.rgb = __import__('docx').shared.RGBColor(0, 0, 0)
    for style_name in ['Title', 'Heading 1', 'Heading 2']:
        st = doc.styles[style_name]
        st.font.name = 'Times New Roman'
        st.font.size = Pt(10)
        st.font.bold = True
        st.font.color.rgb = __import__('docx').shared.RGBColor(0, 0, 0)
    doc.add_paragraph('Analysis 3 — Mutant versus wild-type peptide–HLA binding')
    doc.add_paragraph('Predicted binding affinity (IC50, nM); lower values indicate stronger binding. Fold change = corresponding WT IC50 / best mutant IC50. Ratios >1 indicate stronger predicted mutant binding. Descriptive thresholds are 50 and 500 nM; these predictions are not experimental binding measurements.')
    def table(title, frame):
        doc.add_paragraph(title)
        t = doc.add_table(rows=1, cols=len(frame.columns))
        t.style = 'Table Grid'
        for j, col in enumerate(frame.columns):
            t.rows[0].cells[j].text = str(col)
        for row in frame.itertuples(index=False, name=None):
            cells = t.add_row().cells
            for j, value in enumerate(row):
                cells[j].text = 'NA' if pd.isna(value) else str(value)
        for row in t.rows:
            for cell in row.cells:
                for p in cell.paragraphs:
                    for run in p.runs:
                        run.font.name = 'Times New Roman'
                        run.font.size = Pt(9)
                        run.font.color.rgb = __import__('docx').shared.RGBColor(0, 0, 0)
    table('Table 1. Cohort summary', pd.DataFrame([{'Measure': k, 'Result': v} for k, v in cohort.items()]))
    table('Table 2. Subtype summary', by_subtype[['Subtype', 'Total candidate rows', 'Paired valid rows', 'MT >=2-fold stronger', 'MT-only binder <=500 nM', 'Median WT/MT ratio']])
    table('Table 3. Statistical summary', pd.DataFrame([{'Measure': k, 'Result': v} for k, v in stats.items()]))
    doc.add_paragraph('Quality control: The source is pVACseq filtered.tsv when available, otherwise all_epitopes.tsv. Wild-type comparisons exclude unavailable or nonpositive IC50 values. Cohort percentages use paired evaluable candidate rows unless labeled otherwise. Multiple rows from a patient are not statistically independent; patient-level median log2 ratios are used for the cohort Wilcoxon signed-rank test. Per-allele predictions and transcripts can represent repeated observations of biological variants.')
    doc.add_paragraph('Data sources: ' + '; '.join(sorted(set(sources))))
    doc.save(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', default='.', help='02_80_patients directory')
    ap.add_argument('--manifest', default=None)
    ap.add_argument('--input-root', action='append', default=[])
    ap.add_argument('--output', default=None)
    args = ap.parse_args()
    base = Path(args.base).expanduser().resolve()
    manifest = Path(args.manifest).expanduser().resolve() if args.manifest else base / '02_data/01_TCGA_BRCA_80_Input_Manifest.tsv'
    output = Path(args.output).expanduser().resolve() if args.output else base / '03_results/03_Mutant_vs_Wildtype_Binding'
    output.mkdir(parents=True, exist_ok=True)
    if not manifest.exists():
        raise FileNotFoundError(f'Manifest not found: {manifest}')
    m = clean_columns(pd.read_csv(manifest, sep='\t', dtype='string'))
    roots = [base, manifest.parent, base / '02_data'] + [Path(x).expanduser().resolve() for x in args.input_root]
    records = []
    audit = []
    for _, row in m.iterrows():
        patient = patient_of(row)
        if not patient:
            continue
        subtype = subtype_of(row)
        listed = candidate_files(row, roots)
        if not listed:
            for root in roots:
                for p in root.glob('**/*.filtered.tsv'):
                    if patient in str(p):
                        listed.append(p)
                if not listed:
                    for p in root.glob('**/*.all_epitopes.tsv'):
                        if patient in str(p):
                            listed.append(p)
        selected, source_type = choose_files(listed)
        if not selected:
            audit.append({'Patient': patient, 'Subtype': subtype, 'Status': 'no pVACseq source found', 'Source Type': '', 'Source Files': ''})
            continue
        for p in selected:
            try:
                records.append(extract_file(p, patient, subtype))
                audit.append({'Patient': patient, 'Subtype': subtype, 'Status': 'loaded', 'Source Type': source_type, 'Source Files': str(p)})
            except Exception as e:
                audit.append({'Patient': patient, 'Subtype': subtype, 'Status': str(e), 'Source Type': source_type, 'Source Files': str(p)})
    pd.DataFrame(audit).to_csv(output / '00_Input_Audit.tsv', sep='\t', index=False)
    if not records:
        raise RuntimeError('No compatible pVACseq outputs loaded. Inspect 00_Input_Audit.tsv and pass --input-root /path/to/your/pvacseq/data if needed.')
    data = pd.concat(records, ignore_index=True)
    keys = ['Patient', 'HLA Allele', 'MT Epitope Seq', 'WT Epitope Seq', 'Transcript', 'Mutation', 'Best MT IC50 (nM)', 'Corresponding WT IC50 (nM)']
    data = data.drop_duplicates(subset=keys).reset_index(drop=True)
    data['Paired Valid'] = (data['Best MT IC50 (nM)'] > 0) & (data['Corresponding WT IC50 (nM)'] > 0) & np.isfinite(data['Best MT IC50 (nM)']) & np.isfinite(data['Corresponding WT IC50 (nM)'])
    data['WT/MT ratio'] = np.where(data['Paired Valid'], data['Corresponding WT IC50 (nM)'] / data['Best MT IC50 (nM)'], np.nan)
    data['log2(WT/MT)'] = np.log2(data['WT/MT ratio'])
    data['MT category'] = categorize(data['Best MT IC50 (nM)'])
    data['WT category'] = categorize(data['Corresponding WT IC50 (nM)'])
    data['Binding transition'] = np.select([
        data['Paired Valid'] & (data['Best MT IC50 (nM)'] <= 500) & (data['Corresponding WT IC50 (nM)'] > 500),
        data['Paired Valid'] & (data['Best MT IC50 (nM)'] > 500) & (data['Corresponding WT IC50 (nM)'] <= 500),
        data['Paired Valid'] & (data['Best MT IC50 (nM)'] <= 500) & (data['Corresponding WT IC50 (nM)'] <= 500),
        data['Paired Valid'] & (data['Best MT IC50 (nM)'] > 500) & (data['Corresponding WT IC50 (nM)'] > 500)
    ], ['Mutant-only binder', 'Wild-type-only binder', 'Both bind', 'Neither binds'], default='Not evaluable')
    data['>=2-fold improvement'] = data['WT/MT ratio'] >= 2
    data['>=10-fold improvement'] = data['WT/MT ratio'] >= 10
    data.to_csv(output / '01_All_Candidate_Mutant_WT_Binding.tsv', sep='\t', index=False)
    pairs = data[data['Paired Valid']].copy()
    if pairs.empty:
        raise RuntimeError('No positive paired mutant and wild-type binding scores; detailed input audit and candidate file were saved.')
    patients = pd.DataFrame([{'Patient': p, 'Subtype': g['Subtype'].iloc[0], **summary_counts(g)} for p, g in data.groupby('Patient')]).sort_values('Patient')
    patients.to_csv(output / '02_Patient_Summary.tsv', sep='\t', index=False)
    subtype = pd.DataFrame([{'Subtype': s, 'Patients': g['Patient'].nunique(), **summary_counts(g)} for s, g in data.groupby('Subtype')]).sort_values('Subtype')
    subtype.to_csv(output / '03_Subtype_Summary.tsv', sep='\t', index=False)
    cohort = {'Patients in manifest': len(set(filter(None, (patient_of(r) for _, r in m.iterrows())))), 'Patients loaded': data['Patient'].nunique(), **summary_counts(data)}
    pd.DataFrame([{'Measure': k, 'Value': v} for k, v in cohort.items()]).to_csv(output / '04_Cohort_Summary.tsv', sep='\t', index=False)
    transition = pd.crosstab(pairs['WT category'], pairs['MT category']).reindex(index=['Strong (<=50 nM)', 'Intermediate (50-500 nM)', 'Weak (>500 nM)'], columns=['Strong (<=50 nM)', 'Intermediate (50-500 nM)', 'Weak (>500 nM)'], fill_value=0)
    transition.to_csv(output / '05_Binding_Class_Transition_Matrix.tsv', sep='\t')
    top = pairs.sort_values(['WT/MT ratio', 'Best MT IC50 (nM)'], ascending=[False, True]).drop_duplicates(subset=['Patient', 'HLA Allele', 'MT Epitope Seq'])
    top.head(100).to_csv(output / '06_Top_100_Binding_Gain_Candidates.tsv', sep='\t', index=False)
    top.groupby('Patient', group_keys=False).head(10).to_csv(output / '07_Top_10_Per_Patient.tsv', sep='\t', index=False)
    method = pairs.groupby('Best Method', dropna=False).agg(N=('WT/MT ratio', 'size'), Median_Ratio=('WT/MT ratio', 'median'), Twofold_Gain=('>=2-fold improvement', 'sum')).reset_index()
    method.to_csv(output / '08_Prediction_Method_Summary.tsv', sep='\t', index=False)
    variant = pairs.groupby('Variant Type', dropna=False).agg(N=('WT/MT ratio', 'size'), Median_Ratio=('WT/MT ratio', 'median'), Twofold_Gain=('>=2-fold improvement', 'sum')).reset_index()
    variant.to_csv(output / '09_Variant_Type_Summary.tsv', sep='\t', index=False)
    median_log2 = patients['Median log2(WT/MT)'].dropna()
    if len(median_log2) >= 2 and (median_log2 != 0).any():
        test = wilcoxon(median_log2, alternative='two-sided', zero_method='wilcox', method='auto')
        w_stat = float(test.statistic)
        w_p = float(test.pvalue)
    else:
        w_stat, w_p = np.nan, np.nan
    stats = {'Patient-level Wilcoxon comparison': 'Median per-patient log2(WT IC50 / MT IC50) vs 0', 'N evaluable patients': len(median_log2), 'Wilcoxon statistic': w_stat, 'Two-sided p-value': w_p, 'Median patient log2 ratio': median_log2.median()}
    pd.DataFrame([{'Measure': k, 'Value': v} for k, v in stats.items()]).to_csv(output / '10_Patient_Level_Statistical_Test.tsv', sep='\t', index=False)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'savefig.dpi': 300, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(6, 5.5))
    plot = pairs.sample(min(20000, len(pairs)), random_state=47)
    x = np.log10(plot['Corresponding WT IC50 (nM)'])
    y = np.log10(plot['Best MT IC50 (nM)'])
    sc = ax.scatter(x, y, c=np.clip(plot['log2(WT/MT)'], -8, 8), cmap='viridis', s=9, alpha=.5, rasterized=True)
    bounds = [min(x.min(), y.min()), max(x.max(), y.max())]
    ax.plot(bounds, bounds, '--', color='black', lw=1)
    ax.set(xlabel='log10(wild-type IC50, nM)', ylabel='log10(mutant IC50, nM)', title='Paired predicted peptide–HLA binding')
    fig.colorbar(sc, ax=ax, label='log2(WT / MT)')
    fig.tight_layout()
    fig.savefig(output / '01_MT_vs_WT_Affinity_Scatter.png')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    lim = max(2, np.nanquantile(abs(pairs['log2(WT/MT)']), .99))
    ax.hist(pairs['log2(WT/MT)'].clip(-lim, lim), bins=50, color=PALETTE(.62), edgecolor='white', linewidth=.3)
    ax.axvline(0, color='black', linestyle='--')
    ax.axvline(1, color=PALETTE(.9), linestyle=':')
    ax.set(xlabel='log2(wild-type IC50 / mutant IC50)', ylabel='Candidate peptide–HLA pairs', title='Predicted affinity change (clipped at 99th percentile)')
    fig.tight_layout()
    fig.savefig(output / '02_Log2_Binding_Gain_Distribution.png')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4.8))
    im = ax.imshow(transition.values, cmap='viridis')
    labels = ['Strong', 'Intermediate', 'Weak']
    ax.set_xticks(range(3), labels)
    ax.set_yticks(range(3), labels)
    ax.set(xlabel='Mutant binding class', ylabel='Wild-type binding class', title='Binding class transitions (nM thresholds)')
    for (i, j), val in np.ndenumerate(transition.values):
        ax.text(j, i, f'{val:,}', ha='center', va='center', color='white' if im.norm(val) < .55 else 'black')
    fig.colorbar(im, ax=ax, label='Candidate pairs')
    fig.tight_layout()
    fig.savefig(output / '03_WT_to_MT_Binding_Transitions.png')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.7))
    p = patients.sort_values('Median log2(WT/MT)', ascending=False)
    ax.bar(range(len(p)), p['Median log2(WT/MT)'], color=PALETTE(.62), width=.85)
    ax.axhline(0, color='black', linewidth=.8)
    ax.set(xlabel='Patients ranked by median affinity gain', ylabel='Median log2(WT / MT)', title='Patient-level predicted mutant binding advantage')
    fig.tight_layout()
    fig.savefig(output / '04_Patient_Median_Binding_Gain.png')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.4))
    counts = pairs['Binding transition'].value_counts().reindex(['Mutant-only binder', 'Both bind', 'Wild-type-only binder', 'Neither binds'], fill_value=0)
    ax.bar(range(len(counts)), counts.values, color=[PALETTE(x) for x in [.2, .4, .65, .85]])
    ax.set_xticks(range(len(counts)), counts.index, rotation=18, ha='right')
    ax.set(ylabel='Paired peptide–HLA candidates', title='Binding status at 500 nM')
    fig.tight_layout()
    fig.savefig(output / '05_Binding_Status_at_500nM.png')
    plt.close(fig)
    render_doc(output / '11_Analysis_3_Summary_Tables.docx', cohort, patients, subtype, stats, data['Source File'].unique())
    with open(output / '12_Analysis_Metadata.json', 'w') as handle:
        json.dump({'input_manifest': str(manifest), 'output_directory': str(output), 'source_files': sorted(data['Source File'].unique()), 'thresholds_nM': [50, 500], 'primary_method': 'Best MT IC50 and corresponding WT IC50 from same algorithm', 'statistical_unit': 'patient median log2(WT/MT)', 'data_scope': 'pVACseq filtered if available; otherwise all epitopes', 'note': 'Predicted binding does not establish processing, presentation, T-cell recognition or immunogenicity.'}, handle, indent=2)
    print(f'Completed Analysis 3: {len(data):,} candidate rows, {len(pairs):,} evaluable pairs, {data.Patient.nunique()} patients. Results: {output}')
    print('Source audit saved to 00_Input_Audit.tsv; confirm patient coverage and filtered/all_epitopes provenance before reporting.')


if __name__ == '__main__':
    main()
