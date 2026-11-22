import argparse
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.stats import fisher_exact
from docx import Document
from docx.shared import Inches, Pt, RGBColor

NPG_BLUE = '#3C5488'
NPG_TEAL = '#00A087'
NPG_ORANGE = '#E64B35'


def two_field(x):
    if pd.isna(x):
        return None
    s = str(x).strip().upper().replace('HLA-', '')
    m = re.fullmatch(r'([A-Z0-9]+)\*(\d+):(\d+)(?::\d+)*', s)
    return f'{m.group(1)}*{m.group(2)}:{m.group(3)}' if m else None


def matched_state(row):
    if row['LOH'] != 'Y':
        return 'no_LOH_call'
    allele = row['hla_2field']
    lost = row['lost_2field']
    kept = row['kept_2field']
    if lost is None or kept is None or allele is None:
        return 'unresolved_allele'
    if lost == kept:
        return 'ambiguous_two_field'
    if allele == lost:
        return 'lost_allele'
    if allele == kept:
        return 'retained_allele'
    return 'unmatched_allele'


def doc_table(doc, heading, frame):
    doc.add_paragraph(heading)
    t = doc.add_table(rows=1, cols=len(frame.columns))
    t.style = 'Table Grid'
    for c, v in zip(t.rows[0].cells, frame.columns):
        c.text = str(v)
    for vals in frame.itertuples(index=False, name=None):
        cells = t.add_row().cells
        for c, v in zip(cells, vals):
            c.text = '' if pd.isna(v) else str(v)
    for row in t.rows:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.name = 'Times New Roman'
                    run.font.size = Pt(9)
                    run.font.color.rgb = RGBColor(0, 0, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--loh', required=True)
    ap.add_argument('--candidates', required=True)
    ap.add_argument('--burden-patients', required=True)
    ap.add_argument('--burden-stats', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    loh = pd.read_csv(args.loh, sep='\t', dtype={'patient_id': str, 'HLA': str, 'LOH': str})
    can = pd.read_csv(args.candidates, sep='\t', dtype={'patient_id': str, 'hla_key': str})
    required_loh = {'patient_id', 'HLA', 'LOH', 'LossHLA', 'KeptHLA'}
    required_can = {'patient_id', 'peptide', 'hla_key', 'priority_tier', 'MT_IC50_nM', 'minimum_rank', 'TCR_screening_available'}
    if not required_loh.issubset(loh) or not required_can.issubset(can):
        raise ValueError('Missing required input columns: ' + str(sorted(required_loh - set(loh)) + sorted(required_can - set(can))))
    if loh.duplicated(['patient_id', 'HLA']).any() or can.duplicated(['patient_id', 'peptide', 'hla_key']).any():
        raise ValueError('Duplicate patient/HLA LOH record or patient/peptide/HLA candidate')
    if not loh.LOH.isin(['Y', 'N']).all():
        raise ValueError('Unexpected LOH flag')
    can['hla_2field'] = can.hla_key.map(two_field)
    if can.hla_2field.isna().any():
        raise ValueError('Unparseable pMHC HLA alleles')
    can['HLA'] = can.hla_2field.str.split('*', regex=False).str[0]
    if not can.HLA.isin(['A', 'B', 'C']).all():
        raise ValueError('This analysis expects class I neoantigen candidates only')
    l = loh.loc[loh.HLA.isin(['A', 'B', 'C']), ['patient_id', 'HLA', 'LOH', 'LossHLA', 'KeptHLA']].copy()
    l['lost_2field'] = l.LossHLA.map(two_field)
    l['kept_2field'] = l.KeptHLA.map(two_field)
    can = can.merge(l, on=['patient_id', 'HLA'], how='left', validate='many_to_one', indicator=True)
    can['allele_state'] = can.apply(lambda row: 'no_SpecHLA_result' if row['_merge'] != 'both' else matched_state(row), axis=1)
    can = can.drop(columns=['_merge'])
    can['HLA_binder'] = pd.to_numeric(can.MT_IC50_nM, errors='coerce').le(500)
    can['TCR_evaluated'] = can.TCR_screening_available.astype(str).str.lower().eq('true')
    can['TCR_prioritized'] = can.TCR_evaluated & pd.to_numeric(can.minimum_rank, errors='coerce').le(0.01)
    can['Tier_A_or_B'] = can.priority_tier.astype(str).str.match(r'^[AB]:')
    if 'expression_supported_tpm_ge_1' in can:
        can['expression_supported'] = can.expression_supported_tpm_ge_1.astype(str).str.lower().eq('true')
    else:
        can['expression_supported'] = False
    can.to_csv(out / '16_HLA_LOH_Annotated_pMHC_Candidates.tsv', sep='\t', index=False)
    states = ['lost_allele', 'retained_allele', 'ambiguous_two_field', 'unmatched_allele', 'unresolved_allele', 'no_LOH_call', 'no_SpecHLA_result']
    summary = can.groupby('allele_state').agg(n_pMHC=('peptide', 'size'), n_patients=('patient_id', 'nunique'),
        n_HLA_binders=('HLA_binder', 'sum'), n_TCR_screened=('TCR_evaluated', 'sum'),
        n_TCR_prioritized=('TCR_prioritized', 'sum'), n_tier_A_B=('Tier_A_or_B', 'sum'),
        n_tier_A_B_expression_supported=('expression_supported', lambda x: int((x & can.loc[x.index, 'Tier_A_or_B']).sum()))).reindex(states, fill_value=0).reset_index()
    summary['recognition_fraction_among_screened'] = summary.n_TCR_prioritized.div(summary.n_TCR_screened.replace(0, np.nan))
    summary.to_csv(out / '16_HLA_LOH_Allele_Status_Summary.tsv', sep='\t', index=False)
    per = can.groupby(['patient_id', 'allele_state']).agg(n_pMHC=('peptide', 'size'), n_TCR_prioritized=('TCR_prioritized', 'sum'), n_tier_A_B=('Tier_A_or_B', 'sum')).reset_index()
    per.to_csv(out / '16_HLA_LOH_Allele_Status_Per_Patient.tsv', sep='\t', index=False)
    lost = summary.set_index('allele_state').loc['lost_allele']
    retained = summary.set_index('allele_state').loc['retained_allele']
    cont = [[int(lost.n_TCR_prioritized), int(lost.n_TCR_screened-lost.n_TCR_prioritized)],
            [int(retained.n_TCR_prioritized), int(retained.n_TCR_screened-retained.n_TCR_prioritized)]]
    if all(sum(row) > 0 for row in cont):
        odds, p = fisher_exact(cont, alternative='two-sided')
    else:
        odds, p = np.nan, np.nan
    pd.DataFrame([{'lost_recognized': cont[0][0], 'lost_screened': sum(cont[0]),
                   'retained_recognized': cont[1][0], 'retained_screened': sum(cont[1]),
                   'fisher_odds_ratio': odds, 'fisher_two_sided_p': p,
                   'note': 'Descriptive pooled candidate test; multiple pMHC per patient and different TCR opportunity mean observations are not independent.'}]).to_csv(out / '16_TCR_Recognition_Descriptive_Comparison.tsv', sep='\t', index=False)
    bur = pd.read_csv(args.burden_patients, sep='\t')
    stats = pd.read_csv(args.burden_stats, sep='\t')
    doc = Document()
    sec = doc.sections[0]
    sec.top_margin, sec.bottom_margin = Inches(.7), Inches(.7)
    sec.left_margin, sec.right_margin = Inches(.8), Inches(.8)
    norm = doc.styles['Normal']
    norm.font.name, norm.font.size, norm.font.color.rgb = 'Times New Roman', Pt(10), RGBColor(0, 0, 0)
    doc.add_paragraph('Table 8. HLA-LOH associations with neoantigen burden and predicted T-cell recognition.')
    compact = stats.loc[stats.subset.eq('All eligible'), ['comparison', 'n_positive', 'n_negative', 'median_positive', 'median_negative', 'mann_whitney_two_sided_p']].copy()
    compact.columns = ['LOH group', 'LOH+ N', 'LOH− N', 'LOH+ median burden', 'LOH− median burden', 'P-value']
    compact['P-value'] = compact['P-value'].map(lambda v: f'{v:.3g}' if pd.notna(v) else 'NA')
    doc_table(doc, 'A. Neoantigen burden', compact)
    tab = summary.loc[summary.allele_state.isin(['lost_allele', 'retained_allele']), ['allele_state', 'n_pMHC', 'n_HLA_binders', 'n_TCR_screened', 'n_TCR_prioritized', 'n_tier_A_B', 'n_tier_A_B_expression_supported']].copy()
    tab.columns = ['Allele status', 'pMHC targets', 'HLA binders', 'TCR-screened', 'TCR-prioritized', 'Tier A/B', 'Tier A/B TPM≥1']
    doc_table(doc, 'B. Peptide–HLA targets assigned to lost or retained HLA alleles', tab)
    doc.add_paragraph('Note. Patient burden is based on 5,533 mutation-specific candidates. The allele-status analysis uses distinct patient–peptide–HLA combinations. Only class I loci and HLA alleles matched at two-field resolution are used for lost/retained comparisons. Ambiguous, unmatched, and unavailable calls are retained in separate TSV outputs. TCR recognition uses original pMTnet rank ≤0.01.')
    doc.save(out / 'Table_8_HLA_LOH_Immune_Evasion.docx')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'savefig.dpi': 300, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.7))
    groups = [bur.loc[~bur.class_I_LOH, 'candidates_total'], bur.loc[bur.class_I_LOH, 'candidates_total']]
    bp = ax[0].boxplot(groups, tick_labels=['HLA-I LOH−', 'HLA-I LOH+'], showfliers=True, patch_artist=True)
    for patch, color in zip(bp['boxes'], [NPG_BLUE, NPG_TEAL]):
        patch.set_facecolor(color)
        patch.set_alpha(.75)
    ax[0].set_ylabel('Candidate mutation–peptide–HLA combinations')
    ax[0].set_title('A. Neoantigen burden by HLA-I LOH')
    for i, g in enumerate(groups, 1):
        ax[0].text(i, max(bur.candidates_total) * .89, f'n={len(g)}', ha='center')
    states2 = ['lost_allele', 'retained_allele']
    x = np.arange(2)
    ax[1].bar(x, summary.set_index('allele_state').loc[states2, 'n_pMHC'].to_numpy(), color=[NPG_ORANGE, NPG_TEAL])
    ax[1].set_xticks(x, ['Lost HLA', 'Retained HLA'])
    ax[1].set_ylabel('Distinct patient–peptide–HLA targets')
    ax[1].set_title('B. Targets assigned by HLA status')
    recognized = summary.set_index('allele_state').loc[states2, 'recognition_fraction_among_screened'].to_numpy(float) * 100
    ax[2].bar(x, np.nan_to_num(recognized), color=[NPG_ORANGE, NPG_TEAL])
    ax[2].set_xticks(x, ['Lost HLA', 'Retained HLA'])
    ax[2].set_ylabel('TCR-prioritized pMHCs among screened (%)')
    ax[2].set_ylim(0, max(10, np.nanmax(recognized) * 1.25) if np.isfinite(recognized).any() else 10)
    ax[2].set_title('C. Predicted TCR recognition')
    for i, state in enumerate(states2):
        record = summary.set_index('allele_state').loc[state]
        ax[2].text(i, np.nan_to_num(recognized[i]) + 1, f'{int(record.n_TCR_prioritized)}/{int(record.n_TCR_screened)}', ha='center', fontsize=9)
    fig.tight_layout()
    fig.savefig(out / 'Figure_10_HLA_LOH_Neoantigen_Immune_Evasion.png', dpi=300, bbox_inches='tight')
    fig.savefig(out / 'Figure_10_HLA_LOH_Neoantigen_Immune_Evasion.pdf', bbox_inches='tight')
    plt.close(fig)
    print(summary.to_string(index=False))
    print('Lost/retained Fisher test: P=', p, '; odds ratio=', odds)
    print('Saved:', out)


if __name__ == '__main__':
    main()
