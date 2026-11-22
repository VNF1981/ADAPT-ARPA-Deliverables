import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--loh', required=True)
    ap.add_argument('--burden', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--exclude-patient', default='TCGA-BH-A18G')
    args = ap.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    loh = pd.read_csv(args.loh, sep='\t', dtype={'patient_id': str, 'HLA': str, 'LOH': str})
    burden = pd.read_csv(args.burden, sep='\t', dtype={'patient_id': str})
    if not {'patient_id', 'HLA', 'LOH'}.issubset(loh) or not {'patient_id', 'candidates_total'}.issubset(burden):
        raise ValueError('Expected LOH patient_id/HLA/LOH and burden patient_id/candidates_total columns')
    if loh.duplicated(['patient_id', 'HLA']).any() or burden.patient_id.duplicated().any():
        raise ValueError('Duplicate patient/HLA or patient burden records')
    if not loh.LOH.isin(['Y', 'N']).all():
        raise ValueError('Invalid SpecHLA LOH value')
    loci = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
    if len(loh) != loh.patient_id.nunique() * 8 or not loh.groupby('patient_id').HLA.apply(lambda x: set(x) == set(loci)).all():
        raise ValueError('SpecHLA input must contain exactly eight loci per patient')
    loh['positive'] = loh.LOH.eq('Y')
    per = loh.groupby('patient_id').positive.sum().rename('loh_positive_loci').reset_index()
    ci = loh[loh.HLA.isin(['A', 'B', 'C'])].groupby('patient_id').positive.sum().rename('class_I_positive_loci').reset_index()
    data = per.merge(ci, on='patient_id', validate='one_to_one').merge(burden[['patient_id', 'candidates_total']], on='patient_id', how='left', validate='one_to_one')
    if data.candidates_total.isna().any():
        raise ValueError('Burden input missing eligible patients: ' + ', '.join(data.loc[data.candidates_total.isna(), 'patient_id']))
    data['any_LOH'] = data.loh_positive_loci.gt(0)
    data['class_I_LOH'] = data.class_I_positive_loci.gt(0)
    data.to_csv(out / '15_HLA_LOH_Neoantigen_Burden_Per_Patient.tsv', sep='\t', index=False)
    rows = []
    for label in ['any_LOH', 'class_I_LOH']:
        for sensitivity, subset in [('All eligible', data), ('Without ' + args.exclude_patient, data[data.patient_id.ne(args.exclude_patient)])]:
            pos = subset.loc[subset[label], 'candidates_total'].to_numpy(float)
            neg = subset.loc[~subset[label], 'candidates_total'].to_numpy(float)
            if len(pos) and len(neg):
                p = float(mannwhitneyu(pos, neg, alternative='two-sided', method='auto').pvalue)
                a = np.subtract.outer(pos, neg)
                cliffs_delta = float(((a > 0).sum() - (a < 0).sum()) / a.size)
            else:
                p, cliffs_delta = np.nan, np.nan
            rows.append({'comparison': label, 'subset': sensitivity, 'n_positive': len(pos), 'n_negative': len(neg),
                         'median_positive': float(np.median(pos)) if len(pos) else np.nan,
                         'median_negative': float(np.median(neg)) if len(neg) else np.nan,
                         'mean_positive': float(np.mean(pos)) if len(pos) else np.nan,
                         'mean_negative': float(np.mean(neg)) if len(neg) else np.nan,
                         'mann_whitney_two_sided_p': p, 'cliffs_delta': cliffs_delta})
    summary = pd.DataFrame(rows)
    summary.to_csv(out / '15_HLA_LOH_Neoantigen_Burden_Statistics.tsv', sep='\t', index=False)
    print(summary.to_string(index=False))
    print('Saved:', out)


if __name__ == '__main__':
    main()
