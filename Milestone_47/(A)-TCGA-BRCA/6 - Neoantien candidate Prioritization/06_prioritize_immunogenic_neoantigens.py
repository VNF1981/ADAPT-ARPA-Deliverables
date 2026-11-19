from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from docx import Document
from docx.shared import Pt, Inches, RGBColor

ROOT = Path('/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics')
BASE = ROOT / '05_Downstream_Analysis/02_80_patients'
IN3 = BASE / '03_results/03_Mutant_vs_Wildtype_Binding/01_All_Candidate_Mutant_WT_Binding.tsv'
IN5 = BASE / '03_results/05_pMHC_TCR_Prediction/03_Peptide_HLA_Recognition_Summary.tsv'
IN5P = BASE / '03_results/05_pMHC_TCR_Prediction/02_pMHC_TCR_Patient_Summary.tsv'
IN5RAW = ROOT / '08_pMHC_TCR_screening/01_old_pipeline/03_pMTnet/05_downstream_analysis/02_results/01_TCGA_BRCA_80_pMTnet_Combined_Predictions.tsv'
OUT = BASE / '03_results/06_Immunogenic_Candidate_Prioritization'
OUT.mkdir(parents=True, exist_ok=True)

for f in (IN3, IN5, IN5P, IN5RAW):
    if not f.is_file():
        raise FileNotFoundError(f'Expected input missing: {f}')

a = pd.read_csv(IN3, sep='\t', low_memory=False)
b = pd.read_csv(IN5, sep='\t', low_memory=False)
p = pd.read_csv(IN5P, sep='\t', low_memory=False)
raw = pd.read_csv(IN5RAW, sep='\t', usecols=['patient_id','Antigen','HLA','CDR3','Rank_for_summary'],low_memory=False)
need_a = ['Patient','HLA Allele','MT Epitope Seq','WT Epitope Seq','Gene Name','Variant Type','Mutation','Best MT IC50 (nM)','Corresponding WT IC50 (nM)','WT/MT ratio','Chromosome','Start','Stop']
need_b = ['patient_id','Antigen','HLA','tested_tcrs','minimum_rank','interactions_le_001','interactions_le_005']
need_p = ['patient_id','status','subtype']
for name, df, cols in [('Analysis 3',a,need_a),('Analysis 5 pMHC',b,need_b),('Analysis 5 patient',p,need_p)]:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f'{name}: missing columns {missing}')

if p.patient_id.nunique() != 80 or len(p) != 80:
    raise ValueError('Expected 80 distinct patients in Analysis 5 patient summary')
if not a['Patient'].isin(p.patient_id).all():
    raise ValueError('Analysis 3 includes patients not in Analysis 5 manifest')

def hla(x):
    return str(x).strip().upper().replace('HLA-','')

a['patient_id'] = a['Patient'].astype(str).str.strip()
a['peptide'] = a['MT Epitope Seq'].astype(str).str.strip()
a['hla_key'] = a['HLA Allele'].map(hla)
b['patient_id'] = b['patient_id'].astype(str).str.strip()
b['peptide'] = b['Antigen'].astype(str).str.strip()
b['hla_key'] = b['HLA'].map(hla)
keys = ['patient_id','peptide','hla_key']
if a[keys].isna().any().any() or b[keys].isna().any().any():
    raise ValueError('Null patient, peptide, or HLA key')
if b.duplicated(keys).any():
    raise ValueError('Duplicate pMHC identifiers in Analysis 5 summary')

num_a = ['Best MT IC50 (nM)','Corresponding WT IC50 (nM)','WT/MT ratio']
for c in num_a:
    a[c] = pd.to_numeric(a[c], errors='coerce')
raw['peptide'] = raw['Antigen'].astype(str).str.strip()
raw['hla_key'] = raw['HLA'].map(hla)
raw['Rank_for_summary'] = pd.to_numeric(raw['Rank_for_summary'],errors='coerce')
if raw['Rank_for_summary'].isna().any() or (~raw['Rank_for_summary'].between(0,1)).any():
    raise ValueError('Original pMTnet interaction ranks invalid')
best_tcr = raw.sort_values(['Rank_for_summary','CDR3'],kind='stable').drop_duplicates(keys)
best_tcr = best_tcr[keys + ['CDR3','Rank_for_summary']].rename(columns={'CDR3':'best_predicted_CDR3beta','Rank_for_summary':'best_predicted_TCR_rank'})
num_b = ['tested_tcrs','minimum_rank','interactions_le_001','interactions_le_005']
for c in num_b:
    b[c] = pd.to_numeric(b[c],errors='coerce')
if a['Best MT IC50 (nM)'].isna().any() or (a['Best MT IC50 (nM)'] <= 0).any():
    raise ValueError('Analysis 3 has nonpositive or missing mutant IC50 values')
if b[num_b].isna().any().any() or (b['tested_tcrs'] <= 0).any() or (~b['minimum_rank'].between(0,1)).any():
    raise ValueError('Invalid pMTnet summary values')
if (b['interactions_le_001'] > b['interactions_le_005']).any() or (b['interactions_le_005'] > b['tested_tcrs']).any():
    raise ValueError('Inconsistent rank-threshold counts')

original_rows = len(a)
a = a.sort_values(['patient_id','peptide','hla_key','Best MT IC50 (nM)','Gene Name'], kind='stable')
reps = a.drop_duplicates(keys).copy()
reps['source_prediction_rows'] = reps.set_index(keys).index.map(a.groupby(keys).size())
reps['source_genes'] = reps.set_index(keys).index.map(a.groupby(keys)['Gene Name'].agg(lambda s: '; '.join(sorted(set(s.dropna().astype(str))))))
reps = reps.merge(b[keys + num_b], on=keys, how='left',validate='one_to_one',indicator=True)
reps = reps.merge(p[['patient_id','subtype','status']],on='patient_id',how='left',validate='many_to_one')
reps = reps.merge(best_tcr,on=keys,how='left',validate='one_to_one')
if reps['status'].isna().any():
    raise ValueError('Missing patient screening status')
if (reps.loc[reps['status'] == 'SCREENED','_merge'] != 'both').any():
    bad = reps.loc[(reps['status']=='SCREENED') & (reps['_merge']!='both'),keys]
    raise ValueError(f'Expected all screened patient pMHCs to match; unmatched: {bad.head(10).to_dict("records")}')
if (reps.loc[reps['status'] != 'SCREENED','_merge'] == 'both').any():
    raise ValueError('Unscreened patient unexpectedly matched TCR results')
if (reps.loc[reps['status']=='SCREENED','best_predicted_CDR3beta'].isna()).any():
    raise ValueError('Missing best TCR CDR3beta for screened candidate')
if not np.allclose(reps.loc[reps['status']=='SCREENED','minimum_rank'],reps.loc[reps['status']=='SCREENED','best_predicted_TCR_rank'],rtol=0,atol=1e-10):
    raise ValueError('Best TCR ranks disagree with Analysis 5 summary')
if len(b) != int((reps['_merge']=='both').sum()):
    raise ValueError('Analysis 5 contains pMHCs absent from Analysis 3')

reps['MT_IC50_nM'] = reps['Best MT IC50 (nM)']
reps['WT_IC50_nM'] = reps['Corresponding WT IC50 (nM)']
reps['WT_MT_ratio'] = reps['WT/MT ratio']
reps['HLA_binding_le_500nM'] = reps['MT_IC50_nM'] <= 500
reps['HLA_binding_le_50nM'] = reps['MT_IC50_nM'] <= 50
reps['mutant_only_HLA_binder'] = (reps['MT_IC50_nM'] <= 500) & (reps['WT_IC50_nM'] > 500)
reps['mutant_HLA_gain_ge_2fold'] = reps['WT_MT_ratio'] >= 2
reps['TCR_screening_available'] = reps['status'].eq('SCREENED')
reps['rank_le_001'] = reps['TCR_screening_available'] & reps['interactions_le_001'].fillna(0).gt(0)
reps['rank_le_005'] = reps['TCR_screening_available'] & reps['interactions_le_005'].fillna(0).gt(0)
reps['fraction_TCR_rank_le_001'] = reps['interactions_le_001'] / reps['tested_tcrs']

conditions = [
    ~reps['TCR_screening_available'],
    reps['HLA_binding_le_500nM'] & reps['rank_le_001'] & reps['mutant_only_HLA_binder'],
    reps['HLA_binding_le_500nM'] & reps['rank_le_001'],
    reps['HLA_binding_le_500nM'] & reps['rank_le_005'],
    reps['HLA_binding_le_500nM'],
    reps['rank_le_001']
]
choices = [
    'NOT_SCREENED',
    'A: HLA+TCR+mutant-specific binding gain',
    'B: HLA+TCR prioritized',
    'C: HLA+TCR exploratory',
    'D: HLA binder, no prioritized TCR',
    'E: weak HLA binding, TCR rank <=0.01'
]
reps['priority_tier'] = np.select(conditions, choices, default='F: weak HLA binding, no prioritized TCR')
reps['priority_order'] = reps['priority_tier'].map({name:i for i,name in enumerate(choices + ['F: weak HLA binding, no prioritized TCR'])})
reps['expression_status'] = 'NOT_ASSESSED'
reps['experimental_recognition_status'] = 'NOT_VALIDATED'
reps['TCR_rank_status'] = np.where(reps['TCR_screening_available'],'PREDICTED','NO_USABLE_TCR_DATA')

reps = reps.sort_values(['priority_order','fraction_TCR_rank_le_001','minimum_rank','MT_IC50_nM','patient_id','peptide'],ascending=[True,False,True,True,True,True],na_position='last').reset_index(drop=True)
reps['cohort_priority_position'] = np.arange(1,len(reps)+1)
columns = ['cohort_priority_position','patient_id','subtype','Gene Name','source_genes','Variant Type','Mutation','Chromosome','Start','Stop','peptide','hla_key','WT Epitope Seq','MT_IC50_nM','WT_IC50_nM','WT_MT_ratio','HLA_binding_le_50nM','HLA_binding_le_500nM','mutant_only_HLA_binder','mutant_HLA_gain_ge_2fold','TCR_screening_available','tested_tcrs','minimum_rank','interactions_le_001','interactions_le_005','best_predicted_CDR3beta','best_predicted_TCR_rank','fraction_TCR_rank_le_001','priority_tier','source_prediction_rows','expression_status','experimental_recognition_status']
reps[columns].to_csv(OUT/'01_All_Prioritized_Patient_pMHC_Candidates.tsv',sep='\t',index=False)
reps.loc[reps['priority_tier'].str.startswith(('A:','B:')), columns].head(100).to_csv(OUT/'02_Top_100_High_Priority_Candidates.tsv',sep='\t',index=False)
reps.loc[reps['priority_tier'].str.startswith(('A:','B:')),columns].groupby('patient_id',sort=False).head(10).to_csv(OUT/'03_Top_10_High_Priority_Per_Patient.tsv',sep='\t',index=False)

patients = p[['patient_id','subtype','status']].copy()
all_counts = reps.groupby('patient_id').size().rename('all_candidate_pMHCs')
patients = patients.join(all_counts,on='patient_id')
for tier_code in list('ABCDEF'):
    x = reps.loc[reps.priority_tier.str.startswith(f'{tier_code}:')].groupby('patient_id').size()
    patients[f'tier_{tier_code}'] = patients.patient_id.map(x).fillna(0).astype(int)
patients['high_priority_A_or_B'] = patients['tier_A'] + patients['tier_B']
patients['fraction_high_priority_among_candidate_pMHCs'] = np.where(patients.status=='SCREENED', patients.high_priority_A_or_B/patients.all_candidate_pMHCs,np.nan)
patients['HLA_binders_500nM'] = patients['patient_id'].map(reps.loc[reps.HLA_binding_le_500nM].groupby('patient_id').size()).fillna(0).astype(int)
patients['fraction_high_priority_among_HLA_binders'] = np.where((patients.status=='SCREENED')&(patients.HLA_binders_500nM>0), patients.high_priority_A_or_B/patients.HLA_binders_500nM,np.nan)
patients.to_csv(OUT/'04_Patient_Prioritization_Summary.tsv',sep='\t',index=False)

nA = int(patients.tier_A.sum()); nB=int(patients.tier_B.sum()); total=int(len(reps)); screened=int(p.status.eq('SCREENED').sum()); high=nA+nB
metrics = [
    ('Patients in TCGA-BRCA cohort',len(p)),
    ('Patients with usable pMTnet predictions',screened),
    ('Patients without TCR predictions',len(p)-screened),
    ('Input Analysis 3 prediction rows',original_rows),
    ('Distinct patient-specific mutant peptide-HLA combinations',total),
    ('Patient-peptide-HLA combinations evaluated by pMTnet',len(b)),
    ('Class I mutant HLA binders (IC50 <=500 nM), all 80 patients',int(reps.HLA_binding_le_500nM.sum())),
    ('Tier A: HLA binder + rank <=0.01 + mutant-only HLA binder',nA),
    ('Tier B: HLA binder + rank <=0.01, excluding tier A',nB),
    ('Tier A+B: high-priority candidates',high),
    ('Tier A+B as percentage of TCR-evaluated candidate pMHCs',round(100*high/len(b),2)),
    ('Tier C: HLA binder + rank <=0.05 but not <=0.01',int(patients.tier_C.sum())),
    ('Tier D: HLA binder without rank <=0.05',int(patients.tier_D.sum())),
    ('Tier E: weak HLA binder with rank <=0.01',int(patients.tier_E.sum())),
    ('Tier F: weak HLA binder without rank <=0.01',int(patients.tier_F.sum())),
    ('Unscreened candidate pMHC combinations',int((reps.priority_tier=='NOT_SCREENED').sum())),
    ('Candidate gene expression assessed','No'),
    ('Experimentally confirmed T-cell recognition','No')
]
summary = pd.DataFrame(metrics,columns=['metric','value'])
summary.to_csv(OUT/'05_Cohort_Prioritization_Summary.tsv',sep='\t',index=False)
audit = {
    'analysis': 'Analysis 6: integrated computational candidate prioritization',
    'input_analysis_3':str(IN3),'input_analysis_5':str(IN5),'input_patient_summary':str(IN5P),'input_original_pMTnet_predictions':str(IN5RAW),
    'analysis_3_rows':original_rows,'distinct_pMHC':total,
    'duplicate_analysis_3_rows_collapsed':original_rows-total,
    'matched_pMHC':len(b),'unmatched_pMHC':int((reps._merge!='both').sum()),
    'raw_negative_rank_warning':'Analysis 5 reports 37 negative raw ranks that were transformed to summary ranks; review original transformation before making biological claims about top zero-rank hits',
    'TCR_test_count_warning':'More tested TCRs increases opportunity to find >=1 low percentile rank; tiers are exploratory and not calibrated probabilities',
    'binding_warning':'Predicted HLA binding is not observed tumor-cell presentation',
    'gene_expression':'Not assessed: only 8 of 80 patients had usable NextNEOpi TPM outputs at analysis time',
    'original_pMTnet':'Uses original pMTnet, not pMTnet-Omni',
    'mutation_provenance':'Duplicate patient-peptide-HLA rows are collapsed using best mutant IC50; source_prediction_rows and source_genes show multiplicity',
    'screening_denominator':'Tier A+B fraction uses all matched pMTnet-evaluated pMHC combinations',
    'rank_thresholds':'<=0.01 primary; <=0.05 exploratory; not probabilities of binding',
    'rank_order':'Tier A then B then C then D then E then F; within tier proportion of tested TCRs with rank <=0.01 (descending), then minimum rank (ascending), mutant IC50 (ascending)',
    'missing_TCR_handling':'Unscreened means missing prediction, not no predicted recognition'
}
(OUT/'06_Analysis_6_Metadata_and_QC.json').write_text(json.dumps(audit,indent=2))

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':18,'axes.labelsize':15,'xtick.labelsize':11,'ytick.labelsize':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
colors = {'A':'#E64B35','B':'#4DBBD5','C':'#00A087','D':'#3C5488','E':'#8491B4','F':'#91D1C2'}
fig,ax=plt.subplots(figsize=(11,5.8))
labels = ['A','B','C','D','E','F','Not screened']
values = [int(patients[f'tier_{x}'].sum()) for x in 'ABCDEF'] + [int((reps.priority_tier=='NOT_SCREENED').sum())]
bars = ax.bar(labels,values,color=[colors[x] for x in 'ABCDEF']+['#999999'],edgecolor='black',linewidth=.6)
ax.bar_label(bars,padding=3,fontsize=11)
ax.set_ylabel('Patient-specific peptide–HLA candidates')
ax.set_xlabel('Computational prioritization tier')
ax.set_title('Integrated neoantigen candidate prioritization')
ax.set_ylim(0,max(values)*1.17)
fig.tight_layout();fig.savefig(OUT/'01_Candidate_Priority_Tiers.png');plt.close(fig)

ordered = patients[patients.status=='SCREENED'].sort_values(['high_priority_A_or_B','patient_id'],ascending=[False,True]).reset_index(drop=True)
fig,ax=plt.subplots(figsize=(23,8))
x=np.arange(len(ordered))
ax.bar(x,ordered.tier_A,color=colors['A'],edgecolor='black',linewidth=.4,label='Tier A')
ax.bar(x,ordered.tier_B,bottom=ordered.tier_A,color=colors['B'],edgecolor='black',linewidth=.4,label='Tier B')
ax.set_xticks(x);ax.set_xticklabels(ordered.patient_id,rotation=90,fontsize=9)
ax.set_xlim(-.7,len(x)-.3)
ax.set_ylabel('High-priority candidate pMHCs (A+B)',fontsize=16)
ax.set_xlabel('TCGA-BRCA patients (ordered by A+B candidate count)',fontsize=15)
ax.set_title('Patient-level high-priority neoantigen candidates',fontsize=19)
ax.legend(frameon=False,loc='upper right',fontsize=13)
ax.set_ylim(0,max(ordered.high_priority_A_or_B)*1.2)
fig.tight_layout();fig.savefig(OUT/'02_High_Priority_Candidates_By_Patient.png');plt.close(fig)

visible = reps[reps.TCR_screening_available & (reps.MT_IC50_nM>0)].copy()
fig,ax=plt.subplots(figsize=(9,6.5))
c = np.log2(visible['WT_MT_ratio'].clip(lower=1e-4,upper=1e4))
points = ax.scatter(visible.MT_IC50_nM,visible.minimum_rank,c=c,cmap='viridis',s=13,alpha=.6,vmin=-4,vmax=4,rasterized=True)
ax.set_xscale('log')
ax.axvline(500,color='black',ls='--',lw=1,label='HLA IC50 = 500 nM')
ax.axhline(.01,color='black',ls=':',lw=1,label='pMTnet rank = 0.01')
ax.set_xlabel('Mutant peptide–HLA predicted IC50 (nM)')
ax.set_ylabel('Best predicted TCR recognition rank (lower = stronger)')
ax.set_title('HLA binding versus predicted TCR recognition')
ax.legend(frameon=False,fontsize=10,loc='upper right')
cb=fig.colorbar(points,ax=ax);cb.set_label('log2(WT / MT IC50)')
fig.tight_layout();fig.savefig(OUT/'03_HLA_Binding_vs_TCR_Rank.png');plt.close(fig)

D=Document()
section=D.sections[0]
section.page_width=Inches(11.7);section.page_height=Inches(8.3)
section.left_margin=section.right_margin=Inches(.7)
section.top_margin=section.bottom_margin=Inches(.7)
for sty in ['Normal','Title','Heading 1','Heading 2']:
    s=D.styles[sty];s.font.name='Times New Roman';s.font.color.rgb=RGBColor(0,0,0);s.font.size=Pt(10)
D.add_heading('Analysis 6 — Immunogenic Neoantigen Candidate Prioritization',0)
D.add_paragraph('Original pVACseq HLA-binding predictions from Analysis 3 are integrated with original pMTnet recognition ranks from Analysis 5. These are hypotheses for experimental follow-up, not observed HLA presentation, measured TCR binding, or demonstrated tumor killing.')
D.add_paragraph('Tier A: predicted mutant HLA IC50 <=500 nM, >=1 TCR rank <=0.01, and corresponding wild-type peptide IC50 >500 nM. Tier B: HLA IC50 <=500 nM with >=1 TCR rank <=0.01, excluding Tier A. Tier C: HLA IC50 <=500 nM and TCR rank <=0.05 but not <=0.01. Tiers D-F are lower evidence priority. NOT_SCREENED denotes absent usable TCR predictions and must not be interpreted as biological nonrecognition.')

def add_table(title,df,cols,labels=None):
    q=D.add_paragraph(title);q.runs[0].bold=True
    t=D.add_table(rows=1,cols=len(cols));t.style='Table Grid'
    for cell,s in zip(t.rows[0].cells,labels or cols):cell.text=str(s)
    for record in df[cols].itertuples(index=False,name=None):
        cells=t.add_row().cells
        for cell,value in zip(cells,record):
            cell.text='' if pd.isna(value) else str(value)
    for row in t.rows:
        for cell in row.cells:
            for para in cell.paragraphs:
                for run in para.runs:
                    run.font.name='Times New Roman';run.font.size=Pt(8);run.font.color.rgb=RGBColor(0,0,0)

add_table('Table 1. Cohort-level screening and prioritization',summary,['metric','value'],['Metric','Value'])
D.add_page_break()
short = reps[reps.priority_tier.str.startswith(('A:','B:'))].head(12).copy()
short['MT_IC50_nM']=short.MT_IC50_nM.round(2)
short['minimum_rank']=short.minimum_rank.round(5)
short['tier'] = short.priority_tier.str.slice(0,1)
add_table('Table 2. First 12 high-priority candidate peptide–HLA–TCR combinations',short,['patient_id','Gene Name','peptide','hla_key','MT_IC50_nM','minimum_rank','best_predicted_CDR3beta','tier'],['Patient','Gene','Mutant peptide','HLA','MT IC50 (nM)','Min rank','Best TCR CDR3β','Tier'])
D.add_paragraph('Notes: A patient–peptide–HLA combination can arise from multiple variant rows; duplicate combinations were collapsed, with provenance retained in the complete TSV. A low rank is a relative prioritization score, not a probability of TCR binding. A larger tested TCR repertoire increases the chance of at least one low rank. Gene expression, peptide presentation, paired TCR alpha/beta sequences, functional recognition, and patient treatment response are not established by this analysis. The 37 transformed negative raw ranks in the Analysis 5 input should be checked before interpreting top zero-rank predictions.')
D.save(OUT/'07_Analysis_6_Summary_Tables.docx')

print('===== ANALYSIS 6 — IMMUNOGENIC NEOANTIGEN PRIORITIZATION =====')
print('Patients:',len(p),'screened:',screened)
print('Analysis 3 rows:',original_rows)
print('Distinct patient-peptide-HLA candidates:',total)
print('Matched pMTnet predictions:',len(b))
print('Tier A:',nA,'Tier B:',nB,'Tier A+B:',high)
print('Patients with >=1 Tier A+B:',int((patients.high_priority_A_or_B>0).sum()))
print('Missing TCR predictions:',', '.join(p.loc[p.status!='SCREENED','patient_id']))
print('Expression: NOT ASSESSED; experimental recognition: NOT VALIDATED')
print('Validation: PASS')
print('Outputs:',OUT)
