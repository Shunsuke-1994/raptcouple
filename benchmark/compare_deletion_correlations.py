"""Compare RaptCouple(-dH) vs MPBind correlations with measured deletion effects.
Two dependent correlations sharing the outcome variable (relative binding).
"""
import numpy as np, csv, json
from scipy.stats import spearmanr, rankdata, t as tdist

ISHIDA_MPBIND = {  # notebook exp-10-01, deletion_scan.score.txt
 'GCUGUGUUACGUCCGGAUUGGGGACCUGCACGGCCCAUG':2.64,'GCUGUGUCACGUCCGGAUUGGGGACCUGCACGGCCCAUG':0.81,
 'GCUGUGUCUCGUCCGGAUUGGGGACCUGCACGGCCCAUG':2.52,'GCUGUGUCUAGUCCGGAUUGGGGACCUGCACGGCCCAUG':1.28,
 'GCUGUGUCUACUCCGGAUUGGGGACCUGCACGGCCCAUG':1.28,'GCUGUGUCUACGCCGGAUUGGGGACCUGCACGGCCCAUG':4.93,
 'GCUGUGUCUACGUCGGAUUGGGGACCUGCACGGCCCAUG':-0.89,'GCUGUGUCUACGUCCGAUUGGGGACCUGCACGGCCCAUG':-0.66,
 'GCUGUGUCUACGUCCGGUUGGGGACCUGCACGGCCCAUG':-5.01,'GCUGUGUCUACGUCCGGAUGGGGACCUGCACGGCCCAUG':2.52,
 'GCUGUGUCUACGUCCGGAUUGGGACCUGCACGGCCCAUG':3.22,'GCUGUGUCUACGUCCGGAUUGGGGCCUGCACGGCCCAUG':-0.66,
 'GCUGUGUCUACGUCCGGAUUGGGGACUGCACGGCCCAUG':-2.83,'GCUGUGUCUACGUCCGGAUUGGGGACCGCACGGCCCAUG':-3.07,
 'GCUGUGUCUACGUCCGGAUUGGGGACCUCACGGCCCAUG':-1.13}

def williams(x, y, z, n):
    """Hotelling-Williams t for r(x,y) vs r(x,z), on ranks -> Spearman."""
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    r12 = np.corrcoef(rx, ry)[0,1]; r13 = np.corrcoef(rx, rz)[0,1]; r23 = np.corrcoef(ry, rz)[0,1]
    det = 1 - r12**2 - r13**2 - r23**2 + 2*r12*r13*r23
    num = (r12 - r13) * np.sqrt((n-1)*(1+r23))
    den = np.sqrt(2*((n-1)/(n-3))*det + ((r12+r13)**2/4)*(1-r23)**3)
    tstat = num/den
    return r12, r13, r23, tstat, 2*tdist.sf(abs(tstat), n-3)

def boot(x, y, z, n, reps=20000, seed=0):
    rng = np.random.default_rng(seed); diffs=[]
    for _ in range(reps):
        i = rng.integers(0, n, n)
        if len(set(x[i]))<3 or len(set(y[i]))<3 or len(set(z[i]))<3: continue
        a = spearmanr(x[i], y[i]).statistic; b = spearmanr(x[i], z[i]).statistic
        if np.isnan(a) or np.isnan(b): continue
        diffs.append(a-b)
    d = np.array(diffs)
    lo, hi = np.percentile(d, [2.5, 97.5])
    p = 2*min((d<=0).mean(), (d>=0).mean())
    return len(d), d.mean(), lo, hi, max(p, 1/len(d))

out={}
# --- Ishida2020 -------------------------------------------------------------
rows=[r for r in csv.DictReader(open('/Users/sumishunsuke/Desktop/RNA/VariationalNeuralPotts/'
    'revision_drafts/versions/v2.13/figures/__none__','r'))] if False else None
import zipfile, xml.etree.ElementTree as ET
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
REL='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
z=zipfile.ZipFile('/Users/sumishunsuke/Desktop/RNA/VariationalNeuralPotts/revised_fig/source_data/Source_Data.xlsx')
wb=ET.fromstring(z.read('xl/workbook.xml'))
rels={r.get('Id'):r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
tgt=[rels[s.get(REL+'id')] for s in wb.iter(NS+'sheet') if s.get('name')=='Fig. 6b'][0].lstrip('/')
root=ET.fromstring(z.read(tgt if tgt in z.namelist() else 'xl/'+tgt))
recs=[]
for row in root.iter(NS+'row'):
    cells=[]
    for c in row.iter(NS+'c'):
        if c.get('t')=='inlineStr': cells.append(''.join(x.text or '' for x in c.iter(NS+'t')))
        else:
            v=c.find(NS+'v'); cells.append('' if v is None else v.text)
    if len(cells)==4 and cells[0] in ('Ishida2020','Adachi2024'): recs.append(cells)

for ds, mp in [('Ishida2020', ISHIDA_MPBIND), ('Adachi2024', None)]:
    sub=[r for r in recs if r[0]==ds]
    seqs=[r[1] for r in sub]; bind=np.array([float(r[2]) for r in sub]); ndh=np.array([float(r[3]) for r in sub])
    if mp is None:
        d={r['seq']:float(r['Z_Combined.MetaScore']) for r in csv.DictReader(open(
            '/Users/sumishunsuke/Desktop/RNA/VariationalNeuralPotts/notebooks_sourcedata_couplings/fig5b/adachi_del.csv'))}
        mpv=np.array([d[s] for s in seqs])
    else:
        mpv=np.array([mp[s] for s in seqs])
    n=len(bind)
    rc=spearmanr(bind,ndh); mb=spearmanr(bind,mpv)
    r12,r13,r23,tt,pw = williams(bind,ndh,mpv,n)
    nb,md,lo,hi,pb = boot(bind,ndh,mpv,n)
    out[ds]=dict(n=n, rho_raptcouple=round(rc.statistic,3), p_raptcouple=round(rc.pvalue,4),
                 rho_mpbind=round(mb.statistic,3), p_mpbind=round(mb.pvalue,4),
                 rho_between_predictors=round(r23,3), diff=round(rc.statistic-mb.statistic,3),
                 williams_t=round(tt,3), williams_df=n-3, williams_p=round(pw,4),
                 boot_reps=nb, boot_mean_diff=round(md,3),
                 boot_ci95=[round(lo,3),round(hi,3)], boot_p=round(pb,4))
print(json.dumps(out, indent=2))
