"""Does the T_s collapse survive dropping the diffusive constraint?

Prints the number SI 5.3 quotes.  Equation 2 ties the fast time to one
coefficient, tau_fast(Q) = 1/(D0 Q^2).  The general two-mode form of
Bahadur 2019 frees it in every Q bin instead.  If T_s barely moves between
the two, the collapse does not rest on the constraint.

WHY THE SEED IS PART OF THE ANSWER.  With p_fast free, the two-mode
likelihood has local minima in which the fast term imitates the slow one.
Every fit in this paper therefore starts from the minimum of a p_fast
profile scanned on P_GRID and is then released, and the released fit is
kept only if it does not make chi2 worse.  Seeded that way the general form
moves T_s by -2.2%.  Started from the bare default p_fast = 0.5 with no
profile, three of the 34 groups land in one of those minima at chi2 of 27 to
33, and the move grows to -3.9%.  Both are small against the 5% spread
across Q bins, so the conclusion holds either way -- but the two numbers are
different, and only the seeded one is what the SI quotes.  This script
prints both so the difference cannot be mistaken for a discrepancy.

Reads data_2021-3/ only.  Draws nothing.
"""

import sys, os, csv, glob
from datetime import datetime
import numpy as np, h5py
from scipy.optimize import least_squares
HERE = '/home/8-id-i/QZ/Papers/PubCode-2021-IDPP/analysis/SAXPCS_8id'
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import xpcs_fit
from nexus_read import parse_name, read_start_time, read_g2
DATA = os.path.join(HERE, 'data_2021-3'); TRACE = os.path.join(DATA, 'trace_2021-3.csv')
STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')
QI, Q0, ERR_MAX, CHI2_MAX, AT_TARGET = [0,1,2,3,4], 0, 0.016, 10.0, 0.01
P_GRID = np.round(np.arange(0.20, 1.31, 0.05), 2)
T_REF, TS_SEED, N_BOOT = 28.00, 0.22, 150
RUNS = {'H0402_D8_280C10p':28.00,'H0407_D6_282C10p':28.25,
        'H0404_D7_285C10p':28.50,'H0413_D9_287C10p':28.75}
# xpcs_fit.CONTRAST defaults to the 2022-1 beta.  Overwrite it with this
# beamtime's, exactly as the two figure scripts do: at the default every f here
# comes out high and the collapse returns 0.219 K where Figure 4 returns 0.224.
xpcs_fit.CONTRAST = xpcs_fit.measure_beta(STANDARD, '2021-3 static standard')
def t_zero(run):
    rows=[r for r in csv.DictReader(open(TRACE)) if r['run']==run]
    rows.sort(key=lambda r:int(r['frame']))
    T=np.array([float(r['qnw3_C']) for r in rows])
    i=int(np.argmax(T>=float(rows[0]['qnw3_set_C'])-AT_TARGET))
    return datetime.strptime(rows[i]['start_time'],'%Y-%m-%d %H:%M:%S')
def gen_seeded(tau,g2,g2e):
    prof=[(r['red_chi2'],pf,r) for pf in P_GRID
          for r in [xpcs_fit.fit_g2_global(tau,g2,g2e,QI,p1_fixed=pf)] if r is not None]
    if not prof: return None
    chi0,pf0,best=min(prof,key=lambda x:x[0])
    xpcs_fit.P_EXP_P0[0]=pf0
    free=xpcs_fit.fit_g2_global(tau,g2,g2e,QI)
    xpcs_fit.P_EXP_P0[0]=0.5
    return free if (free is not None and free['red_chi2']<=chi0) else best
def gen_default(tau,g2,g2e):
    xpcs_fit.P_EXP_P0[0]=0.5
    return xpcs_fit.fit_g2_global(tau,g2,g2e,QI)
def diffusive(tau,g2,g2e,q):
    prof=[(r['red_chi2'],pf,r) for pf in P_GRID
          for r in [xpcs_fit.fit_g2_diffusive(tau,g2,g2e,QI,q,p1_fixed=pf)] if r is not None]
    chi0,pf0,best=min(prof,key=lambda x:x[0])
    xpcs_fit.P_EXP_P0[0]=pf0
    free=xpcs_fit.fit_g2_diffusive(tau,g2,g2e,QI,q)
    xpcs_fit.P_EXP_P0[0]=0.5
    return free if (free is not None and free['red_chi2']<=chi0) else best
S={'diff':{}, 'seeded':{}, 'default':{}}
QV=None; bad=[]
for run,T in RUNS.items():
    t0=t_zero(run); acc={k:[[],[],[],[]] for k in S}
    for fp in sorted(glob.glob(os.path.join(DATA,f'Average_{run}_*_results.hdf')),
                     key=lambda p:parse_name(p)[1]):
        with h5py.File(fp,'r') as hf:
            w=(read_start_time(hf)-t0).total_seconds(); tau,g2,g2e,q=read_g2(hf)
        if w<0 or np.median(g2e[(tau>1e-4)&(tau<1e-2),0])>ERR_MAX: continue
        a=diffusive(tau,g2,g2e,q)
        if a is None or a['red_chi2']>=CHI2_MAX: continue
        b=gen_seeded(tau,g2,g2e); c=gen_default(tau,g2,g2e)
        if b is None or c is None: continue
        if c['red_chi2']>10: bad.append((run,w,c['red_chi2'],b['red_chi2']))
        for k,r in (('diff',a),('seeded',b),('default',c)):
            acc[k][0].append(w)
            acc[k][1].append([r['per_q'][i]['f'] for i in QI])
            acc[k][2].append([r['per_q'][i]['f_err'] for i in QI])
            acc[k][3].append([r['per_q'][i]['tau_fast'] for i in QI])
        QV=[q[i] for i in QI]
    for k in S: S[k][run]=(T,np.array(acc[k][0]),np.array(acc[k][1]),np.array(acc[k][2]),np.array(acc[k][3]))
def mres(p,T,t,f,e,tf=None):
    A,lxc,n=p[0],p[1],p[2]; ts=tf if tf is not None else p[3]
    return (A/(1.0+(t*np.exp((T-T_REF)/ts)/np.exp(lxc))**n)-f)/e
def collapse(D):
    def pooled(qi):
        T=np.concatenate([np.full(len(t),Tv) for Tv,t,*_ in D.values()])
        t=np.concatenate([t for _,t,*_ in D.values()])
        f=np.concatenate([F[:,qi] for _,_,F,*_ in D.values()])
        e=np.clip(np.concatenate([E[:,qi] for _,_,_,E,*_ in D.values()]),0.01,None)
        return T,t,f,e
    out=[]
    for qi in range(5):
        T,t,f,e=pooled(qi)
        r=least_squares(mres,[1.0,np.log(np.median(t)*3),2.0,TS_SEED],
                        bounds=([0.5,0.0,0.3,0.05],[1.3,20.0,12.0,3.0]),
                        max_nfev=20000,args=(T,t,f,e,None))
        out.append(r.x[3])
    return np.array(out)
print('\n@@ T_s per Q bin, then the mean')
res={}
for k in ('diff','seeded','default'):
    ts=collapse(S[k]); res[k]=ts
    print(f'@@  {k:8}: ' + ' '.join(f'{v:.3f}' for v in ts) + f'   mean {ts.mean():.4f} K')
print(f"@@ seeded  vs diffusive: {100*(res['seeded'].mean()/res['diff'].mean()-1):+.1f}%")
print(f"@@ default vs diffusive: {100*(res['default'].mean()/res['diff'].mean()-1):+.1f}%")
# SI 5.3 also quotes how far f itself moves, right beside the T_s shift, so it
# is printed from the same fits rather than from a second run.
print('\n@@ |df| at the lowest Q, general form against Equation 2')
for k in ('seeded', 'default'):
    d = np.concatenate([np.abs(S[k][r][2][:, Q0] - S['diff'][r][2][:, Q0]) for r in RUNS])
    g = np.concatenate([np.abs(S[k][r][2][:, Q0] - S['diff'][r][2][:, Q0])[S['diff'][r][2][:, Q0] < 0.15]
                        for r in RUNS])
    print(f'@@  {k:8}: at most {d.max():.3f}, median {np.median(d):.3f}'
          f'   (where the fast mode has nearly gone, f < 0.15: at most {g.max():.3f})')

# SI 5.3 reports the Q exponent the released tau_fast returns.  It has to come
# from the SEEDED fits, the same ones the rest of the paragraph quotes.
print('\n@@ Q exponent of the released tau_fast (seeded), where the fast mode is resolved')
_lq = np.log(np.array(QV))
_sl, _f0 = [], []
for r in RUNS:
    for tf, f in zip(S['seeded'][r][4], S['diff'][r][2][:, Q0]):
        if np.all(np.isfinite(tf)) and np.all(np.asarray(tf) > 0):
            _sl.append(np.polyfit(_lq, np.log(tf), 1)[0]); _f0.append(f)
_sl, _f0 = np.array(_sl), np.array(_f0)
_k = _f0 > 0.15
print(f'@@  median {np.median(_sl[_k]):.2f} over the {_k.sum()} groups with f > 0.15 at the '
      f'lowest Q, 80% between {np.percentile(_sl[_k], 10):.1f} and {np.percentile(_sl[_k], 90):.1f}')

print(f'\n@@ groups where the DEFAULT seed lands above chi2 10 ({len(bad)}):')
for run,w,cd,cs in bad:
    print(f'@@   {run} t_w={w:.0f}s  default chi2 {cd:.1f}  vs seeded {cs:.2f}')
