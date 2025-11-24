#!/usr/bin/env python3
import pandas as pd, numpy as np, os, sys, json
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("data", nargs="?", default="MU_VOB.csv")
parser.add_argument("--spss-mode", action="store_true", default=True,
                    help="Use non-robust SEs and sample SD (ddof=1) for z and ±1SD probes (matches SPSS).")
parser.add_argument("--outdir", default="results")
args = parser.parse_args()

DATA = args.data
MANIFEST = "column_manifest.json"
SPSS = args.spss_mode
OUTDIR = args.outdir
os.makedirs(OUTDIR, exist_ok=True)

with open(MANIFEST) as f:
    col = json.load(f)

df = pd.read_csv(DATA)
to_num = lambda s: pd.to_numeric(s, errors="coerce")

# Coerce numerics
for k, name in col.items():
    if name is None: continue
    df[name] = to_num(df[name])

# --------- Correlations (series-mean imputation) ---------
corr_cols = [col["MU"], col["VOB"], col["age"]]
if col["EPDS"]: corr_cols.append(col["EPDS"])
if "Maternal_Education" in df.columns:
    df["Maternal_Education"] = to_num(df["Maternal_Education"])
    corr_cols.append("Maternal_Education")

corr_df = df[corr_cols].copy()
for c in corr_cols:
    corr_df[c] = corr_df[c].fillna(corr_df[c].mean())

labels = corr_df.columns.tolist()
R = pd.DataFrame(index=labels, columns=labels, dtype=float)
P = pd.DataFrame(index=labels, columns=labels, dtype=float)
for a in labels:
    for b in labels:
        r, p = stats.pearsonr(corr_df[a], corr_df[b])
        R.loc[a,b] = r; P.loc[a,b] = p
R.to_csv(os.path.join(OUTDIR, "correlations_r.csv"))
P.to_csv(os.path.join(OUTDIR, "correlations_p.csv"))

# --------- ANCOVA: VOB ~ age + MU (Type II) ---------
anc = df[[col["VOB"], col["age"], col["MU"]]].dropna()
mod = smf.ols(f'{col["VOB"]} ~ {col["age"]} + {col["MU"]}', data=anc).fit()
a2 = anova_lm(mod, typ=2)
ss_res = a2.loc["Residual","sum_sq"]
a2["MS"] = a2["sum_sq"]/a2["df"]
a2["eta_p2"] = a2["sum_sq"]/(a2["sum_sq"]+ss_res)
a2.to_csv(os.path.join(OUTDIR, "ANCOVA_typeII.csv"))

# --------- Age-only model + residual correlation ---------
age_only = smf.ols(f'{col["VOB"]} ~ {col["age"]}', data=anc).fit()
a_age = anova_lm(age_only, typ=2)
a_age.to_csv(os.path.join(OUTDIR, "AgeOnly_ANOVA.csv"))
resid = age_only.resid
mu_aligned = anc[col["MU"]].loc[resid.index]
r_res, p_res = stats.pearsonr(resid, mu_aligned)
with open(os.path.join(OUTDIR, "ResidualCorr.txt"), "w") as f:
    f.write(f"r(resid VOB, MU) = {r_res:.3f}, p = {p_res:.3g}, df = {len(resid)-2}\n")

# --------- 2×2 mean split chi-square ---------
vres_series = df[col["VOBres"]] if col["VOBres"] and col["VOBres"] in df.columns else None
if vres_series is None or vres_series.isna().all():
    vres_series = age_only.resid.reindex(df.index)
sub = pd.DataFrame({ "MU": df[col["MU"]], "VOB_for_age": vres_series }).dropna()
sub["MU_bin"] = (sub["MU"] > sub["MU"].mean())
sub["VOB_bin"] = (sub["VOB_for_age"] > sub["VOB_for_age"].mean())
ct = pd.crosstab(sub["MU_bin"], sub["VOB_bin"])
chi2, p, dof, exp = stats.chi2_contingency(ct, correction=False)
phi = np.sqrt(chi2 / ct.values.sum())
ct.to_csv(os.path.join(OUTDIR, "ChiSquare_counts.csv"))
pd.DataFrame({"Chi2":[chi2], "df":[dof], "p":[p], "phi":[phi]}).to_csv(os.path.join(OUTDIR, "ChiSquare_results.csv"), index=False)

# --------- Inhibitory Control model (Type III) ---------
ic_cols = [col["IC"], col["sex"], col["EPDS"], col["VOBres"], col["age"], col["MU"]]
ic_dat = df[ic_cols].copy()
ic_dat = ic_dat.dropna(subset=[col["IC"], col["VOBres"], col["age"], col["MU"]])

# series-mean impute covariates to match SPSS table N
for c in [col["sex"], col["EPDS"]]:
    if c and c in ic_dat.columns:
        ic_dat[c] = ic_dat[c].fillna(ic_dat[c].mean())

formula = f'{col["IC"]} ~ {col["sex"]} + {col["EPDS"]} + {col["VOBres"]} + {col["age"]} + {col["MU"]} + {col["age"]}:{col["MU"]} + {col["VOBres"]}:{col["age"]} + {col["VOBres"]}:{col["MU"]} + {col["VOBres"]}:{col["age"]}:{col["MU"]}'

# Non-robust fit (SPSS-style)
mod_ic = smf.ols(formula, data=ic_dat).fit()
a3 = anova_lm(mod_ic, typ=3)
ss_err = a3.loc["Residual","sum_sq"]
a3["eta_p2"] = a3["sum_sq"]/(a3["sum_sq"]+ss_err)
a3.to_csv(os.path.join(OUTDIR, "IC_TypeIII_ANOVA.csv"))

# Simple slopes of MU at ±1 SD age and VOB_for_age
cov = mod_ic.cov_params()  # non-robust
b = mod_ic.params; names = b.index.tolist()
age_mean = ic_dat[col["age"]].mean()
v_mean = ic_dat[col["VOBres"]].mean()
# Sample SD (ddof=1) for SPSS matching
age_sd = ic_dat[col["age"]].std(ddof=1)
v_sd   = ic_dat[col["VOBres"]].std(ddof=1)

def slope(a, v):
    coef = b[col["MU"]] + b.get(f'{col["age"]}:{col["MU"]}',0)*a + b.get(f'{col["VOBres"]}:{col["MU"]}',0)*v + b.get(f'{col["VOBres"]}:{col["age"]}:{col["MU"]}',0)*a*v
    g = np.zeros(len(names))
    for nm, w in [(col["MU"],1.0),(f'{col["age"]}:{col["MU"]}',a),(f'{col["VOBres"]}:{col["MU"]}',v),(f'{col["VOBres"]}:{col["age"]}:{col["MU"]}',a*v)]:
        if nm in names: g[names.index(nm)] = w
    se = float(np.sqrt(g.T @ cov.values @ g))
    from scipy import stats
    t = coef/se if se>0 else np.nan
    p = 2*stats.t.sf(abs(t), df=mod_ic.df_resid)
    return coef, se, p

rows = []
for a_k, a_lab in [(-1,"-1 SD age"), (1,"+1 SD age")]:
    a_val = age_mean + a_k*age_sd
    for v_k, v_lab in [(-1,"-1 SD VOB-for-age"), (1,"+1 SD VOB-for-age")]:
        v_val = v_mean + v_k*v_sd
        coef, se, p = slope(a_val, v_val)
        rows.append({"Age level": a_lab, "VOB level": v_lab, "b(MU→IC)": coef, "SE": se, "p": p})
pd.DataFrame(rows).round(6).to_csv(os.path.join(OUTDIR, "IC_simple_slopes_raw_SPSS.csv"), index=False)

# Standardized DV using sample SD (ddof=1)
ICz = (ic_dat[col["IC"]] - ic_dat[col["IC"]].mean())/ic_dat[col["IC"]].std(ddof=1)
modz = smf.ols(f'ICz ~ {col["sex"]} + {col["EPDS"]} + {col["VOBres"]} + {col["age"]} + {col["MU"]} + {col["age"]}:{col["MU"]} + {col["VOBres"]}:{col["age"]} + {col["VOBres"]}:{col["MU"]} + {col["VOBres"]}:{col["age"]}:{col["MU"]}', data=ic_dat.assign(ICz=ICz)).fit()
bz = modz.params; covz = modz.cov_params(); namesz = bz.index.tolist()

def slopez(a, v):
    coef = bz[col["MU"]] + bz.get(f'{col["age"]}:{col["MU"]}',0)*a + bz.get(f'{col["VOBres"]}:{col["MU"]}',0)*v + bz.get(f'{col["VOBres"]}:{col["age"]}:{col["MU"]}',0)*a*v
    g = np.zeros(len(namesz))
    for nm, w in [(col["MU"],1.0),(f'{col["age"]}:{col["MU"]}',a),(f'{col["VOBres"]}:{col["MU"]}',v),(f'{col["VOBres"]}:{col["age"]}:{col["MU"]}',a*v)]:
        if nm in namesz: g[namesz.index(nm)] = w
    se = float(np.sqrt(g.T @ covz.values @ g))
    from scipy import stats
    t = coef/se if se>0 else np.nan
    p = 2*stats.t.sf(abs(t), df=modz.df_resid)
    return coef, se, p

rows = []
for a_k, a_lab in [(-1,"-1 SD age"), (1,"+1 SD age")]:
    a_val = age_mean + a_k*age_sd
    for v_k, v_lab in [(-1,"-1 SD VOB-for-age"), (1,"+1 SD VOB-for-age")]:
        v_val = v_mean + v_k*v_sd
        coef, se, p = slopez(a_val, v_val)
        rows.append({"Age level": a_lab, "VOB level": v_lab, "b_z(MU→ICz)": coef, "SEz": se, "p": p})
pd.DataFrame(rows).round(6).to_csv(os.path.join(OUTDIR, "IC_simple_slopes_standardized_SPSS.csv"), index=False)

print("Done (SPSS mode =", SPSS, "). Outputs in:", OUTDIR)
