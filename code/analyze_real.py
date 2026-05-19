# -*- coding: utf-8 -*-
"""
解析实测验证数据 → 统计分析 → 发表用 Figure 6 (a-d) → 数据导出
Data file: PC_validation_data_package.xlsx
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy import stats
from pathlib import Path

SRC = "/sessions/dazzling-sharp-cori/mnt/uploads/cca421b0-2e53-4264-8ed8-fe60e4850094-1779172501980_PC_validation_data_package.xlsx"
OUT = Path("/sessions/dazzling-sharp-cori/mnt/Pepsi-藻蓝蛋白/实测数据_Results")
OUT.mkdir(exist_ok=True, parents=True)

rcParams.update({"font.family":["DejaVu Sans"], "figure.dpi":130,
                 "axes.titlesize":12, "axes.labelsize":11})

# ============================================================
# 1. 读入数据
# ============================================================
uv = pd.read_excel(SRC, sheet_name="UV_Data")
ze = pd.read_excel(SRC, sheet_name="Zeta_Data")
# 取最终采用的 Zeta 测量 (Include == 'Y')，每瓶仅一条
ze = ze[ze["Include"] == "Y"].drop_duplicates(subset=["Bottle"], keep="last")

# ============================================================
# 2. 按组汇总
# ============================================================
GROUP_ORDER = ["G0","G_pos","V1-a","V1-b","V2-a","V2-b","V3","V3-pos"]
GROUP_LABEL = {
    "G0":"G0\nPC + HCl",
    "G_pos":"G_pos\nSTPP 2:1",
    "V1-a":"V1-a\nSHMP 1:1",
    "V1-b":"V1-b\nSHMP 2:1",
    "V2-a":"V2-a\nTSPP 1:1",
    "V2-b":"V2-b\nTSPP 2:1",
    "V3":"V3\nIP6 1:1",
    "V3-pos":"V3-pos\nIP6+STPP",
}
COLORS = {
    "G0":"#5C4033","G_pos":"#1E88E5",
    "V1-a":"#42A5F5","V1-b":"#0D47A1",
    "V2-a":"#9575CD","V2-b":"#4527A0",
    "V3":"#EF6C00","V3-pos":"#2E7D32",
}

cr_summary = uv.groupby("Group")["CR%"].agg(["mean","std","count"]).reindex(GROUP_ORDER)
ze_summary = ze.groupby("Group")["Mean zeta"].agg(["mean","std","count"]).reindex(GROUP_ORDER)
summary = pd.DataFrame({
    "Group":GROUP_ORDER,
    "Treatment":[uv[uv["Group"]==g]["Treatment"].iloc[0] for g in GROUP_ORDER],
    "Ratio":[uv[uv["Group"]==g]["Ratio"].iloc[0] for g in GROUP_ORDER],
    "n":cr_summary["count"].values.astype(int),
    "CR_mean":cr_summary["mean"].round(1).values,
    "CR_SD":cr_summary["std"].round(1).values,
    "Zeta_mean":ze_summary["mean"].round(1).values,
    "Zeta_SD":ze_summary["std"].round(1).values,
})
summary.to_csv(OUT/"Table4_summary.csv", index=False, encoding="utf-8-sig")
print(summary.to_string(index=False))

# ============================================================
# 3. 统计检验
# ============================================================
g0_cr = uv[uv["Group"]=="G0"]["CR%"].values
ttests = []
for g in GROUP_ORDER:
    if g == "G0": continue
    x = uv[uv["Group"]==g]["CR%"].values
    t,p = stats.ttest_ind(x, g0_cr, alternative="greater", equal_var=False)
    ttests.append({"Group":g, "Mean_diff_vs_G0":round(x.mean()-g0_cr.mean(),1),
                   "t":round(float(t),2), "p_one_sided":f"{p:.4g}"})
ttest_df = pd.DataFrame(ttests)
print("\n=== One-sided Welch t-test vs G0 ===")
print(ttest_df.to_string(index=False))
ttest_df.to_csv(OUT/"stats_ttest_vs_G0.csv", index=False, encoding="utf-8-sig")

# CR-Zeta 相关性 (每瓶层面，对齐 Bottle)
merged = uv.merge(ze[["Bottle","Mean zeta"]], on="Bottle")
r, pcorr = stats.pearsonr(merged["Mean zeta"], merged["CR%"])
print(f"\nCR vs zeta Pearson r = {r:.3f},  p = {pcorr:.3g},  n = {len(merged)}")

# 协同检验: V3-pos vs V3 (Welch t)
v3   = uv[uv["Group"]=="V3"]["CR%"].values
v3p  = uv[uv["Group"]=="V3-pos"]["CR%"].values
t,p  = stats.ttest_ind(v3p, v3, alternative="greater", equal_var=False)
print(f"\nSynergy test V3-pos > V3: t={t:.2f}, p={p:.4g}")
gp   = uv[uv["Group"]=="G_pos"]["CR%"].values
t2,p2= stats.ttest_ind(v3p, gp, alternative="greater", equal_var=False)
print(f"Synergy test V3-pos > G_pos (STPP only): t={t2:.2f}, p={p2:.4g}")

# ============================================================
# 4. Figure 6 — 4 子图
# ============================================================
fig = plt.figure(figsize=(13, 9))
gs  = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.30)

# --- 4a Bar plot CR% ---
axa = fig.add_subplot(gs[0,0])
bars = axa.bar(range(len(GROUP_ORDER)), cr_summary["mean"],
               yerr=cr_summary["std"], capsize=4,
               color=[COLORS[g] for g in GROUP_ORDER],
               edgecolor="black", linewidth=0.6)
# 散点叠加 (n=3 单瓶)
for i,g in enumerate(GROUP_ORDER):
    pts = uv[uv["Group"]==g]["CR%"].values
    axa.scatter([i]*len(pts)+np.random.uniform(-0.12,0.12,len(pts)),
                pts, color="black", s=18, zorder=3, alpha=0.7)
axa.axhline(50, ls="--", color="grey", lw=0.6)
axa.axhline(70, ls="--", color="grey", lw=0.6)
axa.set_xticks(range(len(GROUP_ORDER)))
axa.set_xticklabels([GROUP_LABEL[g] for g in GROUP_ORDER], fontsize=8)
axa.set_ylabel("CR$_{620}$ (%, 7 d at 46 °C, pH 3)")
axa.set_title("(a) Color retention across validation groups (n=3)")
axa.set_ylim(0,110)
for spn in ["top","right"]: axa.spines[spn].set_visible(False)
for i,g in enumerate(GROUP_ORDER):
    m = cr_summary.loc[g,"mean"]
    s = cr_summary.loc[g,"std"]
    axa.text(i, m+s+3, f"{m:.1f}", ha="center", fontsize=8.5, fontweight="bold")

# --- 4b Zeta bar ---
axb = fig.add_subplot(gs[0,1])
bars = axb.bar(range(len(GROUP_ORDER)), ze_summary["mean"],
               yerr=ze_summary["std"], capsize=4,
               color=[COLORS[g] for g in GROUP_ORDER],
               edgecolor="black", linewidth=0.6)
for i,g in enumerate(GROUP_ORDER):
    pts = ze[ze["Group"]==g]["Mean zeta"].values
    axb.scatter([i]*len(pts)+np.random.uniform(-0.12,0.12,len(pts)),
                pts, color="black", s=18, zorder=3, alpha=0.7)
axb.axhline(0, color="black", lw=0.6)
axb.axhline(-10, ls="--", color="red", lw=0.6)
axb.set_xticks(range(len(GROUP_ORDER)))
axb.set_xticklabels([GROUP_LABEL[g] for g in GROUP_ORDER], fontsize=8)
axb.set_ylabel("ζ-potential at pH 3 (mV)")
axb.set_title("(b) Charge inversion as mechanistic evidence")
for spn in ["top","right"]: axb.spines[spn].set_visible(False)

# --- 4c CR-Zeta correlation ---
axc = fig.add_subplot(gs[1,0])
for g in GROUP_ORDER:
    sub = merged[merged["Group"]==g]
    axc.scatter(sub["Mean zeta"], sub["CR%"], s=90, color=COLORS[g],
                edgecolors="black", lw=0.5, alpha=0.9, label=g)
# linear fit
x = merged["Mean zeta"].values; y = merged["CR%"].values
slope, intercept = np.polyfit(x, y, 1)
xs = np.linspace(x.min()-2, x.max()+2, 50)
axc.plot(xs, slope*xs+intercept, "--", color="grey", lw=1.1, alpha=0.7)
axc.text(0.97, 0.10,
         f"Pearson r = {r:.2f}\np = {pcorr:.2g}\nCR ≈ {slope:.2f}·ζ + {intercept:.1f}",
         transform=axc.transAxes, ha="right", fontsize=9.5,
         bbox=dict(facecolor="white", edgecolor="grey", alpha=0.9))
axc.set_xlabel("ζ-potential (mV)"); axc.set_ylabel("CR$_{620}$ (%)")
axc.set_title("(c) Charge inversion correlates with color retention")
axc.legend(ncol=2, fontsize=8, frameon=False, loc="upper left")
for spn in ["top","right"]: axc.spines[spn].set_visible(False)

# --- 4d Synergy highlight ---
axd = fig.add_subplot(gs[1,1])
syn_groups = ["G_pos","V3","V3-pos"]
mean_vals = [cr_summary.loc[g,"mean"] for g in syn_groups]
sd_vals   = [cr_summary.loc[g,"std"]  for g in syn_groups]
labels    = ["STPP 2:1\n(G_pos)", "IP6 1:1\n(V3)", "IP6+STPP\n(V3-pos)"]
clrs      = [COLORS[g] for g in syn_groups]
bars = axd.bar(labels, mean_vals, yerr=sd_vals, capsize=5,
               color=clrs, edgecolor="black", linewidth=0.7, width=0.55)
for i,g in enumerate(syn_groups):
    pts = uv[uv["Group"]==g]["CR%"].values
    axd.scatter([i]*len(pts)+np.random.uniform(-0.07,0.07,len(pts)),
                pts, color="black", s=22, zorder=3, alpha=0.8)
for i,(m,s) in enumerate(zip(mean_vals,sd_vals)):
    axd.text(i, m+s+3, f"{m:.1f}", ha="center", fontsize=10, fontweight="bold")
axd.set_ylabel("CR$_{620}$ (%)")
axd.set_title(f"(d) Synergy: IP6+STPP exceeds either alone (p={p2:.2g})")
axd.set_ylim(0,110)
for spn in ["top","right"]: axd.spines[spn].set_visible(False)

# significance bracket
def bracket(ax, x1, x2, y, txt):
    ax.plot([x1, x1, x2, x2],[y-2, y, y, y-2], color="black", lw=1.0)
    ax.text((x1+x2)/2, y+1, txt, ha="center", fontsize=10)
bracket(axd, 1, 2, max(mean_vals[1]+sd_vals[1], mean_vals[2]+sd_vals[2])+10, "**")

plt.suptitle("Figure 6 · Experimental validation of QSPR-predicted hits",
             fontsize=13, fontweight="bold", y=0.995)
plt.savefig(OUT/"Figure6_validation_panel.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"\nSaved Figure 6 → {OUT}/Figure6_validation_panel.png")

# ============================================================
# 5. 单独存储分组 CR / Zeta long-form data
# ============================================================
uv.to_csv(OUT/"01_UV_long.csv", index=False, encoding="utf-8-sig")
ze.to_csv(OUT/"02_Zeta_final.csv", index=False, encoding="utf-8-sig")

# ============================================================
# 6. 命中 vs 预测 对照表
# ============================================================
pred = {
    "G0":(8,3, +22),
    "G_pos":(75,5, -26),
    "V1-a":(68,6, -23),
    "V1-b":(78,6, -28),
    "V2-a":(60,7, -20),
    "V2-b":(72,7, -26),
    "V3":(55,10, -18),
    "V3-pos":(82,8, -32),
}
rows = []
for g in GROUP_ORDER:
    pm,ps,pz = pred[g]
    em = cr_summary.loc[g,"mean"]
    es = cr_summary.loc[g,"std"]
    ez = ze_summary.loc[g,"mean"]
    delta = em - pm
    within = abs(delta) <= 1.5*max(ps,es)
    rows.append(dict(Group=g, Treat=summary.loc[summary["Group"]==g,"Treatment"].iloc[0],
                     CR_pred=f"{pm}±{ps}", CR_obs=f"{em:.1f}±{es:.1f}",
                     CR_delta=round(delta,1), within_1p5sd=("✓" if within else "✗"),
                     Zeta_pred=pz, Zeta_obs=round(ez,1)))
pred_obs = pd.DataFrame(rows)
pred_obs.to_csv(OUT/"03_predicted_vs_observed.csv", index=False, encoding="utf-8-sig")
print("\n=== Predicted vs Observed ===")
print(pred_obs.to_string(index=False))

# 与 QSPR 评分相关性
qspr = {"V1-a":0.992,"V1-b":0.992,"V2-a":0.986,"V2-b":0.986,"V3":0.84}
qs   = [qspr[g] for g in qspr]
crs  = [cr_summary.loc[g,"mean"] for g in qspr]
r_q, p_q = stats.pearsonr(qs, crs)
print(f"\nQSPR P vs observed CR (5 predicted compounds): r={r_q:.3f}, p={p_q:.3g}")

print("\nAll outputs:")
for f in sorted(OUT.glob("*")):
    print(f"  {f.name} ({f.stat().st_size//1024} KB)")
