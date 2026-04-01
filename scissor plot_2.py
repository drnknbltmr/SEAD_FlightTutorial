"""
Scissor Plot Generator — AE3211-I Aircraft Tutorial
Group 45 | CRJ-EXX (Part 2) derived from Part 1 aerodynamics
=====================================================================
Modifications applied (per assignment §3.1):
  - Wing effective AR  ×1.25
  - CL_max_wf          ×1.20
  - Nacelle d +20%, l +20% → AC contribution scales ×(1.2²×1.2)=×1.728
  - CL_alpha_wf recalculated via DATCOM with new AR
  - Downwash gradient  recalculated with new AR (geometric formula)
  - Tail geometry / eta_h / CL_h_max: unchanged
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ─────────────────────────────────────────────────────────────────
#  GEOMETRY  (unchanged CRJ-1000 values)
# ─────────────────────────────────────────────────────────────────
S      = 77.39    # [m²]  wing area
b      = 26.17    # [m]   wing span
c_mac  = 3.67     # [m]   MAC
S_h    = 15.91    # [m²]  horizontal tail area  (unchanged — assessed below)
LEMAC  = 19.2084  # [m]   from nose

x_ac_nose_wing = LEMAC + 0.25 * c_mac
x_ac_nose_tail = 37.097
l_h = x_ac_nose_tail - x_ac_nose_wing   # 16.971 m

sweep_qc_rad = np.radians(23.59)

# ─────────────────────────────────────────────────────────────────
#  DATCOM LIFT CURVE SLOPE
# ─────────────────────────────────────────────────────────────────
def datcom_CL_alpha(AR, M, sweep_rad, kappa=1.0):
    beta      = np.sqrt(1.0 - M**2)
    tan_sweep = np.tan(sweep_rad)
    return (2.0 * np.pi * AR /
            (2.0 + np.sqrt(4.0 + (AR * beta / kappa)**2
                           * (1.0 + (tan_sweep / beta)**2))))

# ─────────────────────────────────────────────────────────────────
#  DOWNWASH (geometric, speed-independent)
# ─────────────────────────────────────────────────────────────────
def downwash(AR, sweep_rad):
    c = np.cos(sweep_rad)
    return 4.0 * c / (AR + 2.0 * c)

# ─────────────────────────────────────────────────────────────────
#  PART 1 — CRJ-1000 reference
# ─────────────────────────────────────────────────────────────────
AR_w_ref = 8.85
AR_h_ref = 4.58

CL_alpha_w_cruise_ref   = datcom_CL_alpha(AR_w_ref, 0.78, sweep_qc_rad)
CL_alpha_w_approach_ref = datcom_CL_alpha(AR_w_ref, 0.20, sweep_qc_rad)

# Fuselage correction (additive delta, independent of AR)
# From Part 1: CL_alpha_wf = 6.094 (cruise), 4.710 (approach)
delta_fus_cruise   = 6.094 - CL_alpha_w_cruise_ref    # ≈ -0.055
delta_fus_approach = 4.710 - CL_alpha_w_approach_ref  # ≈ -0.042

CL_alpha_wf_cruise_ref   = CL_alpha_w_cruise_ref   + delta_fus_cruise
CL_alpha_wf_approach_ref = CL_alpha_w_approach_ref + delta_fus_approach

# Tail (AR unchanged → same CLα both variants)
CL_alpha_h_cruise    = 4.480   # /rad
CL_alpha_h_approach  = 3.781   # /rad

# AC breakdown
x_ac_w     =  0.25
x_ac_f     = -0.1695
x_ac_n_ref =  0.0274
x_ac_ref   =  x_ac_w + x_ac_f + x_ac_n_ref   # 0.1079

de_da_ref = downwash(AR_w_ref, sweep_qc_rad)

CL_max_wf_ref        = 2.840
Cm_ac_wf_cruise_ref  =  0.1277
Cm_ac_wf_appro_ref   = -0.4754

# CG range from Part 1 loading diagram
cg_fwd_ref = 0.098
cg_aft_ref = 0.560
oew_cg_ref = 0.4545

# Tail
eta_h    = 1.0
CL_h_max = -0.8
SM       = 0.05

# ─────────────────────────────────────────────────────────────────
#  PART 2 — CRJ-EXX modifications
# ─────────────────────────────────────────────────────────────────

# AR +25%
AR_w_exx = AR_w_ref * 1.25          # = 11.0625

# Recompute wing CLα with new AR
CL_alpha_w_cruise_exx   = datcom_CL_alpha(AR_w_exx, 0.78, sweep_qc_rad)
CL_alpha_w_approach_exx = datcom_CL_alpha(AR_w_exx, 0.20, sweep_qc_rad)

# Apply same fuselage delta
CL_alpha_wf_cruise_exx   = CL_alpha_w_cruise_exx   + delta_fus_cruise
CL_alpha_wf_approach_exx = CL_alpha_w_approach_exx + delta_fus_approach

# Nacelle AC contribution: d +20%, l +20% → scale ×(1.2²×1.2) = ×1.728
x_ac_n_exx = x_ac_n_ref * (1.2**2 * 1.2)
x_ac_exx   = x_ac_w + x_ac_f + x_ac_n_exx

# New downwash (lower, because AR is larger)
de_da_exx = downwash(AR_w_exx, sweep_qc_rad)

# CLmax +20%
CL_max_wf_exx = CL_max_wf_ref * 1.20

# Cm_ac_wf approach: unchanged (conservative assumption —
#   improved HLD increases lift but Cm data not provided)
Cm_ac_wf_appro_exx = Cm_ac_wf_appro_ref

# ─────── UPDATE THESE FROM YOUR PART 2 LOADING DIAGRAM ──────────
cg_fwd_exx = None   # e.g. 0.12 — forward-most CG / MAC
cg_aft_exx = None   # e.g. 0.52 — aft-most CG / MAC
oew_cg_exx = None   # CG@OEW+Batt / MAC
# ────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────
#  SCISSOR PLOT FUNCTIONS
# ─────────────────────────────────────────────────────────────────
def tail_volume(Sh, lh, S, cmac):
    return (Sh * lh) / (S * cmac)

def neutral_point(x_ac, CLa_wf, CLa_h, de_da, eta, Vh):
    return x_ac + (CLa_h / CLa_wf) * eta * (1.0 - de_da) * Vh

def stab_line(Vh_arr, x_ac, CLa_wf, CLa_h, de_da, eta, sm):
    return np.array([neutral_point(x_ac, CLa_wf, CLa_h, de_da, eta, v) - sm
                     for v in Vh_arr])

def neutral_line(Vh_arr, x_ac, CLa_wf, CLa_h, de_da, eta):
    return stab_line(Vh_arr, x_ac, CLa_wf, CLa_h, de_da, eta, sm=0.0)

def ctrl_line(cg_arr, x_ac, Cm_ac, CLmax_wf, CLa_h, CLh_max, eta, de_da):
    """Min Vh to trim at each xcg (approach, full flaps)."""
    out = []
    den = eta * abs(CLh_max) * (1.0 - de_da)
    for xcg in cg_arr:
        num = Cm_ac + CLmax_wf * (xcg - x_ac)
        out.append(np.nan if den == 0 else -num / den)
    return np.array(out)

# ─────────────────────────────────────────────────────────────────
#  COMPUTE
# ─────────────────────────────────────────────────────────────────
Vh_design = tail_volume(S_h, l_h, S, c_mac)
Vh_arr    = np.linspace(0.0, 1.4, 600)
cg_arr    = np.linspace(-0.1, 0.9, 600)

configs = {
    'CRJ-1000': dict(
        x_ac=x_ac_ref,
        CLa_wf_c=CL_alpha_wf_cruise_ref,
        CLa_h_c=CL_alpha_h_cruise,
        de_da=de_da_ref,
        CLa_h_a=CL_alpha_h_approach,
        Cm_appro=Cm_ac_wf_appro_ref,
        CLmax=CL_max_wf_ref,
        cg_fwd=cg_fwd_ref,
        cg_aft=cg_aft_ref,
        oew_cg=oew_cg_ref,
    ),
    'CRJ-EXX': dict(
        x_ac=x_ac_exx,
        CLa_wf_c=CL_alpha_wf_cruise_exx,
        CLa_h_c=CL_alpha_h_cruise,
        de_da=de_da_exx,
        CLa_h_a=CL_alpha_h_approach,
        Cm_appro=Cm_ac_wf_appro_exx,
        CLmax=CL_max_wf_exx,
        cg_fwd=cg_fwd_exx,
        cg_aft=cg_aft_exx,
        oew_cg=oew_cg_exx,
    ),
}

# ─────────────────────────────────────────────────────────────────
#  STATUS REPORT
# ─────────────────────────────────────────────────────────────────
print(f"\n{'='*65}")
print(f"  PARAMETER COMPARISON  CRJ-1000 → CRJ-EXX")
print(f"{'='*65}")
print(f"  {'Parameter':<35} {'CRJ-1000':>10} {'CRJ-EXX':>10}")
print(f"  {'-'*55}")
rows = [
    ("AR_w",                    AR_w_ref,                  AR_w_exx),
    ("CLα_wf cruise  [/rad]",   CL_alpha_wf_cruise_ref,    CL_alpha_wf_cruise_exx),
    ("CLα_wf approach [/rad]",  CL_alpha_wf_approach_ref,  CL_alpha_wf_approach_exx),
    ("x_ac_nacelle  [MAC]",     x_ac_n_ref,                x_ac_n_exx),
    ("x_ac (a/c-less-tail) [MAC]", x_ac_ref,               x_ac_exx),
    ("dε/dα",                   de_da_ref,                  de_da_exx),
    ("CL_max_wf",               CL_max_wf_ref,             CL_max_wf_exx),
    ("Cm_ac_wf approach",       Cm_ac_wf_appro_ref,        Cm_ac_wf_appro_exx),
]
for name, v1, v2 in rows:
    print(f"  {name:<35} {v1:>10.4f} {v2:>10.4f}")

print(f"\n  Design tail volume  V̄_H = {Vh_design:.4f}")
print(f"\n  {'Config':<12} {'x_ac':>8} {'NP [%MAC]':>11} {'Stab aft [%MAC]':>17} {'Ctrl fwd [%MAC]':>17}")
print(f"  {'-'*65}")
for tag, p in configs.items():
    xnp  = neutral_point(p['x_ac'], p['CLa_wf_c'], p['CLa_h_c'],
                          p['de_da'], eta_h, Vh_design)
    den  = eta_h * abs(CL_h_max) * (1.0 - p['de_da'])
    xcf  = ((-Vh_design * den - p['Cm_appro']) / p['CLmax']) + p['x_ac']
    print(f"  {tag:<12} {p['x_ac']:>8.4f} {xnp*100:>11.1f} {(xnp-SM)*100:>17.1f} {xcf*100:>17.1f}")
print(f"{'='*65}")

# ─────────────────────────────────────────────────────────────────
#  PLOT — side-by-side comparison
# ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), sharey=True)
pct = mticker.FuncFormatter(lambda x, _: f'{x*100:.0f}%')

colors = {'neutral': 'gray', 'stab': '#185FA5', 'ctrl': '#1D9E75',
          'Vh': '#888780', 'cg': '#A32D2D', 'oew': 'darkorange'}

for ax, (tag, p) in zip(axes, configs.items()):
    cg_n  = neutral_line(Vh_arr, p['x_ac'], p['CLa_wf_c'], p['CLa_h_c'], p['de_da'], eta_h)
    cg_s  = stab_line   (Vh_arr, p['x_ac'], p['CLa_wf_c'], p['CLa_h_c'], p['de_da'], eta_h, SM)
    Vh_c  = ctrl_line   (cg_arr, p['x_ac'], p['Cm_appro'],  p['CLmax'],
                          p['CLa_h_a'], CL_h_max, eta_h, p['de_da'])

    ax.plot(cg_n, Vh_arr, color=colors['neutral'], lw=1.4, ls='--',
            label='Neutral stability line  (cruise)')
    ax.plot(cg_s, Vh_arr, color=colors['stab'],    lw=2.0,
            label=f'Stability line  (SM=5%, cruise)')
    ax.plot(cg_arr, Vh_c, color=colors['ctrl'],    lw=2.0,
            label='Controllability line  (approach)')

    # Shaded feasible region
    ax.fill_betweenx(Vh_arr,  cg_s, x2=0.9,   alpha=0.06, color=colors['stab'])
    ax.fill_between(cg_arr,  Vh_c,  y2=0.0,   alpha=0.06, color=colors['ctrl'])

    ax.axhline(Vh_design, color=colors['Vh'], lw=1.2, ls=':',
               label=f'Design V̄_H = {Vh_design:.3f}')

    # CG range bar
    if p['cg_fwd'] is not None and p['cg_aft'] is not None:
        ax.plot([p['cg_fwd'], p['cg_aft']], [Vh_design, Vh_design],
                color=colors['cg'], lw=3.5, solid_capstyle='round',
                label=f"CG range  [{p['cg_fwd']*100:.1f}%–{p['cg_aft']*100:.1f}% MAC]")
        ax.plot([p['cg_fwd'], p['cg_aft']], [Vh_design, Vh_design],
                'v', color=colors['cg'], ms=7)
    else:
        ax.annotate('← insert cg_fwd_exx\n   and cg_aft_exx',
                    xy=(0.4, Vh_design), xycoords=('data','data'),
                    color=colors['cg'], fontsize=8.5, ha='center')

    if p['oew_cg'] is not None:
        ax.axvline(p['oew_cg'], color=colors['oew'], lw=1.2, ls=':',
                   label=f"OEW CG = {p['oew_cg']*100:.1f}% MAC")

    ax.set_xlim(-0.1, 0.9)
    ax.set_ylim(0.0,  1.25)
    ax.xaxis.set_major_formatter(pct)
    ax.set_xlabel('CG position  (fraction of MAC from LEMAC)', fontsize=10)
    ax.set_title(f'Scissor plot — {tag}  (Group 45)', fontsize=12, fontweight='bold')
    ax.grid(True, ls='--', lw=0.4, alpha=0.5)
    ax.legend(fontsize=8.0, loc='lower left')

axes[0].set_ylabel('Tail volume coefficient  V̄_H  [–]', fontsize=10)

plt.tight_layout()
plt.savefig('scissor_plot_part2.png', dpi=150, bbox_inches='tight')
plt.savefig('scissor_plot_part2.pdf', dpi=150, bbox_inches='tight')
print("\nSaved: scissor_plot_part2.png / .pdf")