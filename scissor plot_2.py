"""
Scissor Plot — Part 2 with three V̄_H indicators per subplot:
  1. Dark red  — actual CG range bar at design V̄_H = 0.951
  2. Orange    — required V̄_H from Part 1 (CRJ-1000 aero + CRJ-1000 CG range)
  3. Purple    — required V̄_H for CRJ-EXX (EXX aero + EXX CG range)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ─────────────────────────────────────────────────────────────────
#  GEOMETRY
# ─────────────────────────────────────────────────────────────────
S, b, c_mac, S_h = 77.39, 26.17, 3.67, 15.91
LEMAC  = 19.2084
x_ac_nose_wing = LEMAC + 0.25 * c_mac
l_h    = 37.097 - x_ac_nose_wing
sweep_qc_rad = np.radians(26.9)

eta_h, CL_h_max, SM = 1.0, -0.8, 0.05
x_ac_w, x_ac_f = 0.25, -0.1695

# ─────────────────────────────────────────────────────────────────
#  DATCOM & DOWNWASH
# ─────────────────────────────────────────────────────────────────
def datcom_CL_alpha(AR, M, sweep_rad, kappa=1.0):
    beta = np.sqrt(1.0 - M**2)
    return (2.0*np.pi*AR /
            (2.0 + np.sqrt(4.0 + (AR*beta/kappa)**2
                           * (1.0 + (np.tan(sweep_rad)/beta)**2))))

def downwash(AR, sweep_rad):
    c = np.cos(sweep_rad)
    return 4.0*c / (AR + 2.0*c)

# ─────────────────────────────────────────────────────────────────
#  CRJ-1000 (Part 1)
# ─────────────────────────────────────────────────────────────────
AR_w_ref   = 8.85
x_ac_n_ref = 0.0274
x_ac_ref   = x_ac_w + x_ac_f + x_ac_n_ref

CL_alpha_w_c_ref = datcom_CL_alpha(AR_w_ref, 0.78, sweep_qc_rad)
CL_alpha_w_a_ref = datcom_CL_alpha(AR_w_ref, 0.20, sweep_qc_rad)
delta_fus_c = 6.094 - CL_alpha_w_c_ref
delta_fus_a = 4.710 - CL_alpha_w_a_ref
CLa_wf_c_ref = CL_alpha_w_c_ref + delta_fus_c   # 6.094
CLa_wf_a_ref = CL_alpha_w_a_ref + delta_fus_a   # 4.710

CLa_h_c = 4.480;  CLa_h_a = 3.781
de_da_ref    = downwash(AR_w_ref, sweep_qc_rad)
CL_max_ref   = 2.840
Cm_appro_ref = -0.4754
cg_fwd_ref, cg_aft_ref, oew_ref = 0.052, 0.575, 0.4545

# ─────────────────────────────────────────────────────────────────
#  CRJ-EXX (Part 2)
# ─────────────────────────────────────────────────────────────────
AR_w_exx   = AR_w_ref * 1.25
x_ac_n_exx = x_ac_n_ref * (1.2**2 * 1.2)
x_ac_exx   = x_ac_w + x_ac_f + x_ac_n_exx

CLa_wf_c_exx = datcom_CL_alpha(AR_w_exx, 0.78, sweep_qc_rad) + delta_fus_c
CLa_wf_a_exx = datcom_CL_alpha(AR_w_exx, 0.20, sweep_qc_rad) + delta_fus_a
de_da_exx    = downwash(AR_w_exx, sweep_qc_rad)
CL_max_exx   = CL_max_ref * 1.20
Cm_appro_exx = Cm_appro_ref
cg_fwd_exx, cg_aft_exx, oew_exx = 0.099, 0.598, 0.4967

# ─────────────────────────────────────────────────────────────────
#  SCISSOR FUNCTIONS
# ─────────────────────────────────────────────────────────────────
def tail_volume(Sh, lh, S, cmac): return (Sh*lh)/(S*cmac)
def NP(x_ac, CLa_wf, CLa_h, de_da, eta, Vh):
    return x_ac + (CLa_h/CLa_wf)*eta*(1-de_da)*Vh
def stab_line(Vh_arr, x_ac, CLa_wf, CLa_h, de_da, eta, sm):
    return np.array([NP(x_ac,CLa_wf,CLa_h,de_da,eta,v)-sm for v in Vh_arr])
def neutral_line(Vh_arr, x_ac, CLa_wf, CLa_h, de_da, eta):
    return stab_line(Vh_arr, x_ac, CLa_wf, CLa_h, de_da, eta, sm=0.0)
def ctrl_line(cg_arr, x_ac, Cm_ac, CLmax, CLa_h, CLh_max, eta, de_da):
    den = eta*abs(CLh_max)*(1-de_da)
    return np.array([-(Cm_ac + CLmax*(xcg-x_ac))/den for xcg in cg_arr])

def vh_required(x_ac, CLa_wf_c, CLa_h_c, de_da, CLa_h_a, Cm_appro, CLmax, 
                cg_fwd, cg_aft, eta=1.0, CLh_max=-0.8, SM=0.05):
    """Min V̄_H so that [cg_fwd, cg_aft] fits exactly in the feasible region."""
    Vh_s = (cg_aft + SM - x_ac) / ((CLa_h_c/CLa_wf_c)*eta*(1-de_da))
    Vh_c = -(Cm_appro + CLmax*(cg_fwd - x_ac)) / (eta*abs(CLh_max)*(1-de_da))
    return max(Vh_s, Vh_c), Vh_s, Vh_c

# COMPUTE
Vh_design = tail_volume(S_h, l_h, S, c_mac)
Vh_arr    = np.linspace(0.0, 1.6, 600)
cg_arr    = np.linspace(-0.1, 0.9, 600)

Vh_req_ref, Vh_rs_ref, Vh_rc_ref = vh_required(
    x_ac_ref, CLa_wf_c_ref, CLa_h_c, de_da_ref, CLa_h_a,
    Cm_appro_ref, CL_max_ref, cg_fwd_ref, cg_aft_ref)

Vh_req_exx, Vh_rs_exx, Vh_rc_exx = vh_required(
    x_ac_exx, CLa_wf_c_exx, CLa_h_c, de_da_exx, CLa_h_a,
    Cm_appro_exx, CL_max_exx, cg_fwd_exx, cg_aft_exx)

print(f"  CRJ-1000  Vh_req = {Vh_req_ref:.4f}  (stab={Vh_rs_ref:.4f}, ctrl={Vh_rc_ref:.4f})")
print(f"  CRJ-EXX   Vh_req = {Vh_req_exx:.4f}  (stab={Vh_rs_exx:.4f}, ctrl={Vh_rc_exx:.4f})")
print(f"  Design Vh         = {Vh_design:.4f}")

configs = {
    'CRJ-1000': dict(
        x_ac=x_ac_ref, CLa_wf_c=CLa_wf_c_ref, CLa_h_c=CLa_h_c,
        de_da=de_da_ref, CLa_h_a=CLa_h_a,
        Cm_appro=Cm_appro_ref, CLmax=CL_max_ref,
        cg_fwd=cg_fwd_ref, cg_aft=cg_aft_ref, oew=oew_ref,
        Vh_req=Vh_req_ref,
    ),
    'CRJ-EXX': dict(
        x_ac=x_ac_exx, CLa_wf_c=CLa_wf_c_exx, CLa_h_c=CLa_h_c,
        de_da=de_da_exx, CLa_h_a=CLa_h_a,
        Cm_appro=Cm_appro_exx, CLmax=CL_max_exx,
        cg_fwd=cg_fwd_exx, cg_aft=cg_aft_exx, oew=oew_exx,
        Vh_req=Vh_req_exx,
    ),
}

# ─────────────────────────────────────────────────────────────────
#  PLOT
# ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
pct = mticker.FuncFormatter(lambda x, _: f'{x*100:.0f}%')

C_NEUT  = 'gray'
C_STAB  = '#185FA5'
C_CTRL  = '#1D9E75'
C_VH    = '#888780'
C_CG    = '#A32D2D'    # actual CG bar  — dark red
C_REQ1  = '#C94B00'    # Part 1 req     — orange
C_REQ2  = '#7B2D8B'    # EXX req        — purple
C_OEW   = 'darkorange'

for ax, (tag, p) in zip(axes, configs.items()):
    cg_n = neutral_line(Vh_arr, p['x_ac'], p['CLa_wf_c'], p['CLa_h_c'], p['de_da'], eta_h)
    cg_s = stab_line   (Vh_arr, p['x_ac'], p['CLa_wf_c'], p['CLa_h_c'], p['de_da'], eta_h, SM)
    Vh_c = ctrl_line   (cg_arr, p['x_ac'], p['Cm_appro'],  p['CLmax'],
                         p['CLa_h_a'], CL_h_max, eta_h, p['de_da'])

    ax.plot(cg_n, Vh_arr, color=C_NEUT, lw=1.4, ls='--', label='Neutral stability  (cruise)')
    ax.plot(cg_s, Vh_arr, color=C_STAB, lw=2.0,           label='Stability line  (SM=5%, cruise)')
    ax.plot(cg_arr, Vh_c, color=C_CTRL, lw=2.0,           label='Controllability line  (approach)')

    ax.fill_betweenx(Vh_arr, cg_s, x2=0.9,  alpha=0.06, color=C_STAB)
    ax.fill_between(cg_arr, Vh_c,  y2=0.0,  alpha=0.06, color=C_CTRL)

    # ── 1. Design V̄_H dotted reference ─────────────────────────
    ax.axhline(Vh_design, color=C_VH, lw=1.0, ls=':',
               label=f'Design V̄_H = {Vh_design:.3f}')

    # ── 2. Actual CG range bar at design V̄_H (dark red) ────────
    ax.plot([p['cg_fwd'], p['cg_aft']], [Vh_design, Vh_design],
            color=C_CG, lw=3.5, solid_capstyle='round',
            label=f"CG range  [{p['cg_fwd']*100:.1f}%–{p['cg_aft']*100:.1f}% MAC]")
    ax.plot([p['cg_fwd'], p['cg_aft']], [Vh_design, Vh_design],
            'v', color=C_CG, ms=7)

    # ── 3. Required V̄_H from Part 1 (orange, CRJ-1000 CG range) ─
    ax.plot([cg_fwd_ref, cg_aft_ref], [Vh_req_ref, Vh_req_ref],
            color=C_REQ1, lw=3.5, solid_capstyle='round',
            label=f'Min V̄_H (Part 1 ref) = {Vh_req_ref:.3f}')
    ax.plot([cg_fwd_ref, cg_aft_ref], [Vh_req_ref, Vh_req_ref],
            '^', color=C_REQ1, ms=7)

    # ── 4. Required V̄_H for this config (purple, EXX CG range) ──
    #       On CRJ-1000 subplot this equals bar 3, so skip to avoid overlap
    if tag == 'CRJ-EXX':
        ax.plot([p['cg_fwd'], p['cg_aft']], [p['Vh_req'], p['Vh_req']],
                color=C_REQ2, lw=3.5, solid_capstyle='round',
                label=f'Min V̄_H (CRJ-EXX) = {p["Vh_req"]:.3f}')
        ax.plot([p['cg_fwd'], p['cg_aft']], [p['Vh_req'], p['Vh_req']],
                '^', color=C_REQ2, ms=7)

    # OEW CG line
    ax.axvline(p['oew'], color=C_OEW, lw=1.0, ls=':',
               label=f"OEW CG = {p['oew']*100:.1f}% MAC")

    ax.set_xlim(-0.1, 0.9)
    ax.set_ylim(0.0,  1.55)
    ax.xaxis.set_major_formatter(pct)
    ax.set_xlabel('CG position (fraction of MAC)', fontsize=10)
    ax.set_title(f'Scissor plot — {tag}', fontsize=12, fontweight='bold')
    ax.grid(True, ls='--', lw=0.4, alpha=0.5)
    ax.legend(fontsize=8.2, loc='upper left')

axes[0].set_ylabel('Tail volume coefficient  V̄_H  [–]', fontsize=10)
plt.tight_layout()
plt.savefig('scissor_plot_part2_vhreq.png', dpi=150, bbox_inches='tight')
plt.savefig('scissor_plot_part2_vhreq.pdf', dpi=150, bbox_inches='tight')
print("Saved.")
