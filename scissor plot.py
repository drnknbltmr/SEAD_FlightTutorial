import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
#  GEOMETRY & AERO (unchanged from Part 1)
# =============================================================================
S, b, c_mac, S_h = 77.39, 26.17, 3.67, 15.91
AR_w   = 8.85
LEMAC  = 19.2084
x_ac_nose_wing = LEMAC + 0.25 * c_mac
l_h    = 37.097 - x_ac_nose_wing

eta_h, CL_h_max, SM = 1.0, -0.8, 0.05
x_ac_w, x_ac_f, x_ac_n = 0.25, -0.1695, 0.0274
x_ac = x_ac_w + x_ac_f + x_ac_n

sweep_qc_rad      = np.radians(26.9)
d_epsilon_d_alpha = 4*np.cos(sweep_qc_rad) / (AR_w + 2*np.cos(sweep_qc_rad))

CL_alpha_wf_cruise   = 6.094
CL_alpha_wf_approach = 4.710
CL_alpha_h_cruise    = 4.480
CL_alpha_h_approach  = 3.781
Cm_ac_wf_approach    = -0.4754
CL_max_wf            = 2.840

cg_fwd, cg_aft = 0.052, 0.575

# =============================================================================
#  FUNCTIONS
# =============================================================================
def tail_volume(Sh, lh, S, cmac): return (Sh*lh)/(S*cmac)

def neutral_point(x_ac, CLa_wf, CLa_h, de_da, eta, Vh):
    return x_ac + (CLa_h/CLa_wf)*eta*(1-de_da)*Vh

def stability_line(Vh_arr, x_ac, CLa_wf, CLa_h, de_da, eta, SM):
    return np.array([neutral_point(x_ac,CLa_wf,CLa_h,de_da,eta,v)-SM for v in Vh_arr])

def neutral_stability_line(Vh_arr, x_ac, CLa_wf, CLa_h, de_da, eta):
    return stability_line(Vh_arr, x_ac, CLa_wf, CLa_h, de_da, eta, SM=0.0)

def controllability_line(cg_arr, x_ac, Cm_ac, CLmax, CLa_h, CLh_max, eta, de_da):
    den = eta*abs(CLh_max)*(1-de_da)
    return np.array([-(Cm_ac + CLmax*(xcg-x_ac))/den for xcg in cg_arr])

# =============================================================================
#  MINIMUM REQUIRED V_H
# =============================================================================
Vh_req_stab = (cg_aft + SM - x_ac) / (
    (CL_alpha_h_cruise/CL_alpha_wf_cruise)*eta_h*(1-d_epsilon_d_alpha))
Vh_req_ctrl = -(Cm_ac_wf_approach + CL_max_wf*(cg_fwd - x_ac)) / (
    eta_h*abs(CL_h_max)*(1-d_epsilon_d_alpha))
Vh_required = max(Vh_req_stab, Vh_req_ctrl)

# =============================================================================
#  PLOT
# =============================================================================
Vh_design = tail_volume(S_h, l_h, S, c_mac)
Vh_arr    = np.linspace(0.0, 1.4, 500)
cg_arr    = np.linspace(-0.1, 0.9, 500)

cg_neutral = neutral_stability_line(Vh_arr, x_ac, CL_alpha_wf_cruise, CL_alpha_h_cruise, d_epsilon_d_alpha, eta_h)
cg_stab    = stability_line(Vh_arr, x_ac, CL_alpha_wf_cruise, CL_alpha_h_cruise, d_epsilon_d_alpha, eta_h, SM)
Vh_ctrl    = controllability_line(cg_arr, x_ac, Cm_ac_wf_approach, CL_max_wf, CL_alpha_h_approach, CL_h_max, eta_h, d_epsilon_d_alpha)

# Where does each end of the required bar touch the lines?
# Left end (cg_fwd) touches controllability line at Vh_required
# Right end (cg_aft) touches stability line at Vh_required
# (by construction of Vh_required)

fig, ax = plt.subplots(figsize=(9, 6))

ax.plot(cg_neutral, Vh_arr, color='gray',    lw=1.5, ls='--', label='Neutral stability line  (cruise)')
ax.plot(cg_stab,    Vh_arr, color='#185FA5', lw=2,            label='Stability line  (SM=5% MAC, cruise)')
ax.plot(cg_arr,     Vh_ctrl, color='#1D9E75', lw=2,           label='Controllability line  (approach, full flaps)')

ax.fill_betweenx(Vh_arr, cg_stab, x2=0.9,  alpha=0.07, color='#185FA5')
ax.fill_between(cg_arr, Vh_ctrl,  y2=0,    alpha=0.07, color='#1D9E75')

ax.axhline(Vh_design, color='#888780', lw=1.0, ls=':',
           label=f'Design V̄_H = {Vh_design:.3f}')

# ── Actual CG range bar (dark red) ───────────────────────────────
ax.plot([cg_fwd, cg_aft], [Vh_design, Vh_design],
        color='#A32D2D', lw=3, solid_capstyle='round',
        label=f'CG range  [{cg_fwd*100:.1f}%–{cg_aft*100:.1f}% MAC]')
ax.plot([cg_fwd, cg_aft], [Vh_design, Vh_design], 'v', color='#A32D2D', ms=7)

# ── Required CG range bar — translated up to Vh_required ─────────
# Endpoints land exactly on the controllability (left) and stability (right) lines
ax.plot([cg_fwd, cg_aft], [Vh_required, Vh_required],
        color='#C94B00', lw=3, solid_capstyle='round',
        label=f'Min required V̄_H = {Vh_required:.3f}')
ax.plot([cg_fwd, cg_aft], [Vh_required, Vh_required], '^', color='#C94B00', ms=7)

ax.axvline(0.4545, color='orange', lw=1.0, ls=':', label='OEW CG = 45.45% MAC')

ax.set_xlabel('CG position  (fraction of MAC from LEMAC)', fontsize=11)
ax.set_ylabel('Tail volume coefficient  V̄_H  [–]',        fontsize=11)
ax.set_title('Scissor plot — CRJ-1000  (Group 45)',        fontsize=13)
ax.set_xlim(-0.1, 0.9)
ax.set_ylim(0.0,  1.4)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x*100:.0f}%'))
ax.grid(True, ls='--', lw=0.4, alpha=0.5)
ax.legend(fontsize=9, loc='lower left')

plt.tight_layout()
plt.savefig('scissor_plot_vhreq2.png', dpi=150, bbox_inches='tight')
plt.savefig('scissor_plot_vhreq2.pdf', dpi=150, bbox_inches='tight')
print(f"Vh_required = {Vh_required:.4f}  (stab={Vh_req_stab:.4f}, ctrl={Vh_req_ctrl:.4f})")
print("Saved.")