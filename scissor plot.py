"""
Scissor Plot Generator — AE3211-I Aircraft Tutorial
Group 45 | CRJ-1000 reference aircraft
=====================================================
"""

import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
#  GEOMETRY
# =============================================================================
S      = 77.39    # [m²]  wing area
b      = 26.17    # [m]   wing span
c_mac  = 3.67     # [m]   MAC
S_h    = 15.91    # [m²]  horizontal tail area
AR_w   = 8.85
AR_h   = 4.58
LEMAC  = 19.2084  # [m]   from nose

x_ac_nose_wing = LEMAC + 0.25 * c_mac
x_ac_nose_tail = 37.097
l_h = x_ac_nose_tail - x_ac_nose_wing   # 16.971 m

# =============================================================================
#  AERODYNAMIC PARAMETERS
# =============================================================================

# ── From lectures (no calculation needed) ───────────────────────
eta_h    = 1.0     # T-tail  (Lec7 s42)
CL_h_max = -0.8    # adjustable stabiliser  (Lec8 s17)
SM       = 0.05    # 5% MAC  (Lec7 s27)
x_ac_w   = 0.25

# ── Lift curve slopes ────────────────────────────────────────────
CL_alpha_w_cruise   = 6.149   # /rad  (DATCOM, M=0.78)
CL_alpha_w_approach = 4.752   # /rad  (DATCOM, M=0.20)
CL_alpha_wf_cruise   = 6.094  # /rad  (wing+fuselage correction, cruise)
CL_alpha_wf_approach = 4.710  # /rad  (wing+fuselage correction, approach)
CL_alpha_h_cruise    = 4.480  # /rad  (tail DATCOM, cruise)
CL_alpha_h_approach  = 3.781  # /rad  (tail DATCOM, approach)

# ── AC contributions ─────────────────────────────────────────────
x_ac_f = -0.1695  # fuselage (Torenbeek two terms)
x_ac_n =  0.0274  # nacelles (Torenbeek, rear-mounted)
x_ac   = x_ac_w + x_ac_f + x_ac_n   # = 0.1079  (10.79% MAC)

# ── FIX 1: Downwash gradient — geometric formula ─────────────────
# dε/dα = 4·cos(Λ_0.25c) / (AR + 2·cos(Λ_0.25c))
# This is speed-independent (incompressible potential-flow) and is
# used for BOTH the stability line (cruise) and controllability line (approach).
# Previous Slingerland implementation used the M=0.78 compressibility-
# corrected CLα_w, which inflated dε/dα to ~0.80 (physically too high).
sweep_qc_rad = np.radians(23.59)
d_epsilon_d_alpha = 4 * np.cos(sweep_qc_rad) / (AR_w + 2 * np.cos(sweep_qc_rad))
print(f"dε/dα = 4·cos(23.59°)/(8.85+2·cos(23.59°)) = {d_epsilon_d_alpha:.4f}")
print(f"  Sanity check 4/(AR+2) = {4/(AR_w+2):.4f}  ✓ close")

# ── Pitching moments ─────────────────────────────────────────────
Cm_ac_wf_cruise   =  0.1277   # clean wings, cruise       (stability line)
Cm_ac_wf_approach = -0.4754   # full flaps, approach      (controllability line)

# ── Max lift ─────────────────────────────────────────────────────
CL_max_wf = 2.840   # approach, full HLD

# ── CG range from loading diagram ────────────────────────────────
cg_fwd = 0.098    # fill in
cg_aft = 0.56   # fill in

# =============================================================================
#  FUNCTIONS
# =============================================================================

def tail_volume(S_h, l_h, S, c_mac):
    return (S_h * l_h) / (S * c_mac)


def neutral_point(x_ac, CL_alpha_wf, CL_alpha_h,
                  d_epsilon_d_alpha, eta_h, Vh):
    """
    Lec7 s27 — neutral point as fraction of MAC.
    x_ac = combined AC of aircraft-less-tail (wing + fuselage + nacelles)
    """
    tail = (CL_alpha_h / CL_alpha_wf) * eta_h * (1 - d_epsilon_d_alpha) * Vh
    return x_ac + tail


def stability_line(Vh_range, x_ac, CL_alpha_wf, CL_alpha_h,
                   d_epsilon_d_alpha, eta_h, SM):
    """
    Most-aft allowable CG for each Vh. CRUISE aerodynamics (Lec8 s24).
    """
    return np.array([
        neutral_point(x_ac, CL_alpha_wf, CL_alpha_h,
                      d_epsilon_d_alpha, eta_h, Vh) - SM
        for Vh in Vh_range
    ])


def neutral_stability_line(Vh_range, x_ac, CL_alpha_wf, CL_alpha_h,
                           d_epsilon_d_alpha, eta_h):
    return stability_line(Vh_range, x_ac, CL_alpha_wf, CL_alpha_h,
                          d_epsilon_d_alpha, eta_h, SM=0.0)


def controllability_line(cg_range, x_ac, Cm_ac_wf, CL_max_wf,
                         CL_alpha_h, CL_h_max, eta_h,
                         d_epsilon_d_alpha):
    """
    Minimum Vh to trim at each CG position. APPROACH aerodynamics (Lec8 s25).

    FIX 2 — Correct trim equation (Lec8 s9/16):
      Cm_cg = Cm_ac_wf + CL_max_wf*(xcg - x_ac)/c - Vh*eta_h*CLh*(1-de/da) = 0

    Rearranged:
      Vh_min = -(Cm_ac_wf + CL_max_wf*(xcg - x_ac)) / (eta_h*|CLh_max|*(1-de/da))

    Previous version omitted the x_ac term, shifting the line ~22% MAC left.
    """
    Vh_min = []
    for xcg in cg_range:
        num = Cm_ac_wf + CL_max_wf * (xcg - x_ac)   # ← corrected moment arm
        den = eta_h * abs(CL_h_max) * (1 - d_epsilon_d_alpha)
        Vh_min.append(np.nan if den == 0 else -num / den)
    return np.array(Vh_min)


# =============================================================================
#  STATUS REPORT
# =============================================================================

Vh_design = tail_volume(S_h, l_h, S, c_mac)

# Preview the corrected scissor region at design Vh
xnp  = neutral_point(x_ac, CL_alpha_wf_cruise, CL_alpha_h_cruise,
                     d_epsilon_d_alpha, eta_h, Vh_design)
denom_ctrl = eta_h * abs(CL_h_max) * (1 - d_epsilon_d_alpha)
xcg_ctrl   = ((-Vh_design * denom_ctrl - Cm_ac_wf_approach) / CL_max_wf) + x_ac

print(f"\n{'='*55}")
print(f"  CRJ-1000 Scissor Plot — Group 45  (corrected)")
print(f"{'='*55}")
print(f"  Tail volume            V_H       = {Vh_design:.4f}")
print(f"  Aircraft-less-tail AC  x_ac      = {x_ac:.4f}  ({x_ac*100:.2f}% MAC)")
print(f"  Downwash gradient      dε/dα     = {d_epsilon_d_alpha:.4f}  (geometric)")
print(f"  Tail effectiveness     1-dε/dα   = {1-d_epsilon_d_alpha:.4f}")
print(f"\n  At design V_H = {Vh_design:.4f}:")
print(f"    Neutral point                  = {xnp*100:.1f}% MAC")
print(f"    Stability aft-CG limit         = {(xnp-SM)*100:.1f}% MAC")
print(f"    Controllability fwd-CG limit   = {xcg_ctrl*100:.1f}% MAC")
print(f"    Feasible CG band               = {xcg_ctrl*100:.1f}% → {(xnp-SM)*100:.1f}% MAC")
print(f"    OEW CG = 45.45%  {'✓ inside' if xcg_ctrl < 0.4545 < xnp-SM else '⚠ outside'} feasible band")
print(f"{'='*55}")

if cg_fwd is None or cg_aft is None:
    print("\n  Fill in cg_fwd and cg_aft from your loading diagram to generate plot.")
    exit()

# =============================================================================
#  PLOT
# =============================================================================

Vh_range = np.linspace(0.0, 1.4, 500)
cg_range = np.linspace(-0.1, 0.9, 500)

cg_neutral = neutral_stability_line(
    Vh_range, x_ac, CL_alpha_wf_cruise, CL_alpha_h_cruise,
    d_epsilon_d_alpha, eta_h)

cg_stab = stability_line(
    Vh_range, x_ac, CL_alpha_wf_cruise, CL_alpha_h_cruise,
    d_epsilon_d_alpha, eta_h, SM)

Vh_ctrl = controllability_line(
    cg_range, x_ac, Cm_ac_wf_approach, CL_max_wf,
    CL_alpha_h_approach, CL_h_max, eta_h, d_epsilon_d_alpha)

fig, ax = plt.subplots(figsize=(9, 6))

ax.plot(cg_neutral, Vh_range,
        color='gray', linewidth=1.5, linestyle='--',
        label='Neutral stability line  (cruise)')
ax.plot(cg_stab, Vh_range,
        color='#185FA5', linewidth=2,
        label=f'Stability line  (SM={SM*100:.0f}% MAC, cruise)')
ax.plot(cg_range, Vh_ctrl,
        color='#1D9E75', linewidth=2,
        label='Controllability line  (approach, full flaps)')

ax.fill_betweenx(Vh_range, cg_stab, x2=0.9,
                 alpha=0.07, color='#185FA5')
ax.fill_between(cg_range, Vh_ctrl, y2=0,
                alpha=0.07, color='#1D9E75')

ax.axhline(Vh_design, color='#888780', linewidth=1, linestyle=':',
           label=f'Design V̄_H = {Vh_design:.3f}')

ax.plot([cg_fwd, cg_aft], [Vh_design, Vh_design],
        color='#A32D2D', linewidth=3, solid_capstyle='round',
        label=f'CG range  [{cg_fwd*100:.1f}%–{cg_aft*100:.1f}% MAC]')
ax.plot([cg_fwd, cg_aft], [Vh_design, Vh_design],
        'v', color='#A32D2D', markersize=7)

ax.axvline(0.4545, color='orange', linewidth=1, linestyle=':',
           label='OEW CG = 45.45% MAC')

ax.set_xlabel('CG position  (fraction of MAC from LEMAC)', fontsize=11)
ax.set_ylabel('Tail volume coefficient  V̄_H  [–]',        fontsize=11)
ax.set_title('Scissor plot — CRJ-1000  (Group 45)',        fontsize=13)
ax.set_xlim(-0.1, 0.9)
ax.set_ylim(0.0,  1.2)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x*100:.0f}%'))
ax.grid(True, linestyle='--', linewidth=0.4, alpha=0.5)
ax.legend(fontsize=9, loc='lower left')

plt.tight_layout()
plt.savefig('scissor_plot.pdf', dpi=150, bbox_inches='tight')
plt.savefig('scissor_plot.png', dpi=150, bbox_inches='tight')
print("\nSaved:  scissor_plot.pdf  and  scissor_plot.png")
plt.show()