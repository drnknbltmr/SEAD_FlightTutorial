import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

# Shared Geometry
MAC    = 3.48       # m
LEMAC  = 19.2084    # m
MTOW   = 40_823.0   # kg

def to_mac(arm):
    """Convert longitudinal arm [m from nose] to % MAC."""
    return (arm - LEMAC) / MAC * 100.0

def calc_cg(w1, cg1, w2, cg2):
    """Mass-weighted combined centre of gravity."""
    return (w1 * cg1 + w2 * cg2) / (w1 + w2)

# CRJ-1000 Reference Parameters
REF = dict(
    OEW          = 23_178.0,   # kg
    OEW_CG       = 20.79,      # m
    PAX_WEIGHT   = 95.0,       # kg/pax
    CARGO_FWD    = 1_230.0,    # kg
    CARGO_AFT    = 1_245.0,    # kg
    CARGO_FWD_ARM = (9.89 + 12.73) / 2,
    CARGO_AFT_ARM = 29.0,      # m
    FUEL         = 5_670.0,    # kg
    FUEL_ARM     = 19.0,       # m
    rows         = np.linspace(6.0, 27.0, 25),
)

# CRJ-EXX Modified Parameters
BATT1_WEIGHT = 2_025.0     # kg - forward pack
BATT1_ARM    = 9.89        # m
BATT2_WEIGHT = 2_475.0     # kg - aft pack
BATT2_ARM    = 29.0        # m

# New OEW (without batteries)
OEW_EXX     = 22_041.69    # kg
OEW_EXX_CG  = 21.047       # m

# OEW + batteries (fixed starting point)
OEW_BATT       = 26_541.69  # kg
OEW_BATT_CG    = 20.937     # m

# Cargo hold volumes and proportional re-distribution
V_FWD_ORIG, V_AFT_ORIG = 10_600.0, 12_600.0
V_B1,       V_B2       = 1_300.0,  1_600.0

rem_fwd = V_FWD_ORIG - V_B1
rem_aft = V_AFT_ORIG - V_B2
cargo_total = REF["CARGO_FWD"] + REF["CARGO_AFT"]

CARGO_FWD_EXX     = cargo_total * rem_fwd / (rem_fwd + rem_aft)
CARGO_AFT_EXX     = cargo_total * rem_aft / (rem_fwd + rem_aft)
CARGO_FWD_ARM_EXX = 11.6
CARGO_AFT_ARM_EXX = 28.3

# Passengers: 22 rows, 88 pax
PAX_WEIGHT = 95.0
rows_exx   = np.linspace(6.0, 27.0, 25)[:22]

# Fuel adjusted for constant MTOW
total_payload_exx = 88 * PAX_WEIGHT + cargo_total
FUEL_EXX   = MTOW - OEW_BATT - total_payload_exx
FUEL_ARM   = 19.0

def simulate_ref(direction="forward"):
    """CRJ-1000 reference loading: cargo -> pax -> fuel."""
    curr_w, curr_cg = REF["OEW"], REF["OEW_CG"]
    weights, cgs = [curr_w], [curr_cg]

    cargo_seq = [
        (REF["CARGO_FWD"], REF["CARGO_FWD_ARM"]),
        (REF["CARGO_AFT"], REF["CARGO_AFT_ARM"]),
    ]
    if direction == "aft":
        cargo_seq = list(reversed(cargo_seq))

    for w, arm in cargo_seq:
        curr_cg = calc_cg(curr_w, curr_cg, w, arm)
        curr_w += w
        weights.append(curr_w); cgs.append(curr_cg)

    row_list = list(REF["rows"]) if direction == "forward" else list(reversed(REF["rows"]))
    for arm in row_list + row_list:
        curr_cg = calc_cg(curr_w, curr_cg, 2.0 * REF["PAX_WEIGHT"], arm)
        curr_w += 2.0 * REF["PAX_WEIGHT"]
        weights.append(curr_w); cgs.append(curr_cg)

    curr_cg = calc_cg(curr_w, curr_cg, REF["FUEL"], REF["FUEL_ARM"])
    curr_w += REF["FUEL"]
    weights.append(curr_w); cgs.append(curr_cg)

    return weights, [to_mac(c) for c in cgs]

def simulate_exx(direction="forward"):
    """CRJ-EXX loading: cargo -> pax -> fuel (starting from OEW+Batt)."""
    curr_w, curr_cg = OEW_BATT, OEW_BATT_CG
    weights, cgs = [curr_w], [curr_cg]

    cargo_seq = [
        (CARGO_FWD_EXX, CARGO_FWD_ARM_EXX),
        (CARGO_AFT_EXX, CARGO_AFT_ARM_EXX),
    ]
    if direction == "aft":
        cargo_seq = list(reversed(cargo_seq))

    for w, arm in cargo_seq:
        curr_cg = calc_cg(curr_w, curr_cg, w, arm)
        curr_w += w
        weights.append(curr_w); cgs.append(curr_cg)

    row_list = list(rows_exx) if direction == "forward" else list(reversed(rows_exx))
    for arm in row_list + row_list:
        curr_cg = calc_cg(curr_w, curr_cg, 2.0 * PAX_WEIGHT, arm)
        curr_w += 2.0 * PAX_WEIGHT
        weights.append(curr_w); cgs.append(curr_cg)

    curr_cg = calc_cg(curr_w, curr_cg, FUEL_EXX, FUEL_ARM)
    curr_w += FUEL_EXX
    weights.append(curr_w); cgs.append(curr_cg)

    return weights, [to_mac(c) for c in cgs]

# Run Simulations
ref_fwd_w, ref_fwd_cg = simulate_ref("forward")
ref_aft_w, ref_aft_cg = simulate_ref("aft")
ref_all_cg = ref_fwd_cg + ref_aft_cg
ref_cg_fwd, ref_cg_aft = min(ref_all_cg), max(ref_all_cg)
ref_mzfw = REF["OEW"] + REF["CARGO_FWD"] + REF["CARGO_AFT"] + 100 * REF["PAX_WEIGHT"]

exx_fwd_w, exx_fwd_cg = simulate_exx("forward")
exx_aft_w, exx_aft_cg = simulate_exx("aft")
exx_all_cg = exx_fwd_cg + exx_aft_cg
exx_cg_fwd, exx_cg_aft = min(exx_all_cg), max(exx_all_cg)
exx_mzfw = OEW_BATT + CARGO_FWD_EXX + CARGO_AFT_EXX + 88 * PAX_WEIGHT

# Print Summary
print("=" * 60)
print("CRJ-EXX Loading Diagram - Key Values")
print("=" * 60)
print(f"  New OEW (no batt)  : {OEW_EXX:.2f} kg  |  CG = {to_mac(OEW_EXX_CG):.1f} % MAC")
print(f"  OEW + Batteries    : {OEW_BATT:.2f} kg  |  CG = {to_mac(OEW_BATT_CG):.2f} % MAC")
print(f"  Cargo fwd / aft    : {CARGO_FWD_EXX:.1f} / {CARGO_AFT_EXX:.1f} kg")
print(f"  Pax                : 88 x {PAX_WEIGHT:.0f} kg = {88*PAX_WEIGHT:.0f} kg")
print(f"  Fuel               : {FUEL_EXX:.1f} kg")
print(f"  MZFW (EXX)         : {exx_mzfw:.1f} kg")
print(f"  MTOW               : {MTOW:.0f} kg")
print(f"  Most fwd CG (EXX)  : {exx_cg_fwd:.1f} % MAC")
print(f"  Most aft  CG (EXX) : {exx_cg_aft:.1f} % MAC")
print("-" * 60)
print(f"  Most fwd CG (REF)  : {ref_cg_fwd:.1f} % MAC")
print(f"  Most aft  CG (REF) : {ref_cg_aft:.1f} % MAC")
print("=" * 60)

# Plotting
fig, ax = plt.subplots(figsize=(11, 8))

# CRJ-1000 Reference
ax.plot(ref_fwd_cg, ref_fwd_w, color="#7FB3D3", linewidth=1.2, linestyle="--", label="CRJ-1000 Reference")
ax.plot(ref_aft_cg, ref_aft_w, color="#F1948A", linewidth=1.2, linestyle="--")
ax.scatter(to_mac(REF["OEW_CG"]), REF["OEW"], color="#7FB3D3", zorder=5, s=55)
ax.annotate("OEW (CRJ-1000)", xy=(to_mac(REF["OEW_CG"]), REF["OEW"]), xytext=(to_mac(REF["OEW_CG"]) + 5, REF["OEW"] - 1_100),
            fontsize=8, color="#5D8AA8", arrowprops=dict(arrowstyle="->", color="#5D8AA8"))

# CRJ-EXX
ax.plot(exx_fwd_cg, exx_fwd_w, color="#1A5276", linewidth=2.0, marker="o", markersize=3.5, label="CRJ-EXX Loading")
ax.plot(exx_aft_cg, exx_aft_w, color="#C0392B", linewidth=2.0, marker="o", markersize=3.5)
ax.axvline(exx_cg_fwd, color="#1A5276", linestyle="--", linewidth=1.5, alpha=0.8, label=f"Most Fwd ({exx_cg_fwd:.1f}%)")
ax.axvline(exx_cg_aft, color="#C0392B", linestyle="--", linewidth=1.5, alpha=0.8, label=f"Most Aft ({exx_cg_aft:.1f}%)")

# Key Data Points
ax.scatter(to_mac(OEW_EXX_CG), OEW_EXX, color="black", zorder=7, s=70, marker="D")
ax.annotate("OEW (EXX, no batt)", xy=(to_mac(OEW_EXX_CG), OEW_EXX), xytext=(to_mac(OEW_EXX_CG) + 4, OEW_EXX - 1_200),
            fontsize=8, arrowprops=dict(arrowstyle="->", color="black"))

ax.scatter(to_mac(OEW_BATT_CG), OEW_BATT, color="darkorange", zorder=7, s=90, marker="*")
ax.annotate("OEW+Batt (EXX)", xy=(to_mac(OEW_BATT_CG), OEW_BATT), xytext=(to_mac(OEW_BATT_CG) - 18, OEW_BATT + 800),
            fontsize=8, color="darkorange", arrowprops=dict(arrowstyle="->", color="darkorange"))

# Reference Lines
for weight, label in [
    (REF["OEW"],  f"OEW = {REF['OEW']:.0f} kg"),
    (OEW_EXX,     f"OEW EXX = {OEW_EXX:.0f} kg"),
    (OEW_BATT,    f"OEW+Batt = {OEW_BATT:.0f} kg"),
    (ref_mzfw,    f"MZFW (CRJ-1000) = {ref_mzfw:.0f} kg"),
    (exx_mzfw,    f"MZFW (EXX) = {exx_mzfw:.1f} kg"),
    (MTOW,        f"MTOW = {MTOW:.0f} kg"),
]:
    ax.axhline(y=weight, color="gray", linestyle=":", linewidth=0.8, alpha=0.45)
    ax.text(97, weight + 90, label, fontsize=7, color="gray", ha="right", va="bottom")

ax.set_xlim(-5, 100)
ax.set_ylim(20_000, 45_000)
ax.set_xlabel("Centre of Gravity [% MAC]")
ax.set_ylabel("Aircraft Gross Weight [kg]")
ax.grid(True, linestyle=":", alpha=0.40)
ax.legend(fontsize=8, loc="upper left")

plt.tight_layout()
plt.savefig("loading_diagram_part2.png", dpi=300)
print("\nSaved -> loading_diagram_part2.png")
plt.show()
