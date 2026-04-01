import matplotlib.pyplot as plt
import numpy as np

# Aircraft Data
OEW = 23_178.0      # kg - Operational Empty Weight
OEW_CG = 20.59      # m - OEW center of gravity
MAC = 3.48
          # m - Mean Aerodynamic Chord
LEMAC = 19.2084     # m - Leading Edge of MAC

# Payload & Fuel
PAX_WEIGHT = 95.0
CARGO_FRONT_MAX = 1_230.0
CARGO_AFT_MAX = 1_245.0
FUEL_MAX = 5_670.0

# Longitudinal Arms
CARGO_FRONT_ARM = (9.89 + 12.73) / 2
CARGO_AFT_ARM = 29.0
FUEL_ARM = 18.0
rows = np.linspace(6, 27.0, 25)

def to_mac(arm):
    """Convert arm from nose [m] to % MAC."""
    return (arm - LEMAC) / MAC * 100.0

def calc_cg(w1, cg1, w2, cg2):
    """Calculate combined CG of two masses."""
    return (w1 * cg1 + w2 * cg2) / (w1 + w2)

def simulate_loading(direction="forward"):
    """Simulate loading sequence: Cargo -> Passengers (Window then Aisle) -> Fuel."""
    weights, cgs = [OEW], [OEW_CG]
    curr_w, curr_cg = OEW, OEW_CG

    # Cargo
    cargo_order = [(CARGO_FRONT_MAX, CARGO_FRONT_ARM), (CARGO_AFT_MAX, CARGO_AFT_ARM)]
    if direction == "aft":
        cargo_order = list(reversed(cargo_order))

    for w, arm in cargo_order:
        curr_cg = calc_cg(curr_w, curr_cg, w, arm)
        curr_w += w
        weights.append(curr_w)
        cgs.append(curr_cg)

    # Passengers
    row_list = list(rows) if direction == "forward" else list(reversed(rows))
    for arm in row_list + row_list:  # Windows first, then aisles
        curr_cg = calc_cg(curr_w, curr_cg, 2.0 * PAX_WEIGHT, arm)
        curr_w += 2.0 * PAX_WEIGHT
        weights.append(curr_w)
        cgs.append(curr_cg)

    # Fuel
    curr_cg = calc_cg(curr_w, curr_cg, FUEL_MAX, FUEL_ARM)
    curr_w += FUEL_MAX
    weights.append(curr_w)
    cgs.append(curr_cg)

    return weights, [to_mac(c) for c in cgs]

# Generate Loading Envelopes
fwd_w, fwd_cg = simulate_loading("forward")
aft_w, aft_cg = simulate_loading("aft")

all_cg = fwd_cg + aft_cg
cg_fwd, cg_aft = min(all_cg), max(all_cg)
mtow = fwd_w[-1]
mzfw = OEW + CARGO_FRONT_MAX + CARGO_AFT_MAX + 100 * PAX_WEIGHT

# Plotting
fig, ax = plt.subplots(figsize=(10, 7))

ax.plot(fwd_cg, fwd_w, "b-o", markersize=3, linewidth=1.5, label="Forward Loading")
ax.plot(aft_cg, aft_w, "r-o", markersize=3, linewidth=1.5, label="Aft Loading")

# Annotations & Reference Lines
ax.scatter(fwd_cg[0], fwd_w[0], color="black", zorder=6, s=70)
ax.annotate("OEW", xy=(fwd_cg[0], fwd_w[0]), xytext=(fwd_cg[0] + 1.5, fwd_w[0] - 600),
            fontsize=9, arrowprops=dict(arrowstyle="->", color="black"))

ax.axvline(x=cg_fwd, color="green", linestyle="--", linewidth=1.5, alpha=0.7, label=f"Most Fwd ({cg_fwd:.1f}%)")
ax.axvline(x=cg_aft, color="purple", linestyle="--", linewidth=1.5, alpha=0.7, label=f"Most Aft ({cg_aft:.1f}%)")
ax.axvline(x=to_mac(OEW_CG), color="gray", linestyle="-.", alpha=0.6)

for weight in [OEW, mzfw, mtow]:
    ax.axhline(y=weight, color="gray", linestyle=":", alpha=0.35)

ax.set_xlim(0, 100)
ax.set_ylim(20_000, 45_000)
ax.set_xlabel("Center of Gravity [% MAC]")
ax.set_ylabel("Aircraft Gross Weight [kg]")
ax.grid(True, which="both", linestyle=":", alpha=0.45)
ax.legend(fontsize=8, loc="upper right")

plt.tight_layout()
plt.savefig("loading_diagram.png", dpi=300)
plt.show()
