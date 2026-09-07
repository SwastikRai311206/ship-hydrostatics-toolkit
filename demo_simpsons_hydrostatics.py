"""
demo_simpsons_hydrostatics.py
-------------------------------
Uses Simpson's 1/3rd Rule (from simpsons_rule.py) to compute the
waterplane area and moment of inertia from a table of half-ordinates,
and integrate these results into hydrostatic calculations for KB, BM, KM, and GM

Demonstrates how realistic hull forms, including tapering ends rather than a simple
rectangular shape, are handled in practical hydrostatic calculations.

"""

import numpy as np
import matplotlib.pyplot as plt
from simpsons_13_rule import waterplane_area_simpsons, moment_of_inertia_simpsons, lcf_simpsons

RHO_SEAWATER = 1025.0     # density of sea water, kg/m^3

# ---------------------------------------------------------------
# 1. Half-ordinates (half-breadths, m) at 7 equally spaced stations
#    along the waterline -- this barge tapers slightly at bow/stern
#    instead of being a perfect rectangle.
# ---------------------------------------------------------------

half_ordinates = [0.0, 4.5, 5.8, 6.0, 5.9, 4.8, 0.0]  # asymmetric hull
L = 40.0                                              # waterline length
n_intervals = len(half_ordinates) - 1
h = L/n_intervals                                     # station spacing, m

T = 1.8          # draft, m
KG = 2.0         # Center of gravity about keel, m

# ---------------------------------------------------------------
# 2. Waterplane properties via Simpson's Rule
# ---------------------------------------------------------------
A_w = waterplane_area_simpsons(half_ordinates, h)
I_CL = moment_of_inertia_simpsons(half_ordinates, h)
LCF_midship = lcf_simpsons(half_ordinates, h, reference="midship")

# ---------------------------------------------------------------
# 3. Displaced volume -- approximated here as waterplane area x
#    draft (valid for a wall-sided hull with vertical sides down
#    to the keel; a more advanced version would apply Simpson's
#    Rule again over sectional areas at each station).
# ---------------------------------------------------------------

V = A_w * T
mass = RHO_SEAWATER * V

# ------------------------------------------------------------------
# 4. Standard hydrostatics, using values derived from Simpson's Rule
# ------------------------------------------------------------------

KB = T/2                 # Centre of buoyancy
BM = I_CL/V              # Metacentric Radius
KM = KB + BM             # Height of Metacentre above the Keel
GM = KM - KG             # Metacentric Height

print("=" * 55)
print("Hydrostatics of a tapered hull using Simpsons rule")
print("=" * 55)
print(f"Half-ordinates (m)     : {half_ordinates}")
print(f"Station spacing (h)    : {h:.3f} m")
print(f"Waterplane area (A_w)  : {A_w:.3f} m^2")
print(f"Moment of inertia (I)  : {I_CL:.3f} m^4")
print(f"LCF from midship       : {LCF_midship:+.3f} m "
      f"({'fwd of' if LCF_midship >= 0 else 'aft of'} midship)")
print("-"*55)
print(f"Displaced Volume : {V:.3f} m^3")
print(f"Displaced mass   : {mass/1000:.3f} t ")
print(f"Centre of buoyancy (KB) : {KB:.3f} m")
print(f"Metacentric Radius (BM) : {BM:.3f} m")
print(f"Height of Metacentre above the Keel (KM) : {KM:.3f} m")
print(f"Centre of Gravity (KG)  : {KG:.3f} m")
print(f"Metacentric Height (GM) : {GM:.3f} m")
print("-"*55)
print("Result", "STABLE" if GM > 0 else "UUNSTABLE", "(GM", ">0)" if GM > 0 else "<= 0)")

# ---------------------------------------------------------------
# 5. Plotted the waterplane shape (plan view of the hull at the
#    waterline) using the half-ordinates, with LCF marked.
# ---------------------------------------------------------------

stations_x = np.arange(len(half_ordinates)) * h      #station position from AP
y = np.array(half_ordinates)

midship_x = (len(half_ordinates)-1)/2 * h
lcf_x_from_AP = midship_x + LCF_midship

fig, ax = plt.subplots(figsize=(9,4))

# Plot both port and starboard sides to show the full hull outline

ax.plot(stations_x, y, color = "tab:blue", lw=2)
ax.plot(stations_x, -y, color = "tab:blue", lw=2)
ax.fill_between(stations_x, y, -y, color = "tab:blue", alpha =0.15)

# Station markers

ax.scatter(stations_x, np.zeros_like(stations_x), color="brown", s=15, zorder=5)
for i,x in enumerate(stations_x):
    ax.annotate(str(i), (x,0), textcoords="offset points", xytext=(0,-14),
                ha="center", fontsize=8, color="brown")

# Midship and LCF refrence lines
ax.axvline(midship_x, color="black", linestyle=":", lw=1)
ax.text(midship_x, y.max()*1.15, "midship", ha="center", fontsize=8)

ax.axvline(lcf_x_from_AP,color="red", linestyle="--", lw=1.5)
ax.text(lcf_x_from_AP, -y.max()*1.25,
        f"LCF\n({LCF_midship:+.2f} m)", ha="center", fontsize =8)

ax.set_xlabel("Distance of waterline from AP (m)")
ax.set_ylabel("Half-breadth (m)")
ax.set_title("Waterplane Shape (from Half-Ordinates) with LCF")
ax.set_aspect("equal", adjustable="datalim")
ax.grid(alpha=0.3)
 
fig.tight_layout()
fig.savefig("waterplane_lcf.png", dpi=150)
print("Saved plot: waterplane_lcf.png")
plt.show()



