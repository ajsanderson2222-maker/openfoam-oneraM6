"""Top-down Cp map on wing upper surface — matplotlib tricontourf."""
import pyvista as pv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

pv.OFF_SCREEN = True

case_dir = "/home/andrew/openfoam_jobs/openfoam-oneraM6/simulation"
out_dir  = "/home/andrew/openfoam_jobs/openfoam-oneraM6"

p_inf     = 80510.081
rho_inf   = p_inf / (287.0 * 255.556)
U_inf_mag = (268.654**2 + 14.362**2)**0.5
q_inf     = 0.5 * rho_inf * U_inf_mag**2

reader = pv.OpenFOAMReader(os.path.join(case_dir, "case.foam"))
reader.set_active_time_value(2000)
reader.enable_all_patch_arrays()
data   = reader.read()
wing   = data["boundary"]["wing"]

pts    = wing.cell_centers().points
x_w    = pts[:, 0]   # chordwise
y_w    = pts[:, 1]   # thickness
z_w    = pts[:, 2]   # spanwise
p_w    = wing.cell_data["p"]
Cp     = (p_w - p_inf) / q_inf

b_semi = z_w.max()
chord  = x_w.max() - x_w.min()
x_le   = x_w.min()
x_norm = (x_w - x_le) / chord      # x/c
eta    = z_w / b_semi               # 2y/b (spanwise fraction)

# Upper surface only (y > 0)
upper  = y_w >= 0.0

fig, ax = plt.subplots(figsize=(10, 7))
tc = ax.tricontourf(x_norm[upper], eta[upper], Cp[upper],
                    levels=50, cmap="coolwarm_r", vmin=-1.2, vmax=0.3)
cb = fig.colorbar(tc, ax=ax, label="$C_p$", orientation="vertical",
                  fraction=0.03, pad=0.02, extend="neither")
ax.set_xlabel("x/c")
ax.set_ylabel("2y/b  (spanwise fraction)")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_title("ONERA M6 Upper Surface $C_p$  (Ma=0.84, α=3.06°)")
fig.tight_layout()
fig.savefig(os.path.join(out_dir, "wing_cp_topdown.png"), dpi=150)
plt.close(fig)
print("Saved wing_cp_topdown.png")
