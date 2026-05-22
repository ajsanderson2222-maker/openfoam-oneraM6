"""
Post-process hisa ONERA M6 simulation.
Generates:
  - Cp distribution on wing surface vs 2y/b (spanwise stations)
  - Mach contours on symmetry plane (z=0)
  - 3D wing surface colored by Cp
"""
import pyvista as pv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

pv.OFF_SCREEN = True

case_dir = "/home/andrew/openfoam_jobs/openfoam-oneraM6/simulation"
out_dir  = "/home/andrew/openfoam_jobs/openfoam-oneraM6"

foam_file = os.path.join(case_dir, "case.foam")
reader = pv.OpenFOAMReader(foam_file)
reader.set_active_time_value(2000)

# Enable boundary patches
reader.enable_all_patch_arrays()

mesh_data = reader.read()
internal  = mesh_data["internalMesh"]

b = internal.bounds
print(f"Bounds: x=[{b[0]:.3f},{b[1]:.3f}] y=[{b[2]:.3f},{b[3]:.3f}] z=[{b[4]:.3f},{b[5]:.3f}]")

# ── Freestream ────────────────────────────────────────────────────────────────
gamma  = 1.4
R      = 287.0
T_inf  = 255.556
p_inf  = 80510.081
U_inf  = np.array([268.654, 14.362, 0.0])
q_inf  = 0.5 * (p_inf / (R * T_inf)) * np.dot(U_inf, U_inf)  # dynamic pressure

# ── Wing surface patch ─────────────────────────────────────────────────────────
wing = mesh_data["boundary"]["wing"]
print(f"Wing patch: n_cells={wing.n_cells}, n_points={wing.n_points}")

print("Wing fields:", wing.cell_data.keys() if wing.n_cells > 0 else wing.point_data.keys())

# Get pressure on wing
if wing.n_cells > 0:
    p_wing = wing.cell_data["p"]
    pts_w  = wing.cell_centers().points
else:
    p_wing = wing.point_data["p"]
    pts_w  = wing.points

Cp_wing = (p_wing - p_inf) / q_inf

# ONERA M6 mesh orientation: x=chordwise, y=thickness, z=spanwise
x_w = pts_w[:, 0]   # chordwise
y_w = pts_w[:, 1]   # thickness (upper/lower surface)
z_w = pts_w[:, 2]   # spanwise

# Wing semi-span in z
b_semi = z_w.max()
print(f"Wing semi-span b/2 = {b_semi:.4f} m  (full span = {2*b_semi:.4f} m)")

# ── Load Schmitt & Charpin (1979) experimental Cp data ────────────────────────
# Sections 1-7 → eta = 0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99
exp_eta = [0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99]

def load_exp_cp(dat_file):
    """Parse Tecplot POINT format: Section Tap X/L Z/L CP.
    Z/L > 0 → upper surface; Z/L < 0 → lower surface."""
    sections = {}
    cur_sec = None
    with open(dat_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("ZONE"):
                # extract section number from leading data token
                cur_sec = None
            elif line.startswith("VARIABLES") or line.startswith("TITLE"):
                continue
            else:
                vals = line.split()
                if len(vals) == 5:
                    try:
                        sec = int(float(vals[0]))
                        xc  = float(vals[2])
                        zl  = float(vals[3])   # thickness coord — sign → surface
                        cp  = float(vals[4])
                        if sec not in sections:
                            sections[sec] = {"upper": [], "lower": []}
                        if zl >= 0:
                            sections[sec]["upper"].append((xc, cp))
                        else:
                            sections[sec]["lower"].append((xc, cp))
                    except ValueError:
                        pass
    return sections

exp_dat = os.path.join(out_dir, "onera_m6_exp_2308.dat")
exp_sections = load_exp_cp(exp_dat) if os.path.exists(exp_dat) else {}

# ── Cp spanwise stations — match ONERA M6 experiment ─────────────────────────
# TMR stations: 2y/b = 0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99
eta_stations = [0.20, 0.44, 0.65, 0.80, 0.90, 0.96, 0.99]

# ONERA M6 is swept ~30°: x_LE and x_TE vary linearly with z.
# Fit linear LE/TE from upper surface cell data.
upper_all = y_w >= 0
z_bins = np.linspace(0, b_semi, 30)
x_le_fit, x_te_fit, z_mid = [], [], []
for i in range(len(z_bins)-1):
    m = upper_all & (z_w >= z_bins[i]) & (z_w < z_bins[i+1])
    if m.sum() > 0:
        x_le_fit.append(x_w[m].min())
        x_te_fit.append(x_w[m].max())
        z_mid.append((z_bins[i]+z_bins[i+1])/2)
le_coeffs = np.polyfit(z_mid, x_le_fit, 1)
te_coeffs = np.polyfit(z_mid, x_te_fit, 1)

# Per-cell local x/c using swept LE/TE
x_le_cell  = le_coeffs[1] + le_coeffs[0] * z_w
x_te_cell  = te_coeffs[1] + te_coeffs[0] * z_w
chord_cell = x_te_cell - x_le_cell
xc_all     = (x_w - x_le_cell) / chord_cell

# Wing patch has stride=32 in spanwise direction (structured O-grid).
# First 32 cell indices span the full semi-span — use these as the spanwise grid.
# Select exact j-index (no z-band) for clean single-layer sampling.
span_z   = z_w[0:32]           # z coordinate of each spanwise station j=0..31
span_eta = span_z / b_semi

from scipy.signal import savgol_filter

def smooth_cp(x_norm, cp, window=7, poly=3):
    """Sort by x/c then Savitzky-Golay smooth."""
    o = np.argsort(x_norm)
    x_s, c_s = x_norm[o], cp[o]
    w = min(window, len(c_s) if len(c_s) % 2 == 1 else len(c_s) - 1)
    w = max(w, poly + 2) if w > poly else poly + 2
    if w % 2 == 0:
        w += 1
    if len(c_s) >= w:
        c_s = savgol_filter(c_s, w, poly)
    return x_s, c_s

for idx, eta in enumerate(eta_stations):
    # Pick nearest j-index in the structured spanwise grid
    j = int(np.argmin(np.abs(span_eta - eta)))
    # Select all cells at this spanwise index: j, j+32, j+64, ...
    indices = np.arange(j, len(z_w), 32)
    mask = np.zeros(len(z_w), dtype=bool)
    mask[indices] = True

    fig, ax = plt.subplots(figsize=(8, 6))

    if mask.sum() == 0:
        ax.set_title(f"2y/b = {eta:.2f}  (no data)")
    else:
        cp_loc = Cp_wing[mask]
        y_loc  = y_w[mask]
        x_norm = xc_all[mask]   # sweep-corrected x/c
        actual_eta = span_eta[j]

        upper = y_loc >= 0.0
        lower = y_loc < 0.0

        if upper.sum() > 0:
            xb, cb_ = smooth_cp(x_norm[upper], cp_loc[upper])
            ax.plot(xb, cb_, "b-", lw=2, label="CFD upper")
        if lower.sum() > 0:
            xb, cb_ = smooth_cp(x_norm[lower], cp_loc[lower])
            ax.plot(xb, cb_, "r-", lw=2, label="CFD lower")
        if upper.sum() == 0 and lower.sum() == 0:
            xb, cb_ = smooth_cp(x_norm, cp_loc)
            ax.plot(xb, cb_, "b-", lw=2)

        sec_id = idx + 1
        if sec_id in exp_sections:
            edata = exp_sections[sec_id]
            if edata["upper"]:
                ex, ec = zip(*sorted(edata["upper"]))
                ax.plot(ex, ec, "bs", ms=6, mfc="none", lw=1.2, label="Exp upper")
            if edata["lower"]:
                ex, ec = zip(*sorted(edata["lower"]))
                ax.plot(ex, ec, "r^", ms=6, mfc="none", lw=1.2, label="Exp lower")

        ax.set_title(f"2y/b = {eta:.2f}  (mesh η={actual_eta:.3f})", fontsize=12)

    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("x/c", fontsize=11)
    ax.set_ylabel("$C_p$", fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.suptitle("ONERA M6 — $C_p$  (Ma=0.84, α=3.06°)", fontsize=12)
    fig.tight_layout()
    fname = os.path.join(out_dir, f"cp_station_{int(eta*100):02d}.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"Saved {fname}")

# ── Mach on symmetry plane (y=0 or z=0 depending on orientation) ──────────────
# Check mesh orientation: ONERA M6 in hisa example — wing in xz plane, y=span
# Symmetry plane is y=0

centres  = internal.cell_centers().points
x_c = centres[:, 0]
y_c = centres[:, 1]
z_c = centres[:, 2]

T_c  = internal.cell_data["T"]
U_c  = internal.cell_data["U"]
Umag = np.linalg.norm(U_c, axis=1)
a_c  = np.sqrt(gamma * R * T_c)
Ma_c = Umag / a_c

print(f"Mach range: {Ma_c.min():.3f} – {Ma_c.max():.3f}")

# Symmetry plane: z ~ 0 (wing root, symmetry plane of domain)
# Wing is in xz plane, farfield in y, z=spanwise; root is at z=0
z_root_tol = b_semi * 0.04
mask_sym   = np.abs(z_c) < z_root_tol

print(f"Root-plane cells (z~0): {mask_sym.sum()}")

if mask_sym.sum() > 0:
    xp, yp, Mp = x_c[mask_sym], y_c[mask_sym], Ma_c[mask_sym]

    fig, ax = plt.subplots(figsize=(12, 8))
    levels = np.linspace(0, 1.1, 51)
    tc = ax.tricontourf(xp, yp, Mp, levels=levels, cmap="coolwarm", extend="max")
    cb = fig.colorbar(tc, ax=ax, label="Mach", orientation="horizontal",
                      location="top", fraction=0.04, pad=0.12, shrink=0.4, aspect=30,
                      extend="neither")
    cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1])
    ax.set_xlim(-0.3, 1.6)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    fig.suptitle("ONERA M6 — Mach Contours, Root Section z=0 (Ma=0.84, α=3.06°)", y=0.98, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "mach_root_section.png"), dpi=150)
    plt.close(fig)
    print("Saved mach_root_section.png")

# ── Wing surface colored by Cp — 3D view ──────────────────────────────────────
wing_cp = wing.copy()
wing_cp.cell_data["Cp"] = Cp_wing

# Use pyvista plotter for 3D surface
plotter = pv.Plotter(off_screen=True, window_size=[1600, 900])
plotter.add_mesh(wing_cp, scalars="Cp", cmap="coolwarm",
                 clim=[-1.5, 0.5], show_edges=False,
                 scalar_bar_args={"title": "Cp", "vertical": True,
                                  "position_x": 0.88, "position_y": 0.05,
                                  "width": 0.08, "height": 0.7})
plotter.camera_position = "xz"
plotter.camera.elevation = 25
plotter.camera.azimuth   = -30
plotter.camera.zoom(1.2)
plotter.background_color = "white"
plotter.show(screenshot=os.path.join(out_dir, "wing_cp_3d.png"))
print("Saved wing_cp_3d.png")

# ── Cp at root (2y/b≈0) for simple validation check ──────────────────────────
eta_root = 0.05
mask_root = np.abs(z_w - eta_root * b_semi) < b_semi * 0.04
if mask_root.sum() > 0:
    x_loc  = x_w[mask_root]
    cp_loc = Cp_wing[mask_root]
    y_loc  = y_w[mask_root]
    x_le   = x_loc.min(); chord = x_loc.max() - x_le
    x_norm = (x_loc - x_le) / chord
    upper  = y_loc >= 0.0; lower = y_loc < 0.0
    fig, ax = plt.subplots(figsize=(7, 5))
    if upper.sum() > 0:
        o = np.argsort(x_norm[upper])
        ax.plot(x_norm[upper][o], cp_loc[upper][o], "b-", lw=2, label="upper surface")
    if lower.sum() > 0:
        o = np.argsort(x_norm[lower])
        ax.plot(x_norm[lower][o], cp_loc[lower][o], "r-", lw=2, label="lower surface")
    ax.invert_yaxis()
    ax.set_xlim(0, 1); ax.set_xlabel("x/c"); ax.set_ylabel("$C_p$")
    ax.set_title(f"ONERA M6 Root Section Cp  (2y/b≈{eta_root:.2f})")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "cp_root.png"), dpi=150)
    plt.close(fig)
    print("Saved cp_root.png")
