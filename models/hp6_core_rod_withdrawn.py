"""
Compact fast-spectrum microreactor model

Case: 6 sodium heat-pipe channels + central control rod channel (withdrawn/void).

Model features:
- Metallic HEU fuel (fast-spectrum surrogate)
- 6 sodium channels (heat pipes represented as sodium-filled cylinders)
- Central control rod channel present but unfilled (void)
- Full-height BeO reflector
- Eigenvalue calculation and core flux spectrum tally

Outputs:
- Geometry plots (raw + annotated) under results/plots/
- Prints k-effective
- Writes a CSV summary to results/keff_cases.csv (appends one row)
"""

import os
import numpy as np
import openmc
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# -------------------------------------------------
# Annotated plot utility
# -------------------------------------------------

def save_annotated_plot(
    raw_png_path: str,
    out_png_path: str,
    title: str,
    xlabel: str,
    ylabel: str,
    x_span_cm: float,
    y_span_cm: float,
    scalebar_cm: float = 10.0
) -> None:
    """
    Save an annotated copy of an OpenMC geometry plot.

    The annotation adds axis labels, a title, and a scale bar of fixed length.
    The axes are calibrated using the plot width/height in centimeters.
    """
    img = mpimg.imread(raw_png_path)

    fig = plt.figure(figsize=(6, 6), dpi=200)
    ax = fig.add_subplot(1, 1, 1)

    extent = (-x_span_cm / 2.0, x_span_cm / 2.0,
              -y_span_cm / 2.0, y_span_cm / 2.0)

    ax.imshow(img, extent=extent, origin="upper")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    x0 = -x_span_cm / 2.0 + 0.08 * x_span_cm
    y0 = -y_span_cm / 2.0 + 0.08 * y_span_cm
    x1 = x0 + scalebar_cm

    ax.plot([x0, x1], [y0, y0], linewidth=4)
    ax.text((x0 + x1) / 2.0, y0 + 0.04 * y_span_cm,
            f"{scalebar_cm:.0f} cm", ha="center", va="bottom")

    ax.set_xlim(-x_span_cm / 2.0, x_span_cm / 2.0)
    ax.set_ylim(-y_span_cm / 2.0, y_span_cm / 2.0)

    fig.tight_layout()
    fig.savefig(out_png_path)
    plt.close(fig)


# ----------------------------
# Case identifiers
# ----------------------------

CASE_NAME = "hp6_rod_withdrawn_Na_T600_refl10"


# ----------------------------
# Geometry parameters (cm)
# ----------------------------

core_radius = 6.0
core_height = 20.0
reflector_thickness = 9.0

n_pipes = 6
pipe_radius = 1.0
pipe_ring_radius = 5.5

rod_radius = 1.5


# ----------------------------
# Output directories
# ----------------------------

os.makedirs("results", exist_ok=True)
os.makedirs("results/plots", exist_ok=True)


# ----------------------------
# Materials
# ----------------------------

fuel = openmc.Material(name="U_metal_HEU_surrogate")
fuel.set_density("g/cm3", 18.5)
fuel.add_nuclide("U235", 0.93, "ao")
fuel.add_nuclide("U238", 0.07, "ao")
fuel.temperature = 600.0

sodium = openmc.Material(name="Sodium")
sodium.set_density("g/cm3", 0.9)
sodium.add_element("Na", 1.0)

beo = openmc.Material(name="BeO_reflector")
beo.set_density("g/cm3", 3.02)
beo.add_element("Be", 1.0)
beo.add_element("O", 1.0)

materials = openmc.Materials([fuel, sodium, beo])


# ----------------------------
# Surfaces
# ----------------------------

fuel_cyl = openmc.ZCylinder(r=core_radius)
outer_cyl = openmc.ZCylinder(
    r=core_radius + reflector_thickness,
    boundary_type="vacuum"
)

z_bot = openmc.ZPlane(z0=-core_height / 2.0, boundary_type="vacuum")
z_top = openmc.ZPlane(z0=+core_height / 2.0, boundary_type="vacuum")
axial_window = +z_bot & -z_top

rod_cyl = openmc.ZCylinder(r=rod_radius)

pipe_cyls = []
for i in range(n_pipes):
    theta = 2.0 * np.pi * i / n_pipes
    x0 = pipe_ring_radius * np.cos(theta)
    y0 = pipe_ring_radius * np.sin(theta)
    pipe_cyls.append(openmc.ZCylinder(x0=x0, y0=y0, r=pipe_radius))


# ----------------------------
# Cells
# ----------------------------

fuel_region = -fuel_cyl & axial_window
fuel_region &= ~(-rod_cyl)
for pc in pipe_cyls:
    fuel_region &= ~(-pc)

fuel_cell = openmc.Cell(name="Fuel", region=fuel_region, fill=fuel)

rod_region = -rod_cyl & axial_window
rod_cell = openmc.Cell(name="RodChannelVoid", region=rod_region, fill=None)

pipe_cells = []
for j, pc in enumerate(pipe_cyls):
    reg = -pc & axial_window
    pipe_cells.append(openmc.Cell(name=f"HeatPipe_{j+1}", region=reg, fill=sodium))

reflector_region = +fuel_cyl & -outer_cyl & axial_window
reflector_cell = openmc.Cell(name="Reflector", region=reflector_region, fill=beo)

universe = openmc.Universe(cells=[fuel_cell, rod_cell, reflector_cell] + pipe_cells)
geometry = openmc.Geometry(universe)


# ----------------------------
# Geometry plots (raw + annotated)
# ----------------------------

plot_xy = openmc.Plot()
plot_xy.filename = "geometry_xy_raw"
plot_xy.width = (2 * (core_radius + reflector_thickness),
                 2 * (core_radius + reflector_thickness))
plot_xy.pixels = (800, 800)
plot_xy.color_by = "material"
plot_xy.basis = "xy"
plot_xy.colors = {
    fuel: (200, 170, 80),
    sodium: (30, 100, 160),
    beo: (150, 150, 150),
}

plot_xz = openmc.Plot()
plot_xz.filename = "geometry_xz_raw"
plot_xz.width = (2 * (core_radius + reflector_thickness),
                 core_height)
plot_xz.pixels = (800, 800)
plot_xz.color_by = "material"
plot_xz.basis = "xz"
plot_xz.colors = {
    fuel: (200, 170, 80),
    sodium: (30, 100, 160),
    beo: (150, 150, 150),
}

plots = openmc.Plots([plot_xy, plot_xz])
plots.export_to_xml()

openmc.plot_geometry()

xy_raw_local = "geometry_xy_raw.png"
xz_raw_local = "geometry_xz_raw.png"

xy_raw_alt = "geometry_xy.png"
xz_raw_alt = "geometry_xz.png"

if not os.path.exists(xy_raw_local) and os.path.exists(xy_raw_alt):
    xy_raw_local = xy_raw_alt

if not os.path.exists(xz_raw_local) and os.path.exists(xz_raw_alt):
    xz_raw_local = xz_raw_alt

if not os.path.exists(xy_raw_local):
    raise FileNotFoundError(f"Expected plot not found in working directory: {xy_raw_local}")

if not os.path.exists(xz_raw_local):
    raise FileNotFoundError(f"Expected plot not found in working directory: {xz_raw_local}")

xy_raw = "results/plots/geometry_xy_raw.png"
xz_raw = "results/plots/geometry_xz_raw.png"

os.replace(xy_raw_local, xy_raw)
os.replace(xz_raw_local, xz_raw)

xy_span = 2.0 * (core_radius + reflector_thickness)
xz_span_x = 2.0 * (core_radius + reflector_thickness)

save_annotated_plot(
    raw_png_path=xy_raw,
    out_png_path="results/plots/geometry_xy.png",
    title="Microreactor Geometry (XY Midplane)",
    xlabel="x (cm)",
    ylabel="y (cm)",
    x_span_cm=xy_span,
    y_span_cm=xy_span,
    scalebar_cm=10.0
)

save_annotated_plot(
    raw_png_path=xz_raw,
    out_png_path="results/plots/geometry_xz.png",
    title="Microreactor Geometry (XZ Axial Slice)",
    xlabel="x (cm)",
    ylabel="z (cm)",
    x_span_cm=xz_span_x,
    y_span_cm=core_height,
    scalebar_cm=10.0
)


# ----------------------------
# Settings
# ----------------------------

settings = openmc.Settings()
settings.run_mode = "eigenvalue"
settings.batches = 100
settings.inactive = 20
settings.particles = 50000
settings.seed = 42


# ----------------------------
# Tallies (core flux spectrum)
# ----------------------------

energy_bins = np.logspace(-3, 7, 400)
energy_filter = openmc.EnergyFilter(energy_bins)

flux_tally = openmc.Tally(name="core_flux_spectrum")
flux_tally.filters = [energy_filter]
flux_tally.cells = [fuel_cell]
flux_tally.scores = ["flux"]

tallies = openmc.Tallies([flux_tally])


# ----------------------------
# Run
# ----------------------------

model = openmc.Model(
    geometry=geometry,
    materials=materials,
    settings=settings,
    tallies=tallies
)

statepoint_path = model.run()

saved_statepoint = f"results/{CASE_NAME}.h5"
os.replace(statepoint_path, saved_statepoint)
statepoint_path = saved_statepoint

print(f"Saved statepoint to {statepoint_path}")


# ----------------------------
# Extract results
# ----------------------------

sp = openmc.StatePoint(statepoint_path)
keff = sp.keff

try:
    keff_val = float(keff.n)
    keff_unc = float(keff.s)
except Exception:
    keff_val = float(keff)
    keff_unc = float("nan")

print(f"\n{CASE_NAME}")
print(f"k-eff = {keff_val:.5f} +/- {keff_unc:.5f}")

sp.close()


# ----------------------------
# Append to results CSV
# ----------------------------

csv_path = "results/keff_cases.csv"
write_header = not os.path.exists(csv_path)

with open(csv_path, "a", encoding="utf-8") as f:
    if write_header:
        f.write("case,keff,unc\n")
    f.write(f"{CASE_NAME},{keff_val:.8f},{keff_unc:.8f}\n")