"""
Compact fast-spectrum microreactor model (withdrawn case)

Case: 6 sodium heat-pipe channels + central control rod channel (withdrawn/void).

Purpose:
- Run eigenvalue calculations for a specified fuel temperature.
- Save a per-case statepoint for downstream spectrum extraction.
- Append k-effective results to a Doppler sweep CSV.

Inputs:
- --Tfuel: fuel temperature in K
- --plot: enable geometry plot generation (off by default)

Outputs:
- results/statepoints/<CASE_NAME>.h5
- results/doppler_keff.csv (appends one row)
"""

import argparse
import os
import numpy as np
import openmc


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--Tfuel", type=float, required=True, help="Fuel temperature in K.")
    parser.add_argument("--plot", action="store_true", help="Enable geometry plotting.")
    args = parser.parse_args()

    # ----------------------------
    # Case identifiers
    # ----------------------------

    Tfuel = float(args.Tfuel)
    CASE_NAME = f"hp6_rod_withdrawn_Na_T{int(round(Tfuel))}_refl10"

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
    os.makedirs("results/statepoints", exist_ok=True)

    # ----------------------------
    # Materials
    # ----------------------------

    fuel = openmc.Material(name="U_metal_HEU_surrogate")
    fuel.set_density("g/cm3", 18.5)
    fuel.add_nuclide("U235", 0.93, "ao")
    fuel.add_nuclide("U238", 0.07, "ao")
    fuel.temperature = Tfuel

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
    outer_cyl = openmc.ZCylinder(r=core_radius + reflector_thickness, boundary_type="vacuum")

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
    # Settings
    # ----------------------------

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.batches = 200
    settings.inactive = 40
    settings.particles = 150000
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

    model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)

    if args.plot:
        os.makedirs("results/plots", exist_ok=True)

        plot_xy = openmc.Plot()
        plot_xy.filename = "geometry_xy"
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
        plot_xz.filename = "geometry_xz"
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

        if os.path.exists("geometry_xy.png"):
            os.replace("geometry_xy.png", f"results/plots/{CASE_NAME}_xy.png")
        if os.path.exists("geometry_xz.png"):
            os.replace("geometry_xz.png", f"results/plots/{CASE_NAME}_xz.png")

    # ----------------------------
    # Run and save statepoint
    # ----------------------------

    statepoint_path = model.run()

    saved_statepoint = f"results/statepoints/{CASE_NAME}.h5"
    os.replace(statepoint_path, saved_statepoint)
    statepoint_path = saved_statepoint

    # ----------------------------
    # Extract keff
    # ----------------------------

    sp = openmc.StatePoint(statepoint_path)
    keff = sp.keff
    keff_val = float(keff.n)
    keff_unc = float(keff.s)
    sp.close()

    print(f"\n{CASE_NAME}")
    print(f"k-eff = {keff_val:.5f} +/- {keff_unc:.5f}")

    # ----------------------------
    # Append Doppler CSV
    # ----------------------------

    csv_path = "results/doppler_keff.csv"
    write_header = not os.path.exists(csv_path)

    with open(csv_path, "a", encoding="utf-8") as f:
        if write_header:
            f.write("case,Tfuel_K,keff,unc\n")
        f.write(f"{CASE_NAME},{Tfuel:.1f},{keff_val:.8f},{keff_unc:.8f}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())