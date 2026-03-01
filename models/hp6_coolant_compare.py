"""
Compact fast-spectrum microreactor model (coolant comparison, withdrawn/inserted)

Case: 6 heat-pipe channels (coolant-filled) + central control rod channel
      (withdrawn/void OR inserted B4C).

Purpose:
- Run eigenvalue calculations for a set of coolant choices at a fixed reflector thickness.
- Save per-case statepoints for downstream extraction.
- Append k-effective and tally-derived leakage fraction to a CSV.

Inputs:
- --coolant: coolant keywords (one or more): Na, K, NaK, Li
- --rod: control rod state: withdrawn or inserted
- --refl: reflector thickness in cm (single value)
- --Tfuel: fuel temperature in K (default 600)
- --plot: enable geometry plot generation for the first case only (off by default)

Outputs:
- results/statepoints/<CASE_NAME>.h5
- results/coolant_compare.csv (appends one row per case)
"""

import argparse
import os
import numpy as np
import openmc


def build_coolant_material(name: str) -> openmc.Material:
    """
    Create a coolant material by keyword.

    Notes:
    - Densities are placeholders for parametric comparison and should be updated
      if a specific operating state (T/P) is being modeled.
    - NaK is treated as an elemental Na/K mixture placeholder, not a phase-accurate eutectic model.
    """
    key = name.strip().lower()

    if key == "na":
        m = openmc.Material(name="Sodium")
        m.set_density("g/cm3", 0.9)
        m.add_element("Na", 1.0)
        return m

    if key == "k":
        m = openmc.Material(name="Potassium")
        m.set_density("g/cm3", 0.86)
        m.add_element("K", 1.0)
        return m

    if key == "nak":
        m = openmc.Material(name="NaK_eutectic_placeholder")
        m.set_density("g/cm3", 0.86)
        m.add_element("Na", 0.5)
        m.add_element("K", 0.5)
        return m

    if key == "li":
        m = openmc.Material(name="Lithium")
        m.set_density("g/cm3", 0.53)
        m.add_element("Li", 1.0)
        return m

    raise ValueError(f"Unsupported coolant keyword: {name}. Use Na, K, NaK, Li.")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--coolant",
        type=str,
        nargs="+",
        required=True,
        help="Coolant keywords (one or more): Na, K, NaK, Li.",
    )
    parser.add_argument(
        "--rod",
        type=str,
        required=True,
        choices=["withdrawn", "inserted"],
        help="Control rod state: withdrawn or inserted.",
    )
    parser.add_argument("--refl", type=float, required=True, help="Reflector thickness in cm (single value).")
    parser.add_argument("--Tfuel", type=float, default=600.0, help="Fuel temperature in K.")
    parser.add_argument("--plot", action="store_true", help="Enable geometry plotting (first case only).")
    args = parser.parse_args()

    # ----------------------------
    # Geometry parameters (cm)
    # ----------------------------

    core_radius = 6.0
    core_height = 20.0

    n_pipes = 6
    pipe_radius = 1.0
    pipe_ring_radius = 5.5

    rod_radius = 1.5

    # ----------------------------
    # Settings
    # ----------------------------

    batches = 200
    inactive = 40
    particles = 150000
    seed = 42

    # ----------------------------
    # Run inputs
    # ----------------------------

    Tfuel = float(args.Tfuel)
    reflector_thickness = float(args.refl)
    rod_state = str(args.rod).lower()

    coolant_keys = [str(x) for x in args.coolant]
    coolant_keys = [c.strip() for c in coolant_keys]
    coolant_keys = sorted(coolant_keys, key=lambda s: s.lower())

    # ----------------------------
    # Output directories
    # ----------------------------

    os.makedirs("results", exist_ok=True)
    os.makedirs("results/statepoints", exist_ok=True)
    os.makedirs("results/plots", exist_ok=True)

    # ----------------------------
    # CSV output
    # ----------------------------

    csv_path = "results/coolant_compare.csv"
    write_header = not os.path.exists(csv_path)

    # ----------------------------
    # Loop over coolant cases
    # ----------------------------

    for idx, coolant_name in enumerate(coolant_keys):
        coolant_tag = coolant_name.strip()
        CASE_NAME = (
            f"hp6_rod_{rod_state}_{coolant_tag}_T{int(round(Tfuel))}_refl{reflector_thickness:.1f}"
        ).replace(".", "p")

        # ----------------------------
        # Materials
        # ----------------------------

        fuel = openmc.Material(name="U_metal_HEU_surrogate")
        fuel.set_density("g/cm3", 18.5)
        fuel.add_nuclide("U235", 0.93, "ao")
        fuel.add_nuclide("U238", 0.07, "ao")
        fuel.temperature = Tfuel

        coolant = build_coolant_material(coolant_tag)

        beo = openmc.Material(name="BeO_reflector")
        beo.set_density("g/cm3", 3.02)
        beo.add_element("Be", 1.0)
        beo.add_element("O", 1.0)

        b4c = openmc.Material(name="B4C_control_rod")
        b4c.set_density("g/cm3", 2.52)
        b4c.add_nuclide("B10", 0.8)
        b4c.add_nuclide("B11", 0.2)
        b4c.add_element("C", 1.0)

        mats = [fuel, coolant]
        if reflector_thickness > 0.0:
            mats.append(beo)
        if rod_state == "inserted":
            mats.append(b4c)

        materials = openmc.Materials(mats)

        # ----------------------------
        # Surfaces
        # ----------------------------

        if reflector_thickness <= 0.0:
            outer_cyl = openmc.ZCylinder(r=core_radius, boundary_type="vacuum")
            fuel_cyl = outer_cyl
        else:
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

        rod_region = -rod_cyl & -fuel_cyl & axial_window
        if rod_state == "withdrawn":
            rod_cell = openmc.Cell(name="RodChannelVoid", region=rod_region, fill=None)
        else:
            rod_cell = openmc.Cell(name="ControlRod_B4C", region=rod_region, fill=b4c)

        pipe_cells = []
        for j, pc in enumerate(pipe_cyls):
            reg = -pc & -fuel_cyl & axial_window
            pipe_cells.append(openmc.Cell(name=f"HeatPipe_{j+1}", region=reg, fill=coolant))

        cells = [fuel_cell, rod_cell] + pipe_cells

        if reflector_thickness > 0.0:
            reflector_region = +fuel_cyl & -outer_cyl & axial_window
            reflector_cell = openmc.Cell(name="Reflector", region=reflector_region, fill=beo)
            cells.append(reflector_cell)
        else:
            reflector_cell = None

        universe = openmc.Universe(cells=cells)
        geometry = openmc.Geometry(universe)

        # ----------------------------
        # Settings
        # ----------------------------

        settings = openmc.Settings()
        settings.run_mode = "eigenvalue"
        settings.batches = batches
        settings.inactive = inactive
        settings.particles = particles
        settings.seed = seed

        # ----------------------------
        # Tallies
        # ----------------------------

        energy_bins = np.logspace(-3, 7, 400)
        energy_filter = openmc.EnergyFilter(energy_bins)

        flux_tally = openmc.Tally(name="core_flux_spectrum")
        flux_tally.filters = [energy_filter]
        flux_tally.cells = [fuel_cell]
        flux_tally.scores = ["flux"]

        boundary_filter = openmc.SurfaceFilter([outer_cyl, z_top, z_bot])

        leak_tally = openmc.Tally(name="boundary_leakage")
        leak_tally.filters = [boundary_filter]
        leak_tally.scores = ["current"]

        abs_cells = [fuel_cell] + pipe_cells
        if reflector_cell is not None:
            abs_cells.append(reflector_cell)
        if rod_state == "inserted":
            abs_cells.append(rod_cell)

        abs_filter = openmc.CellFilter(abs_cells)

        abs_tally = openmc.Tally(name="system_absorption")
        abs_tally.filters = [abs_filter]
        abs_tally.scores = ["absorption"]

        tallies = openmc.Tallies([flux_tally, leak_tally, abs_tally])

        model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)

        # ----------------------------
        # Optional plotting (first case only)
        # ----------------------------

        if args.plot and idx == 0:
            width_r = core_radius + max(0.0, reflector_thickness)

            plot_xy = openmc.Plot()
            plot_xy.filename = "geometry_xy"
            plot_xy.width = (2 * width_r, 2 * width_r)
            plot_xy.pixels = (800, 800)
            plot_xy.color_by = "material"
            plot_xy.basis = "xy"
            plot_xy.colors = {
                fuel: (200, 170, 80),
                coolant: (30, 100, 160),
                beo: (150, 150, 150),
                b4c: (40, 40, 40),
            }

            plot_xz = openmc.Plot()
            plot_xz.filename = "geometry_xz"
            plot_xz.width = (2 * width_r, core_height)
            plot_xz.pixels = (800, 800)
            plot_xz.color_by = "material"
            plot_xz.basis = "xz"
            plot_xz.colors = {
                fuel: (200, 170, 80),
                coolant: (30, 100, 160),
                beo: (150, 150, 150),
                b4c: (40, 40, 40),
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
        # Extract keff and leakage fraction (tally-derived)
        # ----------------------------

        sp = openmc.StatePoint(statepoint_path)

        keff = sp.keff
        keff_val = float(keff.n)
        keff_unc = float(keff.s)

        leak = sp.get_tally(name="boundary_leakage")
        absr = sp.get_tally(name="system_absorption")

        leak_mean_bins = leak.mean.ravel()
        leak_sd_bins = leak.std_dev.ravel()

        abs_mean_bins = absr.mean.ravel()
        abs_sd_bins = absr.std_dev.ravel()

        leak_rate = float(np.sum(np.abs(leak_mean_bins)))
        abs_rate = float(np.sum(abs_mean_bins))

        leak_rate_sd = float(np.sqrt(np.sum(leak_sd_bins ** 2)))
        abs_rate_sd = float(np.sqrt(np.sum(abs_sd_bins ** 2)))

        den = leak_rate + abs_rate
        leak_frac = leak_rate / den

        dfrac_dleak = abs_rate / (den ** 2)
        dfrac_dabs = -leak_rate / (den ** 2)
        leak_frac_sd = float(np.sqrt((dfrac_dleak * leak_rate_sd) ** 2 + (dfrac_dabs * abs_rate_sd) ** 2))

        sp.close()

        print(f"\n{CASE_NAME}")
        print(f"k-eff = {keff_val:.5f} +/- {keff_unc:.5f}")
        print(f"leakage fraction (tally) = {leak_frac:.5f} +/- {leak_frac_sd:.5f}")

        # ----------------------------
        # Append coolant comparison CSV
        # ----------------------------

        with open(csv_path, "a", encoding="utf-8") as f:
            if write_header:
                f.write(
                    "case,Tfuel_K,reflector_thickness_cm,rod_state,coolant,keff,unc,leak_frac,leak_frac_unc\n"
                )
                write_header = False
            f.write(
                f"{CASE_NAME},{Tfuel:.1f},{reflector_thickness:.3f},{rod_state},{coolant_tag},"
                f"{keff_val:.8f},{keff_unc:.8f},{leak_frac:.8f},{leak_frac_sd:.8f}\n"
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())