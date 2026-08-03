# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A design tool for multiplex modular fiber-positioner telescopes (MOS instruments). It computes, for a
given telescope focal surface and a triangular "module" of fiber positioner robots:

1. **Coverage** — fraction of the vignetting area reachable by the summed workspaces of all robots.
2. **Positions** — 3D positions (x, y, z, theta, phi) of every robot and module across the curved
   focal surface.

Approach is inspired by Joe Silber's (LBNL) [raft-design](https://github.com/joesilber/raft-design).

## Running

```
python main.py
```

There is no CLI/argparse for `main.py` — all run parameters (project, robot pitch, HR fiber indices,
inner/global gaps, GFA layout, output toggles) are set by editing constants near the top of `main.py`
directly. The `#%%` markers throughout mark it as a Jupyter-cell-style script (run interactively as
notebook cells, e.g. in VS Code / Spyder, or straight through as a plain script).

Python environment: a `.venv/` is checked out locally (gitignored). `requirements.txt` is UTF-16
encoded (re-save as UTF-8 if editing it with tools that assume UTF-8). Key third-party deps: `shapely`,
`geopandas`, `numpy`, `scipy`, `pandas`, `matplotlib`, `ezdxf` (DXF export), `progress` (CLI bar).

Each module file (`Robot.py`, `Module.py`, `Focal_Suface.py`, `GFAs.py`) also has an
`if __name__ == "__main__":` block with a minimal standalone demo/plot — run any of them directly
(e.g. `python Module.py`) to sanity-check that piece in isolation.

There is no test suite, linter, or build step configured in this repo.

## Architecture

`main.py` is the entry point and orchestrates everything; it does not contain reusable logic itself.
The pipeline it drives:

1. **`projects.json`** — per-project optical parameters (curvature, vignetting diameter `vigD`,
   f-number, BFS, FoV, WST "donut hole" diameter, etc.). Selected via the `PROJECT` string constant
   in `main.py` (e.g. `"WST25"`, `"MUST"`, `"Spec-S5"`, `"VLT_2030"`).
2. **`Focal_Suface.py` → `FocalSurf`** — loads the corresponding Zemax-exported CSV from
   `Data_focal_planes/{project}.csv` and builds interpolated transfer functions (`R2Z`, `R2CRD`,
   `R2NORM`, `R2NUT`, `S2R`) mapping radial position on the flat layout to height/angle on the real
   curved focal surface. Also owns the vignetting disk / trimming polygon geometry (shapely `Polygon`)
   and unit conversions (mm ↔ arcsec/arcmin/deg, used for WST's angular axes).
3. **`Robot.py` → `Robot`** — one fiber positioner: two arm lengths (`l_alpha`, `l_beta`) define an
   annular workspace `Polygon` (lazily cached). Tracks both its original flat-grid position (`x0,y0,z0`)
   and its updated 3D position (`x1,y1,z1`) after projection onto the curved surface.
4. **`Module.py` → `Module`** — a triangular block of `nb_robots` `Robot`s arranged on a hex-packed
   grid with configurable pitch, wall/no-wall coverage clipping, chamfered corners, and per-robot
   HR/LR fiber-type assignment (`HR_fibers` = list of robot indices getting the high-resolution arm
   lengths). Computes `LR_coverage` / `HR_coverage` / `module_coverage` as unioned shapely geometries.
   `mod0` in `main.py` is a single reference module built at the origin (used only for the isolated
   module plot); one more `Module` instance is then built per grid cell.
5. **`updown_tri.py`** — pure-geometry helpers for laying out up/down-pointing triangles on a
   triangular lattice; used by `Grid.py`.
6. **`GFAs.py` → `GFA`** — places Guide/Focus/Alignment sensor rectangles evenly around `vigR`; their
   footprints are hard exclusion zones (no module may overlap a GFA).
7. **`Grid.py` → `Grid`** — the layout engine: tiles modules (grouped 4-at-a-time into "intermediate
   triangles" per `inner_gap`/`global_gap`) across the flat plane (`flat_grid`), projects each module
   centroid onto the curved surface via `FocalSurf`'s transfer functions (`grid_3d`), computes fiducial
   positions, the concave hull of the resulting layout, and trims modules outside the field.
8. **`main.py` main loop** — for every grid cell, builds a real `Module`, checks whether it sticks out
   past the vignetting/limiting polygon or overlaps a GFA (drop it, or clip its coverage if it only
   partially overlaps), and accumulates per-module HR/LR coverage areas plus a flat `robots_workspaces`
   table (one row per robot, all positions/geometries).
9. **`SavingResults.py` → `SavingResults`** — all file output funnels through this: PNG plots
   (`save_plots`), DXF export of module/robot boundaries (`save_dxf`, via `ezdxf`), grid dumps as
   `.txt`/`.csv` (`save_txt`/`save_csv`). Auto-creates `Results/{project_name}/` and timestamps every
   filename. Toggled by the `save = SavingResults({...})` dict at the top of `main.py`; nothing is
   written unless the corresponding flag is `True`.
10. **`CustomLegends.py`** — pure helper functions building matplotlib legend handles and plot titles
    (HR/LR/GFA/fiducial markers, the summary title strings) shared by `Module.plot_module` and
    `main.py`'s final layout plot.

Coordinate/units convention used throughout: origin at focal surface center, z toward the fiber tips,
`theta` = polar angle from +z, `phi` = azimuthal angle in the xy-plane. Positions carry both a `0`
suffix (initial/flat-grid placement) and a `1` suffix (final placement after projection onto the
curved surface) on `Robot`, `Module`, and the grid dataframes.

For WST projects specifically: there's a circular "donut hole" cut out of the center (IFU mode
exclusion, `FocalSurf.donut_hole`), and the final plot's y-axis is relabeled from mm to degrees via
`FocalSurf.mm2deg`.

## Legacy / parallel code — not part of the `main.py` pipeline

- **`parameters.py`** and **`focal_plane_coverage.py`** are an older, self-contained, monolithic
  predecessor of the current `FocalSurf`/`Module`/`Grid`/`SavingResults` split. They are not imported
  by `main.py` or any of its dependencies. `optimize.py`, `read_Ansys_results.py`, and several scripts
  under `Additional_scripts/` (`draw_stairs.py`, `focus_differences_within_module.py`,
  `wks_for_astrobots_tp.py`, `compare_focal_surfaces.py`) still depend on `parameters.py`, so it can't
  be deleted outright, but new coverage/layout work should go through the `main.py` pipeline instead.
- **`Additional_scripts/`** holds one-off analyses/experiments (throughput comparisons, focus-error
  studies, angle-arrangement sketches) rather than reusable library code.
- **`SolidworksAutomationTool/`** is a separate sub-tool (see its own README) for driving SolidWorks
  automation from exported layouts; unrelated to the Python geometry pipeline.

## Data

- `Data_focal_planes/` holds the per-project Zemax-exported optics CSVs (columns typically `R`, `Z`,
  optionally `CRD`) that `FocalSurf.optics_data()` reads by convention from
  `./Data_focal_planes/{PROJECT}.csv`. `MUST` uses `;` as separator; everything else uses `,`.
- `Results/` (gitignored) is created on demand by `SavingResults`, one subfolder per project name.
- `Results_examples/` has committed sample outputs referenced from `README.md`.
