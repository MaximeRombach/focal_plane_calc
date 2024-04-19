# Focal Plane Calculations

The current repo aims to be a design tool for multiplex modular telescopes. It was inspired by the late work of Joe Silber from LBNL in [raft-design](https://github.com/joesilber/raft-design).

It produces 2 main information:

1. Coverage: fraction of the vignetting area covered by the summed up workspaces of all the robots
2. Positions: 3D positions (x,y,z, $\theta$, $\phi$) of robots & modules within the focal plane

# Results examples

Typical resulting files can be found in *Results_examples* folder:

## Robots workspace in individual modules

The first output shows the smallest building blovk of this program. It shows the shape of a module of positioners, along with each individual workspace. As we are interested in the net coverages of each modules are summed up to calculate the net area covered by the positioners with a tirangular module.

Individual robots workspaces             |  Summed robots workspace
:-------------------------:|:-------------------------:
![2023-10-20-11-39-57_Module_coverage_raw_63_robots_per_module](https://github.com/MaximeRombach/focal_plane_calc/assets/73754365/3b3756ea-3b74-4ccb-b611-998703b1c3b1) |  ![Net_robots_wks](Results_examples/2023-10-20-11-39-59_Module_cov_w_walls__63_robots_per_mod.png)

## Robots workspace in intermediate modules

The goal of this program is to study the impact of locally changing the gap between modules on the total coverage. They are thus packed by groups of 4 ideally without physical walls between them. The said gap is the *inner gap* variable in *focal_plane_covergae.py*.

![](Results_examples/2023-10-20-11-40-00_Intermediate_plot_63_robots_per_mod.png)

## Total coverage

The modules are then placed 4 by 4 (in green) in the focal plane bounded by its vignetting radius (*vigR*). Individual modules (in blue) are used as we get closer to the edge of the vigR, more can be inserted as we increase the *out_allowance* parameter, which sets how much the workspace of a module is allowed to stick out of vigR.

![](Results_examples/2023-10-20-11-40-06_Coverage_global_63_rob__Inner_1_mm__Global_3_mm.png)
