# SolidworksAutomationTool

# Brief

The Solidworks Automation Tool project aims to automate complex and repetitive modelling tasks. 

Though there are parameter-based CAD software available, we chose to try automating Solidworks, as it is one of the most widely used CAD modelling software. 

The original idea was to have Python talking to Solidworks and give modelling commands. However, as of July 2023, Solidworks does not officially support Python. Instead, Solidworks has official APIs for C#, VBA, and C++.

We chose C# because of its high performance, high-level syntax and built-in garbage collector. 

# Development Environment Setup

## Pre-requisites

- Visual Studio (the author is using version 2022)

- Solidworks

- The `focal_plane_calc` repo, which is the parent directory of the automation tool.

# Current Progress
![image](https://github.com/MaximeRombach/focal_plane_calc/assets/85515041/f7f39f5d-9507-4501-9bc2-6f7f0cd8d46a)

Can repeatedly create sketches of modules(both chamfered and unchamfered) and extrude using them.

# TODOs

## On the logic inside the Automation Tool

- [x] Complete the routine to import and place point clouds

- [x] Complete the routine to make extrusion axes

- [x] Complete the routine for making a "full pizza slice"

- [x] Complete the routine for making a full-triangle block

- [x] Complete the routine for making a chamfered-triangle block

- [ ] Reliably position the (full and chamfered) triangles on the "full pizza slice"

- [ ] Extrude the modules

- [ ] Repeat the "extruded pizza slice" to make a full focol plane

## On the high-level planning

- [ ] Make GUI to help the user adjust various parameters easily

## On optimizations

- [ ] Make wrapper functions to reduce boilerplat code
