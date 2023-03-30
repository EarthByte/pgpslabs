# pgpslabs

Code to make simple reconstructions of slab geometries using gplates-format topological plate reconstructions.

The idea is principally inspired by various papers of Derek Thorkelson to map the extent of slab windows through geological time. The windows form where mid-ocean ridges intersect with subduction zones. This code maps points in the slabs themselves, so that slab-windows are visualised as gaps between different slab segments.

Assumptions:

* all calculations assume a single dip angle throughout (workflow can handle plate reconstructions inputs with different dip angles, but will treat them as if they all have the same dip angle)
* all slabs convergence rates are taken from the plate reconstruction
* In current version, the slabs are represented as points along lines of equal subduction age (currently no method is implemented to turn these into a surface). The attributes at each point are the age of subuction, and optionally the age of the seafloor (requires seafloor age grids from which values would be interpolated).

### Notes
The [initial version](https://github.com/siwill22/pgpslabs) of this workflow was developed by Simon Williams, and was used within [McGirr et al. (2021)](https://doi.org/10.1130/B35595.1). This respository takes the inital workflow and impliments a co-ordinate handling system to plot slab onto a map with continental landmasses, mid-ocean ridges, and trenches. It also adds velocity information and allows an animation to be created so that changes in slab position can be observed across time.

