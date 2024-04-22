# wlp-dual-energy

## Project Goals
- Quantify the amount of water and lipid in the simulated coronary arteries using [dual energy material decomposition](https://github.com/MolloiLab/MaterialDecomposition.jl) and [single energy volume fraction](https://github.com/MolloiLab/CalciumScoring.jl)
- Compare the accuracy, sensitivity, robustness, etc. of the two techniques and publish a paper

## TODO
- [ ] Update the `predictions.jl` notebook to directly utilize density (mg/cm^3) of lipid, instead of lipid/water percentage (%)
- [ ] Convert lipid density calculations to lipid mass, utilizing the volume of the segmented coronary artery
- [ ] Include volume fraction in the `predictions.jl` script for calculating mass, from the CalciumScoring.jl package, and compare against material decomposition


## DrWatson Stuff...
This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> wlp-dual-energy

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "wlp-dual-energy"
```
which auto-activate the project and enable local path handling from DrWatson.
