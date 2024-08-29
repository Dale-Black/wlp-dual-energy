# wlp-dual-energy

## Project Goals
1. Quantify the amount of water and lipid in the simulated coronary arteries using [dual energy material decomposition](https://github.com/MolloiLab/MaterialDecomposition.jl) and [single energy volume fraction](https://github.com/MolloiLab/CalciumScoring.jl)
2. Compare the accuracy, sensitivity, robustness, etc. of the two techniques and publish a paper

## TODO
- [ ] Update the `predictions.jl` notebook to directly utilize density (mg/cm^3) of lipid, instead of lipid/water percentage (%)
- [ ] Convert lipid density calculations to lipid mass, utilizing the volume of the segmented coronary artery
- [ ] Include volume fraction in the `predictions.jl` script for calculating mass, from the CalciumScoring.jl package, and compare against material decomposition