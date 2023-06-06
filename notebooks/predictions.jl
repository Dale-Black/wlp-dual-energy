### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 794bad23-f9f6-457b-8460-07526b3e2ac8
using DrWatson

# ╔═╡ c1c3df8c-06ed-47ca-98fe-79c092516a67
# ╠═╡ show_logs = false
@quickactivate "wlp-dual-energy"

# ╔═╡ 1689f739-21c8-41ff-87ed-21d1a8f4f630
using PlutoUI, CairoMakie, Statistics, CSV, DataFrames, DICOM, CSVFiles

# ╔═╡ 94e15038-a37e-4743-87de-8ddc2a40dca7
using StatsBase: quantile!, rmsd

# ╔═╡ dfaf0959-524e-4109-ba8e-868f5784f91c
# ╠═╡ show_logs = false
using DICOMUtils, CalciumScoring

# ╔═╡ 49374a5f-6d60-42c8-901d-5a284d7af8fd
include(srcdir("masks.jl"))

# ╔═╡ d17a81d2-eab0-441e-8236-f34fde8b814c
TableOfContents()

# ╔═╡ d682c750-5843-457b-b52e-c4986e6a48ca
md"""
# Calibration
"""

# ╔═╡ 7149a771-1db8-4571-a691-a9d62e1cfdfa
begin
	densities_cal = vcat(collect(0:10:90), 95) # percentage lipid
	sizes_cal = [30]
	energies = [80, 135]
end

# ╔═╡ f2191ea7-1555-4845-839d-1c6d06e290a7
md"""
## Low Energy
"""

# ╔═╡ 6ce19529-ab2a-4e4c-8c60-07e9e6054b44
begin
	means_80 = Dict(:density => densities_cal, :means => zeros(length(densities_cal)))
	for (i, density) in enumerate(densities_cal)
		for _size in sizes_cal
			path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[1]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
		
			center_insert1, center_insert2 = 187, 318
			
			calibration_rod = zeros(25, 25, size(dcm_array, 3))
			
			for z in axes(dcm_array, 3)
				rows, cols, depth = size(dcm_array)
				half_row, half_col = center_insert1, center_insert2
				offset = 12
				row_range = half_row-offset:half_row+offset
				col_range = half_col-offset:half_col+offset	
				calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
			end
		
			means_80[:means][i] = mean(calibration_rod)
		end
	end
end

# ╔═╡ fc8ce2dc-a7fa-4831-971c-1c23df2ba75a
md"""
## High Energy
"""

# ╔═╡ fc894a43-3b49-4eae-8f39-a115e12c1545
begin
	means_135 = Dict(:density => densities_cal, :means => zeros(length(densities_cal)))
	for (i, density) in enumerate(densities_cal)
		for _size in sizes_cal
			path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[2]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
		
			center_insert1, center_insert2 = 187, 318
			
			calibration_rod = zeros(25, 25, size(dcm_array, 3))
			
			for z in axes(dcm_array, 3)
				rows, cols, depth = size(dcm_array)
				half_row, half_col = center_insert1, center_insert2
				offset = 12
				row_range = half_row-offset:half_row+offset
				col_range = half_col-offset:half_col+offset	
				calibration_rod[:, :, z] .= dcm_array[row_range, col_range, z];
			end
		
			means_135[:means][i] = mean(calibration_rod)
		end
	end
end

# ╔═╡ e3df324f-a1ff-40a8-b24c-2ca3fe8e1981
md"""
## Fit Parameters
"""

# ╔═╡ 991ac293-1499-45fc-889c-c06fe2d1bb7a
calculated_intensities = hcat(means_80[:means], means_135[:means]) # low energy, high energy

# ╔═╡ 5642af46-e30c-48d1-8b97-bb09fc9ca7ac
ps = fit_calibration(calculated_intensities, densities_cal)

# ╔═╡ a59cd31d-22eb-443f-9d23-f12026f09650
md"""
## Check Results
"""

# ╔═╡ c461cd60-01c9-47bf-8a6d-8f0d512f759f
begin
	predicted_densities = []
	
	for i in 1:length(densities_cal)
		append!(
			predicted_densities, 
			score(calculated_intensities[i, 1], calculated_intensities[i, 2], ps, MaterialDecomposition()
			)
		)
	end
end

# ╔═╡ c0fdd285-af1d-4e59-8f95-6bc1e1a757b7
df = DataFrame(
	densities = densities_cal,
	predicted_densities = predicted_densities,
	mean_intensities_low = means_80[:means],
	mean_intensities_high = means_135[:means],
)

# ╔═╡ 2cc2c40b-e1c5-46ae-bbce-9953b117e84c
md"""
# Validation
"""

# ╔═╡ 1fb32e9e-e61e-4641-b8df-1e16fc882945
begin
	densities_val = [
		"61_66_71"
		"73_78_82"
		# "85_89_94"
	] # percentage lipid

	sizes_val = ["small", "medium", "large"]
end;

# ╔═╡ e57a0c47-8037-482a-979e-58f814766c8d
md"""
## Load masks
"""

# ╔═╡ 3d4c1088-4e25-4f78-95ce-e2ea46075679
begin
	# Load masks
	masks_large = Dict(
		:mask_L_HD => Array(CSV.read(datadir("julia_arrays", "large", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_M_HD => Array(CSV.read(datadir("julia_arrays", "large", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_S_HD => Array(CSV.read(datadir("julia_arrays", "large", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_L_MD => Array(CSV.read(datadir("julia_arrays", "large", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_M_MD => Array(CSV.read(datadir("julia_arrays", "large", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_S_MD => Array(CSV.read(datadir("julia_arrays", "large", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_M_LD => Array(CSV.read(datadir("julia_arrays", "large", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_L_LD => Array(CSV.read(datadir("julia_arrays", "large", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_S_LD => Array(CSV.read(datadir("julia_arrays", "large", "mask_S_LD.csv"), DataFrame; header=false)),
	)
	
	masks_medium = Dict(
		:mask_L_HD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_M_HD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_S_HD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_L_MD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_M_MD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_S_MD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_M_LD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_L_LD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_S_LD => Array(CSV.read(datadir("julia_arrays", "medium", "mask_S_LD.csv"), DataFrame; header=false)),
	)

	masks_small = Dict(
		:mask_L_HD => Array(CSV.read(datadir("julia_arrays", "small", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_M_HD => Array(CSV.read(datadir("julia_arrays", "small", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_S_HD => Array(CSV.read(datadir("julia_arrays", "small", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_L_MD => Array(CSV.read(datadir("julia_arrays", "small", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_M_MD => Array(CSV.read(datadir("julia_arrays", "small", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_S_MD => Array(CSV.read(datadir("julia_arrays", "small", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_M_LD => Array(CSV.read(datadir("julia_arrays", "small", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_L_LD => Array(CSV.read(datadir("julia_arrays", "small", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_S_LD => Array(CSV.read(datadir("julia_arrays", "small", "mask_S_LD.csv"), DataFrame; header=false)),
	)
end;

# ╔═╡ d5cdd336-ba50-497a-9397-b659d75c3029
md"""
## Run Predictions
"""

# ╔═╡ 5d6df93e-cf98-4eec-a3e2-6c7e8c7f2421
begin
	dfs = []
	for density in densities_val
		for _size in sizes_val
			if _size == "small"
				mask_L_HD = masks_small[:mask_L_HD]
				mask_M_HD = masks_small[:mask_M_HD]
				mask_S_HD = masks_small[:mask_S_HD]
				mask_L_MD = masks_small[:mask_L_MD]
				mask_M_MD = masks_small[:mask_M_MD]
				mask_S_MD = masks_small[:mask_S_MD]
				mask_L_LD = masks_small[:mask_L_LD]
				mask_M_LD = masks_small[:mask_M_LD]
				mask_S_LD = masks_small[:mask_S_LD]
			elseif _size == "medium"
				mask_L_HD = masks_medium[:mask_L_HD]
				mask_M_HD = masks_medium[:mask_M_HD]
				mask_S_HD = masks_medium[:mask_S_HD]
				mask_L_MD = masks_medium[:mask_L_MD]
				mask_M_MD = masks_medium[:mask_M_MD]
				mask_S_MD = masks_medium[:mask_S_MD]
				mask_L_LD = masks_medium[:mask_L_LD]
				mask_M_LD = masks_medium[:mask_M_LD]
				mask_S_LD = masks_medium[:mask_S_LD]
			else
				mask_L_HD = masks_large[:mask_L_HD]
				mask_M_HD = masks_large[:mask_M_HD]
				mask_S_HD = masks_large[:mask_S_HD]
				mask_L_MD = masks_large[:mask_L_MD]
				mask_M_MD = masks_large[:mask_M_MD]
				mask_S_MD = masks_large[:mask_S_MD]
				mask_L_LD = masks_large[:mask_L_LD]
				mask_M_LD = masks_large[:mask_M_LD]
				mask_S_LD = masks_large[:mask_S_LD]
			end

			

			# Large Inserts
			large_erosions = 3
			eroded_mask_L_HD = erode_recursively(mask_L_HD, large_erosions)
			eroded_mask_L_HD_3D = cat(eroded_mask_L_HD, eroded_mask_L_HD, eroded_mask_L_HD, dims=3)

			eroded_mask_L_MD = erode_recursively(mask_L_MD, large_erosions)
			eroded_mask_L_MD_3D = cat(eroded_mask_L_MD, eroded_mask_L_MD, eroded_mask_L_MD, dims=3)

			eroded_mask_L_LD = erode_recursively(mask_L_LD, large_erosions)
			eroded_mask_L_LD_3D = cat(eroded_mask_L_LD, eroded_mask_L_LD, eroded_mask_L_LD, dims=3)

			# Medium Inserts
			medium_erosions = 1
			eroded_mask_M_HD = erode_recursively(mask_M_HD, medium_erosions)
			eroded_mask_M_HD_3D = cat(eroded_mask_M_HD, eroded_mask_M_HD, eroded_mask_M_HD, dims=3)
		
			eroded_mask_M_MD = erode_recursively(mask_M_MD, medium_erosions)
			eroded_mask_M_MD_3D = cat(eroded_mask_M_MD, eroded_mask_M_MD, eroded_mask_M_MD, dims=3)
		
			eroded_mask_M_LD = erode_recursively(mask_M_LD, medium_erosions)
			eroded_mask_M_LD_3D = cat(eroded_mask_M_LD, eroded_mask_M_LD, eroded_mask_M_LD, dims=3)
			

			# Low Energy
			path_80 = datadir("dcms", "val", density, _size, string(energies[1]))
			dcm_80 = dcmdir_parse(path_80)
			dcm_array_80 = load_dcm_array(dcm_80)
			pixel_size = DICOMUtils.get_pixel_size(dcm_80[1].meta)
			
			means_80 = [
				mean(dcm_array_80[eroded_mask_L_HD_3D]), mean(dcm_array_80[eroded_mask_L_MD_3D]), mean(dcm_array_80[eroded_mask_L_LD_3D]),
				mean(dcm_array_80[eroded_mask_M_HD_3D]), mean(dcm_array_80[eroded_mask_M_MD_3D]), mean(dcm_array_80[eroded_mask_M_LD_3D]),
			]

			# High Energy
			path_135 = datadir("dcms", "val", density, _size, string(energies[2]))
			dcm_135 = dcmdir_parse(path_135)
			dcm_array_135 = load_dcm_array(dcm_135)
			
			means_135 = [
				mean(dcm_array_135[eroded_mask_L_HD_3D]), mean(dcm_array_135[eroded_mask_L_MD_3D]), mean(dcm_array_135[eroded_mask_L_LD_3D]),
				mean(dcm_array_135[eroded_mask_M_HD_3D]), mean(dcm_array_135[eroded_mask_M_MD_3D]), mean(dcm_array_135[eroded_mask_M_LD_3D]),
			]

			calculated_intensities = hcat(means_80, means_135)
			predicted_densities = zeros(length(means_80))
			for i in eachindex(predicted_densities)
				predicted_densities[i] = score(means_80[i], means_80[i], ps, MaterialDecomposition())
			end

			
			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				ground_truth_density_large_inserts = parse.(Int, split(density, "_")),
				predicted_densities_large_inserts = predicted_densities[1:3],
				ground_truth_density_medium_inserts = parse.(Int, split(density, "_")),
				predicted_densities_medium_inserts = predicted_densities[4:end]
			)
			push!(dfs, df_results)
		end
	end
end

# ╔═╡ 653ab0d9-8dda-41a3-a4cf-59a5b7beaa03
dfs

# ╔═╡ Cell order:
# ╠═794bad23-f9f6-457b-8460-07526b3e2ac8
# ╠═c1c3df8c-06ed-47ca-98fe-79c092516a67
# ╠═1689f739-21c8-41ff-87ed-21d1a8f4f630
# ╠═94e15038-a37e-4743-87de-8ddc2a40dca7
# ╠═dfaf0959-524e-4109-ba8e-868f5784f91c
# ╠═d17a81d2-eab0-441e-8236-f34fde8b814c
# ╟─d682c750-5843-457b-b52e-c4986e6a48ca
# ╠═7149a771-1db8-4571-a691-a9d62e1cfdfa
# ╟─f2191ea7-1555-4845-839d-1c6d06e290a7
# ╠═6ce19529-ab2a-4e4c-8c60-07e9e6054b44
# ╟─fc8ce2dc-a7fa-4831-971c-1c23df2ba75a
# ╠═fc894a43-3b49-4eae-8f39-a115e12c1545
# ╟─e3df324f-a1ff-40a8-b24c-2ca3fe8e1981
# ╠═991ac293-1499-45fc-889c-c06fe2d1bb7a
# ╠═5642af46-e30c-48d1-8b97-bb09fc9ca7ac
# ╟─a59cd31d-22eb-443f-9d23-f12026f09650
# ╠═c461cd60-01c9-47bf-8a6d-8f0d512f759f
# ╠═c0fdd285-af1d-4e59-8f95-6bc1e1a757b7
# ╟─2cc2c40b-e1c5-46ae-bbce-9953b117e84c
# ╠═1fb32e9e-e61e-4641-b8df-1e16fc882945
# ╟─e57a0c47-8037-482a-979e-58f814766c8d
# ╠═3d4c1088-4e25-4f78-95ce-e2ea46075679
# ╠═49374a5f-6d60-42c8-901d-5a284d7af8fd
# ╟─d5cdd336-ba50-497a-9397-b659d75c3029
# ╠═5d6df93e-cf98-4eec-a3e2-6c7e8c7f2421
# ╠═653ab0d9-8dda-41a3-a4cf-59a5b7beaa03
