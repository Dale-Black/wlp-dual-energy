### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d7077f76-dff8-44bb-b646-a8b691e4c435
# ╠═╡ show_logs = false
using Pkg; Pkg.activate(".."); Pkg.instantiate()

# ╔═╡ 8083b53c-2c5a-4a75-98f9-41cc7d060ec9
using DrWatson

# ╔═╡ 7d4d0074-6954-4f0c-98ef-f5f9dacfd7e7
using DICOM: dcmdir_parse

# ╔═╡ 355a986b-a188-458d-a2a2-c1936d4eda7e
using CairoMakie: Figure, Axis, heatmap!, scatter!

# ╔═╡ f4960eae-b3c8-4f87-ad1d-189ccd630cbb
using DataFrames: DataFrame

# ╔═╡ d8d7437e-3adc-4803-82da-2ebb95f80305
using Statistics: mean

# ╔═╡ 39fece51-1e13-4555-9649-6ded938e35da
using PlutoUI: TableOfContents, Slider

# ╔═╡ a5cb80f6-4db8-48a1-9524-095b8b13c787
using StatsBase: quantile!, rmsd

# ╔═╡ 2425710e-c4be-4b8a-8d4f-694e9892eb5a
using MaterialDecomposition: fit_calibration, quantify

# ╔═╡ da809d52-98fc-4ce8-aa55-19d1f97c5273
include(srcdir("dicom_utils.jl"));

# ╔═╡ 71e251f6-75c0-459d-b7a7-4f01674ac817
include(srcdir("masks.jl"));

# ╔═╡ d9cf4c69-91c1-44c5-a8fc-d3ee173e3cb8
import CSV

# ╔═╡ 1a2c0525-30cc-4c5d-9631-c2d3f04df3f5
TableOfContents()

# ╔═╡ fec4a658-ddc5-43e5-b798-7ed584c550ac
md"""
# Calibration
"""

# ╔═╡ 2b83b4a6-a46d-42f0-9a4e-e2efb61e9017
begin
	densities_cal = vcat(collect(0:10:90), 95) # percentage lipid
	sizes_cal = [30]
	energies = [80, 135]
end

# ╔═╡ 20d50bc7-cf28-4ad4-9a0d-f29bef38758f
md"""
## Low Energy
"""

# ╔═╡ e9e79430-ef65-41cb-aaf9-89e39a685a2e
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

# ╔═╡ 1c63ea37-46b4-4194-a9df-e839076ec68a
md"""
## High Energy
"""

# ╔═╡ 2077e981-7c66-4d04-85cb-0ca4b2261d7d
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

# ╔═╡ 957dfa92-2b27-4706-9c89-137acba46550
md"""
## Fit Parameters
"""

# ╔═╡ a51cc08f-9774-428e-bd5b-c9a748d51a11
calculated_intensities = hcat(means_80[:means], means_135[:means]) # low energy, high energy

# ╔═╡ a8ecae34-afe4-4bbf-a95c-daf1bd1c3a14
ps = fit_calibration(calculated_intensities, densities_cal)

# ╔═╡ 6743cc7c-9b47-4df5-9335-ae673af7a017
md"""
## Check Results
"""

# ╔═╡ b8366af7-2219-4147-b172-02d1b3ba2319
begin
	predicted_densities = []
	
	for i in 1:length(densities_cal)
		append!(
			predicted_densities, 
			quantify(calculated_intensities[i, 1], calculated_intensities[i, 2], ps)
		)
	end
end

# ╔═╡ 551db71c-bd5e-4da5-aac5-4e87cc885d85
df = DataFrame(
	densities = densities_cal,
	predicted_densities = predicted_densities,
	mean_intensities_low = means_80[:means],
	mean_intensities_high = means_135[:means],
)

# ╔═╡ ef09b04f-8f5b-41ff-9f44-2559ff432bd8
md"""
# Validation
"""

# ╔═╡ b90d765e-128f-4635-a68b-1baa728ce9b7
begin
	densities_val = [
		"61_66_71"
		"73_78_82"
		"85_89_94"
	] # percentage lipid

	sizes_val = ["small", "medium", "large"]
end;

# ╔═╡ 4739d184-2484-4e7a-9cc8-f8632d9543e8
md"""
## Load masks
"""

# ╔═╡ efe5352f-b221-41c0-adf2-118c99719b26
density = densities_val[2]

# ╔═╡ 1c16f854-f66a-40b1-a450-9d6d5d815a02
_size = sizes_val[2]

# ╔═╡ 8aed5131-d475-4d25-9181-a7b837d9bc83
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

# ╔═╡ a9ede6a4-c40e-4f9a-b3ae-82897928701b
begin
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
end;

# ╔═╡ 98e8db0e-5617-422e-8ffd-1e2bc3edd06b
begin
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

	# # Small Inserts
	# small_erosions = 0
	# eroded_mask_S_HD = erode_recursively(mask_S_HD, small_erosions)
	# eroded_mask_S_HD_3D = cat(eroded_mask_S_HD, eroded_mask_S_HD, eroded_mask_S_HD, dims=3)

	# eroded_mask_S_MD = erode_recursively(mask_S_MD, small_erosions)
	# eroded_mask_S_MD_3D = cat(eroded_mask_S_MD, eroded_mask_S_MD, eroded_mask_S_MD, dims=3)

	# eroded_mask_S_LD = erode_recursively(mask_S_LD, small_erosions)
	# eroded_mask_S_LD_3D = cat(eroded_mask_S_LD, eroded_mask_S_LD, eroded_mask_S_LD, dims=3)
end;

# ╔═╡ 8079214c-089d-4b67-bb0d-16d8edc1d45e
md"""
## Low Energy
"""

# ╔═╡ 80257e95-938e-4265-a243-991591cfa5f3
begin
	path_80 = datadir("dcms", "val", density, _size, string(energies[1]))
	dcm_80 = dcmdir_parse(path_80)
	dcm_array_80 = load_dcm_array(dcm_80)
	pixel_size = get_pixel_size(dcm_80[1].meta)
end;

# ╔═╡ eec7cba6-68b1-4acb-b118-991c17d45dca
@bind z1 Slider(axes(dcm_array_80, 3); show_value = true, default = 2)

# ╔═╡ 8eeca471-6794-43e6-a646-cfc5042cb057
let
	idxs_large = getindex.(findall(isone, eroded_mask_L_HD_3D[:, :, z1]), [1 2])
	idxs_mediun = getindex.(findall(isone, eroded_mask_M_HD_3D[:, :, z1]), [1 2])
	# idxs_small = getindex.(findall(isone, eroded_mask_S_HD_3D[:, :, z1]), [1 2])
	
	f = Figure()
	ax = Axis(
		f[1, 1]
	)
	heatmap!(transpose(dcm_array_80[:, :, z1]); colormap = :grays)
	scatter!(idxs_large[:, 2], idxs_large[:, 1]; color = :red, markersize = 1)
	scatter!(idxs_mediun[:, 2], idxs_mediun[:, 1]; color = :blue, markersize = 1)
	f
end

# ╔═╡ Cell order:
# ╠═8083b53c-2c5a-4a75-98f9-41cc7d060ec9
# ╠═d7077f76-dff8-44bb-b646-a8b691e4c435
# ╠═7d4d0074-6954-4f0c-98ef-f5f9dacfd7e7
# ╠═355a986b-a188-458d-a2a2-c1936d4eda7e
# ╠═f4960eae-b3c8-4f87-ad1d-189ccd630cbb
# ╠═d9cf4c69-91c1-44c5-a8fc-d3ee173e3cb8
# ╠═d8d7437e-3adc-4803-82da-2ebb95f80305
# ╠═39fece51-1e13-4555-9649-6ded938e35da
# ╠═a5cb80f6-4db8-48a1-9524-095b8b13c787
# ╠═2425710e-c4be-4b8a-8d4f-694e9892eb5a
# ╠═da809d52-98fc-4ce8-aa55-19d1f97c5273
# ╠═71e251f6-75c0-459d-b7a7-4f01674ac817
# ╠═1a2c0525-30cc-4c5d-9631-c2d3f04df3f5
# ╟─fec4a658-ddc5-43e5-b798-7ed584c550ac
# ╠═2b83b4a6-a46d-42f0-9a4e-e2efb61e9017
# ╟─20d50bc7-cf28-4ad4-9a0d-f29bef38758f
# ╠═e9e79430-ef65-41cb-aaf9-89e39a685a2e
# ╟─1c63ea37-46b4-4194-a9df-e839076ec68a
# ╠═2077e981-7c66-4d04-85cb-0ca4b2261d7d
# ╟─957dfa92-2b27-4706-9c89-137acba46550
# ╠═a51cc08f-9774-428e-bd5b-c9a748d51a11
# ╠═a8ecae34-afe4-4bbf-a95c-daf1bd1c3a14
# ╟─6743cc7c-9b47-4df5-9335-ae673af7a017
# ╠═b8366af7-2219-4147-b172-02d1b3ba2319
# ╠═551db71c-bd5e-4da5-aac5-4e87cc885d85
# ╟─ef09b04f-8f5b-41ff-9f44-2559ff432bd8
# ╠═b90d765e-128f-4635-a68b-1baa728ce9b7
# ╟─4739d184-2484-4e7a-9cc8-f8632d9543e8
# ╠═efe5352f-b221-41c0-adf2-118c99719b26
# ╠═1c16f854-f66a-40b1-a450-9d6d5d815a02
# ╠═8aed5131-d475-4d25-9181-a7b837d9bc83
# ╠═a9ede6a4-c40e-4f9a-b3ae-82897928701b
# ╠═98e8db0e-5617-422e-8ffd-1e2bc3edd06b
# ╟─8079214c-089d-4b67-bb0d-16d8edc1d45e
# ╠═80257e95-938e-4265-a243-991591cfa5f3
# ╟─eec7cba6-68b1-4acb-b118-991c17d45dca
# ╟─8eeca471-6794-43e6-a646-cfc5042cb057
