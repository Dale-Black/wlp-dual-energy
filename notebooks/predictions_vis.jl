### A Pluto.jl notebook ###
# v0.19.45

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

# ╔═╡ 7d4d0074-6954-4f0c-98ef-f5f9dacfd7e7
using DICOM: dcmdir_parse, DICOMData, rescale!

# ╔═╡ 355a986b-a188-458d-a2a2-c1936d4eda7e
using CairoMakie: Figure, Axis, heatmap!, scatter!

# ╔═╡ c4ab9bfa-3a6f-440f-8a9d-adf1636299e2
using ImageMorphology: erode, dilate

# ╔═╡ f4960eae-b3c8-4f87-ad1d-189ccd630cbb
using DataFrames: DataFrame

# ╔═╡ d9cf4c69-91c1-44c5-a8fc-d3ee173e3cb8
using CSV: read

# ╔═╡ d8d7437e-3adc-4803-82da-2ebb95f80305
using Statistics: mean

# ╔═╡ 39fece51-1e13-4555-9649-6ded938e35da
using PlutoUI: TableOfContents, Slider

# ╔═╡ a5cb80f6-4db8-48a1-9524-095b8b13c787
using StatsBase: quantile!, rmsd

# ╔═╡ 2425710e-c4be-4b8a-8d4f-694e9892eb5a
using MaterialDecomposition: fit_calibration, quantify

# ╔═╡ 5da4e3b7-d6f8-42c1-9144-42bf6b39693e
using CalciumScoring: score, VolumeFraction

# ╔═╡ 08276075-156a-4e49-8b06-6c317b6b9541
using ImageMorphology: centered

# ╔═╡ 1a2c0525-30cc-4c5d-9631-c2d3f04df3f5
TableOfContents()

# ╔═╡ 900c728d-e89c-42a4-9ac8-34b2d9b473cd
md"""
## Helper functions
"""

# ╔═╡ ac679d1e-7bc7-4f93-af9c-0945ee304b11
"""
	load_dcm_array(dcm_data::Vector{DICOMData})
Given some DICOM.DICOMData, `load_dcm_array` loads the pixels
of each slice into a 3D array and returns the array
"""
function load_dcm_array(dcm_data::Vector{DICOMData})
	return array = cat(
		[rescale!(dcm_data[i])[(0x7fe0, 0x0010)] for i in 1:length(dcm_data)]...; dims=3
	)
end

# ╔═╡ 32bd1a54-b2f5-463f-8db6-a0802d1ab108
"""
    get_pixel_size(header)

Get the pixel information of the DICOM image given the `header` info.
Returns the x, y, and z values, where `z` corresponds to slice thickness

"""
function get_pixel_size(header)
    head = copy(header)
    pixel_size =
        try
            pixel_size = (head[(0x0028, 0x0030)])
            push!(pixel_size, head[(0x0018, 0x0050)])
        catch
            FOV = (head[(0x0018, 0x1100)])
            matrix_size = head[(0x0028, 0x0010)]

            pixel_size = FOV / matrix_size
            push!(pixel_size, head[(0x0018, 0x0050)])
        end
    return pixel_size
end

# ╔═╡ b0dbf47e-da61-4f8e-aa7b-23712fbb5184
begin
	function create_mask(array, mask)
	    @assert size(array) == size(mask)
	    idxs = findall(x -> x == true, mask)
	    overlayed_mask = zeros(size(array))
	    for idx in idxs
	        overlayed_mask[idx] = array[idx]
	    end
	    return overlayed_mask
	end
	
	function create_circular_mask(h, w, center_circle, radius_circle)
	    Y, X = collect(1:h), collect(1:w)'
	    dist_from_center = sqrt.((X .- center_circle[1]) .^ 2 .+ (Y .- center_circle[2]) .^ 2)
	    mask = dist_from_center .<= radius_circle
	    return mask
	end
	
	"""
	    erode_recursively(mask, n)
	
	Recursively erode a `mask`
	
	#### Inputs
	- `mask`: boolean array to erode
	- `n`: number of erosions
	
	#### Returns
	- `eroded_mask`: eroded mask with `n` erosions
	"""
	function erode_recursively(mask, n)
	    eroded_mask = copy(mask)
	    for _ in 1:n
	        eroded_mask = erode(eroded_mask)
	    end
	    return eroded_mask
	end
	
	"""
	    dilate_recursively(mask, n)
	
	Recursively erode a `mask`
	
	#### Inputs
	- `mask`: boolean array to erode
	- `n`: number of erosions
	
	#### Returns
	- `dilated_mask`: dilated mask with `n` erosions
	"""
	function dilate_recursively(mask, n)
	    dilated_mask = copy(mask)
	    for _ in 1:n
	        dilated_mask = dilate(dilated_mask)
	    end
	    return dilated_mask
	end
	
	function dilate_mask_large(mask)
	    return dilate_recursively(mask, 1)
	end
	
	function ring_mask_large(dilated_mask)
	    return Bool.(dilate_recursively(dilated_mask, 6) - dilated_mask)
	end
	
	function dilate_mask_medium(mask)
	    return dilate_recursively(mask, 3)
	end
	
	function ring_mask_medium(dilated_mask)
	    return Bool.(dilate_recursively(dilated_mask, 6) - dilated_mask)
	end
	
	function dilate_mask_small(mask)
	    return dilate_recursively(mask, 3)
	end
end

# ╔═╡ bc546e51-c327-4624-80b2-1c012a960d6f
function ring_mask_small(dilated_mask)
    return Bool.(dilate_recursively(dilated_mask, 6) - dilated_mask)
end

# ╔═╡ 206b4b46-69c4-4cc2-af0e-71caba005b51
begin
	function dilate_mask_large_bkg(mask)
	    return dilate_recursively(mask, 2)
	end
	
	function dilate_mask_medium_bkg(mask)
	    return dilate_recursively(mask, 1)
	end
	
	function dilate_mask_small_bkg(mask)
	    return mask
	end
end

# ╔═╡ fec4a658-ddc5-43e5-b798-7ed584c550ac
md"""
# Calibration Phantom
"""

# ╔═╡ 7b03cadf-af98-4850-96c4-19f458319e0f
data_dir = joinpath(dirname(pwd()), "data")

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
			path = joinpath(data_dir, "dcms", "cal", string(density), string(_size), string(energies[1]))
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
			path = joinpath(data_dir, "dcms", "cal", string(density), string(_size), string(energies[2]))
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
# Measurement Phantom
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
phantom_size = sizes_val[1] # size of the actual phantom (not the inserts)

# ╔═╡ 8aed5131-d475-4d25-9181-a7b837d9bc83
begin
	masks_large = Dict(
		:mask_large_low_lipid => Array(read(joinpath(data_dir, "julia_arrays", "large", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_medium_low_lipid => Array(read(joinpath(data_dir, "julia_arrays", "large", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_small_low_lipid => Array(read(joinpath(data_dir, "julia_arrays", "large", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_large_high_lipid => Array(read(joinpath(data_dir, "julia_arrays", "large", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_medium_high_lipid => Array(read(joinpath(data_dir, "julia_arrays", "large", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_small_high_lipid => Array(read(joinpath(data_dir, "julia_arrays", "large", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_large_medium_lipid => Array(read(joinpath(data_dir, "julia_arrays", "large", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_medium_medium_lipid => Array(read(joinpath(data_dir, "julia_arrays", "large", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_small_medium_lipid => Array(read(joinpath(data_dir, "julia_arrays", "large", "mask_S_LD.csv"), DataFrame; header=false)),
	)
	
	masks_medium = Dict(
		:mask_large_low_lipid => Array(read(joinpath(data_dir, "julia_arrays", "medium", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_medium_low_lipid => Array(read(joinpath(data_dir, "julia_arrays", "medium", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_small_low_lipid => Array(read(joinpath(data_dir, "julia_arrays", "medium", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_large_high_lipid => Array(read(joinpath(data_dir, "julia_arrays", "medium", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_medium_high_lipid => Array(read(joinpath(data_dir, "julia_arrays", "medium", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_small_high_lipid => Array(read(joinpath(data_dir, "julia_arrays", "medium", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_large_medium_lipid => Array(read(joinpath(data_dir, "julia_arrays", "medium", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_medium_medium_lipid => Array(read(joinpath(data_dir, "julia_arrays", "medium", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_small_medium_lipid => Array(read(joinpath(data_dir, "julia_arrays", "medium", "mask_S_LD.csv"), DataFrame; header=false)),
	)

	masks_small = Dict(
		:mask_large_low_lipid => Array(read(joinpath(data_dir, "julia_arrays", "small", "mask_L_HD.csv"), DataFrame; header=false)),
		:mask_medium_low_lipid => Array(read(joinpath(data_dir, "julia_arrays", "small", "mask_M_HD.csv"), DataFrame; header=false)),
		:mask_small_low_lipid => Array(read(joinpath(data_dir, "julia_arrays", "small", "mask_S_HD.csv"), DataFrame; header=false)),
		:mask_large_high_lipid => Array(read(joinpath(data_dir, "julia_arrays", "small", "mask_L_MD.csv"), DataFrame; header=false)),
		:mask_medium_high_lipid => Array(read(joinpath(data_dir, "julia_arrays", "small", "mask_M_MD.csv"), DataFrame; header=false)),
		:mask_small_high_lipid => Array(read(joinpath(data_dir, "julia_arrays", "small", "mask_S_MD.csv"), DataFrame; header=false)),
		:mask_large_medium_lipid => Array(read(joinpath(data_dir, "julia_arrays", "small", "mask_L_LD.csv"), DataFrame; header=false)),
		:mask_medium_medium_lipid => Array(read(joinpath(data_dir, "julia_arrays", "small", "mask_M_LD.csv"), DataFrame; header=false)),
		:mask_small_medium_lipid => Array(read(joinpath(data_dir, "julia_arrays", "small", "mask_S_LD.csv"), DataFrame; header=false)),
	)
end;

# ╔═╡ 4bfb72e7-d8d6-4182-8a01-6bcc76c376f7
function create_circular_se(radius)
	# Ensure radius is at least 1
	radius = max(1, radius)
	
	# Create a grid of coordinates
	dim = 2 * radius + 1
	center = (radius + 1, radius + 1)
	
	# Create the structuring element
	se = [((i - center[1])^2 + (j - center[2])^2) <= radius^2 for i in 1:dim, j in 1:dim]
	
	# Ensure the center pixel is always included
	se[center...] = true
	
	return centered(se)
end

# ╔═╡ 3c55665b-ff97-44d8-b7b4-be9b3df74927
begin
	# Create structuring elements
	circular_se_middle = create_circular_se(1)  # for middle masks
	circular_se_inner = create_circular_se(4)   # larger SE for inner masks
	
    if phantom_size == "small"
        masks_dict = masks_small
    elseif phantom_size == "medium"
        masks_dict = masks_medium
    else # phantom size == large
        masks_dict = masks_large
    end

    # Outer masks (unchanged)
    insert_mask_outer_low_lipid_percentage = cat(
        masks_dict[:mask_large_low_lipid],
        masks_dict[:mask_large_low_lipid],
        masks_dict[:mask_large_low_lipid];
        dims=3
    )
    insert_mask_outer_high_lipid_percentage = cat(
        masks_dict[:mask_large_high_lipid],
        masks_dict[:mask_large_high_lipid],
        masks_dict[:mask_large_high_lipid];
        dims=3
    )
    insert_mask_outer_medium_lipid_percentage = cat(
        masks_dict[:mask_large_medium_lipid],
        masks_dict[:mask_large_medium_lipid],
        masks_dict[:mask_large_medium_lipid];
        dims=3
    )

    # Middle masks (dilated with circular_se_middle)
    insert_mask_middle_low_lipid_percentage = cat(
        dilate(masks_dict[:mask_medium_low_lipid], circular_se_middle),
        dilate(masks_dict[:mask_medium_low_lipid], circular_se_middle),
        dilate(masks_dict[:mask_medium_low_lipid], circular_se_middle);
        dims=3
    )
    insert_mask_middle_high_lipid_percentage = cat(
        dilate(masks_dict[:mask_medium_high_lipid], circular_se_middle),
        dilate(masks_dict[:mask_medium_high_lipid], circular_se_middle),
        dilate(masks_dict[:mask_medium_high_lipid], circular_se_middle);
        dims=3
    )
    insert_mask_middle_medium_lipid_percentage = cat(
        dilate(masks_dict[:mask_medium_medium_lipid], circular_se_middle),
        dilate(masks_dict[:mask_medium_medium_lipid], circular_se_middle),
        dilate(masks_dict[:mask_medium_medium_lipid], circular_se_middle);
        dims=3
    )

    # Inner masks (dilated with circular_se_inner)
    insert_mask_inner_low_lipid_percentage = cat(
        dilate(masks_dict[:mask_small_low_lipid], circular_se_inner),
        dilate(masks_dict[:mask_small_low_lipid], circular_se_inner),
        dilate(masks_dict[:mask_small_low_lipid], circular_se_inner);
        dims=3
    )
    insert_mask_inner_high_lipid_percentage = cat(
        dilate(masks_dict[:mask_small_high_lipid], circular_se_inner),
        dilate(masks_dict[:mask_small_high_lipid], circular_se_inner),
        dilate(masks_dict[:mask_small_high_lipid], circular_se_inner);
        dims=3
    )
    insert_mask_inner_medium_lipid_percentage = cat(
        dilate(masks_dict[:mask_small_medium_lipid], circular_se_inner),
        dilate(masks_dict[:mask_small_medium_lipid], circular_se_inner),
        dilate(masks_dict[:mask_small_medium_lipid], circular_se_inner);
        dims=3
    )
end;

# ╔═╡ de09b176-3bff-4e49-9bfa-c60df13116e1
md"""
## Ground Truth Lipid Mass
"""

# ╔═╡ e7f988f9-dac1-4e0a-a58c-ccf5aefd2588
lipid_percentages = parse.(Int, split(density, "_"))

# ╔═╡ 35d8c23d-1cc1-470c-9ffa-c96265dfcd72
begin
	lipid_percentage_high = lipid_percentages[3]
	lipid_percentage_medium = lipid_percentages[2]
	lipid_percentage_low = lipid_percentages[1]
end;

# ╔═╡ 0d479ec4-be83-4c2d-b392-6616fa9bd4b4
begin
	const insert_diameter = 0.5  # cm
	const insert_length = 3 * 3  # cm (3 inserts * 3.0 mm slice thickness)
	const insert_volume = π * (insert_diameter/2)^2 * insert_length  # cm³
end;

# ╔═╡ c666fcd2-6316-4a71-b54f-bf1c449ae05c
function calculate_lipid_mass(insert_volume, lipid_percentage; lipid_density = 0.92)
	if lipid_percentage > 1
		lipid_percentage *= 1e-2
	end
	lipid_volume = insert_volume * (lipid_percentage)
	lipid_mass = lipid_volume * lipid_density
	return lipid_mass
end;

# ╔═╡ 5f452c4a-7e69-4717-a5dd-f01f4bdcfa14
lipid_mass_high_percentage_ground_truth = calculate_lipid_mass(insert_volume, lipid_percentage_high) # grams

# ╔═╡ dd6eba57-0a41-4069-90b6-fb708c11dc4a
lipid_mass_medium_percentage_ground_truth = calculate_lipid_mass(insert_volume, lipid_percentage_medium) # grams

# ╔═╡ 8622fe63-b052-47cd-aeda-a84adab2811e
lipid_mass_low_percentage_ground_truth = calculate_lipid_mass(insert_volume, lipid_percentage_low) # grams

# ╔═╡ 8079214c-089d-4b67-bb0d-16d8edc1d45e
md"""
## Low Energy
"""

# ╔═╡ 80257e95-938e-4265-a243-991591cfa5f3
begin
	path_80 = joinpath(data_dir, "dcms", "val", density, phantom_size, string(energies[1]))
	dcm_80 = dcmdir_parse(path_80)
	dcm_array_80 = load_dcm_array(dcm_80)
	pixel_size = get_pixel_size(dcm_80[1].meta)
end;

# ╔═╡ 8033437b-abc6-4b32-aed2-cddef2a16d3b
md"""
### Visualize
"""

# ╔═╡ eec7cba6-68b1-4acb-b118-991c17d45dca
@bind z1 Slider(axes(dcm_array_80, 3); show_value = true, default = 2)

# ╔═╡ 8eeca471-6794-43e6-a646-cfc5042cb057
let
    alpha = 0.4
    
    f = Figure(size = (1000, 1500))
    ax = Axis(
        f[1, 1],
        title = "Outer insert"
    )
    heatmap!(dcm_array_80[:, :, z1]; colormap = :grays)
    heatmap!(
		insert_mask_outer_low_lipid_percentage[:, :, z1];
		colormap = (:blues, alpha)
	)

    ax = Axis(
        f[1, 2],
        title = "Middle insert"
    )
    heatmap!(dcm_array_80[:, :, z1]; colormap = :grays)
    heatmap!(
		insert_mask_middle_low_lipid_percentage[:, :, z1];
		colormap = (:blues, alpha)
	)
	
    ax = Axis(
        f[2, 1],
        title = "Inner insert"
    )
    heatmap!(dcm_array_80[:, :, z1]; colormap = :grays)
    heatmap!(
		insert_mask_inner_low_lipid_percentage[:, :, z1];
		colormap = (:blues, alpha)
	)
	
    f
end

# ╔═╡ 2eec7b8f-bd4c-4aa5-bc9c-f9dedb69341e
md"""
### Volume Fraction Calculations
"""

# ╔═╡ 563cc9d1-5eab-439b-90c7-99403249eef2
begin
	idx = 9
	calibration_density_lipid = df[9, :densities] * 1e-4
	calibration_signal_low_energy = df[9, :mean_intensities_low]
	calibration_signal_high_energy = df[9, :mean_intensities_high]
end

# ╔═╡ a87294ca-ab86-4561-994c-cff73bbdb57f
voxel_size = prod(pixel_size)

# ╔═╡ 03f8c810-3e1b-475d-86af-f81d7470a01c
md"""
**High Density Outer Insert**
"""

# ╔═╡ 499a8c37-0b59-4616-9fa6-f89579c3b46b
begin
	calibration_background_mask_high_lipid_outer = dilate_recursively(
		insert_mask_outer_high_lipid_percentage, 2
	)
	calibration_background_mask_high_lipid_outer = Bool.(
		calibration_background_mask_high_lipid_outer - insert_mask_outer_high_lipid_percentage
	)
	calibration_background_low_energy_high_lipid_outer = mean(
		dcm_array_80[calibration_background_mask_high_lipid_outer]
	)
end

# ╔═╡ 0627ae02-319f-45a5-8d03-b21c02190819
@bind z3 Slider(axes(dcm_array_80, 3); show_value = true, default = 2)

# ╔═╡ efd60d28-4b31-43af-9af6-1d627dbe9f00
let
	f = Figure()
	ax = Axis(
		f[1, 1]
	)
	heatmap!(dcm_array_80[:, :, z3]; colormap = :grays)
	heatmap!(insert_mask_outer_high_lipid_percentage[:, :, z3]; colormap = (:jet, 0.3))
	heatmap!(calibration_background_mask_high_lipid_outer[:, :, z3]; colormap = (:viridis, 0.3))
	f
end

# ╔═╡ 54d3cab0-d797-4742-a3e0-8b5fd63bdd3c
lipid_mass_vf_high_lipid_percentage_outer = score(
    dcm_array_80[insert_mask_outer_high_lipid_percentage],
    calibration_signal_low_energy,
    calibration_background_low_energy_high_lipid_outer,
    voxel_size,
    calibration_density_lipid,
    VolumeFraction()
)

# ╔═╡ 192b3805-94f3-4797-a1af-2f1ff023783f
md"""
**Medium Density Outer Insert**
"""

# ╔═╡ 39eb363c-18f9-46c5-96d5-695b781acc94
begin
	calibration_background_mask_medium_lipid_outer = dilate_recursively(
		insert_mask_outer_medium_lipid_percentage, 2
	)
	calibration_background_mask_medium_lipid_outer = Bool.(
		calibration_background_mask_medium_lipid_outer - insert_mask_outer_medium_lipid_percentage
	)
	calibration_background_low_energy_medium_lipid_outer = mean(
		dcm_array_80[calibration_background_mask_medium_lipid_outer]
	)
end

# ╔═╡ 3c7d6058-a065-4ca8-9dd6-f5818b69ee0a
@bind z4 Slider(axes(dcm_array_80, 3); show_value = true, default = 2)

# ╔═╡ 4e649a2f-4797-4b52-984f-d318227ecd61
let
	f = Figure()
	ax = Axis(
		f[1, 1]
	)
	heatmap!(dcm_array_80[:, :, z4]; colormap = :grays)
	heatmap!(insert_mask_outer_medium_lipid_percentage[:, :, z4]; colormap = (:jet, 0.3))
	heatmap!(calibration_background_mask_medium_lipid_outer[:, :, z4]; colormap = (:viridis, 0.3))
	f
end

# ╔═╡ 7b809715-a972-4414-a936-b2bde03a335b
lipid_mass_vf_medium_lipid_percentage_outer = score(
    dcm_array_80[insert_mask_outer_medium_lipid_percentage],
    calibration_signal_low_energy,
    calibration_background_low_energy_medium_lipid_outer,
    voxel_size,
    calibration_density_lipid,
    VolumeFraction()
)

# ╔═╡ e3af684d-07af-42b3-a4b2-680db74147c7
md"""
**Low Density Outer Insert**
"""

# ╔═╡ e70aeb34-eb69-4b26-bb35-5b41740f1c10
begin
	calibration_background_mask_low_lipid_outer = dilate_recursively(
		insert_mask_outer_low_lipid_percentage, 2
	)
	calibration_background_mask_low_lipid_outer = Bool.(
		calibration_background_mask_low_lipid_outer - insert_mask_outer_low_lipid_percentage
	)
	calibration_background_low_energy_low_lipid_outer = mean(
		dcm_array_80[calibration_background_mask_low_lipid_outer]
	)
end

# ╔═╡ 5408c332-080c-4e0d-8d74-3ecc5c6bad60
@bind z5 Slider(axes(dcm_array_80, 3); show_value = true, default = 2)

# ╔═╡ 1cd9faa5-e791-4a12-9efd-9956a2e77820
let
	f = Figure()
	ax = Axis(
		f[1, 1]
	)
	heatmap!(dcm_array_80[:, :, z5]; colormap = :grays)
	heatmap!(insert_mask_outer_low_lipid_percentage[:, :, z5]; colormap = (:jet, 0.3))
	heatmap!(calibration_background_mask_low_lipid_outer[:, :, z5]; colormap = (:viridis, 0.3))
	f
end

# ╔═╡ 51418e01-1be7-4dc7-826b-d8980d5b7d54
lipid_mass_vf_low_lipid_percentage_outer = score(
    dcm_array_80[insert_mask_outer_low_lipid_percentage],
    calibration_signal_low_energy,
    calibration_background_low_energy_low_lipid_outer,
    voxel_size,
    calibration_density_lipid,
    VolumeFraction()
)

# ╔═╡ 840219c7-c85a-4a34-8aa9-9455d698c50c
md"""
**High Density Middle Insert**
"""

# ╔═╡ 932ed040-96f0-4578-b1a4-9c953c907c94
begin
	calibration_background_mask_high_lipid_middle = dilate_recursively(
		insert_mask_middle_high_lipid_percentage, 1
	)
	calibration_background_mask_high_lipid_middle = Bool.(
		calibration_background_mask_high_lipid_middle - insert_mask_middle_high_lipid_percentage
	)
	calibration_background_low_energy_high_lipid_middle = mean(
		dcm_array_80[calibration_background_mask_high_lipid_middle]
	)
end

# ╔═╡ c49f2b92-0a28-4721-9d3f-9befec50ce43
@bind z6 Slider(axes(dcm_array_80, 3); show_value = true, default = 2)

# ╔═╡ 91d3b329-163d-43db-a7f8-11216dab1de0
let
	f = Figure()
	ax = Axis(
		f[1, 1]
	)
	heatmap!(dcm_array_80[:, :, z6]; colormap = :grays)
	heatmap!(insert_mask_middle_high_lipid_percentage[:, :, z6]; colormap = (:jet, 0.3))
	heatmap!(calibration_background_mask_high_lipid_middle[:, :, z6]; colormap = (:viridis, 0.3))
	f
end

# ╔═╡ bbe0855d-7deb-49c6-9719-0d5d24d74058
lipid_mass_vf_high_lipid_percentage_middle = score(
    dcm_array_80[insert_mask_middle_high_lipid_percentage],
    calibration_signal_low_energy,
    calibration_background_low_energy_high_lipid_middle,
    voxel_size,
    calibration_density_lipid,
    VolumeFraction()
)

# ╔═╡ Cell order:
# ╠═d7077f76-dff8-44bb-b646-a8b691e4c435
# ╠═7d4d0074-6954-4f0c-98ef-f5f9dacfd7e7
# ╠═355a986b-a188-458d-a2a2-c1936d4eda7e
# ╠═c4ab9bfa-3a6f-440f-8a9d-adf1636299e2
# ╠═f4960eae-b3c8-4f87-ad1d-189ccd630cbb
# ╠═d9cf4c69-91c1-44c5-a8fc-d3ee173e3cb8
# ╠═d8d7437e-3adc-4803-82da-2ebb95f80305
# ╠═39fece51-1e13-4555-9649-6ded938e35da
# ╠═a5cb80f6-4db8-48a1-9524-095b8b13c787
# ╠═2425710e-c4be-4b8a-8d4f-694e9892eb5a
# ╠═5da4e3b7-d6f8-42c1-9144-42bf6b39693e
# ╠═1a2c0525-30cc-4c5d-9631-c2d3f04df3f5
# ╟─900c728d-e89c-42a4-9ac8-34b2d9b473cd
# ╠═ac679d1e-7bc7-4f93-af9c-0945ee304b11
# ╠═32bd1a54-b2f5-463f-8db6-a0802d1ab108
# ╠═b0dbf47e-da61-4f8e-aa7b-23712fbb5184
# ╠═bc546e51-c327-4624-80b2-1c012a960d6f
# ╠═206b4b46-69c4-4cc2-af0e-71caba005b51
# ╟─fec4a658-ddc5-43e5-b798-7ed584c550ac
# ╠═7b03cadf-af98-4850-96c4-19f458319e0f
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
# ╠═4bfb72e7-d8d6-4182-8a01-6bcc76c376f7
# ╠═3c55665b-ff97-44d8-b7b4-be9b3df74927
# ╠═08276075-156a-4e49-8b06-6c317b6b9541
# ╟─de09b176-3bff-4e49-9bfa-c60df13116e1
# ╠═e7f988f9-dac1-4e0a-a58c-ccf5aefd2588
# ╠═35d8c23d-1cc1-470c-9ffa-c96265dfcd72
# ╠═0d479ec4-be83-4c2d-b392-6616fa9bd4b4
# ╠═c666fcd2-6316-4a71-b54f-bf1c449ae05c
# ╠═5f452c4a-7e69-4717-a5dd-f01f4bdcfa14
# ╠═dd6eba57-0a41-4069-90b6-fb708c11dc4a
# ╠═8622fe63-b052-47cd-aeda-a84adab2811e
# ╟─8079214c-089d-4b67-bb0d-16d8edc1d45e
# ╠═80257e95-938e-4265-a243-991591cfa5f3
# ╟─8033437b-abc6-4b32-aed2-cddef2a16d3b
# ╟─eec7cba6-68b1-4acb-b118-991c17d45dca
# ╟─8eeca471-6794-43e6-a646-cfc5042cb057
# ╟─2eec7b8f-bd4c-4aa5-bc9c-f9dedb69341e
# ╠═563cc9d1-5eab-439b-90c7-99403249eef2
# ╠═a87294ca-ab86-4561-994c-cff73bbdb57f
# ╟─03f8c810-3e1b-475d-86af-f81d7470a01c
# ╠═499a8c37-0b59-4616-9fa6-f89579c3b46b
# ╟─0627ae02-319f-45a5-8d03-b21c02190819
# ╟─efd60d28-4b31-43af-9af6-1d627dbe9f00
# ╠═54d3cab0-d797-4742-a3e0-8b5fd63bdd3c
# ╟─192b3805-94f3-4797-a1af-2f1ff023783f
# ╠═39eb363c-18f9-46c5-96d5-695b781acc94
# ╟─3c7d6058-a065-4ca8-9dd6-f5818b69ee0a
# ╟─4e649a2f-4797-4b52-984f-d318227ecd61
# ╠═7b809715-a972-4414-a936-b2bde03a335b
# ╟─e3af684d-07af-42b3-a4b2-680db74147c7
# ╠═e70aeb34-eb69-4b26-bb35-5b41740f1c10
# ╟─5408c332-080c-4e0d-8d74-3ecc5c6bad60
# ╟─1cd9faa5-e791-4a12-9efd-9956a2e77820
# ╠═51418e01-1be7-4dc7-826b-d8980d5b7d54
# ╟─840219c7-c85a-4a34-8aa9-9455d698c50c
# ╠═932ed040-96f0-4578-b1a4-9c953c907c94
# ╟─c49f2b92-0a28-4721-9d3f-9befec50ce43
# ╟─91d3b329-163d-43db-a7f8-11216dab1de0
# ╠═bbe0855d-7deb-49c6-9719-0d5d24d74058
