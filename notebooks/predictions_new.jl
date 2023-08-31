### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 51258732-faf5-4705-82b6-64cdafe55340
using DrWatson

# ╔═╡ 482eb9fb-0254-40a3-8ccb-82c3fecb2080
# ╠═╡ show_logs = false
@quickactivate "wlp-dual-energy"

# ╔═╡ f3404432-63b2-4655-b4d9-1ae40afc68ab
using PlutoUI, CairoMakie, MAT, DICOM, Statistics, CSV, DataFrames

# ╔═╡ 28b8a6aa-d229-47b7-b103-d45074e13efb
using StatsBase: quantile!, rmsd

# ╔═╡ 9ce934b8-aae9-491d-927e-f71aab15e592
# ╠═╡ show_logs = false
using  CalciumScoring

# ╔═╡ ea3afefa-01e9-4297-8559-bb2edef19ada
using MaterialDecomposition

# ╔═╡ dc427017-af81-405c-bc92-5ac4ccd2b41f
using GLM, MLJBase

# ╔═╡ dae2e015-0a39-4ce3-845a-b9d87989712b
include(srcdir("dicom_utils.jl"))

# ╔═╡ eec36f96-4290-4d07-9df1-45419f8bc195
include(srcdir("masks.jl"))

# ╔═╡ 52d6b7bb-225d-4144-8565-48165c9d0f36
include(srcdir("helper_functions.jl"));

# ╔═╡ 9725eee4-6cd8-4c98-ab1c-e9f4f1b05b48
include(srcdir("plot_utils.jl"));

# ╔═╡ 75872b5b-97f2-4755-ae88-701c91eb129e
TableOfContents()

# ╔═╡ 4618af4a-6fa0-4db4-b460-9cb701e19ecf
md"""
# Calibration
"""

# ╔═╡ b1955f82-70aa-4fe9-aa95-ec93615f2368
begin
	densities_cal = vcat(collect(0:10:90), 95) # percentage lipid
	sizes_cal = [30]
	energies = [80, 135]
end

# ╔═╡ a992c94d-42e3-4732-be30-e54c3f7fe595
begin	
	hu_lipid_array_high_energy = []
	hu_lipid_array_low_energy = []
end

# ╔═╡ 3e736a87-a9be-45e3-bff6-ed720ad7f184
begin
	percentages = [0, 20, 40, 70]

	for percentage in percentages
	
		path = joinpath(datadir("dcms", "cal"),
						string(percentage), 
						"30", 
						string(energies[2]) )
	
		dcm_cal = dcmdir_parse(path)
		dcm_array_cal = load_dcm_array(dcm_cal)
	
		center_insert1, center_insert2 = 187, 318
		offset = 5
		calibration_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array_cal, 3))
		
		for z in axes(dcm_array_cal, 3)
			rows, cols, depth = size(dcm_array_cal)
			half_row, half_col = center_insert1, center_insert2
			row_range = half_row-offset:half_row+offset
			col_range = half_col-offset:half_col+offset	
			calibration_rod[:, :, z] .= dcm_array_cal[row_range, col_range, z];
		end
	
		hu_lipid = mean(calibration_rod)
		push!(hu_lipid_array_high_energy,hu_lipid)
		# ρ_calcium_100 = 0.100 # mg/cm^3
	end
end

# ╔═╡ 7f85ff08-0048-4d01-9e80-dee769b7ec64
begin
	for percentage in percentages
	
		path = joinpath(datadir("dcms", "cal"),
						string(percentage), 
						"30", 
						string(energies[1]) )
	
		dcm_cal = dcmdir_parse(path)
		dcm_array_cal = load_dcm_array(dcm_cal)
	
		center_insert1, center_insert2 = 187, 318
		offset = 5
		calibration_rod = zeros(offset*2 + 1, offset*2 + 1, size(dcm_array_cal, 3))
		
		for z in axes(dcm_array_cal, 3)
			rows, cols, depth = size(dcm_array_cal)
			half_row, half_col = center_insert1, center_insert2
			row_range = half_row-offset:half_row+offset
			col_range = half_col-offset:half_col+offset	
			calibration_rod[:, :, z] .= dcm_array_cal[row_range, col_range, z];
		end
	
		hu_lipid = mean(calibration_rod)
		push!(hu_lipid_array_low_energy,hu_lipid)
		# ρ_calcium_100 = 0.100 # mg/cm^3
	end
end

# ╔═╡ 9bc5706a-44aa-4d9f-a770-386602220692
begin
	hu_background_he = mean(hu_lipid_array_high_energy[1])
	hu_lipid_lp_he = mean(hu_lipid_array_high_energy[2])
	hu_lipid_mp_he = mean(hu_lipid_array_high_energy[3])
	hu_lipid_hp_he = mean(hu_lipid_array_high_energy[4])

	hu_background_le = mean(hu_lipid_array_low_energy[1])
	hu_lipid_lp_le = mean(hu_lipid_array_low_energy[2])
	hu_lipid_mp_le = mean(hu_lipid_array_low_energy[3])
	hu_lipid_hp_le = mean(hu_lipid_array_low_energy[4])
end;

# ╔═╡ f0dc0f2c-b2bd-4170-82a4-6acb49beaac2
begin 
	# Low Energy
	means_80 = Dict(:density => densities_cal, :means => zeros(length(densities_cal)))
	for (i, density) in enumerate(densities_cal)
		for _size in sizes_cal
			path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[1]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
		
			center_insert1, center_insert2 = 187, 318
			
			calibration_rod = zeros(25, 25, size(dcm_array, 3))
			push!(hu_lipid_array_low_energy, calibration_rod)
			
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
	
	# High Energy
	means_135 = Dict(:density => densities_cal, :means => zeros(length(densities_cal)))
	for (i, density) in enumerate(densities_cal)
		for _size in sizes_cal
			path = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energies[2]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
		
			center_insert1, center_insert2 = 187, 318
			
			calibration_rod = zeros(25, 25, size(dcm_array, 3))
			push!(hu_lipid_array_high_energy, calibration_rod)
			
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

	# # hu_lipid_75 = mean(calibration_rod)
	# hu_lipid_lp = mean
	# # mixture_limid_75 = 0.75

	# Fit Parameters
	calculated_intensities = hcat(means_80[:means], means_135[:means]) # low energy, high energy
	ps = fit_calibration(calculated_intensities, densities_cal)
end

# ╔═╡ f2dda856-7f0c-4f6a-ad91-f3dc8a9ed7a7
md"""
## Check Results
"""

# ╔═╡ c81998f3-44cd-48fa-9ff5-f47284f7e41d
begin
	predicted_densities = []
	
	for i in 1:length(densities_cal)
		append!(
			predicted_densities, 
			quantify(calculated_intensities[i, 1], calculated_intensities[i, 2], ps,
			)
		)
	end
end

# ╔═╡ 99c6f606-200d-44a7-bcc4-a3d73b729815
df = DataFrame(
	densities = densities_cal,
	predicted_densities = predicted_densities,
)

# ╔═╡ 27d11ad6-ed5c-41e6-a097-a342a0f52d50
md"""
# Validation
"""

# ╔═╡ 5eb0eef2-1bd3-42c8-aa6d-b264d72000dd
begin
	densities_val = [
		"61_66_71"
		"73_78_82"
		# "85_89_94"
	] # percentage lipid

	sizes_val = ["small", "medium", "large"]
end;

# ╔═╡ 99530752-6e29-44dd-93df-c44b768b7553
energies_val = [80, 135]

# ╔═╡ 07d9df84-c369-444d-b165-ca81da1f1cef
md"""
### Load masks
"""

# ╔═╡ 910228ab-9ac7-4e53-b3e2-eb8f6f95e0b8
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

# ╔═╡ 03272ae3-c8de-47d2-afb6-f7e4291c9515
md"""
## Run Predictions
"""

# ╔═╡ ca3ee165-af7f-49a0-8c8f-7e92004a7b7a
begin
	dfs = []
	dfs_low_vf = []
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

			# Large Insert Masks
			large_erosions = 3
			eroded_mask_L_HD = erode_recursively(mask_L_HD, large_erosions)
			eroded_mask_L_HD_3D = cat(eroded_mask_L_HD, eroded_mask_L_HD, eroded_mask_L_HD, dims=3)

			eroded_mask_L_MD = erode_recursively(mask_L_MD, large_erosions)
			eroded_mask_L_MD_3D = cat(eroded_mask_L_MD, eroded_mask_L_MD, eroded_mask_L_MD, dims=3)

			eroded_mask_L_LD = erode_recursively(mask_L_LD, large_erosions)
			eroded_mask_L_LD_3D = cat(eroded_mask_L_LD, eroded_mask_L_LD, eroded_mask_L_LD, dims=3)

			# Medium Insert Masks
			medium_erosions = 1
			eroded_mask_M_HD = erode_recursively(mask_M_HD, medium_erosions)
			eroded_mask_M_HD_3D = cat(eroded_mask_M_HD, eroded_mask_M_HD, eroded_mask_M_HD, dims=3)
		
			eroded_mask_M_MD = erode_recursively(mask_M_MD, medium_erosions)
			eroded_mask_M_MD_3D = cat(eroded_mask_M_MD, eroded_mask_M_MD, eroded_mask_M_MD, dims=3)
		
			eroded_mask_M_LD = erode_recursively(mask_M_LD, medium_erosions)
			eroded_mask_M_LD_3D = cat(eroded_mask_M_LD, eroded_mask_M_LD, eroded_mask_M_LD, dims=3)

			#------- Material Decomposition -------#

			# Low Energy
			path_80 = datadir("dcms", "val", density, _size, string(energies[1]))
			dcm_80 = dcmdir_parse(path_80)
			dcm_array_80 = load_dcm_array(dcm_80)
			pixel_size = get_pixel_size(dcm_80[1].meta)
			
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
				predicted_densities[i] = quantify(means_80[i], means_80[i], ps)
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


			
			#------- Volume Fraction -------#
			# 80 KeV
			path = datadir("dcms", "val", density, _size, string(energies[1]))
			dcm = dcmdir_parse(path)
			dcm_array = load_dcm_array(dcm)
			voxel_size = pixel_size[1] * pixel_size[2] * pixel_size[3]

			# Score Inserts
			# hp stands for high percentage
			hu_background = hu_background_le
			hu_lipid_hp = hu_lipid_hp_le
			hu_lipid_mp = hu_lipid_mp_le
			hu_lipid_lp = hu_lipid_hp_le

			#finding percentage of large size, high percentage
			percentage_lipid_large_hp = score(dcm_array[eroded_mask_L_HD_3D], hu_lipid_hp, hu_background, VolumeFraction())
			
			#finding percentage of large size, medium percentage
			percentage_lipid_large_mp = score(dcm_array[eroded_mask_L_MD_3D], hu_lipid_mp, hu_background, VolumeFraction())

			#finding percentage of large size, low percentage
			percentage_lipid_large_lp = score(dcm_array[eroded_mask_L_LD_3D], hu_lipid_lp, hu_background, VolumeFraction())

			#finding percentage of medium size, high percentage
			percentage_lipid_medium_hp = score(dcm_array[eroded_mask_M_HD_3D], hu_lipid_hp, hu_background, VolumeFraction())

			#finding percentage of medium size, medium percentage
			percentage_lipid_medium_mp = score(dcm_array[eroded_mask_M_MD_3D], hu_lipid_mp, hu_background, VolumeFraction())

			#finding percentage of medium size, low percentage
			percentage_lipid_medium_lp = score(dcm_array[eroded_mask_M_LD_3D], hu_lipid_lp, hu_background, VolumeFraction())

			

			predicted_percentage_large_inserts = [
				percentage_lipid_large_hp, 
				percentage_lipid_large_mp, 
				percentage_lipid_large_lp
			]
	
			predicted_percentage_medium_inserts = [
				percentage_lipid_medium_hp, 
				percentage_lipid_medium_mp, 
				percentage_lipid_medium_lp
			]

			df_results = DataFrame(
				phantom_size = _size,
				density = density,
				insert_densities = [:low_percentage, :medium_percentage, :high_percentage],
				ground_truth_density_large_inserts = parse.(Int, split(density, "_")),
				predicted_percentage_large_inserts = predicted_percentage_large_inserts,
				
				ground_truth_density_medium_inserts = parse.(Int, split(density, "_")),
				predicted_mass_medium_inserts = predicted_percentage_medium_inserts,
			)
			push!(dfs_low_vf, df_results)

			
		end
	end
end

# ╔═╡ cc661c6e-3deb-49fe-a119-49c21b507e98
dfs

# ╔═╡ 771d3743-e1b0-4924-af9f-f63e2c63b3ff
dfs_low_vf

# ╔═╡ 12fa26dc-d3ab-4789-a608-2f72b56d9f21
md"""
# Analyze
"""

# ╔═╡ 60febdca-124e-4455-a3ae-9afb93189566
medphys_theme = Theme(
    Axis = (
        backgroundcolor = :white,
		xgridcolor = :gray,
		xgridwidth = 0.1,
		xlabelsize = 15,
		xticklabelsize = 15,
		ygridcolor = :gray,
		ygridwidth = 0.1,
		ylabelsize = 15,
		yticklabelsize = 15,
		bottomsplinecolor = :black,
		leftspinecolor = :black,
		titlesize = 25
	)
);

# ╔═╡ 87f4a121-e385-4d20-92df-82bba4eff328
md"""
### Accuracy
"""

# ╔═╡ 54140655-4dc8-4e53-96e0-821c2ca1e515
new_df_low_md = vcat(dfs[1:length(dfs)]...);

# ╔═╡ 1be32f5e-5c85-4118-bf5c-e53dc12809a1
new_df_low_vf = vcat(dfs_low_vf[1:length(dfs_low_vf)]...);

# ╔═╡ 32c1a063-94c5-4ff1-8057-122d24dd1923
new_df_low_md

# ╔═╡ ad06c138-05db-4226-85f2-2391c4916b6b
new_df_low_vf

# ╔═╡ 53115e78-9bcf-48ec-b76a-1a2d5bcf97cd
co_1_low_md, r_squared_1_low_md, rms_values_1_low_md, pred_1_low_md = calculate_coefficients(new_df_low_md);

# ╔═╡ d87751cd-f639-4b67-94a0-0590710c05b5
co_1_low_vf, r_squared_1_low_vf, rms_values_1_low_vf, pred_1_low_vf = calculate_coefficients_vf(new_df_low_vf);

# ╔═╡ 70951b1e-f0e7-40af-bb71-dc1366c33601
function accuracy()
	f = Figure()

	##-- A --##
	ax = Axis(
		f[1, 1],
		xticks = collect(-5:10:150),
		yticks = collect(-5:10:150),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Material Decomposition (Low Density)",
	)
	
	df = new_df_low_md
	sc1 = scatter!(
		df[!, :ground_truth_density_large_inserts], df[!, :predicted_densities_large_inserts]
	)
	sc2 = scatter!(
		df[!, :ground_truth_density_medium_inserts], df[!, :predicted_densities_medium_inserts]
	)
	# sc3 = scatter!(
	# 	df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	# )
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_low_md, linestyle=:dashdot)
	create_textbox(f[1, 1], co_1_low_md, r_squared_1_low_md, rms_values_1_low_md)
	
	xlims!(ax, low=-5, high=150)
	ylims!(ax, low=-5, high=150)

	ax = Axis(
		f[2, 1],
		xticks = collect(-5:10:200),
		yticks = collect(-5:10:200),
		xlabel = "Known Mass (mg)",
		ylabel = "Calculated Mass (mg)",
		title = "Volume Fraction (Low Density)",
	)
	
	df = new_df_low_vf
	sc1 = scatter!(
		df[!, :ground_truth_density_large_inserts], df[!, :predicted_percentage_large_inserts]
	)
	sc2 = scatter!(
		df[!, :ground_truth_density_medium_inserts], df[!, :predicted_mass_medium_inserts]
	)
	# sc3 = scatter!(
	# 	df[!, :gt_mass_small_inserts], df[!, :predicted_mass_small_inserts], color=:red
	# )
	ln1 = lines!([-1000, 1000], [-1000, 1000])
	ln2 = lines!(collect(1:1000), pred_1_low_vf, linestyle=:dashdot)
	create_textbox(f[2, 1], co_1_low_vf, r_squared_1_low_vf, rms_values_1_low_vf)
	
	xlims!(ax, low=-5, high=150)
	ylims!(ax, low=-5, high=200)


	#-- LABELS --##
	f[2, 2] = Legend(f, [sc1, sc2, ln1, ln2], ["Large Inserts", "Medium Inserts", "Unity", "Fitted Line"], framevisible = false)

	
	for (label, layout) in zip([["A"], ["B"]], [f[1,1], f[2, 1]])
	    Label(layout[1, 1, TopLeft()], label,
	        fontsize = 25,
	        padding = (0, 60, 25, 0),
	        halign = :right)
	end
	f
end

# ╔═╡ 128155d9-34c0-4f59-ada3-53e3d670c274
with_theme(accuracy, medphys_theme)

# ╔═╡ Cell order:
# ╠═51258732-faf5-4705-82b6-64cdafe55340
# ╠═482eb9fb-0254-40a3-8ccb-82c3fecb2080
# ╠═f3404432-63b2-4655-b4d9-1ae40afc68ab
# ╠═28b8a6aa-d229-47b7-b103-d45074e13efb
# ╠═9ce934b8-aae9-491d-927e-f71aab15e592
# ╠═ea3afefa-01e9-4297-8559-bb2edef19ada
# ╠═dae2e015-0a39-4ce3-845a-b9d87989712b
# ╠═75872b5b-97f2-4755-ae88-701c91eb129e
# ╠═eec36f96-4290-4d07-9df1-45419f8bc195
# ╟─4618af4a-6fa0-4db4-b460-9cb701e19ecf
# ╠═b1955f82-70aa-4fe9-aa95-ec93615f2368
# ╠═a992c94d-42e3-4732-be30-e54c3f7fe595
# ╠═3e736a87-a9be-45e3-bff6-ed720ad7f184
# ╠═7f85ff08-0048-4d01-9e80-dee769b7ec64
# ╠═9bc5706a-44aa-4d9f-a770-386602220692
# ╠═f0dc0f2c-b2bd-4170-82a4-6acb49beaac2
# ╟─f2dda856-7f0c-4f6a-ad91-f3dc8a9ed7a7
# ╠═c81998f3-44cd-48fa-9ff5-f47284f7e41d
# ╠═99c6f606-200d-44a7-bcc4-a3d73b729815
# ╟─27d11ad6-ed5c-41e6-a097-a342a0f52d50
# ╠═5eb0eef2-1bd3-42c8-aa6d-b264d72000dd
# ╠═99530752-6e29-44dd-93df-c44b768b7553
# ╟─07d9df84-c369-444d-b165-ca81da1f1cef
# ╠═910228ab-9ac7-4e53-b3e2-eb8f6f95e0b8
# ╟─03272ae3-c8de-47d2-afb6-f7e4291c9515
# ╠═ca3ee165-af7f-49a0-8c8f-7e92004a7b7a
# ╠═cc661c6e-3deb-49fe-a119-49c21b507e98
# ╠═771d3743-e1b0-4924-af9f-f63e2c63b3ff
# ╟─12fa26dc-d3ab-4789-a608-2f72b56d9f21
# ╠═dc427017-af81-405c-bc92-5ac4ccd2b41f
# ╠═52d6b7bb-225d-4144-8565-48165c9d0f36
# ╠═9725eee4-6cd8-4c98-ab1c-e9f4f1b05b48
# ╠═60febdca-124e-4455-a3ae-9afb93189566
# ╟─87f4a121-e385-4d20-92df-82bba4eff328
# ╠═54140655-4dc8-4e53-96e0-821c2ca1e515
# ╠═1be32f5e-5c85-4118-bf5c-e53dc12809a1
# ╠═32c1a063-94c5-4ff1-8057-122d24dd1923
# ╠═ad06c138-05db-4226-85f2-2391c4916b6b
# ╠═53115e78-9bcf-48ec-b76a-1a2d5bcf97cd
# ╠═d87751cd-f639-4b67-94a0-0590710c05b5
# ╠═70951b1e-f0e7-40af-bb71-dc1366c33601
# ╠═128155d9-34c0-4f59-ada3-53e3d670c274
