### A Pluto.jl notebook ###
# v0.19.26

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

# ╔═╡ dfe2e586-c08e-4abe-8f2c-d0b80c56de2d
using DrWatson

# ╔═╡ 5eee0e35-490f-4436-8c82-e0785a983b39
# ╠═╡ show_logs = false
@quickactivate "wlp-dual-energy"

# ╔═╡ 1bf5f2a3-e0e0-4e90-b665-e7ca4cffb1aa
using PlutoUI, CairoMakie, MAT, DICOM, DICOMUtils

# ╔═╡ b51a3623-a24a-45f8-a544-5ba17bb4fde3
TableOfContents()

# ╔═╡ 09012a97-7c3a-46a4-ae42-8f6315894e4c
md"""
# Calibration
"""

# ╔═╡ ebcd93f4-adff-42d8-8272-1fe758d0e47e
densities_cal = vcat(collect(0:10:90), 95) # percentage lipid

# ╔═╡ e42245dc-0620-427b-9d03-936e8171ef2f
energies_cal = [80, 135]

# ╔═╡ 7e006f7d-37cd-4001-8c46-3dcee8f4f18b
sizes_cal = [30]

# ╔═╡ ff597f41-7c8e-4edb-989b-30bcc046e732
for density in densities_cal
	for energy in energies_cal
		for _size in sizes_cal
			# Convert first slice into dcm
			path1 = joinpath(datadir("mats", "cal"), string(density, "rod", energy, "kV", _size, "_1.mat"))
			var1 = matread(path1)
			img1 = var1[string("I")]
			img1 = Int16.(round.(img1))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img1
			dcm[tag"Instance Number"] = 1
			dcm[tag"Rows"] = size(img1, 1)
			dcm[tag"Columns"] = size(img1, 2)

			output_dir = joinpath(datadir("dcms", "cal"), string(density), string(_size), string(energy))
			if !isdir(output_dir)
				mkpath(output_dir)
			end

			output_path1 = joinpath(output_dir, "1.dcm") 
			dcm_write(output_path1, dcm)

			# Convert second slice into dcm
			path2 = joinpath(datadir("mats", "cal"), string(density, "rod", energy, "kV", _size, "_2.mat"))
			var2 = matread(path2)
			img2 = var2[string("I")]
			img2 = Int16.(round.(img2))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img2
			dcm[tag"Instance Number"] = 2
			dcm[tag"Rows"] = size(img2, 1)
			dcm[tag"Columns"] = size(img2, 2)

			output_path2 = joinpath(output_dir, "2.dcm") 
			dcm_write(output_path2, dcm)

			# Convert third slice into dcm
			path3 = joinpath(datadir("mats", "cal"), string(density, "rod", energy, "kV", _size, "_3.mat"))
			var3 = matread(path3)
			img3 = var3[string("I")]
			img3 = Int16.(round.(img3))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img3
			dcm[tag"Instance Number"] = 2
			dcm[tag"Rows"] = size(img3, 1)
			dcm[tag"Columns"] = size(img3, 2)

			global output_path3_cal = joinpath(output_dir, "3.dcm") 
			dcm_write(output_path3_cal, dcm)
		end
	end
end

# ╔═╡ 34be1b53-e6d5-468b-9efc-8b5d678fd8f3
md"""
## Check DICOM Calibration
"""

# ╔═╡ 879fe54d-8d21-41dc-8661-6873092eb8f6
dcmdir_combined = dcmdir_parse(dirname(output_path3_cal));

# ╔═╡ aeedbc91-c64f-4088-9fb7-5bd1cf5eb80a
vol_combined = load_dcm_array(dcmdir_combined);

# ╔═╡ af612106-012f-4121-a769-0ca476164bc6
@bind z PlutoUI.Slider(1:size(vol_combined, 3); default=1, show_value=true)

# ╔═╡ e376957c-1e0b-4e5e-94d9-9c8c828b951c
heatmap(transpose(vol_combined[:, :, z]); colormap=:grays)

# ╔═╡ f43157c9-7403-46da-b9a2-6be1b42bfdba
md"""
# Validation
"""

# ╔═╡ 014ce041-e042-4925-9d0b-2198711b909a
densities_val = [
	"61_66_71",
	"73_78_82",
	"85_89_94"
] # percentage lipid

# ╔═╡ 0a4eebc2-c4bb-4744-9964-eecf7f5f9e2c
energies_val = [80, 135]

# ╔═╡ 6bfa5d93-5607-424d-81ab-3b3b082d4a2d
sizes_val = ["small", "medium", "large"]

# ╔═╡ abccce08-cd6d-4f5c-af77-2412dcfa173b
for (i, density) in enumerate(densities_val)
	for energy in energies_val
		for _size in sizes_val
			# Convert first slice into dcm
			path1 = joinpath(datadir("mats", "val"), string(density, "energy", energy, _size, "_1.mat"))
			var1 = matread(path1)
			img1 = var1[string("I")]
			img1 = Int16.(round.(img1))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img1
			dcm[tag"Instance Number"] = 1
			dcm[tag"Rows"] = size(img1, 1)
			dcm[tag"Columns"] = size(img1, 2)
			
			output_dir = joinpath(datadir("dcms", "val"), string(density), string(_size), string(energy))
			if !isdir(output_dir)
				mkpath(output_dir)
			end

			output_path1 = joinpath(output_dir, "1.dcm") 
			dcm_write(output_path1, dcm)

			# Convert second slice into dcm
			path2 = joinpath(datadir("mats", "val"), string(density, "energy", energy, _size, "_2.mat"))
			var2 = matread(path2)
			img2 = var2[string("I")]
			img2 = Int16.(round.(img2))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img2
			dcm[tag"Instance Number"] = 1
			dcm[tag"Rows"] = size(img2, 1)
			dcm[tag"Columns"] = size(img2, 2)

			output_path2 = joinpath(output_dir, "2.dcm") 
			dcm_write(output_path2, dcm)

			# Convert second slice into dcm
			path3 = joinpath(datadir("mats", "val"), string(density, "energy", energy, _size, "_3.mat"))
			var3 = matread(path3)
			img3 = var3[string("I")]
			img3 = Int16.(round.(img3))

			dcm_path = datadir("sample.dcm")
			dcm = dcm_parse(dcm_path)
			dcm[tag"Pixel Data"] = img3
			dcm[tag"Instance Number"] = 1
			dcm[tag"Rows"] = size(img3, 1)
			dcm[tag"Columns"] = size(img3, 2)

			global output_path3_val = joinpath(output_dir, "3.dcm") 
			dcm_write(output_path3_val, dcm)
		end
	end
end

# ╔═╡ 36e3a167-00bd-4aa8-a3a1-1aae18bb63a9
md"""
## Check DICOM Validation
"""

# ╔═╡ 085216eb-4501-4c59-9a36-49318146a751
dcmdir_combined_val = dcmdir_parse(dirname(output_path3_val));

# ╔═╡ a4a4d861-e9f4-491d-a0ce-36fb5d96afdb
vol_combined_val = load_dcm_array(dcmdir_combined_val);

# ╔═╡ 817fe66d-f9ab-4191-afc0-aa2bc685152a
@bind z1 PlutoUI.Slider(1:size(vol_combined_val, 3); default=1, show_value=true)

# ╔═╡ 72f4c923-0814-4a49-ace1-c891cb9df12a
heatmap(transpose(vol_combined_val[:, :, z1]); colormap=:grays)

# ╔═╡ Cell order:
# ╠═dfe2e586-c08e-4abe-8f2c-d0b80c56de2d
# ╠═5eee0e35-490f-4436-8c82-e0785a983b39
# ╠═1bf5f2a3-e0e0-4e90-b665-e7ca4cffb1aa
# ╠═b51a3623-a24a-45f8-a544-5ba17bb4fde3
# ╟─09012a97-7c3a-46a4-ae42-8f6315894e4c
# ╠═ebcd93f4-adff-42d8-8272-1fe758d0e47e
# ╠═e42245dc-0620-427b-9d03-936e8171ef2f
# ╠═7e006f7d-37cd-4001-8c46-3dcee8f4f18b
# ╠═ff597f41-7c8e-4edb-989b-30bcc046e732
# ╟─34be1b53-e6d5-468b-9efc-8b5d678fd8f3
# ╠═879fe54d-8d21-41dc-8661-6873092eb8f6
# ╠═aeedbc91-c64f-4088-9fb7-5bd1cf5eb80a
# ╟─af612106-012f-4121-a769-0ca476164bc6
# ╟─e376957c-1e0b-4e5e-94d9-9c8c828b951c
# ╟─f43157c9-7403-46da-b9a2-6be1b42bfdba
# ╠═014ce041-e042-4925-9d0b-2198711b909a
# ╠═0a4eebc2-c4bb-4744-9964-eecf7f5f9e2c
# ╠═6bfa5d93-5607-424d-81ab-3b3b082d4a2d
# ╠═abccce08-cd6d-4f5c-af77-2412dcfa173b
# ╟─36e3a167-00bd-4aa8-a3a1-1aae18bb63a9
# ╠═085216eb-4501-4c59-9a36-49318146a751
# ╠═a4a4d861-e9f4-491d-a0ce-36fb5d96afdb
# ╟─817fe66d-f9ab-4191-afc0-aa2bc685152a
# ╟─72f4c923-0814-4a49-ace1-c891cb9df12a
