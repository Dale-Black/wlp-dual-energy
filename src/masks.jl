using ImageMorphology: erode, dilate

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

# ╔═╡ 603f5390-bb65-4854-9af0-abaf8e980f9b
function ring_mask_small(dilated_mask)
    return Bool.(dilate_recursively(dilated_mask, 6) - dilated_mask)
end

# ╔═╡ 78a5a6b9-8b13-43e7-9672-5459afd0f0b9
function dilate_mask_large_bkg(mask)
    return dilate_recursively(mask, 2)
end

function dilate_mask_medium_bkg(mask)
    return dilate_recursively(mask, 1)
end

function dilate_mask_small_bkg(mask)
    return mask
end