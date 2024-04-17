import DICOM
"""
    load_dcm_array(dcm_data::Vector{DICOM.DICOMData})
Given some DICOM.DICOMData, `load_dcm_array` loads the pixels
of each slice into a 3D array and returns the array
"""
function load_dcm_array(dcm_data::Vector{DICOM.DICOMData})
    return array = cat(
        [dcm_data[i][(0x7fe0, 0x0010)] for i in 1:length(dcm_data)]...; dims=3
    )
end

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