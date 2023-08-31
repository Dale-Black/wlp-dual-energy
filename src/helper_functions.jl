function collect_tuple(tuple_array)
    row_num = size(tuple_array)
    col_num = length(tuple_array[1])
    container = zeros(Int64, row_num..., col_num)
    for i in 1:length(tuple_array)
        container[i, :] = collect(tuple_array[i])
    end
    return container
end

function calculate_coefficients(df;
    label1=:ground_truth_density_large_inserts,
    label2=:ground_truth_density_medium_inserts,
    label4=:predicted_densities_large_inserts,
    label5=:predicted_densities_medium_inserts,
)
    gt_array = vec(hcat(df[!, label1], df[!, label2]))
    calc_array = vec(hcat(df[!, label4], df[!, label5]))
    data = DataFrame(X=gt_array, Y=calc_array)
    model = lm(@formula(Y ~ X), data)
    r_squared = GLM.r2(model)
    rms_values = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model))
    ]
    pred = GLM.predict(model, DataFrame(X=collect(1:1000)))

    return coef(model), r_squared, rms_values, pred
end

function calculate_coefficients_vf(df;
    label1=:ground_truth_density_large_inserts,
    label2=:ground_truth_density_medium_inserts,
    label4=:predicted_percentage_large_inserts,
    label5=:predicted_mass_medium_inserts,
)
    gt_array = vec(hcat(df[!, label1], df[!, label2]))
    calc_array = vec(hcat(df[!, label4], df[!, label5]))
    data = DataFrame(X=gt_array, Y=calc_array)
    model = lm(@formula(Y ~ X), data)
    r_squared = GLM.r2(model)
    rms_values = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model))
    ]
    pred = GLM.predict(model, DataFrame(X=collect(1:1000)))

    return coef(model), r_squared, rms_values, pred
end

function calculate_coefficients_ld(df;  
    label1=:ground_truth_mass_ld,
    label2=:predicted_mass_ld
)
    gt_array = vec(hcat(df[!, label1]))
    calc_array = vec(hcat(df[!, label2]))
    data = DataFrame(X=gt_array, Y=calc_array)
    model = lm(@formula(Y ~ X), data)
    r_squared = GLM.r2(model)
    rms_values = [
        rms(data[!, :X], data[!, :Y]),
        rmsd(data[!, :Y], GLM.predict(model))
    ]
    pred = GLM.predict(model, DataFrame(X=collect(1:1000)))

    return coef(model), r_squared, rms_values, pred
end

function overlay_mask_bind(mask)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    return PlutoUI.Slider(1:length(zs); default=3, show_value=true)
end

function overlay_mask_plot(array, mask, var, title::AbstractString)
    indices = findall(x -> x == 1, mask)
    indices = Tuple.(indices)
    label_array = collect_tuple(indices)
    zs = unique(label_array[:, 3])
    indices_lbl = findall(x -> x == zs[var], label_array[:, 3])

    fig = Figure()
    ax = Makie.Axis(fig[1, 1])
    ax.title = title
    heatmap!(array[:, :, zs[var]]; colormap=:grays)
    scatter!(
        label_array[:, 1][indices_lbl],
        label_array[:, 2][indices_lbl];
        markersize=1,
        color=:red
    )
    return fig
end

function show_matrix(A::Matrix, red::Union{BitMatrix,Matrix{Bool}}=zeros(Bool, size(A)))
    base = RGB.(Gray.(A))

    base[red] .= RGB(1.0, 0.1, 0.1)

    # some tricks to show the pixels more clearly:
    s = max(size(A)...)
    if s >= 20
        min_size = 1200
        factor = min(5, min_size ÷ s)

        kron(base, ones(factor, factor))
    else
        base
    end
end

function create_mask(array, mask)
    @assert size(array) == size(mask)
    idxs = findall(x -> x == true, mask)
    overlayed_mask = zeros(size(array))
    for idx in idxs
        overlayed_mask[idx] = array[idx]
    end
    return overlayed_mask
end

function predict_concentration(x, y, p)
    A = p[1] + (p[2] * x) + (p[3] * y) + (p[4] * x^2) + (p[5] * x * y) + (p[6] * y^2)
    B = 1 + (p[7] * x) + (p[8] * y)
    F = A / B
end