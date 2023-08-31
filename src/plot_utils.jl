using Printf
"""
    create_textbox(loc, coef, r, rms)

Create a textbox, inside a Makie.jl plot that accounts for the sign of the coefficient

#### Inputs
- `loc`: figure location, e.g. `loc=f[1, 1]`
- `coef`: coefficients for the best fit line
- `r`: r squared values
- `rms`: rms values
"""
function create_textbox(loc, coef, r, rms; p=5)
    local tb
    if coef[1] > 0
        tb = Label(
            loc;
            text="y = $(@sprintf "%.2f" coef[2])x+$(@sprintf "%.2f" coef[1]) \nr = $(@sprintf "%.2f" r) \nRMSE: $(@sprintf "%.2f" rms[1]) \nRMSD: $(@sprintf "%.2f" rms[2])",
            padding=(p, p, p, p),
            tellheight=false,
            tellwidth=false,
            halign=:left,
            justification=:left,
            valign=:top,
            fontsize=12
        )
    else
        tb = Label(
            loc,
            text="y = $(@sprintf "%.2f" coef[2])x$(@sprintf "%.2f" coef[1]) \nr = $(@sprintf "%.2f" r) \nRMSE: $(@sprintf "%.2f" rms[1]) \nRMSD: $(@sprintf "%.2f" rms[2])",
            padding=(p, p, p, p),
            tellheight=false,
            tellwidth=false,
            halign=:left,
            justification=:left,
            valign=:top,
            fontsize=12
        )
    end
    return tb
end