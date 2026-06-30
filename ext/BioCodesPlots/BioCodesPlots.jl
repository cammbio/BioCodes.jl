module BioCodesPlots

import Plots
using BioCodes
using BioSymbols
using Colors: RGB

# Same color assignments as BioCodes.GeneticCodes:
const _aa_colors = BioCodes.amino_acid_colors

_bg_rgb(label) = begin
    c = get(_aa_colors, label, (200, 200, 200))
    RGB(c[1] / 255, c[2] / 255, c[3] / 255)
end

# Choose black or white text based on perceived luminance of the cell background.
_fg_color(label) = begin
    c = get(_aa_colors, label, (200, 200, 200))
    L = 0.299 * c[1] + 0.587 * c[2] + 0.114 * c[3]
    L < 140 ? :white : :black
end

_rect(x0, y0, x1, y1) = Plots.Shape([x0, x1, x1, x0], [y0, y0, y1, y1])

"""
    Plots.plot(gc::GeneticCode; title=nothing, kwargs...) -> Plots.Plot

Plot a `GeneticCode` as a colored table in the standard textbook layout:

- **Top axis** – second base of the codon/tuple.
- **Left axis** – first base, one label per group of N rows (N = alphabet size).
- **Right axis** – third base, one label per row (only for tuple length 3).
- **Cells** – amino acid (or generic label) with a distinctive background color.

The layout reproduces the appearance of classical genetic code tables in the
literature (e.g. T/U rows × C/A/G columns × T/U/C/A/G inner rows). Only tuple
lengths 2 and 3 are supported.

# Arguments
- `gc`    – a `BioCodes.GeneticCode` object.
- `title` – optional plot title; a descriptive default is generated if omitted.
- Any remaining keyword arguments are forwarded to `Plots.plot`.

# Example
```julia
using BioCodes, BioCodesPlots, Plots
plot(StandardGeneticCode())
plot(StandardGeneticCode(; alphabet=RNAAlphabet))
plot(GeneticCode(2, [DNA_T, DNA_C, DNA_A, DNA_G]))
```
"""
function Plots.plot(gc::GeneticCode; title::Union{String,Nothing}=nothing, kwargs...)
    l = tuple_length(gc)
    if l < 2 || l > 3
        throw(ArgumentError("BioCodesPlots.plot: unsupported tuple length $l " *
                            "(only 2 and 3 are supported)"))
    end

    M = Matrix(gc)          # Matrix{GeneticCodeCell}, row 1 = top
    nrows, ncols = size(M)
    N = length(alphabet(gc))
    Σ = alphabet(gc)        # alphabet in display order

    # ── Coordinate system ────────────────────────────────────────────────────
    # x: [0, label_w]           left label column  (1st base)
    #    [label_w, label_w+ncols] data columns (one unit each)
    #    [label_w+ncols, total_w] right label column (3rd base, l==3 only)
    # y: [0, nrows]             data rows (row i of M at y ∈ [nrows-i, nrows-i+1])
    #    [nrows, total_h]       top header row (2nd base)
    # ─────────────────────────────────────────────────────────────────────────
    label_w = 1.0
    header_h = 1.0
    total_w = label_w + ncols + label_w
    total_h = header_h + nrows

    cell_cx(j) = label_w + j - 0.5          # x-center of column j
    cell_cy(i) = nrows - i + 0.5            # y-center of row i

    Σ_str = join(string.(Σ), ", ")
    plot_title = isnothing(title) ? "Genetic Code" : title

    # Pixels per coordinate unit - drives the output raster size.
    px_x = 70
    px_y = 30
    p = Plots.plot(;
        xlims=(0, total_w),
        ylims=(0, total_h),
        framestyle=:none,
        legend=false,
        ticks=nothing,
        size=(round(Int, total_w * px_x), round(Int, total_h * px_y)),
        title=plot_title,
        titlefontsize=11,
        kwargs...
    )

    # Data cells
    for i in 1:nrows, j in 1:ncols
        cell = M[i, j]
        x0 = label_w + j - 1
        y0 = nrows - i
        Plots.plot!(p, _rect(x0, y0, x0 + 1, y0 + 1);
            fillcolor=_bg_rgb(cell.label),
            linecolor=:white,
            linewidth=0.5,
            label=nothing)
        Plots.annotate!(p, cell_cx(j), cell_cy(i),
            Plots.text(string(cell.label), 11, :center, _fg_color(cell.label)))
    end

    # Column headers: 2nd base
    for j in 1:ncols
        Plots.annotate!(p, cell_cx(j), nrows + header_h / 2,
            Plots.text(string(Σ[j]), 11, :center, :black))
    end

    if l == 3
        # Left row headers: 1st base (one label centred on each group of N rows)
        for k in 1:N
            # Group k covers rows (k-1)*N+1 … k*N.
            # y-range of that group: [nrows - k*N,  nrows - (k-1)*N]
            cy = nrows - (k - 1) * N - N / 2
            Plots.annotate!(p, label_w / 2, cy,
                Plots.text(string(Σ[k]), 11, :center, :black))
        end

        # Right row headers: 3rd base (one label per row)
        for i in 1:nrows
            i3 = ((i - 1) % N) + 1
            Plots.annotate!(p, label_w + ncols + label_w / 2, cell_cy(i),
                Plots.text(string(Σ[i3]), 11, :center, :black))
        end

        # Horizontal separator lines between 1st-base groups
        for k in 1:(N-1)
            y = Float64(nrows - k * N)
            Plots.plot!(p, [label_w, label_w + ncols], [y, y];
                color=:black,
                linewidth=1.5,
                label=nothing)
        end

    elseif l == 2
        # Left row headers: 1st base (one label per row)
        for i in 1:nrows
            i1 = ((i - 1) % N) + 1
            Plots.annotate!(p, label_w / 2, cell_cy(i),
                Plots.text(string(Σ[i1]), 11, :center, :black))
        end
    end

    return p
end

end # module
