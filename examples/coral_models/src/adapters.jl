#=
adapters.jl -- one function per model: the model's OWN output files  ->  the "standard layout" ADRIAIndicators needs.

    read_cscape(dir, reef)     C~scape: Array_scenario_*.rds (R `saveRDS` of the 6-D out_array) + size_classes.rds
    read_coralpop(dir, reef)   CoralPop: scenario_*.npz (the engine's returned arrays, `np.savez`)

This is the only model-specific code in the example. Both return the same NamedTuple (see `check_model` in
pipeline.jl); README.md explains every step.
=#

using RData, NPZ, DelimitedFiles

# ===== The reef both models were run on =====

"""
    read_reef(path) -> NamedTuple

The one-reef domain the models ran on: location id, total reef area (m2) and habitable ("k") area (m2). Neither
model's output file contains these, and ADRIAIndicators needs both, so they are read from the domain description.
"""
function read_reef(path::AbstractString)
    data, header = readdlm(path, ','; header=true)
    col(name) = data[:, findfirst(==(name), vec(header))]
    return (; ids=String.(col("location_id")), reef_area=Float64.(col("reef_area_m2")), k_area=Float64.(col("k_area_m2")))
end

# ===== Planar-area parameters (shelter volume) =====

## log(S) = b + a log(x), (b, a) per ADRIA functional group -- Aston et al. (2022), as ADRIA.planar_area_params().
## Rows: tabular Acropora, corymbose Acropora, corymbose non-Acropora, small massives, large massives.
const PLANAR_AREA = [
    -8.95 2.80
    -9.13 2.94
    -8.90 2.94
    -8.87 2.30
    -8.87 2.30
]

"[groups, sizes, 2] planar-area parameters when group `g` uses ADRIA row `rows[g]`, the same for every size class."
planar_area_for(rows, n_sizes) = repeat(reshape(PLANAR_AREA[rows, :], length(rows), 1, 2), 1, n_sizes, 1)

# ===== C~scape =====

const CSCAPE_TYPES = ["acro_table", "acro_corym", "corym_non_acro", "small_massive", "large_massive", "brooder"]
const CSCAPE_SIZE_SLOTS = 7:106        # value slots 7..106 of out_array hold the 100 size classes (colony COUNTS)

"""
    read_cscape(dir, reef; years) -> NamedTuple

Read every `Array_scenario_*.rds` in `dir` (one file per run, sorted by scenario id).
"""
function read_cscape(dir::AbstractString, reef; years=2025:2099)
    files = filter(f -> occursin(r"^Array_scenario_.*\.rds$", f), readdir(dir))
    sort!(files; by=f -> parse(Int, match(r"scenario_(\d+)", f).captures[1]))
    sc = RData.load(joinpath(dir, "size_classes.rds"))
    d_cm = permutedims(reduce(hcat, [Float64.(sc[ft]["meshpoints_diam"]) for ft in CSCAPE_TYPES]))   # [6, 100]
    juvenile = permutedims(reduce(hcat, [Bool.(sc[ft]["juvenile"]) for ft in CSCAPE_TYPES]))         # [6, 100]
    G, S = size(d_cm)
    L, N = length(reef.ids), length(files)

    colony_area_m2 = π .* (d_cm ./ 200) .^ 2                       # diameter in cm -> radius in m -> area in m2
    cover = zeros(length(years), G, S, L, N)
    for (n, f) in enumerate(files)
        out = RData.load(joinpath(dir, f))                         # [year, site, intervened, type, heat class, value]
        size(out) == (length(years) + 1, L, 2, G, 14, 107) || error("$f has unexpected shape $(size(out))")
        ## counts: skip slot 1 (2024), keep the 100 size classes, add the two `intervened` sub-areas and the 14 heat classes
        counts = dropdims(sum(out[2:end, :, :, :, :, CSCAPE_SIZE_SLOTS]; dims=(3, 5)); dims=(3, 5))   # [T, L, G, S]
        counts = permutedims(counts, (1, 3, 4, 2))                                                    # [T, G, S, L]
        cover[:, :, :, :, n] = clamp.(counts .* reshape(colony_area_m2, 1, G, S, 1) ./ reshape(reef.k_area, 1, 1, 1, L), 0, 1)
    end
    return (; label="C~scape", years=collect(years), groups=CSCAPE_TYPES, locations=reef.ids, cover,
        k_area=reef.k_area, reef_area=reef.reef_area, diam_cm=repeat(d_cm, 1, 1, N),
        pap=planar_area_for([1, 2, 3, 4, 5, 3], S), is_juvenile=juvenile)          # brooders use the corymbose non-Acropora row
end

# ===== CoralPop =====

const CORALPOP_MORPHS = ["Branching", "Foliose", "Other"]
const CORALPOP_BIN_CM = 5.0            # 20 size bins of 5 cm

"Mean colony diameter (cm) of CoralPop's size bin `b` (1-based), as the engine computes it: sqrt(lo^2/2 + hi^2/2)."
coralpop_bin_diameter(b; bin=CORALPOP_BIN_CM) = sqrt(((b - 1) * bin)^2 / 2 + (b * bin)^2 / 2)

"""
    read_coralpop(dir, reef) -> NamedTuple

Read every `scenario_*.npz` in `dir`. The 75 rows of the engine's time axis are the years 2025-2099.
"""
function read_coralpop(dir::AbstractString, reef; years=2025:2099)
    files = filter(f -> occursin(r"^scenario_\d+\.npz$", f), readdir(dir))
    sort!(files; by=f -> parse(Int, match(r"scenario_(\d+)", f).captures[1]))
    B, M = 20, length(CORALPOP_MORPHS)
    L, N = length(reef.ids), length(files)

    cover = zeros(length(years), M, B, L, N)
    for (n, f) in enumerate(files)
        ## NPZ.jl cannot read the file's text array (benthic_columns), so ask for the array we need
        area = npzread(joinpath(dir, f), ["surface_area_m2"])["surface_area_m2"]    # [T, bin, morph, site], m2
        size(area) == (length(years), B, M, L) || error("$f has unexpected shape $(size(area))")
        cover[:, :, :, :, n] = clamp.(permutedims(area, (1, 3, 2, 4)) ./ reshape(reef.k_area, 1, 1, 1, L), 0, 1)
    end
    d_cm = repeat(reshape([coralpop_bin_diameter(b) for b in 1:B], 1, B), M, 1)       # [morph, bin]
    juvenile = falses(M, B)
    juvenile[:, 1:2] .= true                                                          # bins 1-2: colonies < 10 cm
    return (; label="CoralPop", years=collect(years), groups=CORALPOP_MORPHS, locations=reef.ids, cover,
        k_area=reef.k_area, reef_area=reef.reef_area, diam_cm=repeat(d_cm, 1, 1, N),
        pap=planar_area_for([2, 1, 5], B), is_juvenile=juvenile)                      # Branching, Foliose, Other -> ADRIA rows 2, 1, 5
end

# ===== Each model's own cover number, for checking the adapters =====

"""
    native_cover_cscape(dir, reef; years) -> Array{Float64,3}

C~scape's own total cover (value slot 1 of `out_array`), as relative cover `[T, L, N]`. It is stored in PERCENT of
the site area, per type and heat class, so it is summed and converted to a fraction of the habitable area.
"""
function native_cover_cscape(dir::AbstractString, reef; years=2025:2099)
    files = filter(f -> occursin(r"^Array_scenario_.*\.rds$", f), readdir(dir))
    sort!(files; by=f -> parse(Int, match(r"scenario_(\d+)", f).captures[1]))
    out = zeros(length(years), length(reef.ids), length(files))
    for (n, f) in enumerate(files)
        a = RData.load(joinpath(dir, f))
        pct = dropdims(sum(a[2:end, :, :, :, :, 1]; dims=(3, 4, 5)); dims=(3, 4, 5))       # [T, L], % of site area
        out[:, :, n] = pct ./ 100 .* reshape(reef.reef_area, 1, :) ./ reshape(reef.k_area, 1, :)
    end
    return out
end

"""
    native_cover_coralpop(dir, reef; years) -> Array{Float64,3}

CoralPop's own `cover_pct` (percent of REEF area per morphology), as relative cover `[T, L, N]`.
"""
function native_cover_coralpop(dir::AbstractString, reef; years=2025:2099)
    files = filter(f -> occursin(r"^scenario_\d+\.npz$", f), readdir(dir))
    sort!(files; by=f -> parse(Int, match(r"scenario_(\d+)", f).captures[1]))
    out = zeros(length(years), length(reef.ids), length(files))
    for (n, f) in enumerate(files)
        pct = npzread(joinpath(dir, f), ["cover_pct"])["cover_pct"]                          # [T, morph, L]
        out[:, :, n] = dropdims(sum(pct; dims=2); dims=2) ./ 100 .* reshape(reef.reef_area, 1, :) ./ reshape(reef.k_area, 1, :)
    end
    return out
end
