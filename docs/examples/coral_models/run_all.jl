#=
Calculate ADRIAIndicators metrics from the saved output of C~scape and CoralPOP, and plot them.

    julia --project=. run_all.jl

The code is explained in docs/src/examples/cscape_coralpop.md
=#

using CSV, DataFrames
using RData, NPZ
using CairoMakie

import ADRIAIndicators as AI

cd(@__DIR__)  # the data and output folders are relative to this file

## 1. The reef and its heat stress

domain = CSV.read("data/domain/singlereef_domain.csv", DataFrame)
reef_area_m2 = domain.reef_area_m2       # [locations], m², includes non-habitable area
habitable_area_m2 = domain.k_area_m2     # [locations], m²

years = 2025:2099
dhw = CSV.read("data/domain/dhw_member17.csv", DataFrame).dhw

## 2. C~scape: counts to cover

cscape_types = ["acro_table", "acro_corym", "corym_non_acro", "small_massive", "large_massive", "brooder"]
size_classes = RData.load("data/cscape/size_classes.rds")

# Colony diameter (cm) and juvenile classes [groups ⋅ sizes]
cscape_diam_cm = permutedims(reduce(hcat, [Float64.(size_classes[t]["meshpoints_diam"]) for t in cscape_types]))
cscape_juvenile = permutedims(reduce(hcat, [Bool.(size_classes[t]["juvenile"]) for t in cscape_types]))

cscape_out = RData.load("data/cscape/Array_scenario_101000100001_rep_1_draw_NA.rds")  # [year ⋅ site ⋅ intervened ⋅ type ⋅ enhancement ⋅ value]

# Drop the initial state (first year), keep the 100 size classes (value slots 7-106)
# and add up the `intervened` sub-areas and enhancement classes. [timesteps ⋅ locations ⋅ groups ⋅ sizes]
counts = dropdims(sum(cscape_out[2:end, :, :, :, :, 7:106]; dims=(3, 5)); dims=(3, 5))
counts = permutedims(counts, (1, 3, 4, 2))  # [timesteps ⋅ groups ⋅ sizes ⋅ locations]

colony_area_m2 = π .* (cscape_diam_cm ./ 200) .^ 2
cscape_cover = counts .* reshape(colony_area_m2, 1, size(colony_area_m2)...) ./
               reshape(habitable_area_m2, 1, 1, 1, :)

## 3. CoralPOP: area to cover

# NPZ.jl does not support NumPy's unicode string arrays (the `benthic_columns` array in the
# file), so name the arrays to read.
coralpop_output = npzread("data/coralpop/scenario_1.npz", ["surface_area_m2", "cover_pct"])

# The surface area of all colonies [time ⋅ size bin ⋅ morphology ⋅ site], m²
coralpop_area_m2 = coralpop_output["surface_area_m2"]
coralpop_cover = permutedims(coralpop_area_m2, (1, 3, 2, 4)) ./ reshape(habitable_area_m2, 1, 1, 1, :)

# Bins are 5 cm wide. The first two (colonies under 10 cm) are juveniles
bin_width_cm = 5.0
bin_diam_cm = [sqrt(((b - 1) * bin_width_cm)^2 / 2 + (b * bin_width_cm)^2 / 2) for b in 1:20]
coralpop_diam_cm = repeat(reshape(bin_diam_cm, 1, :), 3, 1)  # [groups ⋅ sizes]
coralpop_juvenile = falses(3, 20)
coralpop_juvenile[:, 1:2] .= true

## 4. Indicators

# Parameters (intercept b, coefficient a) of the log-log planar area model log(S) = b + a log(x)
# for each of ADRIA's functional groups, as in ADRIA.jl's `planar_area_params()`. They come from
# Fig 2B of Aston et al. (2022), https://doi.org/10.3389/fmars.2022.854395 (corymbose non-Acropora
# uses the branching Pocillopora values). Rows: tabular Acropora, corymbose Acropora,
# corymbose non-Acropora, small massives, large massives
adria_planar_area = [
    -8.95 2.80
    -9.13 2.94
    -8.90 2.94
    -8.87 2.30
    -8.87 2.30
]

# Each model's groups take one of the rows above, the same for every size class [groups ⋅ sizes ⋅ 2]
planar_area_params(rows, n_sizes) = repeat(reshape(adria_planar_area[rows, :], length(rows), 1, 2), 1, n_sizes, 1)

cscape_planar_area = planar_area_params([1, 2, 3, 4, 5, 3], 100)  # brooders use corymbose non-Acropora
coralpop_planar_area = planar_area_params([2, 1, 5], 20)          # Branching, Foliose, Other

# Reference colony defining the maximum shelter volume: (diameter cm, intercept, coefficient)
reference = (95.0, -8.95, 2.80)

max_juvenile_density = 15.0  # individuals/m²

function model_indicators(rel_cover, diam_cm, planar_area, is_juvenile)
    # Cover [timesteps ⋅ locations]
    relative_cover = AI.relative_cover(rel_cover)
    ltmp_cover = AI.ltmp_cover(rel_cover, habitable_area_m2, reef_area_m2)

    # Shelter volume, summed over groups and sizes [timesteps ⋅ locations]
    rsv = AI.relative_shelter_volume(rel_cover, diam_cm, planar_area, habitable_area_m2, reference)
    relative_shelter_volume = clamp.(dropdims(sum(rsv; dims=(2, 3)); dims=(2, 3)), 0.0, 1.0)

    asv = AI.absolute_shelter_volume(rel_cover, diam_cm, planar_area, habitable_area_m2)
    absolute_shelter_volume = dropdims(sum(asv; dims=(2, 3)); dims=(2, 3))

    # Diversity, from the cover of each group [timesteps ⋅ groups ⋅ locations]
    relative_loc_taxa_cover = AI.relative_loc_taxa_cover(rel_cover)
    coral_diversity = AI.coral_diversity(relative_loc_taxa_cover)
    coral_evenness = AI.coral_evenness(relative_loc_taxa_cover)

    # Juveniles. Mean colony diameters are expected in m here, whereas the shelter volume
    # functions take cm.
    relative_juveniles = AI.relative_juveniles(rel_cover, is_juvenile)
    juvenile_indicator = AI.juvenile_indicator(
        rel_cover, is_juvenile, habitable_area_m2, diam_cm ./ 100, max_juvenile_density
    )

    # Neither model has rubble, so the indices that do not need it are used
    reef_condition_index = AI.reef_condition_index(ltmp_cover, relative_shelter_volume, juvenile_indicator)
    reef_tourism_index = AI.reef_tourism_index_no_rubble(
        ltmp_cover, coral_evenness, relative_shelter_volume, relative_juveniles
    )
    reef_fish_index = AI.reef_fish_index(relative_cover)
    reef_biodiversity_condition_index = AI.reef_biodiversity_condition_index(
        relative_cover, coral_diversity, relative_shelter_volume
    )

    return (;
        relative_cover, relative_juveniles, juvenile_indicator, relative_shelter_volume,
        absolute_shelter_volume, coral_evenness, coral_diversity, reef_condition_index,
        reef_tourism_index, reef_fish_index, reef_biodiversity_condition_index
    )
end

cscape = model_indicators(cscape_cover, cscape_diam_cm, cscape_planar_area, cscape_juvenile)
coralpop = model_indicators(coralpop_cover, coralpop_diam_cm, coralpop_planar_area, coralpop_juvenile)

metric_titles = [
    :relative_cover => "Relative cover",
    :relative_juveniles => "Relative juveniles",
    :juvenile_indicator => "Juvenile indicator",
    :relative_shelter_volume => "Relative shelter volume",
    :absolute_shelter_volume => "Absolute shelter volume (m³)",
    :coral_evenness => "Coral evenness",
    :coral_diversity => "Coral diversity",
    :reef_condition_index => "Reef Condition Index",
    :reef_tourism_index => "Reef Tourism Index",
    :reef_fish_index => "Reef Fish Index",
    :reef_biodiversity_condition_index => "Reef Biodiversity Condition Index"
]

println("Values in $(last(years)):")
println(rpad("", 36), rpad("C~scape", 10), "CoralPOP")
for (key, title) in metric_titles
    println(rpad(title, 36), rpad(round(getproperty(cscape, key)[end, 1]; sigdigits=3), 10),
        round(getproperty(coralpop, key)[end, 1]; sigdigits=3))
end

## 5. Checking the conversions against each model's own cover

# C~scape: value slot 1, in % of the site area. CoralPOP: `cover_pct`, in % of the reef area.
# Both are converted to relative cover (of the habitable area).
cscape_native = dropdims(sum(cscape_out[2:end, :, :, :, :, 1]; dims=(3, 4, 5)); dims=(3, 4, 5)) ./ 100 .*
                reef_area_m2' ./ habitable_area_m2'
coralpop_pct = coralpop_output["cover_pct"]  # [time ⋅ morphology ⋅ site]
coralpop_native = dropdims(sum(coralpop_pct; dims=2); dims=2) ./ 100 .* reef_area_m2' ./ habitable_area_m2'

println("\nRelative cover from the colony data vs the model's own cover (max |difference|):")
for (name, model, native) in (("C~scape", cscape, cscape_native), ("CoralPOP", coralpop, coralpop_native))
    difference = abs.(model.relative_cover .- native)
    println("  $name: all years $(maximum(difference)), 2026-2099 $(maximum(difference[2:end, :]))")
end

## 6. Plots

mkpath("output")

styles = Dict(
    "C~scape" => (color="#D55E00", linestyle=:dash, marker=:rect),
    "CoralPOP" => (color="#009E73", linestyle=:dot, marker=:utriangle)
)
models = ["C~scape" => cscape, "CoralPOP" => coralpop]

fig_metrics = Figure(size=(1400, 1140))
ax_dhw = Axis(fig_metrics[1, 1]; title="Heat stress forcing (DHW)", ylabel="DHW (°C-weeks)")
barplot!(ax_dhw, years, dhw; color=:grey60, gap=0.1)
hlines!(ax_dhw, [4, 8]; color=:grey30, linestyle=:dash, linewidth=1)
for (i, (key, title)) in enumerate(metric_titles)
    row, col = fld1(i + 1, 3), mod1(i + 1, 3)
    ax = Axis(fig_metrics[row, col]; title, xlabel=row == 4 ? "Year" : "")
    for (name, model) in models
        s = styles[name]
        y = getproperty(model, key)[:, 1]
        lines!(ax, years, y; color=s.color, linestyle=s.linestyle, linewidth=2)
        scatter!(ax, years[1:10:end], y[1:10:end]; color=s.color, marker=s.marker, markersize=8)
    end
end
Legend(
    fig_metrics[5, :],
    [[LineElement(; color=s.color, linestyle=s.linestyle, linewidth=2),
        MarkerElement(; color=s.color, marker=s.marker, markersize=9)]
     for s in getindex.(Ref(styles), first.(models))],
    first.(models); orientation=:horizontal, framevisible=false
)
save("output/domain_metrics.png", fig_metrics; px_per_unit=2)

# Cover of each model's own groups, which helps to read the diversity and evenness differences
group_colors = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
fig_groups = Figure(size=(1200, 450))
for (j, (name, rel_cover, groups)) in enumerate((
    ("C~scape", cscape_cover, replace.(cscape_types, "_" => " ")),
    ("CoralPOP", coralpop_cover, ["Branching", "Foliose", "Other"])
))
    ax = Axis(fig_groups[1, j]; title="$name: cover by its own groups", xlabel="Year",
        ylabel=j == 1 ? "Relative cover" : "")
    for (g, group) in enumerate(groups)
        lines!(ax, years, vec(sum(rel_cover[:, g, :, 1]; dims=2)); color=group_colors[g], linewidth=2, label=group)
    end
    axislegend(ax; position=:lt, labelsize=10, framevisible=false)
end
save("output/cover_by_native_group.png", fig_groups; px_per_unit=2)
