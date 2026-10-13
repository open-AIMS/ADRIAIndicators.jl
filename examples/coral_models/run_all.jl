#=
run_all.jl -- the whole example, from the two models' saved output files to the figures.

    julia --project=. run_all.jl          (run from this folder; see README.md)

    data/cscape/*.rds    --read_cscape-->                                         ┐
    data/coralpop/*.npz  --read_coralpop--> standard layout --ADRIAIndicators--> indicators --> figures + CSV
    data/domain/*.csv    (reef areas, heat stress)                                ┘

One run of each model: C~scape and CoralPop with their default parameters (no sampling).

The model-specific loading is in src/adapters.jl, the indicator pipeline in src/pipeline.jl.
=#

include(joinpath(@__DIR__, "src", "adapters.jl"))
include(joinpath(@__DIR__, "src", "pipeline.jl"))
import CairoMakie as Mk
using Printf

const DATA_DIR = joinpath(@__DIR__, "data")
const OUT_DIR = mkpath(joinpath(@__DIR__, "output"))

# ===== STEP 1: the reef and its heat stress =====

reef = read_reef(joinpath(DATA_DIR, "domain", "singlereef_domain.csv"))
dhw_table, dhw_header = readdlm(joinpath(DATA_DIR, "domain", "dhw_member17.csv"), ','; header=true)
YEARS = collect(2025:2099)
dhw = Float64.(dhw_table[:, findfirst(==("dhw"), vec(dhw_header))])      # one reef: one DHW value per year
@assert Int.(dhw_table[:, 1]) == YEARS

# ===== STEP 2: read each model's own output files =====

cscape = read_cscape(joinpath(DATA_DIR, "cscape"), reef; years=YEARS)
coralpop = read_coralpop(joinpath(DATA_DIR, "coralpop"), reef; years=YEARS)

# ===== STEP 3: check the standard layout, then compute the indicators =====

models = (cscape, coralpop)
results = map(models) do m
    check_model(m)
    ind = compute_indicators(m)                       # [timesteps, locations, runs]
    dom = domain_metrics(m, ind)                      # [timesteps, runs]
    @info "$(m.label): $(size(m.cover, 5)) run, $(size(m.cover, 4)) location(s), " *
          "$(length(m.years)) years ($(first(m.years))-$(last(m.years)))"
    (; m, ind, dom)
end

# ===== STEP 4: check the adapters against each model's own cover =====

native = (native_cover_cscape(joinpath(DATA_DIR, "cscape"), reef; years=YEARS),
    native_cover_coralpop(joinpath(DATA_DIR, "coralpop"), reef; years=YEARS))
println("\nRelative cover recomputed from the colony data vs the model's own cover (max |difference|):")
for (r, nat) in zip(results, native)
    mine = r.ind[:relative_cover]
    ## CoralPop's first row is the cover it was GIVEN as input, the others come from its size classes: check apart
    @printf("  %-9s all years %.2e   years 2026-2099 %.2e\n", r.m.label, maximum(abs.(mine .- nat)),
        maximum(abs.(mine[2:end, :, :] .- nat[2:end, :, :])))
end

# ===== STEP 5: write the indicators to a CSV =====

open(joinpath(OUT_DIR, "domain_metrics.csv"), "w") do io
    println(io, "model,metric,year,value")
    for r in results, (key, _) in METRICS
        v = vec(r.dom[key])                           # one run: [timesteps, 1] -> [timesteps]
        for (t, yr) in enumerate(r.m.years)
            @printf(io, "%s,%s,%d,%.6g
", r.m.label, key, yr, v[t])
        end
    end
end

println("
Indicators in $(last(YEARS)):")
@printf("  %-36s %12s %12s
", "", cscape.label, coralpop.label)
for (key, title) in METRICS
    v = [only(r.dom[key][end, :]) for r in results]
    @printf("  %-36s %12.4g %12.4g
", title, v[1], v[2])
end

# ===== STEP 6: figures =====

## Okabe-Ito colours, line style and marker per model (as in the six-engine comparison figure).
const STYLE = Dict(
    "C~scape" => (color=Mk.Makie.to_color("#D55E00"), linestyle=:dash, marker=:rect),
    "CoralPop" => (color=Mk.Makie.to_color("#009E73"), linestyle=:dot, marker=:utriangle),
)

"One model's line on `ax`; a marker every 10 years so that the legend's marker can be found on the line."
function model_line!(ax, years, y, label)
    st = STYLE[label]
    Mk.lines!(ax, years, vec(y); color=st.color, linewidth=2, linestyle=st.linestyle, label=label)
    Mk.scatter!(ax, years[1:10:end], vec(y)[1:10:end]; color=st.color, marker=st.marker, markersize=8)
end

"A legend with colour + line style + marker per model."
function model_legend!(pos, labels; kw...)
    elems = [[Mk.LineElement(color=STYLE[l].color, linestyle=STYLE[l].linestyle, linewidth=2),
        Mk.MarkerElement(color=STYLE[l].color, marker=STYLE[l].marker, markersize=9)] for l in labels]
    Mk.Legend(pos, elems, labels; framevisible=false, kw...)
end

function domain_metrics_figure(results, years, dhw)
    ncols = 3
    nrows = cld(length(METRICS) + 1, ncols)
    fig = Mk.Figure(; size=(1400, 380 * nrows))
    ax = Mk.Axis(fig[1, 1]; title="Heat stress forcing (DHW)", titlesize=15, ylabel="DHW (°C-weeks)")
    Mk.barplot!(ax, years, dhw; color=:grey60, gap=0.1)
    Mk.hlines!(ax, [4, 8]; color=:grey30, linestyle=:dash, linewidth=1)
    for (i, (key, title)) in enumerate(METRICS)
        row, col = fld1(i + 1, ncols), mod1(i + 1, ncols)
        ax = Mk.Axis(fig[row, col]; title, titlesize=15, xlabel=row == nrows ? "Year" : "")
        for r in results
            model_line!(ax, r.m.years, r.dom[key], r.m.label)
        end
    end
    model_legend!(fig[nrows + 1, :], [r.m.label for r in results]; orientation=:horizontal)
    Mk.Label(fig[0, :], "Domain-mean ADRIAIndicators metrics (one run per model; SingleReef, DHW member 17)"; fontsize=18)
    return fig
end

fig = domain_metrics_figure(results, YEARS, dhw)
Mk.save(joinpath(OUT_DIR, "domain_metrics.png"), fig; px_per_unit=2)

## Cover of each model's OWN functional groups: why diversity and evenness differ between the models
const GROUP_COLORS = Mk.Makie.to_color.(["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"])

function native_groups_figure(results, years)
    fig = Mk.Figure(; size=(1200, 450))
    for (j, r) in enumerate(results)
        ax = Mk.Axis(fig[1, j]; title="$(r.m.label): cover by its own groups", xlabel="Year",
            ylabel=j == 1 ? "Relative cover" : "")
        for g in eachindex(r.m.groups)
            ## cover of group g: sum over size classes, as a fraction of habitable area; one location, one run
            y = vec(sum(r.m.cover[:, g, :, 1, 1]; dims=2))
            Mk.lines!(ax, years, y; color=GROUP_COLORS[mod1(g, length(GROUP_COLORS))], linewidth=2,
                label=replace(r.m.groups[g], "_" => " "))
        end
        Mk.axislegend(ax; position=:lt, labelsize=10, framevisible=false)
    end
    return fig
end

Mk.save(joinpath(OUT_DIR, "cover_by_native_group.png"), native_groups_figure(results, YEARS); px_per_unit=2)
@info "Wrote domain_metrics.png, cover_by_native_group.png and domain_metrics.csv to $OUT_DIR"
