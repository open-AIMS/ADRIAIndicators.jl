#=
coralpop_figures.jl -- makes the figures of the CoralPOP page of the cookbook
(docs/src/cookbook/coralpop.md).

The cookbook code is not run by the documentation build (NPZ.jl and CairoMakie.jl are not
docs dependencies), so the figures are generated here and committed to
docs/src/assets/cookbook/. This script runs every ```julia block of the page, in order,
exactly as written, and then saves the two figures that the page's last section builds.
The page reads the CoralPOP output file from its own folder, so the blocks are run from
docs/examples/coral_models/data/coralpop.

From the repository root:

    julia --project=docs/figures -e 'using Pkg; Pkg.instantiate()'
    julia --project=docs/figures docs/figures/coralpop_figures.jl
=#

const PAGE = joinpath(@__DIR__, "..", "src", "cookbook", "coralpop.md")
const DATA = joinpath(@__DIR__, "..", "examples", "coral_models", "data", "coralpop")
const OUT = mkpath(joinpath(@__DIR__, "..", "src", "assets", "cookbook"))

text = replace(read(PAGE, String), "\r\n" => "\n")
blocks = [m.captures[1] for m in eachmatch(r"```julia\n(.*?)```"s, text)]
@info "Running $(length(blocks)) code blocks of $(basename(PAGE))"

mod = Module(:CoralPOPCookbook)
cd(DATA) do
    for (i, code) in enumerate(blocks)
        Core.eval(mod, Meta.parseall(code; filename="$(basename(PAGE)) block $i"))
    end
end

for name in (:fig_groups, :fig_metrics)
    isdefined(mod, name) || error("The page's code did not define `$name`.")
end
mod.save(joinpath(OUT, "coralpop_groups.png"), mod.fig_groups; px_per_unit=2)
mod.save(joinpath(OUT, "coralpop_metrics.png"), mod.fig_metrics; px_per_unit=2)
@info "Wrote the figures to $OUT"
