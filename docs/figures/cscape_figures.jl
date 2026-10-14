#=
cscape_figures.jl -- makes the figures of the C~scape page of the cookbook
(docs/src/cookbook/cscape.md).

The cookbook code is not run by the documentation build (RData.jl and CairoMakie.jl are not
docs dependencies), so the figures are generated here and committed to
docs/src/assets/cookbook/. This script runs every ```julia block of the page, in order,
exactly as written, and then saves the two figures that the page's last section builds.
The page reads the C~scape output files from its own folder, so the blocks are run from
docs/examples/coral_models/data/cscape.

From the repository root:

    julia --project=docs/figures -e 'using Pkg; Pkg.instantiate()'
    julia --project=docs/figures docs/figures/cscape_figures.jl
=#

const PAGE = joinpath(@__DIR__, "..", "src", "cookbook", "cscape.md")
const DATA = joinpath(@__DIR__, "..", "examples", "coral_models", "data", "cscape")
const OUT = mkpath(joinpath(@__DIR__, "..", "src", "assets", "cookbook"))

text = replace(read(PAGE, String), "\r\n" => "\n")
blocks = [m.captures[1] for m in eachmatch(r"```julia\n(.*?)```"s, text)]
@info "Running $(length(blocks)) code blocks of $(basename(PAGE))"

mod = Module(:CScapeCookbook)
cd(DATA) do
    for (i, code) in enumerate(blocks)
        Core.eval(mod, Meta.parseall(code; filename="$(basename(PAGE)) block $i"))
    end
end

for name in (:fig_groups, :fig_metrics)
    isdefined(mod, name) || error("The page's code did not define `$name`.")
end
mod.save(joinpath(OUT, "cscape_groups.png"), mod.fig_groups; px_per_unit=2)
mod.save(joinpath(OUT, "cscape_metrics.png"), mod.fig_metrics; px_per_unit=2)
@info "Wrote the figures to $OUT"
