using PackageCompiler

target_dir = "$(@__DIR__)/ReefMetricsCompiled"
target_dir = replace(target_dir, "\\"=>"/")

println("Creating library in $target_dir")
PackageCompiler.create_library("$(@__DIR__)/..", target_dir;
                                lib_name="ReefMetrics",
                                precompile_execution_file=["$(@__DIR__)/precompile.jl"],
                                incremental=false,
                                filter_stdlibs=true,
                                force=true, # Overwrite target_dir.
                            )
