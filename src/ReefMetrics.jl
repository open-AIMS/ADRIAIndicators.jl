module ReefMetrics

include("aggregation.jl")
include("conversion.jl")
include("metrics.jl")

export coral_diversity, coral_evenness, absolute_shelter_volume, 
    relative_shelter_volume

end
