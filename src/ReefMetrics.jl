module ReefMetrics

include("aggregation.jl")
include("conversions.jl")
include("metrics.jl")

export coral_diversity, coral_evenness, absolute_shelter_volume, 
    relative_shelter_volume

export rhc_to_rrc, rrc_to_rhc

end
