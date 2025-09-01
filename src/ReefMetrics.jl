module ReefMetrics

include("conversions.jl")

export rhc_to_rrc, rrc_to_rhc

include("metrics.jl")

export coral_diversity, coral_evenness, absolute_shelter_volume,
    relative_shelter_volume

include("cover_metrics.jl")

export relative_cover, relative_loc_cover, relative_taxa_cover,
    relative_loc_taxa_cover

include("juvenile_metrics.jl")

export relative_juveniles, relative_loc_juveniles, relative_taxa_juveniles, 
    relative_loc_taxa_juveniles

export absolute_juveniles, absolute_loc_juveniles, absolute_taxa_juveniles, 
    absolute_loc_taxa_juveniles

export juvenile_indicator

include("indices.jl")

export reef_fish_index, reef_tourism_index, reef_biodiversity_condition_index,
    reef_condition_index

end
