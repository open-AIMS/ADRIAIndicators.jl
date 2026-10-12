module ADRIAIndicators

include("conversions.jl")

include("metrics.jl")

include("cover_metrics.jl")

include("juvenile_metrics.jl")

include("indices.jl")

# Public API. These names are not exported: call them as `ADRIAIndicators.metric(...)`
# or bring them into scope with `import ADRIAIndicators: metric`.
# Conversions between cover representations
public ltmp_cover_to_relative_cover, ltmp_cover_to_relative_cover!,
    relative_cover_to_ltmp_cover, relative_cover_to_ltmp_cover!

# Shelter volume, diversity and other general metrics
public absolute_shelter_volume, absolute_shelter_volume!, coral_diversity,
    coral_diversity!, coral_evenness, coral_evenness!, maximum_colony_volume,
    relative_shelter_volume, relative_shelter_volume!, scenario_metric

# Cover aggregations
public absolute_coral_area, absolute_coral_area!, coral_area_saved, coral_area_saved!,
    ltmp_cover, ltmp_cover!, ltmp_loc_taxa_cover, ltmp_loc_taxa_cover!,
    ltmp_taxa_cover, ltmp_taxa_cover!, relative_cover, relative_cover!,
    relative_loc_taxa_cover, relative_loc_taxa_cover!, relative_taxa_cover,
    relative_taxa_cover!

# Juvenile aggregations and indicators
public absolute_juveniles, absolute_juveniles!, absolute_loc_taxa_juveniles,
    absolute_loc_taxa_juveniles!, absolute_taxa_juveniles, absolute_taxa_juveniles!,
    juvenile_indicator, juvenile_indicator!, relative_juveniles, relative_juveniles!,
    relative_loc_taxa_juveniles, relative_loc_taxa_juveniles!, relative_taxa_juveniles,
    relative_taxa_juveniles!

# Composite reef indices
public reef_biodiversity_condition_index, reef_biodiversity_condition_index!,
    reef_condition_index, reef_condition_index!, reef_fish_index, reef_fish_index!,
    reef_tourism_index, reef_tourism_index!, reef_tourism_index_no_rubble,
    reef_tourism_index_no_rubble!

end
