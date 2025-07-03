#ifndef REEFMETRICS_H
#define REEFMETRICS_H

void coral_diversity(
    int32_t n_tsteps,
    int32_t n_groups,
    int32_t n_locs,
    const double* relative_taxa_cover,
    double* output_taxa_cover
);

#endif
