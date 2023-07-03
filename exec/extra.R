# how many genus?
length(unique(plant_df$genus))
# 75
# how many species, besides sp./spp.?
length(unique(plant_df$epithet)) - 2
# 86
# how many species, including sp. but excluding spp.?
length(unique(plant_df$species_fixed))
# 115
# abundance of each genus?
genus_time_table <- table(plant_df$genus, plant_df$yyyy_mm)
genus_time_table <- matrix(genus_time_table, ncol = ncol(genus_time_table), dimnames = dimnames(genus_time_table))
max(genus_time_table)
# most abundant genus had 98 records across all plot-years
