
## Model lists for selecting scale
ag_formulae <- list(n.compounds.T ~ totalag_15 + (1|RegionalID),
                    n.compounds.T ~ totalag_30 + (1|RegionalID),
                    n.compounds.T ~ totalag_60 + (1|RegionalID),
                    n.compounds.T ~ crops_15 + (1|RegionalID),
                    n.compounds.T ~ crops_30 + (1|RegionalID),
                    n.compounds.T ~ crops_60 + (1|RegionalID),
                    n.compounds.T ~ pasture_15 + (1|RegionalID),
                    n.compounds.T ~ pasture_30 + (1|RegionalID),
                    n.compounds.T ~ pasture_60 + (1|RegionalID))
names(ag_formulae) <- as.character(ag_formulae)

wui_formulae <- list(n.compounds.T ~ mix_15_100 + (1|RegionalID),
                     n.compounds.T ~ mix_30_100 + (1|RegionalID),
                     n.compounds.T ~ mix_60_100 + (1|RegionalID),
                     n.compounds.T ~ mix_15_250 + (1|RegionalID),
                     n.compounds.T ~ mix_30_250 + (1|RegionalID),
                     n.compounds.T ~ mix_60_250 + (1|RegionalID),
                     n.compounds.T ~ mix_15_500 + (1|RegionalID),
                     n.compounds.T ~ mix_30_500 + (1|RegionalID),
                     n.compounds.T ~ mix_60_500 + (1|RegionalID),
                     n.compounds.T ~ face_15_100 + (1|RegionalID),
                     n.compounds.T ~ face_30_100 + (1|RegionalID),
                     n.compounds.T ~ face_60_100 + (1|RegionalID),
                     n.compounds.T ~ face_15_250 + (1|RegionalID),
                     n.compounds.T ~ face_30_250 + (1|RegionalID),
                     n.compounds.T ~ face_60_250 + (1|RegionalID),
                     n.compounds.T ~ face_15_500 + (1|RegionalID),
                     n.compounds.T ~ face_30_500 + (1|RegionalID),
                     n.compounds.T ~ face_60_500 + (1|RegionalID),
                     n.compounds.T ~ wui_15_100 + (1|RegionalID),
                     n.compounds.T ~ wui_30_100 + (1|RegionalID),
                     n.compounds.T ~ wui_60_100 + (1|RegionalID),
                     n.compounds.T ~ wui_15_250 + (1|RegionalID),
                     n.compounds.T ~ wui_30_250 + (1|RegionalID),
                     n.compounds.T ~ wui_60_250 + (1|RegionalID),
                     n.compounds.T ~ wui_15_500 + (1|RegionalID),
                     n.compounds.T ~ wui_30_500 + (1|RegionalID),
                     n.compounds.T ~ wui_60_500 + (1|RegionalID))
names(wui_formulae) <- as.character(wui_formulae)

beech_formulae <- list(n.compounds.T ~ BBA_15 + (1|RegionalID),
                       n.compounds.T ~ BBA_30 + (1|RegionalID),
                       n.compounds.T ~ BBA_60 + (1|RegionalID)
                       # n.compounds.T ~ mast + (1|RegionalID),
                       # n.compounds.T ~ mast + (1|RegionalID),
                       # n.compounds.T ~ mast + (1|RegionalID),
                       # n.compounds.T ~ beechnuts + (1|RegionalID),
                       # n.compounds.T ~ beechnuts + (1|RegionalID),
                       # n.compounds.T ~ beechnuts + (1|RegionalID),
                       # n.compounds.T ~ lag_beechnuts + (1|RegionalID),
                       # n.compounds.T ~ lag_beechnuts + (1|RegionalID),
                       # n.compounds.T ~ lag_beechnuts + (1|RegionalID),
                       # n.compounds.T ~ BBA_15 + mast + (1|RegionalID),
                       # n.compounds.T ~ BBA_30 + mast + (1|RegionalID),
                       # n.compounds.T ~ BBA_60 + mast + (1|RegionalID),
                       # n.compounds.T ~ BBA_15 * mast + (1|RegionalID),
                       # n.compounds.T ~ BBA_30 * mast + (1|RegionalID),
                       # n.compounds.T ~ BBA_60 * mast + (1|RegionalID),
                       # n.compounds.T ~ BBA_15 + year + (1|RegionalID),
                       # n.compounds.T ~ BBA_30 + year + (1|RegionalID),
                       # n.compounds.T ~ BBA_60 + year + (1|RegionalID),
                       # n.compounds.T ~ BBA_15 * year + (1|RegionalID),
                       # n.compounds.T ~ BBA_30 * year + (1|RegionalID),
                       # n.compounds.T ~ BBA_60 * year + (1|RegionalID),
                       # n.compounds.T ~ BBA_15 * beechnuts + (1|RegionalID),
                       # n.compounds.T ~ BBA_30 * beechnuts + (1|RegionalID),
                       # n.compounds.T ~ BBA_60 * beechnuts + (1|RegionalID),
                       # n.compounds.T ~ BBA_15 * lag_beechnuts + (1|RegionalID),
                       # n.compounds.T ~ BBA_30 * lag_beechnuts + (1|RegionalID),
                       # n.compounds.T ~ BBA_60 * lag_beechnuts + (1|RegionalID)
                       )
names(beech_formulae) <- as.character(beech_formulae)

lsm_formulae <- list(n.compounds.T ~ ed_15 + (1|RegionalID),
                     n.compounds.T ~ ed_30 + (1|RegionalID),
                     n.compounds.T ~ ed_60 + (1|RegionalID),
                     n.compounds.T ~ dcad_15 + (1|RegionalID),
                     n.compounds.T ~ dcad_30 + (1|RegionalID),
                     n.compounds.T ~ dcad_60 + (1|RegionalID),
                     n.compounds.T ~ cohesion_15 + (1|RegionalID),
                     n.compounds.T ~ cohesion_30 + (1|RegionalID),
                     n.compounds.T ~ cohesion_60 + (1|RegionalID),
                     n.compounds.T ~ contig_mn_15 + (1|RegionalID),
                     n.compounds.T ~ contig_mn_30 + (1|RegionalID),
                     n.compounds.T ~ contig_mn_60 + (1|RegionalID),
                     n.compounds.T ~ ai_15 + (1|RegionalID),
                     n.compounds.T ~ ai_30 + (1|RegionalID),
                     n.compounds.T ~ ai_60 + (1|RegionalID),
                     n.compounds.T ~ mesh_15 + (1|RegionalID),
                     n.compounds.T ~ mesh_30 + (1|RegionalID),
                     n.compounds.T ~ mesh_60 + (1|RegionalID),
                     n.compounds.T ~ pd_15 + (1|RegionalID),
                     n.compounds.T ~ pd_30 + (1|RegionalID),
                     n.compounds.T ~ pd_60 + (1|RegionalID),
                     n.compounds.T ~ shape_mn_15 + (1|RegionalID),
                     n.compounds.T ~ shape_mn_30 + (1|RegionalID),
                     n.compounds.T ~ shape_mn_60 + (1|RegionalID),
                     n.compounds.T ~ clumpy_15 + (1|RegionalID),
                     n.compounds.T ~ clumpy_30 + (1|RegionalID),
                     n.compounds.T ~ clumpy_60 + (1|RegionalID))
names(lsm_formulae) <- as.character(lsm_formulae)

forest_formulae <- list(n.compounds.T ~ totalforest_15 + (1|RegionalID),
                        n.compounds.T ~ totalforest_30 + (1|RegionalID),
                        n.compounds.T ~ totalforest_60 + (1|RegionalID),
                        n.compounds.T ~ deciduous_15 + (1|RegionalID),
                        n.compounds.T ~ deciduous_30 + (1|RegionalID),
                        n.compounds.T ~ deciduous_60 + (1|RegionalID),
                        n.compounds.T ~ evergreen_15 + (1|RegionalID),
                        n.compounds.T ~ evergreen_30 + (1|RegionalID),
                        n.compounds.T ~ evergreen_60 + (1|RegionalID),
                        n.compounds.T ~ mixed_15 + (1|RegionalID),
                        n.compounds.T ~ mixed_30 + (1|RegionalID),
                        n.compounds.T ~ mixed_60 + (1|RegionalID))
names(forest_formulae) <- as.character(forest_formulae)

build_formulae <- list(n.compounds.T ~ nbuildings_15 + (1|RegionalID),
                       n.compounds.T ~ nbuildings_30 + (1|RegionalID),
                       n.compounds.T ~ nbuildings_60 + (1|RegionalID),
                       n.compounds.T ~ build_cat_15 + (1|RegionalID),
                       n.compounds.T ~ build_cat_30 + (1|RegionalID),
                       n.compounds.T ~ build_cat_60 + (1|RegionalID))
names(build_formulae) <- as.character(build_formulae)


agesex_formulae <- list(n.compounds.T ~ Age + Sex + (1|RegionalID),
                        n.compounds.T ~ Age * Sex + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + (1|RegionalID))
names(agesex_formulae) <- as.character(agesex_formulae)

global_formulae <- list(n.compounds.T ~ Sex + Age + I(Age^2) + wui_60_100 + pasture_60 + BBA_60 * year + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + wui_60_100 + pasture_60 + BBA_60 * lag_beechnuts + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + nbuildings_60*ai_60 + pasture_60 + BBA_60 * year + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + nbuildings_60*mixed_30 + pasture_60 + BBA_60 * year + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + nbuildings_60*ai_60 + pasture_60 + BBA_60 + year + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + nbuildings_60*mixed_30 + pasture_60 + BBA_60 + year + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + wui_60_100 + pasture_60 + BBA_60 * year + ai_60 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + wui_60_100 + pasture_60 + BBA_60 * lag_beechnuts + ai_60 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + wui_60_100 + pasture_60 + BBA_60 * lag_beechnuts + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + wui_60_100 + pasture_60 + BBA_60 * year + (1|RegionalID)
)
names(global_formulae) <- as.character(global_formulae)



# ag_formulae <- list(n.compounds.T ~ totalag_15 + (1|WMU) + (1|year),
#                     n.compounds.T ~ totalag_30 + (1|WMU) + (1|year),
#                     n.compounds.T ~ totalag_60 + (1|WMU) + (1|year),
#                     n.compounds.T ~ crops_15 + (1|WMU) + (1|year),
#                     n.compounds.T ~ crops_30 + (1|WMU) + (1|year),
#                     n.compounds.T ~ crops_60 + (1|WMU) + (1|year),
#                     n.compounds.T ~ pasture_15 + (1|WMU) + (1|year),
#                     n.compounds.T ~ pasture_30 + (1|WMU) + (1|year),
#                     n.compounds.T ~ pasture_60 + (1|WMU) + (1|year))
# names(ag_formulae) <- as.character(ag_formulae)

# wui_formulae <- list(n.compounds.T ~ mix_15_100 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mix_30_100 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mix_60_100 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mix_15_250 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mix_30_250 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mix_60_250 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mix_15_500 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mix_30_500 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mix_60_500 + (1|WMU) + (1|year),
#                      n.compounds.T ~ face_15_100 + (1|WMU) + (1|year),
#                      n.compounds.T ~ face_30_100 + (1|WMU) + (1|year),
#                      n.compounds.T ~ face_60_100 + (1|WMU) + (1|year),
#                      n.compounds.T ~ face_15_250 + (1|WMU) + (1|year),
#                      n.compounds.T ~ face_30_250 + (1|WMU) + (1|year),
#                      n.compounds.T ~ face_60_250 + (1|WMU) + (1|year),
#                      n.compounds.T ~ face_15_500 + (1|WMU) + (1|year),
#                      n.compounds.T ~ face_30_500 + (1|WMU) + (1|year),
#                      n.compounds.T ~ face_60_500 + (1|WMU) + (1|year),
#                      n.compounds.T ~ wui_15_100 + (1|WMU) + (1|year),
#                      n.compounds.T ~ wui_30_100 + (1|WMU) + (1|year),
#                      n.compounds.T ~ wui_60_100 + (1|WMU) + (1|year),
#                      n.compounds.T ~ wui_15_250 + (1|WMU) + (1|year),
#                      n.compounds.T ~ wui_30_250 + (1|WMU) + (1|year),
#                      n.compounds.T ~ wui_60_250 + (1|WMU) + (1|year),
#                      n.compounds.T ~ wui_15_500 + (1|WMU) + (1|year),
#                      n.compounds.T ~ wui_30_500 + (1|WMU) + (1|year),
#                      n.compounds.T ~ wui_60_500 + (1|WMU) + (1|year))
# names(wui_formulae) <- as.character(wui_formulae)
# 
# beech_formulae <- list(n.compounds.T ~ laggedBMI_15 + (1|WMU),
#                        n.compounds.T ~ laggedBMI_30 + (1|WMU),
#                        n.compounds.T ~ laggedBMI_60 + (1|WMU),
#                        n.compounds.T ~ BMI_15 + (1|WMU),
#                        n.compounds.T ~ BMI_30 + (1|WMU),
#                        n.compounds.T ~ BMI_60 + (1|WMU),
#                        n.compounds.T ~ BBA_15 + (1|WMU),
#                        n.compounds.T ~ BBA_30 + (1|WMU),
#                        n.compounds.T ~ BBA_60 + (1|WMU),
#                        n.compounds.T ~ BBA_15 * year + (1|WMU),
#                        n.compounds.T ~ BBA_30 * year + (1|WMU),
#                        n.compounds.T ~ BBA_60 * year + (1|WMU),
#                        n.compounds.T ~ BBA_15 * beechnuts + (1|WMU),
#                        n.compounds.T ~ BBA_30 * beechnuts + (1|WMU),
#                        n.compounds.T ~ BBA_60 * beechnuts + (1|WMU),
#                        n.compounds.T ~ BBA_15 * lag_beechnuts + (1|WMU),
#                        n.compounds.T ~ BBA_30 * lag_beechnuts + (1|WMU),
#                        n.compounds.T ~ BBA_60 * lag_beechnuts + (1|WMU)
#                        )
# names(beech_formulae) <- as.character(beech_formulae)
# 
# lsm_formulae <- list(n.compounds.T ~ ed_15 + (1|WMU) + (1|year),
#                      n.compounds.T ~ ed_30 + (1|WMU) + (1|year),
#                      n.compounds.T ~ ed_60 + (1|WMU) + (1|year),
#                      n.compounds.T ~ dcad_15 + (1|WMU) + (1|year),
#                      n.compounds.T ~ dcad_30 + (1|WMU) + (1|year),
#                      n.compounds.T ~ dcad_60 + (1|WMU) + (1|year),                     
#                      n.compounds.T ~ cohesion_15 + (1|WMU) + (1|year),
#                      n.compounds.T ~ cohesion_30 + (1|WMU) + (1|year),
#                      n.compounds.T ~ cohesion_60 + (1|WMU) + (1|year),
#                      n.compounds.T ~ contig_mn_15 + (1|WMU) + (1|year),
#                      n.compounds.T ~ contig_mn_30 + (1|WMU) + (1|year),
#                      n.compounds.T ~ contig_mn_60 + (1|WMU) + (1|year),
#                      n.compounds.T ~ ai_15 + (1|WMU) + (1|year),
#                      n.compounds.T ~ ai_30 + (1|WMU) + (1|year),
#                      n.compounds.T ~ ai_60 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mesh_15 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mesh_30 + (1|WMU) + (1|year),
#                      n.compounds.T ~ mesh_60 + (1|WMU) + (1|year),
#                      n.compounds.T ~ pd_15 + (1|WMU) + (1|year),
#                      n.compounds.T ~ pd_30 + (1|WMU) + (1|year),
#                      n.compounds.T ~ pd_60 + (1|WMU) + (1|year),
#                      n.compounds.T ~ shape_mn_15 + (1|WMU) + (1|year),
#                      n.compounds.T ~ shape_mn_30 + (1|WMU) + (1|year),
#                      n.compounds.T ~ shape_mn_60 + (1|WMU) + (1|year),
#                      n.compounds.T ~ clumpy_15 + (1|WMU) + (1|year),
#                      n.compounds.T ~ clumpy_30 + (1|WMU) + (1|year),
#                      n.compounds.T ~ clumpy_60 + (1|WMU) + (1|year))
# names(lsm_formulae) <- as.character(lsm_formulae)
# 
# forest_formulae <- list(n.compounds.T ~ totalforest_15 + (1|WMU) + (1|year),
#                         n.compounds.T ~ totalforest_30 + (1|WMU) + (1|year),
#                         n.compounds.T ~ totalforest_60 + (1|WMU) + (1|year),
#                         n.compounds.T ~ deciduous_15 + (1|WMU) + (1|year),
#                         n.compounds.T ~ deciduous_30 + (1|WMU) + (1|year),
#                         n.compounds.T ~ deciduous_60 + (1|WMU) + (1|year),
#                         n.compounds.T ~ evergreen_15 + (1|WMU) + (1|year),
#                         n.compounds.T ~ evergreen_30 + (1|WMU) + (1|year),
#                         n.compounds.T ~ evergreen_60 + (1|WMU) + (1|year),
#                         n.compounds.T ~ mixed_15 + (1|WMU) + (1|year),
#                         n.compounds.T ~ mixed_30 + (1|WMU) + (1|year),
#                         n.compounds.T ~ mixed_60 + (1|WMU) + (1|year))
# names(forest_formulae) <- as.character(forest_formulae)
# 
# build_formulae <- list(n.compounds.T ~ nbuildings_15 + (1|WMU) + (1|year),
#                        n.compounds.T ~ nbuildings_30 + (1|WMU) + (1|year),
#                        n.compounds.T ~ nbuildings_60 + (1|WMU) + (1|year),
#                        n.compounds.T ~ build_cat_15 + (1|WMU) + (1|year),
#                        n.compounds.T ~ build_cat_30 + (1|WMU) + (1|year),
#                        n.compounds.T ~ build_cat_60 + (1|WMU) + (1|year))
# names(build_formulae) <- as.character(build_formulae)



