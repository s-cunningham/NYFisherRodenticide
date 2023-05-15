
## Model lists for selecting scale
ag_formulae <- list(n.compounds.T ~ Sex + Age + I(Age^2) + totalag_15 + (1|RegionalID),
                    n.compounds.T ~ Sex + Age + I(Age^2) + totalag_30 + (1|RegionalID),
                    n.compounds.T ~ Sex + Age + I(Age^2) + totalag_60 + (1|RegionalID),
                    n.compounds.T ~ Sex + Age + I(Age^2) + crops_15 + (1|RegionalID),
                    n.compounds.T ~ Sex + Age + I(Age^2) + crops_30 + (1|RegionalID),
                    n.compounds.T ~ Sex + Age + I(Age^2) + crops_60 + (1|RegionalID),
                    n.compounds.T ~ Sex + Age + I(Age^2) + pasture_15 + (1|RegionalID),
                    n.compounds.T ~ Sex + Age + I(Age^2) + pasture_30 + (1|RegionalID),
                    n.compounds.T ~ Sex + Age + I(Age^2) + pasture_60 + (1|RegionalID))
names(ag_formulae) <- as.character(ag_formulae)

wui_formulae <- list(n.compounds.T ~ Sex + Age + I(Age^2) + wui_15_100 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + wui_30_100 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + wui_60_100 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + wui_15_250 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + wui_30_250 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + wui_60_250 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + wui_15_500 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + wui_30_500 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + wui_60_500 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + mix_15_100 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + mix_30_100 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + mix_60_100 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + mix_15_250 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + mix_30_250 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + mix_60_250 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + mix_15_500 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + mix_30_500 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + mix_60_500 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + face_15_100 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + face_30_100 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + face_60_100 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + face_15_250 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + face_30_250 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + face_60_250 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + face_15_500 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + face_30_500 + (1|RegionalID),
                     n.compounds.T ~ Sex + Age + I(Age^2) + face_60_500 + (1|RegionalID))
names(wui_formulae) <- as.character(wui_formulae)

beech_formulae <- list(n.compounds.T ~ Sex + Age + I(Age^2) + BBA_15 + (1|RegionalID),
                       n.compounds.T ~ Sex + Age + I(Age^2) + BBA_30 + (1|RegionalID),
                       n.compounds.T ~ Sex + Age + I(Age^2) + BBA_60 + (1|RegionalID))
names(beech_formulae) <- as.character(beech_formulae)

forest_formulae <- list(n.compounds.T ~ Sex + Age + I(Age^2) + totalforest_15 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + totalforest_30 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + totalforest_60 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + deciduous_15 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + deciduous_30 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + deciduous_60 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + evergreen_15 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + evergreen_30 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + evergreen_60 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + mixed_15 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + mixed_30 + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + mixed_60 + (1|RegionalID))
names(forest_formulae) <- as.character(forest_formulae)

build_formulae <- list(n.compounds.T ~ Sex + Age + I(Age^2) + nbuildings_15 + (1|RegionalID),
                       n.compounds.T ~ Sex + Age + I(Age^2) + nbuildings_30 + (1|RegionalID),
                       n.compounds.T ~ Sex + Age + I(Age^2) + nbuildings_60 + (1|RegionalID),
                       n.compounds.T ~ Sex + Age + I(Age^2) + build_cat_15 + (1|RegionalID),
                       n.compounds.T ~ Sex + Age + I(Age^2) + build_cat_30 + (1|RegionalID),
                       n.compounds.T ~ Sex + Age + I(Age^2) + build_cat_60 + (1|RegionalID))
names(build_formulae) <- as.character(build_formulae)

agesex_formulae <- list(n.compounds.T ~ Age + Sex + (1|RegionalID),
                        n.compounds.T ~ Age * Sex + (1|RegionalID),
                        n.compounds.T ~ Sex + (1|RegionalID),
                        n.compounds.T ~ Age + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + (1|RegionalID))
names(agesex_formulae) <- as.character(agesex_formulae)

global_formulae <- list(n.compounds.T ~ Sex + Age + I(Age^2) + wui_30_500 + pasture_15 + BBA_15 * lag_beechnuts + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + wui_30_500 + pasture_15 + BBA_15 * mast + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + wui_30_500 + pasture_15 + BBA_15 * beechnuts + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * lag_beechnuts + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * mast + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + mix_15_100 + pasture_15 + BBA_15 * beechnuts + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + build_cat_15*deciduous_15 + pasture_15 + lag_beechnuts + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + build_cat_15*deciduous_15 + pasture_15 + mast + (1|RegionalID),
                        n.compounds.T ~ Sex + Age + I(Age^2) + build_cat_15*deciduous_15 + pasture_15 + beechnuts + (1|RegionalID)
                        )
names(global_formulae) <- as.character(global_formulae)



