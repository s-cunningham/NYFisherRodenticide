agesex_forms <- list(bin.exp ~ Age + Sex + (1|RegionalID),
                     bin.exp ~ Age * Sex + (1|RegionalID),
                     bin.exp ~ Sex + (1|RegionalID),
                     bin.exp ~ Age + (1|RegionalID),
                     bin.exp ~ Sex + Age + I(Age^2) + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + (1|RegionalID))
names(agesex_forms) <- as.character(agesex_forms)

ag_formulae <- list(bin.exp ~ totalag_15 + (1|RegionalID),
                    bin.exp ~ totalag_30 + (1|RegionalID),
                    bin.exp ~ totalag_60 + (1|RegionalID),
                    bin.exp ~ crops_15 + (1|RegionalID),
                    bin.exp ~ crops_30 + (1|RegionalID),
                    bin.exp ~ crops_60 + (1|RegionalID),
                    bin.exp ~ pasture_15 + (1|RegionalID),
                    bin.exp ~ pasture_30 + (1|RegionalID),
                    bin.exp ~ pasture_60 + (1|RegionalID))
names(ag_formulae) <- as.character(ag_formulae)

wui_formulae <- list(bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mix_15_100 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mix_30_100 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mix_60_100 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mix_15_250 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mix_30_250 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mix_60_250 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mix_15_500 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mix_30_500 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mix_60_500 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + face_15_100 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + face_30_100 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + face_60_100 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + face_15_250 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + face_30_250 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + face_60_250 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + face_15_500 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + face_30_500 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + face_60_500 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + wui_15_100 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + wui_30_100 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + wui_60_100 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + wui_15_250 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + wui_30_250 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + wui_60_250 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + wui_15_500 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + wui_30_500 + (1|RegionalID),
                     bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + wui_60_500 + (1|RegionalID))
names(wui_formulae) <- as.character(wui_formulae)

beech_formulae <- list(bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + BBA_15 + (1|RegionalID),
                       bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + BBA_30 + (1|RegionalID),
                       bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + BBA_60 + (1|RegionalID))
names(beech_formulae) <- as.character(beech_formulae)


forest_formulae <- list(bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + totalforest_15 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + totalforest_30 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + totalforest_60 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + deciduous_15 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + deciduous_30 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + deciduous_60 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + evergreen_15 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + evergreen_30 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + evergreen_60 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mixed_15 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mixed_30 + (1|RegionalID),
                        bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + mixed_60 + (1|RegionalID))
names(forest_formulae) <- as.character(forest_formulae)

build_formulae <- list(bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + nbuildings_15 + (1|RegionalID),
                       bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + nbuildings_30 + (1|RegionalID),
                       bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + nbuildings_60 + (1|RegionalID),
                       bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + build_cat_15 + (1|RegionalID),
                       bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + build_cat_30 + (1|RegionalID),
                       bin.exp ~ Sex + Age + Sex:Age + I(Age^2) + Sex:I(Age^2) + build_cat_60 + (1|RegionalID))
names(build_formulae) <- as.character(build_formulae)


# Full models for 

full_models <- list(bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + totalag_15 + laggedBMI_15 + ed_15 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + totalag_30 + laggedBMI_30 + ed_30 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + totalag_60 + laggedBMI_60 + ed_60 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + crops_15 + laggedBMI_15 + ed_15 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + crops_30 + laggedBMI_30 + ed_30 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + crops_60 + laggedBMI_60 + ed_60 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + pasture_15 + laggedBMI_15 + ed_15 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + pasture_30 + laggedBMI_30 + ed_30 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + pasture_60 + laggedBMI_60 + ed_60 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + totalag_15 + BMI_15 + ed_15 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + totalag_30 + BMI_30 + ed_30 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + totalag_60 + BMI_60 + ed_60 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + crops_15 + BMI_15 + ed_15 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + crops_30 + BMI_30 + ed_30 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + crops_60 + BMI_60 + ed_60 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + pasture_15 + BMI_15 + ed_15 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + pasture_30 + BMI_30 + ed_30 + (1|RegionalID),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + pasture_60 + BMI_60 + ed_60 + (1|RegionalID))
names(full_models) <- as.character(full_models)