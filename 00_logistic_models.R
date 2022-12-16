

ag_formulae <- list(bin.exp ~ totalag_15 + (1|WMU) + (1|year),
                    bin.exp ~ totalag_30 + (1|WMU) + (1|year),
                    bin.exp ~ totalag_60 + (1|WMU) + (1|year),
                    bin.exp ~ crops_15 + (1|WMU) + (1|year),
                    bin.exp ~ crops_30 + (1|WMU) + (1|year),
                    bin.exp ~ crops_60 + (1|WMU) + (1|year),
                    bin.exp ~ pasture_15 + (1|WMU) + (1|year),
                    bin.exp ~ pasture_30 + (1|WMU) + (1|year),
                    bin.exp ~ pasture_60 + (1|WMU) + (1|year))
names(ag_formulae) <- as.character(ag_formulae)

wui_formulae <- list(bin.exp ~ mix_15_100 + (1|WMU) + (1|year),
                     bin.exp ~ mix_30_100 + (1|WMU) + (1|year),
                     bin.exp ~ mix_60_100 + (1|WMU) + (1|year),
                     bin.exp ~ mix_15_250 + (1|WMU) + (1|year),
                     bin.exp ~ mix_30_250 + (1|WMU) + (1|year),
                     bin.exp ~ mix_60_250 + (1|WMU) + (1|year),
                     bin.exp ~ mix_15_500 + (1|WMU) + (1|year),
                     bin.exp ~ mix_30_500 + (1|WMU) + (1|year),
                     bin.exp ~ mix_60_500 + (1|WMU) + (1|year),
                     bin.exp ~ face_15_100 + (1|WMU) + (1|year),
                     bin.exp ~ face_30_100 + (1|WMU) + (1|year),
                     bin.exp ~ face_60_100 + (1|WMU) + (1|year),
                     bin.exp ~ face_15_250 + (1|WMU) + (1|year),
                     bin.exp ~ face_30_250 + (1|WMU) + (1|year),
                     bin.exp ~ face_60_250 + (1|WMU) + (1|year),
                     bin.exp ~ face_15_500 + (1|WMU) + (1|year),
                     bin.exp ~ face_30_500 + (1|WMU) + (1|year),
                     bin.exp ~ face_60_500 + (1|WMU) + (1|year),
                     bin.exp ~ wui_15_100 + (1|WMU) + (1|year),
                     bin.exp ~ wui_30_100 + (1|WMU) + (1|year),
                     bin.exp ~ wui_60_100 + (1|WMU) + (1|year),
                     bin.exp ~ wui_15_250 + (1|WMU) + (1|year),
                     bin.exp ~ wui_30_250 + (1|WMU) + (1|year),
                     bin.exp ~ wui_60_250 + (1|WMU) + (1|year),
                     bin.exp ~ wui_15_500 + (1|WMU) + (1|year),
                     bin.exp ~ wui_30_500 + (1|WMU) + (1|year),
                     bin.exp ~ wui_60_500 + (1|WMU) + (1|year))
names(wui_formulae) <- as.character(wui_formulae)

beech_formulae <- list(bin.exp ~ laggedBMI_15 + (1|WMU) + (1|year),
                       bin.exp ~ laggedBMI_30 + (1|WMU) + (1|year),
                       bin.exp ~ laggedBMI_60 + (1|WMU) + (1|year),
                       bin.exp ~ BMI_15 + (1|WMU) + (1|year),
                       bin.exp ~ BMI_30 + (1|WMU) + (1|year),
                       bin.exp ~ BMI_60 + (1|WMU) + (1|year))
names(beech_formulae) <- as.character(beech_formulae)

lsm_formulae <- list(bin.exp ~ ed_15 + (1|WMU) + (1|year),
                     bin.exp ~ ed_30 + (1|WMU) + (1|year),
                     bin.exp ~ ed_60 + (1|WMU) + (1|year),
                     bin.exp ~ dcad_15 + (1|WMU) + (1|year),
                     bin.exp ~ dcad_30 + (1|WMU) + (1|year),
                     bin.exp ~ dcad_60 + (1|WMU) + (1|year),                     
                     bin.exp ~ cohesion_15 + (1|WMU) + (1|year),
                     bin.exp ~ cohesion_30 + (1|WMU) + (1|year),
                     bin.exp ~ cohesion_60 + (1|WMU) + (1|year),
                     bin.exp ~ contig_mn_15 + (1|WMU) + (1|year),
                     bin.exp ~ contig_mn_30 + (1|WMU) + (1|year),
                     bin.exp ~ contig_mn_60 + (1|WMU) + (1|year),
                     bin.exp ~ ai_15 + (1|WMU) + (1|year),
                     bin.exp ~ ai_30 + (1|WMU) + (1|year),
                     bin.exp ~ ai_60 + (1|WMU) + (1|year),
                     bin.exp ~ mesh_15 + (1|WMU) + (1|year),
                     bin.exp ~ mesh_30 + (1|WMU) + (1|year),
                     bin.exp ~ mesh_60 + (1|WMU) + (1|year),
                     bin.exp ~ pd_15 + (1|WMU) + (1|year),
                     bin.exp ~ pd_30 + (1|WMU) + (1|year),
                     bin.exp ~ pd_60 + (1|WMU) + (1|year),
                     bin.exp ~ shape_mn_15 + (1|WMU) + (1|year),
                     bin.exp ~ shape_mn_30 + (1|WMU) + (1|year),
                     bin.exp ~ shape_mn_60 + (1|WMU) + (1|year),
                     bin.exp ~ clumpy_15 + (1|WMU) + (1|year),
                     bin.exp ~ clumpy_30 + (1|WMU) + (1|year),
                     bin.exp ~ clumpy_60 + (1|WMU) + (1|year))
names(lsm_formulae) <- as.character(lsm_formulae)

forest_formulae <- list(bin.exp ~ totalforest_15 + (1|WMU) + (1|year),
                        bin.exp ~ totalforest_30 + (1|WMU) + (1|year),
                        bin.exp ~ totalforest_60 + (1|WMU) + (1|year),
                        bin.exp ~ deciduous_15 + (1|WMU) + (1|year),
                        bin.exp ~ deciduous_30 + (1|WMU) + (1|year),
                        bin.exp ~ deciduous_60 + (1|WMU) + (1|year),
                        bin.exp ~ evergreen_15 + (1|WMU) + (1|year),
                        bin.exp ~ evergreen_30 + (1|WMU) + (1|year),
                        bin.exp ~ evergreen_60 + (1|WMU) + (1|year),
                        bin.exp ~ mixed_15 + (1|WMU) + (1|year),
                        bin.exp ~ mixed_30 + (1|WMU) + (1|year),
                        bin.exp ~ mixed_60 + (1|WMU) + (1|year))
names(forest_formulae) <- as.character(forest_formulae)

build_formulae <- list(bin.exp ~ nbuildings_15 + (1|WMU) + (1|year),
                       bin.exp ~ nbuildings_30 + (1|WMU) + (1|year),
                       bin.exp ~ nbuildings_60 + (1|WMU) + (1|year),
                       bin.exp ~ build_cat_15 + (1|WMU) + (1|year),
                       bin.exp ~ build_cat_30 + (1|WMU) + (1|year),
                       bin.exp ~ build_cat_60 + (1|WMU) + (1|year))
names(build_formulae) <- as.character(build_formulae)


# Full models for 

full_models <- list(bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + totalag_15 + laggedBMI_15 + ed_15 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + totalag_30 + laggedBMI_30 + ed_30 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + totalag_60 + laggedBMI_60 + ed_60 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + crops_15 + laggedBMI_15 + ed_15 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + crops_30 + laggedBMI_30 + ed_30 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + crops_60 + laggedBMI_60 + ed_60 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + pasture_15 + laggedBMI_15 + ed_15 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + pasture_30 + laggedBMI_30 + ed_30 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + pasture_60 + laggedBMI_60 + ed_60 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + totalag_15 + BMI_15 + ed_15 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + totalag_30 + BMI_30 + ed_30 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + totalag_60 + BMI_60 + ed_60 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + crops_15 + BMI_15 + ed_15 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + crops_30 + BMI_30 + ed_30 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + crops_60 + BMI_60 + ed_60 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_15*totalforest_15 + pasture_15 + BMI_15 + ed_15 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_30*totalforest_30 + pasture_30 + BMI_30 + ed_30 + (1|WMU) + (1|year),
                    bin.exp ~ Sex + Age + build_cat_60*totalforest_60 + pasture_60 + BMI_60 + ed_60 + (1|WMU) + (1|year))
names(full_models) <- as.character(full_models)