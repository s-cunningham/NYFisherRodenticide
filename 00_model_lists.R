
## Model lists for selecting scale


ag_formulae <- list(n.compounds.T ~ Sex*Age + totalag_15 + (1|RegionalID),
                    n.compounds.T ~ Sex*Age + totalag_30 + (1|RegionalID),
                    n.compounds.T ~ Sex*Age + totalag_60 + (1|RegionalID),
                    n.compounds.T ~ Sex*Age + crops_15 + (1|RegionalID),
                    n.compounds.T ~ Sex*Age + crops_30 + (1|RegionalID),
                    n.compounds.T ~ Sex*Age + crops_60 + (1|RegionalID),
                    n.compounds.T ~ Sex*Age + pasture_15 + (1|RegionalID),
                    n.compounds.T ~ Sex*Age + pasture_30 + (1|RegionalID),
                    n.compounds.T ~ Sex*Age + pasture_60 + (1|RegionalID))
names(ag_formulae) <- as.character(ag_formulae)

wui_formulae <- list(n.compounds.T ~ Sex*Age + mix_15_100 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + mix_30_100 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + mix_60_100 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + mix_15_250 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + mix_30_250 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + mix_60_250 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + mix_15_500 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + mix_30_500 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + mix_60_500 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + face_15_100 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + face_30_100 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + face_60_100 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + face_15_250 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + face_30_250 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + face_60_250 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + face_15_500 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + face_30_500 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + face_60_500 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + wui_15_100 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + wui_30_100 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + wui_60_100 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + wui_15_250 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + wui_30_250 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + wui_60_250 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + wui_15_500 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + wui_30_500 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + wui_60_500 + (1|RegionalID))
names(wui_formulae) <- as.character(wui_formulae)

beech_formulae <- list(n.compounds.T ~ Sex*Age + laggedBMI_15 + (1|RegionalID),
                       n.compounds.T ~ Sex*Age + laggedBMI_30 + (1|RegionalID),
                       n.compounds.T ~ Sex*Age + laggedBMI_60 + (1|RegionalID),
                       n.compounds.T ~ Sex*Age + BMI_15 + (1|RegionalID),
                       n.compounds.T ~ Sex*Age + BMI_30 + (1|RegionalID),
                       n.compounds.T ~ Sex*Age + BMI_60 + (1|RegionalID))
names(beech_formulae) <- as.character(beech_formulae)

lsm_formulae <- list(n.compounds.T ~ Sex*Age + ed_15 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + ed_30 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + ed_60 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + dcad_15 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + dcad_30 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + dcad_60 + (1|RegionalID),                     
                     n.compounds.T ~ Sex*Age + cohesion_15 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + cohesion_30 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + cohesion_60 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + cpland_15 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + cpland_30 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + cpland_60 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + contig_mn_15 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + contig_mn_30 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + contig_mn_60 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + ai_15 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + ai_30 + (1|RegionalID),
                     n.compounds.T ~ Sex*Age + ai_60 + (1|RegionalID))
vnames(lsm_formulae) <- as.character(lsm_formulae)


