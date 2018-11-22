library(ggplot2)
library(ggthemes)
library(reshape2)
library(survminer)

load("BRCA.RData")


fml<-as.formula("new_death ~ . -death_event - new_death - radiation_therapy - death_event -vital_status")

brcaSurvival<- causalTree(fml,
                    BRCA , treatment = BRCA$radiation_therapy, completeCase = BRCA$death_event,
                    split.Honest=F, cv.option="CT", minsize = 80,
                    split.Rule="survival2",
                    cv.alpha = 1, xval=0, cp=0,
                    propensity = rep(0.5, length(BRCA$radiation_therapy)))


brcaCT<- causalTree(fml,
                    BRCA , treatment = BRCA$radiation_therapy,
                     split.Rule="CT", split.Honest=T, cv.option="CT", minsize = 80,
                     cv.alpha = 1, xval=0, cp=0, split.alpha = 0.5,
                     propensity = rep(0.5, length(BRCA$radiation_therapy)))




