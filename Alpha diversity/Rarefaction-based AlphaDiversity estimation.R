### abundance

bac.melt <- bac.clean.ss %>% psmelt() %>% group_by(PlantSpecies)%>% summarise(SumAbund = sum(Abundance))

bac.melt.total <- bac.clean.ss %>% psmelt() %>% group_by(PlantSpecies)%>% summarise(TotalAbund = sum(Abundance))
head(bac.melt)


write.xlsx(bac.melt, "bacterial abundance in each plantspecies.xlsx")



bac.sample <- bac.clean.ss %>% psmelt() %>% group_by(SampleID)%>% summarise(SumAbund = sum(Abundance))

write.xlsx(bac.sample, "bacterial abundance in each sample.xlsx")




### Diversity - Shannon and Chao1 or ACE
# bac.clean.ss.group <- merge_phyloseq(bac.clean.ss, "Group")
### Rarefaction-based diversity estimation
### For calculating alpha diversity, rarefy step has to be performed. (To remove the effect of different sequencing depth on alpha diversity)

bac.rarefy<- rarefy_even_depth(bac.clean.ss, sample.size = min(sample_sums(bac.clean.ss)),
                  rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#231OTUs were removed because they are no longer present in any sample after random subsampling
tab.bac.alpha<-microbiome::alpha(bac.rarefy, index = "all")
b.meta <-map
b.meta$Richness <- tab.bac.alpha$chao1
b.meta$Shannon <- tab.bac.alpha$diversity_shannon

#To work with species names
b.meta$PS <- gsub(" ", "", b.meta$PlantSpecies)


## Richness of PlantSpecies
## Depth
chao1.AC <- mean(b.meta$Richness[which(b.meta$PlantSpecies == "Araliacordata")])
chao1.AG <- mean(b.meta$Richness[which(b.meta$PlantSpecies == "Angelicagigas")])
chao1.AP <- mean(b.meta$Richness[which(b.meta$PlantSpecies == "Astragaluspropinquus")])
chao1.CA <- mean(b.meta$Richness[which(b.meta$PlantSpecies == "Cynanchumauriculatum")])
chao1.CW <- mean(b.meta$Richness[which(b.meta$PlantSpecies == "Cynanchumwilfordii")])
chao1.GU <- mean(b.meta$Richness[which(b.meta$PlantSpecies == "Glycyrrhizauralensis")])
chao1.PG <- mean(b.meta$Richness[which(b.meta$PlantSpecies == "Panaxginseng")])
chao1.PJ <- mean(b.meta$Richness[which(b.meta$PlantSpecies == "Peucedanumjaponicum")])




##Sample
chao1.AC_1_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AC_1_1")])
chao1.AC_1_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AC_1_2")])
chao1.AC_2_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AC_2_1")])
chao1.AC_2_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AC_2_2")])
chao1.AC_3_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AC_3_1")])
chao1.AC_3_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AC_3_2")])
chao1.AC_4_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AC_4_1")])
chao1.AC_4_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AC_4_2")])
chao1.AG_1_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AG_1_1")])
chao1.AG_1_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AG_1_2")])
chao1.AG_2_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AG_2_1")])
chao1.AG_2_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AG_2_2")])
chao1.AG_3_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AG_3_1")])
chao1.AG_3_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AG_3_2")])
chao1.AG_4_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AG_4_1")])
chao1.AG_4_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AG_4_2")])
chao1.AP_1_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AP_1_1")])
chao1.AP_1_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AP_1_2")])
chao1.AP_2_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AP_2_1")])
chao1.AP_2_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AP_2_2")])
chao1.AP_3_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AP_3_1")])
chao1.AP_3_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AP_3_2")])
chao1.AP_4_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "AP_4_1")])
chao1.AP_4_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "AP_4_2")])
chao1.CA_1_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_1_1")])
chao1.CA_1_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_1_2")])
chao1.CA_2_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_2_1")])
chao1.CA_2_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_2_2")])
chao1.CA_3_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_3_1")])
chao1.CA_3_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_3_2")])
chao1.CA_4_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_4_1")])
chao1.CA_4_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_4_2")])
chao1.CW_1_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_1_1")])
chao1.CW_1_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_1_2")])
chao1.CW_2_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_2_1")])
chao1.CW_2_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_2_2")])
chao1.CW_3_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_3_1")])
chao1.CW_3_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_3_2")])
chao1.CW_4_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_4_1")])
chao1.CW_4_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "CA_4_2")])
chao1.GU_1_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "GU_1_1")])
chao1.GU_1_2<- mean(b.meta$Richness[which(b.meta$SampleID == "GU_1_2")])
chao1.GU_2_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "GU_2_1")])
chao1.GU_2_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "GU_2_2")])
chao1.GU_3_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "GU_3_1")])
chao1.GU_3_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "GU_3_2")])
chao1.GU_4_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "GU_4_1")])
chao1.GU_4_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "GU_4_2")])
chao1.PG_1_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "PG_1_1")])
chao1.PG_1_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "PG_1_2")])
chao1.PG_2_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "PG_2_1")])
chao1.PG_2_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "PG_2_2")])
chao1.PG_3_1<- mean(b.meta$Richness[which(b.meta$SampleID == "PG_3_1")])
chao1.PG_3_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "PG_3_2")])
chao1.PG_4_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "PG_4_1")])
chao1.PG_4_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "PG_4_2")])
chao1.PJ1_1_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "PJ1_1_1")])
chao1.PJ1_1_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "PJ1_1_2")])
chao1.PJ1_2_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "PJ1_2_1")])
chao1.PJ1_2_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "PJ1_2_2")])
chao1.PJ1_3_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "PJ1_3_1")])
chao1.PJ1_3_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "PJ1_3_2")])
chao1.PJ1_4_1 <- mean(b.meta$Richness[which(b.meta$SampleID == "PJ1_4_1")])
chao1.PJ1_4_2 <- mean(b.meta$Richness[which(b.meta$SampleID == "PJ1_4_2")])


## Shannon of depth compartments
## Depth
shannon1.AC <- mean(b.meta$Shannon[which(b.meta$PlantSpecies == "Araliacordata")])
shannon1.AG <- mean(b.meta$Shannon[which(b.meta$PlantSpecies == "Angelicagigas")])
shannon1.AP <- mean(b.meta$Shannon[which(b.meta$PlantSpecies == "Astragaluspropinquus")])
shannon1.CA <- mean(b.meta$Shannon[which(b.meta$PlantSpecies == "Cynanchumauriculatum")])
shannon1.CW <- mean(b.meta$Shannon[which(b.meta$PlantSpecies == "Cynanchumwilfordii")])
shannon1.GU <- mean(b.meta$Shannon[which(b.meta$PlantSpecies == "Glycyrrhizauralensis")])
shannon1.PG <- mean(b.meta$Shannon[which(b.meta$PlantSpecies == "Panaxginseng")])
shannon1.PJ <- mean(b.meta$Shannon[which(b.meta$PlantSpecies == "Peucedanumjaponicum")])


##SampleID
shannon1.AC_1_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AC_1_1")])
shannon1.AC_1_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AC_1_2")])
shannon1.AC_2_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AC_2_1")])
shannon1.AC_2_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AC_2_2")])
shannon1.AC_3_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AC_3_1")])
shannon1.AC_3_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AC_3_2")])
shannon1.AC_4_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AC_4_1")])
shannon1.AC_4_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AC_4_2")])
shannon1.AG_1_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AG_1_1")])
shannon1.AG_1_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AG_1_2")])
shannon1.AG_2_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AG_2_1")])
shannon1.AG_2_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AG_2_2")])
shannon1.AG_3_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AG_3_1")])
shannon1.AG_3_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AG_3_2")])
shannon1.AG_4_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AG_4_1")])
shannon1.AG_4_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AG_4_2")])
shannon1.AP_1_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AP_1_1")])
shannon1.AP_1_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AP_1_2")])
shannon1.AP_2_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AP_2_1")])
shannon1.AP_2_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AP_2_2")])
shannon1.AP_3_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AP_3_1")])
shannon1.AP_3_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AP_3_2")])
shannon1.AP_4_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AP_4_1")])
shannon1.AP_4_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "AP_4_2")])
shannon1.CA_1_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_1_1")])
shannon1.CA_1_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_1_2")])
shannon1.CA_2_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_2_1")])
shannon1.CA_2_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_2_2")])
shannon1.CA_3_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_3_1")])
shannon1.CA_3_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_3_2")])
shannon1.CA_4_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_4_1")])
shannon1.CA_4_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_4_2")])
shannon1.CW_1_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_1_1")])
shannon1.CW_1_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_1_2")])
shannon1.CW_2_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_2_1")])
shannon1.CW_2_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_2_2")])
shannon1.CW_3_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_3_1")])
shannon1.CW_3_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_3_2")])
shannon1.CW_4_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_4_1")])
shannon1.CW_4_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "CA_4_2")])
shannon1.GU_1_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "GU_1_1")])
shannon1.GU_1_2<- mean(b.meta$Shannon[which(b.meta$SampleID == "GU_1_2")])
shannon1.GU_2_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "GU_2_1")])
shannon1.GU_2_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "GU_2_2")])
shannon1.GU_3_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "GU_3_1")])
shannon1.GU_3_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "GU_3_2")])
shannon1.GU_4_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "GU_4_1")])
shannon1.GU_4_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "GU_4_2")])
shannon1.PG_1_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PG_1_1")])
shannon1.PG_1_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PG_1_2")])
shannon1.PG_2_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PG_2_1")])
shannon1.PG_2_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PG_2_2")])
shannon1.PG_3_1<- mean(b.meta$Shannon[which(b.meta$SampleID == "PG_3_1")])
shannon1.PG_3_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PG_3_2")])
shannon1.PG_4_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PG_4_1")])
shannon1.PG_4_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PG_4_2")])
shannon1.PJ1_1_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PJ1_1_1")])
shannon1.PJ1_1_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PJ1_1_2")])
shannon1.PJ1_2_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PJ1_2_1")])
shannon1.PJ1_2_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PJ1_2_2")])
shannon1.PJ1_3_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PJ1_3_1")])
shannon1.PJ1_3_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PJ1_3_2")])
shannon1.PJ1_4_1 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PJ1_4_1")])
shannon1.PJ1_4_2 <- mean(b.meta$Shannon[which(b.meta$SampleID == "PJ1_4_2")])



# wilcoxon test
##Shannon
###PlantSpecies
x <- b.meta$Shannon[which(b.meta$PlantSpecies == "Araliacordata")]
y <- b.meta$Shannon[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 379, p-value = 0.001711


x <- b.meta$Shannon[which(b.meta$PlantSpecies == "Angelicagigas")]
y <- b.meta$Shannon[which(b.meta$PlantSpecies =="Araliacordata" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 341, p-value = 0.01803


x <- b.meta$Shannon[which(b.meta$PlantSpecies == "Astragaluspropinquus")]
y <- b.meta$Shannon[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Araliacordata" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 331, p-value = 0.03062


x <- b.meta$Shannon[which(b.meta$PlantSpecies == "Cynanchumauriculatum")]
y <- b.meta$Shannon[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Araliacordata"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 247, p-value = 0.6479

x <- b.meta$Shannon[which(b.meta$PlantSpecies == "Cynanchumwilfordii")]
y <- b.meta$Shannon[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Araliacordata"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 58, p-value = 0.0007804

x <- b.meta$Shannon[which(b.meta$PlantSpecies == "Glycyrrhizauralensis")]
y <- b.meta$Shannon[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Araliacordata"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 144, p-value = 0.1066

x <- b.meta$Shannon[which(b.meta$PlantSpecies == "Panaxginseng")]
y <- b.meta$Shannon[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Araliacordata"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 86, p-value = 0.005251

x <- b.meta$Shannon[which(b.meta$PlantSpecies == "Peucedanumjaponicum")]
y <- b.meta$Shannon[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Araliacordata" )]
wilcox.test(x, y, conf.int = TRUE) # W = 206, p-value = 0.7224





##Richness
###PlantSpecies
x <- b.meta$Richness[which(b.meta$PlantSpecies == "Araliacordata")]
y <- b.meta$Richness[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 388, p-value = 0.0009031


x <- b.meta$Richness[which(b.meta$PlantSpecies == "Angelicagigas")]
y <- b.meta$Richness[which(b.meta$PlantSpecies =="Araliacordata" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 243, p-value = 0.7072


x <- b.meta$Richness[which(b.meta$PlantSpecies == "Astragaluspropinquus")]
y <- b.meta$Richness[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Araliacordata" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 293.5, p-value = 0.1613


x <- b.meta$Richness[which(b.meta$PlantSpecies == "Cynanchumauriculatum")]
y <- b.meta$Richness[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Araliacordata"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 295, p-value = 0.1524

x <- b.meta$Richness[which(b.meta$PlantSpecies == "Cynanchumwilfordii")]
y <- b.meta$Richness[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Araliacordata"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 14, p-value = 2.11e-05

x <- b.meta$Richness[which(b.meta$PlantSpecies == "Glycyrrhizauralensis")]
y <- b.meta$Richness[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Araliacordata"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 177.5, p-value = 0.3504

x <- b.meta$Richness[which(b.meta$PlantSpecies == "Panaxginseng")]
y <- b.meta$Richness[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Araliacordata"| b.meta$PlantSpecies =="Peucedanumjaponicum" )]
wilcox.test(x, y, conf.int = TRUE) # W = 186, p-value = 0.4465

x <- b.meta$Richness[which(b.meta$PlantSpecies == "Peucedanumjaponicum")]
y <- b.meta$Richness[which(b.meta$PlantSpecies =="Angelicagigas" | b.meta$PlantSpecies =="Astragaluspropinquus" | b.meta$PlantSpecies =="Cynanchumauriculatum"| b.meta$PlantSpecies =="Cynanchumwilfordii"| b.meta$PlantSpecies =="Glycyrrhizauralensis"| b.meta$PlantSpecies =="Panaxginseng"| b.meta$PlantSpecies =="Araliacordata" )]
wilcox.test(x, y, conf.int = TRUE) # W = 195, p-value = 0.5629

p <- ggplot(data = b.meta, aes(x=PlantSpecies, y=Richness)) + geom_boxplot(fill="plum2", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('ANOVA, P-value = 2.43e-07 '))+ 
  ylab("Chao1 Indice\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 1)

p


#Fungi chao
fun.melt <- fun.clean.ss %>% psmelt() %>% group_by(PlantSpecies)%>% summarise(SumAbund = sum(Abundance))

fun.melt.total <- fun.clean.ss %>% psmelt() %>% group_by(PlantSpecies)%>% summarise(TotalAbund = sum(Abundance))
head(fun.melt)


write.xlsx(fun.melt, "fungi abundance in each plantspecies.xlsx")



fun.sample <- fun.clean.ss %>% psmelt() %>% group_by(SampleID)%>% summarise(SumAbund = sum(Abundance))

write.xlsx(fun.sample, "fungi abundance in each sample.xlsx")




### Diversity - Shannon and Chao1 or ACE
# fun.clean.ss.group <- merge_phyloseq(fun.clean.ss, "Group")
### Rarefaction-based diversity estimation
### For calculating alpha diversity, rarefy step has to be performed. (To remove the effect of different sequencing depth on alpha diversity)

fun.rarefy<- rarefy_even_depth(fun.clean.ss, sample.size = min(sample_sums(fun.clean.ss)),
                               rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#35 OTUs were removed because they are no longer present in any sample after random subsampling
tab.fun.alpha<-microbiome::alpha(fun.rarefy, index = "all")
b.meta <-f.map
b.meta$Richness <- tab.fun.alpha$chao1
b.meta$Shannon <- tab.fun.alpha$diversity_shannon

#To work with species names
b.meta$PS <- gsub(" ", "", b.meta$PlantSpecies)


