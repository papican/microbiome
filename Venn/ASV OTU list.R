#venn

##bacteria
### Divide the phyloseq to two different files for each species
bac.CW.ss <- subset_samples(bac.clean.ss, PlantSpecies == 'Cynanchumwilfordii')
bac.CA.ss <- subset_samples(bac.clean.ss, PlantSpecies == 'Cynanchumauriculatum')
bac.PG.ss <- subset_samples(bac.clean.ss, PlantSpecies == 'Panaxginseng')                                 
bac.AC.ss <- subset_samples(bac.clean.ss, PlantSpecies == 'Araliacordata')
bac.AG.ss <- subset_samples(bac.clean.ss, PlantSpecies == 'Angelicagigas')
bac.PJ.ss <- subset_samples(bac.clean.ss, PlantSpecies == 'Peucedanumjaponicum')
bac.AP.ss <- subset_samples(bac.clean.ss, PlantSpecies == 'Astragaluspropinquus')
bac.GU.ss <- subset_samples(bac.clean.ss, PlantSpecies == 'Glycyrrhizauralensis')
bac.AR.ss <- subset_samples(bac.AR.AP.ss, new == 'Araliaceae')
bac.APi.ss <- subset_samples(bac.AR.AP.ss, new == 'Apiaceae')

bac.CW.filt <- phyloseq::filter_taxa(bac.CW.ss, function(x) sum(x) != 0, TRUE)
bac.CA.filt <- phyloseq::filter_taxa(bac.CA.ss, function(x) sum(x) != 0, TRUE)
bac.PG.filt <- phyloseq::filter_taxa(bac.PG.ss, function(x) sum(x) != 0, TRUE)
bac.AC.filt <- phyloseq::filter_taxa(bac.AC.ss, function(x) sum(x) != 0, TRUE)
bac.AG.filt <- phyloseq::filter_taxa(bac.AG.ss, function(x) sum(x) != 0, TRUE)
bac.PJ.filt <- phyloseq::filter_taxa(bac.PJ.ss, function(x) sum(x) != 0, TRUE)
bac.AP.filt <- phyloseq::filter_taxa(bac.AP.ss, function(x) sum(x) != 0, TRUE)
bac.GU.filt <- phyloseq::filter_taxa(bac.GU.ss, function(x) sum(x) != 0, TRUE)
bac.AR.filt <- phyloseq::filter_taxa(bac.AR.ss, function(x) sum(x) != 0, TRUE)
bac.APi.filt <- phyloseq::filter_taxa(bac.APi.ss, function(x) sum(x) != 0, TRUE)

### Extract ASV lists from each phyloseq object
CW.asv.bac<-taxa_names(bac.CW.filt)
CA.asv.bac<-taxa_names(bac.CA.filt)
PG.asv.bac<-taxa_names(bac.PG.filt)
AC.asv.bac<-taxa_names(bac.AC.filt)
AG.asv.bac<-taxa_names(bac.AG.filt)
PJ.asv.bac<-taxa_names(bac.PJ.filt)
AP.asv.bac<-taxa_names(bac.AP.filt)
GU.asv.bac<-taxa_names(bac.GU.filt)
AR.asv.bac<-taxa_names(bac.AR.filt)
APi.asv.bac<-taxa_names(bac.APi.filt)



### Check co-detected ASVs
codetected <-intersect(AR.asv.bac, APi.asv.bac)

##codetected asv
codetected <- as.data.frame(codetected)
colnames(codetected)[1] <- "OTU" 

library(dplyr)
codetected <- bac.list %>%
  filter(OTU %in% codetected$OTU)
write.csv(codetected, "co-detected_AR-AP_bac.csv")


##CA asv
AR.asv.bac <- as.data.frame(AR.asv.bac)
colnames(AR.asv.bac)[1] <- "OTU" 

library(dplyr)
AR.asv.bac <- bac.list %>%
  filter(OTU %in% AR.asv.bac$OTU)

## CW asv
APi.asv.bac <- as.data.frame(APi.asv.bac)
colnames(APi.asv.bac)[1] <- "OTU" 

library(dplyr)
APi.asv.bac <- bac.list %>%
  filter(OTU %in% APi.asv.bac$OTU)

### Save the lists as excel-executable files
write.csv(APi.asv.bac,"ASV list of APiaceae bac.csv")
write.csv(AR.asv.bac,"ASV list of ARaliacea bac.csv")



##fungi
#venn
### Divide the phyloseq to two different files for each species
fun.CW.ss <- subset_samples(fun.clean.ss, PlantSpecies == 'Cynanchumwilfordii')
fun.CA.ss <- subset_samples(fun.clean.ss, PlantSpecies == 'Cynanchumauriculatum')
fun.PG.ss <- subset_samples(fun.clean.ss, PlantSpecies == 'Panaxginseng')                                 
fun.AC.ss <- subset_samples(fun.clean.ss, PlantSpecies == 'Araliacordata')
fun.AG.ss <- subset_samples(fun.clean.ss, PlantSpecies == 'Angelicagigas')
fun.PJ.ss <- subset_samples(fun.clean.ss, PlantSpecies == 'Peucedanumjaponicum')
fun.AP.ss <- subset_samples(fun.clean.ss, PlantSpecies == 'Astragaluspropinquus')
fun.GU.ss <- subset_samples(fun.clean.ss, PlantSpecies == 'Glycyrrhizauralensis')
fun.AR.ss <- subset_samples(fun.AR.AP.ss, new == 'Araliaceae')
fun.APi.ss <- subset_samples(fun.AR.AP.ss, new == 'Apiaceae')

fun.CW.filt <- phyloseq::filter_taxa(fun.CW.ss, function(x) sum(x) != 0, TRUE)
fun.CA.filt <- phyloseq::filter_taxa(fun.CA.ss, function(x) sum(x) != 0, TRUE)
fun.PG.filt <- phyloseq::filter_taxa(fun.PG.ss, function(x) sum(x) != 0, TRUE)
fun.AC.filt <- phyloseq::filter_taxa(fun.AC.ss, function(x) sum(x) != 0, TRUE)
fun.AG.filt <- phyloseq::filter_taxa(fun.AG.ss, function(x) sum(x) != 0, TRUE)
fun.PJ.filt <- phyloseq::filter_taxa(fun.PJ.ss, function(x) sum(x) != 0, TRUE)
fun.AP.filt <- phyloseq::filter_taxa(fun.AP.ss, function(x) sum(x) != 0, TRUE)
fun.GU.filt <- phyloseq::filter_taxa(fun.GU.ss, function(x) sum(x) != 0, TRUE)
fun.AR.filt <- phyloseq::filter_taxa(fun.AR.ss, function(x) sum(x) != 0, TRUE)
fun.APi.filt <- phyloseq::filter_taxa(fun.APi.ss, function(x) sum(x) != 0, TRUE)

### Extract ASV lists from each phyloseq object
CW.asv.fun<-taxa_names(fun.CW.filt)
CA.asv.fun<-taxa_names(fun.CA.filt)
PG.asv.fun<-taxa_names(fun.PG.filt)
AC.asv.fun<-taxa_names(fun.AC.filt)
AG.asv.fun<-taxa_names(fun.AG.filt)
PJ.asv.fun<-taxa_names(fun.PJ.filt)
AP.asv.fun<-taxa_names(fun.AP.filt)
GU.asv.fun<-taxa_names(fun.GU.filt)
AR.asv.fun<-taxa_names(fun.AR.filt)
APi.asv.fun<-taxa_names(fun.APi.filt)



### Check co-detected ASVs
codetected <-intersect(AR.asv.fun, APi.asv.fun)

##codetected asv
codetected <- as.data.frame(codetected)
colnames(codetected)[1] <- "OTU" 

library(dplyr)
codetected <- fun.list %>%
  filter(OTU %in% codetected$OTU)
write.csv(codetected, "co-detected_AR-AP_fun.csv")


##CA asv
AR.asv.fun <- as.data.frame(AR.asv.fun)
colnames(AR.asv.fun)[1] <- "OTU" 

library(dplyr)
AR.asv.fun <- fun.list %>%
  filter(OTU %in% AR.asv.fun$OTU)

## CW asv
APi.asv.fun <- as.data.frame(APi.asv.fun)
colnames(APi.asv.fun)[1] <- "OTU" 

library(dplyr)
APi.asv.fun <- fun.list %>%
  filter(OTU %in% APi.asv.fun$OTU)

### Save the lists as excel-executable files
write.csv(APi.asv.fun,"ASV list of Apiaceae fun.csv")
write.csv(AR.asv.fun,"ASV list of ARaliaceae fun.csv")
