## Defining daOTUs
#Phyloseq files
bac.clean.ss
bac.clean.ss@sam_data[["PlantSpecies"]] <- gsub(" ", "", bac.clean.ss@sam_data[["PlantSpecies"]])
fun.clean.ss
fun.clean.ss@sam_data[["PlantSpecies"]] <- gsub(" ", "", fun.clean.ss@sam_data[["PlantSpecies"]])



fun.AC.PG.ss <- subset_samples(fun.clean.ss, PlantSpecies %in% c('Araliacordata','Panaxginseng'))
fun.CW.CA.ss <- subset_samples(fun.clean.ss, PlantSpecies %in% c('Cynanchumwilfordii','Cynanchumauriculatum'))
fun.AG.PJ.ss <- subset_samples(fun.clean.ss, PlantSpecies %in% c('Angelicagigas','Peucedanumjaponicum'))
fun.AS.GU.ss <- subset_samples(fun.clean.ss, PlantSpecies %in% c('Astragaluspropinquus','Glycyrrhizauralensis'))

zero.clean.ac.pg.filt <- phyloseq::filter_taxa(fun.AC.PG.ss, function(x) sum(x) != 0, TRUE)
zero.clean.cw.ca.filt <- phyloseq::filter_taxa(fun.CW.CA.ss, function(x) sum(x) != 0, TRUE)
zero.clean.ag.pj.filt <- phyloseq::filter_taxa(fun.AG.PJ.ss, function(x) sum(x) != 0, TRUE)
zero.clean.as.gu.filt <- phyloseq::filter_taxa(fun.AS.GU.ss, function(x) sum(x) != 0, TRUE)

sum(taxa_sums(zero.clean.ac.pg.filt) == 0)
sum(taxa_sums(zero.clean.cw.ca.filt) == 0)
sum(taxa_sums(zero.clean.ag.pj.filt) == 0)
sum(taxa_sums(zero.clean.as.gu.filt) == 0)
## CODE for CSS normalization using preloaded data
sort(sample_sums(zero.clean.ac.pg.filt))
sort(sample_sums(zero.clean.cw.ca.filt))
sort(sample_sums(zero.clean.ag.pj.filt))
sort(sample_sums(zero.clean.as.gu.filt))



obj <- phyloseq_to_metagenomeSeq(zero.clean.cw.ca.filt)
obj 

## fitZig sample
obj <-  cumNorm(obj, p = cumNormStatFast(obj))
normFactor <-  normFactors(obj)
normFactor <-  log2(normFactor/median(normFactor) + 1)
settings <-  zigControl(maxit = 30, verbose = TRUE)

Type  <-  pData(obj)$PlantSpecies
Type <- gsub(" ", "", Type)
mod  <-  model.matrix(~Type)
mod
colnames(mod)  <-  unique(Type)
colnames(mod)

res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)


zigFit = res@fit
finalMod= res@fit[["design"]]

contrast.matrix.1 =makeContrasts(Panaxginseng - Araliacordata, levels = finalMod)
contrast.matrix.1 =makeContrasts(Cynanchumwilfordii - Cynanchumauriculatum, levels = finalMod)
contrast.matrix.1 =makeContrasts(Angelicagigas - Peucedanumjaponicum, levels = finalMod)
contrast.matrix.1 =makeContrasts(Astragaluspropinquus - Glycyrrhizauralensis, levels = finalMod)


fit2_diagnosis = contrasts.fit(zigFit, contrasts=contrast.matrix.1) 

fit3.1 = eBayes(fit2_diagnosis)

topTable(fit3.1, coef="Cynanchumwilfordii - Cynanchumauriculatum")

res.1 <- topTable(fit3.1,coef=1,adjust="fdr",n=Inf, p.value = 0.05, lfc =1)
head(res.1)
dim(res.1)

## Taxonomy
tax_fun <- tax_table(fun.CW.CA.ss) 

res.tax_fun <- merge(res.1,tax_fun, by = "row.names")

write.xlsx(res.tax_fun, 'fun_3_OTU_diagnosis.xlsx')


## MA plot
log2AverageAbundance <- psmelt(fun.CW.CA.ss) %>% group_by(OTU) %>% summarise(log2AverageAbundance=log2(mean(Abundance)))
log2AverageAbundance
Ta <- psmelt(fun.CW.CA.ss) %>% group_by(OTU) %>% select(Phylum,Class,Order,Family, Genus)
Ta <- unique(Ta)
Ta
resSig = res.1[!is.na(res.1$adj.P.Val), ]
resSig = data.frame(resSig)
head(resSig)
resSig <- tibble::rownames_to_column(resSig, 'OTU')


## Construct volcano plot
h = -log10(max(resSig[resSig$adj.P.Val<0.01,]$adj.P.Val))
# FDR Q<0.05 & logFC>|2|, log2FC<|2|
mydata <- resSig %>% mutate(volcano_y = -log10(adj.P.Val))
mydata <- mydata %>% mutate(threshold = ifelse(logFC >= 2 & volcano_y >= h, 'Angelicagigas', ifelse(logFC <= -2 & volcano_y >= h,"Peucedanumjaponicum", ifelse(volcano_y < h , "Bottom", 'Middle'))))

theme_set(theme_bw())
ggplot(mydata, aes(x=logFC, y=volcano_y)) +
  xlab('\n Fold Change (Log2)')+
  ylab("-log10 adjusted P\n") +
  geom_point(aes(color = threshold),  alpha=0.6) +
  scale_color_manual(labels = c('Angelicagigas'="Angelica gigas enriched\n(FDR Q<0.01 & log2FC>2)     ",'Peucedanumjaponicum'= "Peucedanum japonicum enriched\n(FDR Q<0.01 & log2FC>2)     ",'Middle'='log2FC<|2|     ','Bottom'='P-value > 0.01     '), values = c("Peucedanumjaponicum"= "deeppink4", 'Angelicagigas'='darkseagreen',"Middle"="lightblue","Bottom"= "grey"))+
  theme(legend.text=element_text(size=13)) + 
  theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
  geom_hline(yintercept=h, color="maroon4")+
  geom_vline(xintercept=2, color="maroon4")+
  geom_vline(xintercept=0, color="black")+
  geom_vline(xintercept=-2, color="maroon4")+
  theme(legend.position="top") +
  theme(legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=8),reverse = TRUE))+
  guides(size=FALSE) +
  scale_x_continuous(breaks=seq(-100,100,5))+
  scale_y_continuous(breaks=seq(0,10,1))+
  theme(axis.title.x = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 18,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())






