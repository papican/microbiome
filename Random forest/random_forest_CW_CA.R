
## Ridge plot


## testing


## get imp.tax
get_imp_tax <- function(phy.clean.ss.5, classifier){
  imp <- data.frame(importance(classifier))
  imp <- tibble::rownames_to_column(imp, 'OTU')
  imp.desc <- imp %>% arrange(desc(MeanDecreaseGini))
  tax <- tax_table(phy.clean.ss.5)
  tax <- as.data.frame(tax)
  tax
  tax <- tibble::rownames_to_column(tax,'OTU')
  
  df.clean.ss.5 <- psmelt(phy.clean.ss.5)
  head(df.clean.ss.5)
  abun <- df.clean.ss.5 %>% group_by(OTU) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
  imp.tax <- left_join(imp.desc, tax, by=c('OTU')) %>% left_join(abun, by=c('OTU'))
  return(imp.tax)
}

#Subset of Conventional and no fertilization 
imp_tax_bac <- get_imp_tax(bac.clean.log, classifier_diagnosis)





#head(mydata)
imp_tax_bac.edit <- subset(imp_tax_bac, OTU != "result1")
head(imp_tax_bac.edit)
imp.arrange.b <- imp_tax_bac.edit %>% arrange(desc(MeanDecreaseGini))
imp.arrange.b$rank <- paste0('RF',1:dim(imp_tax_bac.edit)[1])
another.b <- imp.arrange.b[c('OTU','rank')]






##Construct the plot object

imp.total.arrange <- imp_tax_bac %>% arrange(desc(total))
imp.total.arrange$number <- 1:dim(imp_tax_bac)[1]
imp.total.arrange
num <- imp.total.arrange[c('OTU','number')]
num 

##Construct the plot object
#bac
Ab <- psmelt(bac.CW.CA.ss) %>% group_by(OTU) %>% summarise(Abun=log(sum(Abundance)))
Ta <- psmelt(bac.CW.CA.ss) %>% group_by(OTU) %>% select(OTU, Phylum,Class,Order,Family, Genus,Species)
Ta <- unique(Ta)








## indexing the enrichment classification

zero.clean.filt <- bac.clean.log


sum(taxa_sums(zero.clean.filt) == 0) #2346

sample_names(zero.clean.filt)

obj <- phyloseq_to_metagenomeSeq(zero.clean.filt)

obj <- filterData(obj, depth = 100)
obj   ## removed 0 samples below 100 reads
##
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
contrast.matrix = makeContrasts(Cynanchumwilfordii - Cynanchumauriculatum, levels = finalMod)

fit2=  contrasts.fit(zigFit, contrasts=contrast.matrix) #Cynanchumwilfordii-Cynanchumauriculatum

fit3 = eBayes(fit2)
fit3
topTable(fit3, coef="Cynanchumwilfordii - Cynanchumauriculatum")

res <- topTable(fit3,coef=1,adjust="fdr",n=Inf)

head(res)
dim(res)


imp.tax <- imp_tax_bac


resSig = res[!is.na(res$adj.P.Val), ]
resSig = data.frame(resSig)
head(resSig)
dim(resSig)

resSig <- tibble::rownames_to_column(resSig, 'OTU')
resSig <- dplyr::left_join(resSig, Ab,by= c('OTU'))
resSig <- dplyr::left_join(resSig,Ta,by=c('OTU')) %>% dplyr::left_join(num,by='OTU')
imp.3 <- imp.tax[c('OTU','MeanDecreaseGini','total')]
resSig <- merge(resSig,imp.3,by='OTU')
resSig <- dplyr::left_join(resSig, another.b, by=('OTU'='OTU'))
head(resSig)
unique(resSig)
resSig.arrange <- resSig %>% arrange(desc(total))

resSig$id <- ifelse(is.na(resSig$Genus),paste0(resSig$number,'_f_',resSig$Family),paste0(resSig$number,'_',resSig$Genus))


write.xlsx(resSig, "Bac_CW-CA_daOTUs_with RF.xlsx")




# resSig$padj<0.05
# max(resSig[resSig$padj<0.05,]$padj)


h = -log10(max(resSig[resSig$adj.P.Val<0.05,]$adj.P.Val))
resSig[resSig$adj.P.Val<0.05,]$adj.P.Val
max(resSig[resSig$adj.P.Val<0.05,]$adj.P.Val)
resSig[resSig$adj.P.Val>0.05,]$adj.P.Val
h #1.344532



mydata <- resSig %>% mutate(volcano_y = -log10(adj.P.Val))
mydata <- mydata %>% mutate(threshold = ifelse(logFC > 0 & volcano_y >= h, 'Cynanchumwilfordii', ifelse(logFC < 0 & volcano_y >= h,"Cynanchumauriculatum", "Non-differential")))
# A FDR Q<0.05 & logFC>|2|, C log2FC>|2|, D NS, B FDR Q<0.05
# mydata <- mydata %>% filter(log2FoldChange <20)
size <- ifelse((res$adj.P.Val < 0.05 & abs(res$logFC) > 0), 4, 2)
sub.mydata <- mydata %>% arrange(desc(MeanDecreaseGini)) %>% head(70) #b:head(70), f:head(80)
unique(sub.mydata$threshold)
# }

## get df.ridge

  df.otu <- bac.clean.log %>% psmelt()
  head(df.otu)
  # we need to group by samples
  df.otu.rel <- df.otu %>%  
    group_by(SampleID) %>%                         # Filter out at absolute read of 20       
    mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance
  
  df.otu.rel
  imp.top20 <- imp.tax %>% arrange(desc(MeanDecreaseGini)) %>% head(20) #b: 70, f: 30
  
  imp20 <- imp.top20$OTU
  imp20
  
  imp.top20
  
  df.selected.rel <- df.otu.rel %>% filter(OTU %in% imp20)
  df.ridge <- df.selected.rel %>% select(OTU, Sample,RelAbundance)
  str(df.ridge)
  b.thresh <- sub.mydata[c('OTU','threshold')]
  b.thresh
  df.ridge_2 <- left_join(df.ridge, b.thresh, by=c('OTU'))
  
  unique(df.ridge$OTU)
  unique(b.thresh$OTU)
  
  df.ridge_2
  str(df.ridge_2)
  dim(df.ridge_2)
  





bac.list.OTU_id <- bac.list[c("OTU", "OTU_id")]

sub.mydata <- merge(sub.mydata, bac.list.OTU_id, by  = "OTU")
df.ridge_2 <- merge(df.ridge_2, bac.list.OTU_id, by  = "OTU")

b.manipulate <- sub.mydata[c('rank','threshold','MeanDecreaseGini','OTU_id')]
ord.mani <- b.manipulate %>% group_by(OTU_id) %>%arrange(MeanDecreaseGini)

b.id <- ord.mani$OTU_id
df.ridge_2$OTU_id <- factor(df.ridge_2$OTU_id, levels=b.id)



library(ggplot2)
library(ggridges)

plot_ridge <- function(df_ridge){
  ggplot(df.ridge_2, aes(x=Sample, y=OTU_id, height=RelAbundance, group = OTU_id, fill=threshold)) + 
     #geom_density_ridges_gradient(stat = "identity",scale = 3, rel_min_height = 0.01)+
    geom_density_ridges2(stat = "identity", scale=5, color='white',size=0.001)+
    xlab('\n PlantSpecies')+
    ylab("Relative abundance (‰) \n") +
    theme(legend.text=element_text(size=12)) + 
    # theme(plot.title = element_text(size = 20,hjust = 0.5)) + 
    theme(legend.position="top", legend.spacing.x = unit(0.4, 'cm')) +
    theme(legend.title=element_blank()) +
    guides(colour = guide_legend(override.aes = list(size=8)))+
    guides(size=FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
    theme(axis.title.y = element_text(size = 12,hjust = 0.5, face='bold'))+
    theme(axis.title.x = element_text(size = 12,hjust = 0.5, face='bold')) + 
    theme(axis.text.y = element_text(size=12,face='bold',color='black'))+
    scale_fill_manual(labels = c('Cynanchumwilfordii'="Cynanchum wilfordii",'Cynanchumauriculatum'= "Cynanchum auriculatum",'Non-differential'='Non-differential'), values = c("Cynanchumwilfordii"= "#6699CC", 'Cynanchumauriculatum'='#CC9900',"Non-differential"= "light grey"))+
    theme(panel.grid.major = element_blank())+
    theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())+
    
    geom_vline(xintercept=8.5, color="darkslategray", linetype='dashed',size=1)
  
}

df.ridge_2$Sample <- factor(df.ridge_2$Sample)

plot_ridge(df.ridge_2)


df.ridge_2$Sample <- factor(df.ridge_2$Sample)




daOTU_bac <- read.xlsx("Bac_daOTUs_and_RF.xlsx",1)

daOTU_bac.id <- merge(daOTU_bac, bac.list.OTU_id, by = "OTU")
write.xlsx(daOTU_bac.id, "Bac_daOTUs_and_RF_with id.xlsx")



daOTU_fun <- read.xlsx("Fun_daOTUs_with RF.xlsx",1)

daOTU_fun.id <- merge(daOTU_fun, fun.list.OTU_id, by = "OTU")
write.xlsx(daOTU_fun.id, "Fun_daOTUs_and_RF_with id.xlsx")