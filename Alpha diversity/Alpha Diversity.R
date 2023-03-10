library(egg)
### Alpha diversity
###Rarefying
bac.clean.ss

bac.rarefy <- rarefy_even_depth(bac.clean.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
fun.rarefy <- rarefy_even_depth(fun.clean.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
bac.rarefy #3029 taxa #217OTUs were removed because they are no longer present in any sample after random subsampling
fun.rarefy #854 taxa #20OTUs were removed because they are no longer present in any sample after random subsampling

bac.rarefy.f <- prune_taxa(taxa_sums(bac.rarefy) > 0, bac.rarefy)
#tab_all <- microbiome::alpha(phy.rarefy.2017.t.oligo.t, index = "all")
#kable(head(tab_all))
#write.table(tab_all, "Alpha diversity_all.txt", sep = "\t", row.names = TRUE,  quote = TRUE, na = "NA")
ps1.meta <- data.frame(sample_data(bac.clean.ss))
ps1.meta

tab_shannon <- microbiome::alpha(bac.rarefy.f, index = "shannon")
tab_simpson <- microbiome::alpha(bac.rarefy.f, index = "evenness_simpson")
tab_observed <- microbiome::alpha(bac.rarefy.f, index = "observed")

ps1.meta$observed <- tab_observed$observed
ps1.meta$shannon <- tab_shannon$diversity_shannon
ps1.meta$simpson <- tab_simpson$evenness_simpson


ps1.meta$Group <- factor(ps1.meta$PlantSpecies, levels = order.sample)

max.diversity <- aggregate(ps1.meta$observed, by = list(ps1.meta$PlantSpecies), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

#Normality Test
shapiro.test(ps1.meta$observed) #W = 0.93806, p-value = 0.003084 
##p value lower than 0.05 data does not follow normal distribution, therefore;

##Kruskal-Wallis test
kw<-kruskal.test(ps1.meta$observed ~ ps1.meta$PlantSpecies)
kw$p.value #p-value=3.817474e-05 rounded up to 3.82e-05

library(FSA)
DT = dunnTest(observed ~ PlantSpecies,
              data=ps1.meta,
              method="bh")
PT = DT$res
PT
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn$Group <- max.diversity$Group
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=PlantSpecies, y=observed)) + geom_boxplot(fill="indianred", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = 3.82e-05 '))+ 
  ylab("Observed OTUs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 1)

p




write.csv(PT, "Statistical analysis_richness_Bacteria.csv")

#Normality test 
shapiro.test(ps1.meta$shannon) #W = 0.97532, p-value = 0.2272 => Data follows normal distribution -> ANOVA

##ANOVA
max.diversity <- aggregate(ps1.meta$shannon, by = list(ps1.meta$PlantSpecies), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")
aov<-aov(ps1.meta$observed ~ ps1.meta$PlantSpecies)
summary(aov)#3.27e-06
#library(FSA)
DT = dunnTest(shannon ~ PlantSpecies,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn$Group <- max.diversity$Group
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=PlantSpecies, y=shannon)) + geom_boxplot(fill="darkolivegreen", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('ANOVA, P-value = 3.27e-06'))+
  ylab("Shannon Index\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 1)

p

write.csv(PT, "Statistical analysis_Shannon_Bacteria.csv")


#Normality test
shapiro.test(ps1.meta$simpson)# W = 0.91697, p-value = 0.0003723, p<0.05 so kruskalwallis

##ANOVA
max.diversity <- aggregate(ps1.meta$simpson, by = list(ps1.meta$PlantSpecies), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simpson ~ PlantSpecies, data = ps1.meta)
kw$p.value # 8.517741e-06 rounded up to 8.52e-06


#library(FSA)
DT = dunnTest(simpson ~ PlantSpecies,
              data=ps1.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn$Group <- max.diversity$Group
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps1.meta, aes(x=PlantSpecies, y=simpson)) + geom_boxplot(fill="darkorchid", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = 8.52e-06'))+
  ylab("Inverse Simpson Index\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
theme(aspect.ratio = 1)
p
write.csv(PT, "Statistical analysis_simpson_bacteria.csv")


###Fungi
fun.rarefy.f <- prune_taxa(taxa_sums(fun.rarefy) > 0, fun.rarefy)
#tab_all <- microbiome::alpha(phy.rarefy.2017.t.oligo.t, index = "all")
#kable(head(tab_all))
#write.table(tab_all, "Alpha diversity_all.txt", sep = "\t", row.names = TRUE,  quote = TRUE, na = "NA")
ps2.meta <- data.frame(sample_data(fun.clean.ss))
ps2.meta

tab_shannon <- microbiome::alpha(fun.rarefy.f, index = "shannon")
tab_simpson <- microbiome::alpha(fun.rarefy.f, index = "evenness_simpson")
tab_observed <- microbiome::alpha(fun.rarefy.f, index = "observed")

ps2.meta$observed <- tab_observed$observed
ps2.meta$shannon <- tab_shannon$diversity_shannon
ps2.meta$simpson <- tab_simpson$evenness_simpson

ps2.meta$Group <-  factor(ps2.meta$PlantSpecies, levels = order.sample)

max.diversity <- aggregate(ps2.meta$observed, by = list(ps2.meta$PlantSpecies), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

#Normality test
shapiro.test(ps2.meta$observed)  #  W = 0.97695, p-value = 0.2742 -> Normal Distribution -> ANOVA
##ANOVA
max.diversity <- aggregate(ps2.meta$shannon, by = list(ps2.meta$PlantSpecies), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")
aov<-aov(ps2.meta$observed ~ ps2.meta$PlantSpecies)
summary(aov)#5.08e-07


#library(FSA)
DT = dunnTest(observed ~ PlantSpecies,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn$Group <- max.diversity$Group
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=PlantSpecies, y=observed)) + geom_boxplot(fill="indianred", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('ANOVA, P-value = 5.08e-07'))+
  ylab("Observed OTUs\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 1)

p

write.csv(PT, "Statistical analysis_richness_Fungi.csv")

#Normality test
shapiro.test(ps2.meta$shannon) #W = 0.95046, p-value = 0.01205 No normal distribution


max.diversity <- aggregate(ps2.meta$shannon, by = list(ps2.meta$PlantSpecies), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(shannon ~ PlantSpecies, data = ps2.meta)
kw$p.value #0.002308936

#library(FSA)
DT = dunnTest(shannon ~ PlantSpecies,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn$Group <- max.diversity$Group
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=PlantSpecies, y=shannon)) + geom_boxplot(fill="darkolivegreen", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = 2.30e-03'))+
  ylab("Shannon Index\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 1)

p



write.csv(PT, "Statistical analysis_Shannon_Fungi.csv")

#Normality test
shapiro.test(ps2.meta$simpson) #W = 0.89938, p-value = 7.584e-05 -> kruskalwallis
max.diversity <- aggregate(ps2.meta$simpson, by = list(ps2.meta$PlantSpecies), max)
colnames(max.diversity) <- c("Group", "MaxDiversity")

##Kruskal-Wallis test
kw<-kruskal.test(simpson ~ PlantSpecies, data = ps2.meta)
kw$p.value #0.04368092


#library(FSA)
DT = dunnTest(simpson ~ PlantSpecies,
              data=ps2.meta,
              method="bh")
PT = DT$res
#library(rcompanion)
dunn<-cldList(P.adj ~ Comparison,
              data = PT,
              threshold = 0.05)
dunn$Group <- max.diversity$Group
hsd1 <- merge(dunn,max.diversity, by = 'Group')
p <- ggplot(data = ps2.meta, aes(x=PlantSpecies, y=simpson)) + geom_boxplot(fill="darkorchid", width = 0.8) +
  theme_bw() +
  geom_point(position='jitter',shape=1, alpha=.5)+
  geom_text(data=hsd1,aes(x=Group,y=MaxDiversity, label=Letter), vjust=-1) +
  xlab(paste('Kruskal-Wallis, P-value = 0.04'))+
  ylab("Inverse Simpson Index\n") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.4)) + theme(legend.position = "none")+
  theme(aspect.ratio = 1)

p

write.csv(PT, "Statistical analysis_Simpson_Fungi.csv")
