##Calculate contribution of each factor to compositional varations estimated by beta diversity (PCoA)
bray1.bac <-  ordinate(bac.clean.log, 'PCoA', 'bray')
bray1.fun <-  ordinate(fun.clean.log, 'PCoA', 'bray')
b.meta.numeric <- b.meta %>% select(SampleID, PlantSpecies)
b.meta.numeric
f.meta.numeric <- f.meta %>% select(SampleID, PlantSpecies)
f.meta.numeric

corr.all<-cbind(b.meta.numeric, bray1.bac$vector)
class(corr.all)
corr.all$PlantSpecies <- factor(corr.all$PlantSpecies)

corr.all<-as.matrix(corr.all)



cor_5 <- Hmisc::rcorr(as.matrix(corr.all), type="spearman")


Error in Hmisc::rcorr(as.matrix(corr.all), type = "spearman") : 
  NA/NaN/Inf in foreign function call (arg 1)
In addition: Warning message:
  In storage.mode(x) <- "double" : NAs introduced by coercion


M <- cor_5$r
p_mat <- cor_5$P

cormat <- reorder_cormat(M)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat<-subset(melted_cormat, Var1 %in% colnames(b.meta.numeric))
melted_cormat<-subset(melted_cormat, Var2 %in% colnames(corr.all[,-c(1:4)]))


##p-value table
p_mat[is.na(p_mat)] <-1
cormat.p <- reorder_cormat(p_mat)
upper_tri.p <- get_upper_tri(cormat.p)
# Melt the correlation matrix
melted_cormat.p <- melt(upper_tri.p, na.rm = TRUE)

melted_cormat.p<-subset(melted_cormat.p, Var2 %in% colnames(corr.all[,-c(1:4)]))
melted_cormat.p<-subset(melted_cormat.p, Var1 %in% colnames(b.meta.numeric))
head(melted_cormat.p)
head(melted_cormat)

names(melted_cormat.p)[3] <- "P_value"
melted_cormat.r.p <-merge(melted_cormat, melted_cormat.p, by = c("Var1" = "Var1", "Var2" = "Var2"))
melted_cormat.r.p
melted_cormat.r.p$sig <- gtools::stars.pval(melted_cormat.r.p$P_value)

###Ordering subject factors
melted_cormat.r.p$Var1 <- factor(melted_cormat.r.p$Var1, levels = c("Age","BMI","anti_CCP","RA_factor"))

melted_cormat.r.p.sig <- subset(melted_cormat.r.p, P_value < 0.05)


melted_cormat.r.p.sig.bac <- melted_cormat.r.p.sig




#Fungi
corr.all<-cbind(f.meta.numeric, bray1.fun$vector)
class(corr.all)
corr.all<-as.matrix(corr.all)

cor_5 <- Hmisc::rcorr(corr.all, type="spearman")
M <- cor_5$r
p_mat <- cor_5$P

cormat <- reorder_cormat(M)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat.1<-subset(melted_cormat, Var2 %in% colnames(f.meta.numeric))
melted_cormat.2<-subset(melted_cormat, Var1 %in% colnames(f.meta.numeric))

melted_cormat.1<-melted_cormat.1[,c(2,1,3)]
names(melted_cormat.1)[1] <- "Var1"
names(melted_cormat.1)[2] <- "Var2"

melted_cormat <- rbind(melted_cormat.1, melted_cormat.2)

melted_cormat<-subset(melted_cormat, Var2 %in% colnames(corr.all[,-c(1:6)]))


##p-value table
p_mat[is.na(p_mat)] <-1
cormat.p <- reorder_cormat(p_mat)
upper_tri.p <- get_upper_tri(cormat.p)
# Melt the correlation matrix
melted_cormat.p <- melt(upper_tri.p, na.rm = TRUE)

melted_cormat.p.1<-subset(melted_cormat.p, Var2 %in% colnames(f.meta.numeric))
melted_cormat.p.2<-subset(melted_cormat.p, Var1 %in% colnames(f.meta.numeric))

melted_cormat.p.1<-melted_cormat.p.1[,c(2,1,3)]
names(melted_cormat.p.1)[1] <- "Var1"
names(melted_cormat.p.1)[2] <- "Var2"

melted_cormat.p <- rbind(melted_cormat.p.1, melted_cormat.p.2)

melted_cormat.p<-subset(melted_cormat.p, Var2 %in% colnames(corr.all[,-c(1:6)]))


names(melted_cormat.p)[3] <- "P_value"

melted_cormat.r.p <-merge(melted_cormat, melted_cormat.p, by = c("Var1" = "Var1", "Var2" = "Var2"))
melted_cormat.r.p
melted_cormat.r.p$sig <- gtools::stars.pval(melted_cormat.r.p$P_value)

###Ordering subject factors
melted_cormat.r.p$Var1 <- factor(melted_cormat.r.p$Var1, levels = c("Age","BMI","anti_CCP","RA_factor","Candida","Aspergillus"))

melted_cormat.r.p.sig <- subset(melted_cormat.r.p, P_value < 0.05)
melted_cormat.r.p.sig.bac

##Unconstrained PCoA
plot_ordination(bac.clean.log, bray1.bac, type = "samples", color='PlantSpecies', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())


plot_ordination(fun.clean.log, bray1.fun, type = "samples", color='RA_factor', axes = c(2,3))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "darkgreen",
                        high = "gold3",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())




### Add Candida and Aspergillus abundance in meta data (fungi)

f.meta
data.frame(t(otu.fun.genus.log.candida))
t(otu.fun.genus.log.asper)

f.meta.its1 <- f.meta

f.meta.its1$Candida <- data.frame(t(otu.fun.genus.log.its1.candida))$Candida
f.meta.its1$Aspergillus <- data.frame(t(otu.fun.genus.log.its1.asper))$Aspergillus

fun.clean.log.2.its1 <- fun.its1.clean.log
sample_data(fun.clean.log.2.its1) <-sample_data(f.meta.its1)

f.meta.numeric.2.its1 <- f.meta %>% select(Age, BMI, anti_CCP, RA_factor, Candida, Aspergillus)


bray1.fun <-  ordinate(fun.clean.log.2.its1, 'PCoA', 'bray') #fun.clean.log.2 (its2), fun.clean.log.2.its1 (its1) 
#Fungi
corr.all<-cbind(f.meta.numeric.2.its1, bray1.fun$vector)
class(corr.all)
corr.all<-as.matrix(corr.all)

cor_5 <- Hmisc::rcorr(corr.all, type="spearman")
M <- cor_5$r
p_mat <- cor_5$P

cormat <- reorder_cormat(M)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat.1<-subset(melted_cormat, Var2 %in% colnames(f.meta.numeric.2))
melted_cormat.2<-subset(melted_cormat, Var1 %in% colnames(f.meta.numeric.2))

melted_cormat.1<-melted_cormat.1[,c(2,1,3)]
names(melted_cormat.1)[1] <- "Var1"
names(melted_cormat.1)[2] <- "Var2"

melted_cormat <- rbind(melted_cormat.1, melted_cormat.2)

melted_cormat<-subset(melted_cormat, Var2 %in% colnames(corr.all[,-c(1:6)]))


##p-value table
p_mat[is.na(p_mat)] <-1
cormat.p <- reorder_cormat(p_mat)
upper_tri.p <- get_upper_tri(cormat.p)
# Melt the correlation matrix
melted_cormat.p <- melt(upper_tri.p, na.rm = TRUE)

melted_cormat.p.1<-subset(melted_cormat.p, Var2 %in% colnames(f.meta.numeric.2))
melted_cormat.p.2<-subset(melted_cormat.p, Var1 %in% colnames(f.meta.numeric.2))

melted_cormat.p.1<-melted_cormat.p.1[,c(2,1,3)]
names(melted_cormat.p.1)[1] <- "Var1"
names(melted_cormat.p.1)[2] <- "Var2"

melted_cormat.p <- rbind(melted_cormat.p.1, melted_cormat.p.2)

melted_cormat.p<-subset(melted_cormat.p, Var2 %in% colnames(corr.all[,-c(1:6)]))


names(melted_cormat.p)[3] <- "P_value"

melted_cormat.r.p <-merge(melted_cormat, melted_cormat.p, by = c("Var1" = "Var1", "Var2" = "Var2"))
melted_cormat.r.p
melted_cormat.r.p$sig <- gtools::stars.pval(melted_cormat.r.p$P_value)

###Ordering subject factors
melted_cormat.r.p$Var1 <- factor(melted_cormat.r.p$Var1, levels = c("Age","BMI","anti_CCP","RA_factor","Candida","Aspergillus"))

melted_cormat.r.p.sig <- subset(melted_cormat.r.p, P_value < 0.05)

write.xlsx(melted_cormat.r.p.sig,"PCo and variables_ITS1.xlsx")

plot_ordination(fun.clean.log.2.its1, bray1.fun, type = "samples", color='Diagnosis', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  # scale_colour_gradient(low = "darkgreen",
  #                       high = "gold3",
  #                       space = "Lab",
  #                       na.value = "grey50",
  #                       guide = "colourbar",
  #                       aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot_ordination(fun.clean.log.2, bray1.fun, type = "samples", color='Candida', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "#CCCCFF",
                        high = "#333366",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot_ordination(fun.clean.log.2, bray1.fun, type = "samples", color='Aspergillus', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "#CCFFFF",
                        high = "#003333",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

#### Bacterial genus
b.meta$Blautia <- data.frame(t(otu.bac.genus.log.blautia))$Blautia
b.meta$Streptococcus <- data.frame(t(otu.bac.genus.log.strep))$Streptococcus
b.meta$Bifidobacterium <- data.frame(t(otu.bac.genus.log.bifido))$Bifidobacterium

bac.clean.log.2 <- bac.clean.log
sample_data(bac.clean.log.2) <-sample_data(b.meta)

b.meta.numeric.2 <- b.meta %>% select(Age, BMI, anti_CCP, RA_factor, Blautia, Streptococcus, Bifidobacterium)


bray1.bac <-  ordinate(bac.clean.log.2, 'PCoA', 'bray')
#bacteria
corr.all<-cbind(b.meta.numeric.2, bray1.bac$vector)
class(corr.all)
corr.all<-as.matrix(corr.all)

cor_5 <- Hmisc::rcorr(corr.all, type="spearman")
M <- cor_5$r
p_mat <- cor_5$P

cormat <- reorder_cormat(M)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
melted_cormat.1<-subset(melted_cormat, Var2 %in% colnames(b.meta.numeric.2))
melted_cormat.2<-subset(melted_cormat, Var1 %in% colnames(b.meta.numeric.2))

melted_cormat.1<-melted_cormat.1[,c(2,1,3)]
names(melted_cormat.1)[1] <- "Var1"
names(melted_cormat.1)[2] <- "Var2"

melted_cormat <- rbind(melted_cormat.1, melted_cormat.2)

melted_cormat<-subset(melted_cormat, Var2 %in% colnames(corr.all[,-c(1:7)]))


##p-value table
p_mat[is.na(p_mat)] <-1
cormat.p <- reorder_cormat(p_mat)
upper_tri.p <- get_upper_tri(cormat.p)
# Melt the correlation matrix
melted_cormat.p <- melt(upper_tri.p, na.rm = TRUE)

melted_cormat.p.1<-subset(melted_cormat.p, Var2 %in% colnames(b.meta.numeric.2))
melted_cormat.p.2<-subset(melted_cormat.p, Var1 %in% colnames(b.meta.numeric.2))

melted_cormat.p.1<-melted_cormat.p.1[,c(2,1,3)]
names(melted_cormat.p.1)[1] <- "Var1"
names(melted_cormat.p.1)[2] <- "Var2"

melted_cormat.p <- rbind(melted_cormat.p.1, melted_cormat.p.2)

melted_cormat.p<-subset(melted_cormat.p, Var2 %in% colnames(corr.all[,-c(1:7)]))


names(melted_cormat.p)[3] <- "P_value"

melted_cormat.r.p <-merge(melted_cormat, melted_cormat.p, by = c("Var1" = "Var1", "Var2" = "Var2"))
melted_cormat.r.p
melted_cormat.r.p$sig <- gtools::stars.pval(melted_cormat.r.p$P_value)

###Ordering subject factors
melted_cormat.r.p$Var1 <- factor(melted_cormat.r.p$Var1, levels = c("Age","BMI","anti_CCP","RA_factor","Blautia","Streptococcus", "Bifidobacterium"))

melted_cormat.r.p.sig <- subset(melted_cormat.r.p, P_value < 0.05)

write.xlsx(melted_cormat.r.p.sig,"PCo and variables_bacteria.xlsx")


plot_ordination(bac.clean.log.2, bray1.bac, type = "samples", color='Streptococcus', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "#CCFF99",
                        high = "#336600",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())



plot_ordination(bac.clean.log.2, bray1.bac, type = "samples", color='Diagnosis', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  # scale_colour_gradient(low = "#CCFFFF",
  #                       high = "#003333",
  #                       space = "Lab",
  #                       na.value = "grey50",
  #                       guide = "colourbar",
  #                       aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

plot_ordination(bac.clean.log.2, bray1.bac, type = "samples", color='Bifidobacterium', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+
  scale_colour_gradient(low = "#FFCC99",
                        high = "#660000",
                        space = "Lab",
                        na.value = "grey50",
                        guide = "colourbar",
                        aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

dev.off()


Hmisc::rcorr(f.meta.numeric.2$Candida, f.meta.numeric.2$Aspergillus, type = "spearman")





f.cap.diagnosis <- ordinate(fun.clean.log.2, "CAP", "bray", ~ Diagnosis)

perm_anova.ord <- anova.cca(f.cap.diagnosis)
perm_anova.ord2 <- permutest(f.cap.diagnosis)

## Plotting
# Diagnosis
plot.f.cap.type     <- plot_ordination(fun.clean.log.2, f.cap.diagnosis, type = "samples", color='Candida', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ #scale_color_manual(values=c("control" = "#6699CC", "RA" = "#CC9900"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type

plot.f.cap.type     <- plot_ordination(fun.clean.log.2, f.cap.diagnosis, type = "samples", color='Candida', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+scale_colour_gradient(low = "#CCCCFF",
                                             high = "#333366",
                                             space = "Lab",
                                             na.value = "grey50",
                                             guide = "colourbar",
                                             aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type

plot.f.cap.type     <- plot_ordination(fun.clean.log.2, f.cap.diagnosis, type = "samples", color='Aspergillus', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+scale_colour_gradient(low = "#CCFFFF",
                                             high = "#003333",
                                             space = "Lab",
                                             na.value = "grey50",
                                             guide = "colourbar",
                                             aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type




### Bacteria CAP
plot.b.cap.type     <- plot_ordination(bac.clean.log.2, b.cap.diagnosis, type = "samples", color='Bifidobacterium', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+scale_colour_gradient(low = "#FFCC99",
                                             high = "#660000",
                                             space = "Lab",
                                             na.value = "grey50",
                                             guide = "colourbar",
                                             aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type


plot.b.cap.type     <- plot_ordination(bac.clean.log.2, b.cap.diagnosis, type = "samples", color='Streptococcus', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+  scale_colour_gradient(low = "#CCFF99",
                                               high = "#336600",
                                               space = "Lab",
                                               na.value = "grey50",
                                               guide = "colourbar",
                                               aesthetics = "colour")+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type