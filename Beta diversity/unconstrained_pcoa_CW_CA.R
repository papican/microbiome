##Calculate contribution of each factor to compositional varations estimated by beta diversity (PCoA)
order.plantspecies <- c("Cynanchum auriculatum", "Cynanchum wilfordii")
#control
sample_data(bac.clean.ss)$PlantSpecies


b.meta <- sample_data(bac.clean.log)
b.meta <- data.frame(b.meta)

b.meta$PlantSpecies <- factor(b.meta$PlantSpecies)
# Calculate compositional version of the data ####### from core otu analysis script I have to get bac.clean.ss.rel data
# (relative abundances)
bac.clean.ss.rel <- microbiome::transform(bac.CW.CA.ss, "compositional")
t(otu_table(bac.clean.ss.rel))

fun.clean.ss.rel <- microbiome::transform(bac.CW.CA.ss, "compositional")
t(otu_table(fun.clean.ss.rel))
###

sample_data(bac.clean.log) <- sample_data(b.meta)
sample_data(bac.CW.CA.ss) <- sample_data(b.meta)
sample_data(bac.clean.nolog) <-sample_data(b.meta)

bac.clean.log
bac.clean.nolog
bac.clean.ss.rel

bray1.bac <-  ordinate(bac.clean.log, 'PCoA', 'bray', ~ PlantSpecies)
sample_data(bac.clean.log) <- sample_data(b.meta)


write.csv(bray1.bac$vectors, "Supplementary Fig. A_PCoA.csv")
write.csv(b.meta, "Source data for PCoA and CAP_bac.csv")
write.csv(f.meta, "Source data for PCoA and CAP_fun.csv") ###f.meta dont forget


##Unconstrained PCoA
#Control vs RA
plot_ordination(bac.CW.CA.ss, bray1.bac, type = "samples", color='PlantSpecies', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c( 
    "Cynanchumauriculatum" = "indianred",
    "Cynanchumwilfordii" = "darkseagreen4"
  ))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### CAP (Canonical analysis of principals)
##Relative abundance
b.cap.PlantSpecies <- ordinate(bac.clean.ss.rel, "CAP", "bray", ~ PlantSpecies)

perm_anova.ord <- anova.cca(b.cap.PlantSpecies)
perm_anova.ord2 <- permutest(b.cap.PlantSpecies)

## Plotting
# PlantSpecies
plot.b.cap.type     <- plot_ordination(bac.clean.ss.rel, b.cap.PlantSpecies, type = "samples", color='PlantSpecies', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c(
    "Cynanchumauriculatum" = "indianred",
    "Cynanchumwilfordii" = "darkseagreen4"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type


##Normalized abundance and log transformed
b.cap.PlantSpecies <- ordinate(bac.clean.log, "CAP", "bray", ~ PlantSpecies)

b.cap.PlantSpecies$CCA$Xbar

perm_anova.ord <- anova.cca(b.cap.PlantSpecies)
perm_anova.ord2 <- permutest(b.cap.PlantSpecies)

## Plotting
# PlantSpecies
plot.b.cap.type     <- plot_ordination(bac.clean.log, b.cap.PlantSpecies, type = "samples", color='PlantSpecies', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c( 
    "Cynanchumauriculatum" = "indianred",
    "Cynanchumwilfordii" = "darkseagreen4"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type



###PERMANOVA
## ALL
bac.clean.ss.rel
bac.clean.log
bac.clean.nolog

b.otu <- otu_table(bac.clean.log)
b.meta <- sample_data(bac.clean.log)
b.meta <- data.frame(b.meta)
b.otu
##new variable for active disease
b.meta$active_disease_2 <- b.meta$active_disease
b.meta$active_disease_2<-as.character(b.meta$active_disease_2)
b.meta$active_disease_2[which(b.meta$active_disease_2 %in% c("0", "1"))] <- "low"
b.meta$active_disease_2[which(b.meta$active_disease_2 %in% c("2", "3"))] <- "high"

b.permanova <- adonis(formula = t(b.otu) ~ (PlantSpecies), data = b.meta, permutations=9999, method = "bray")
b.permanova

b.permanova <- adonis(formula = t(b.otu) ~ (TNF_inhibitor+TX), data = b.meta, permutations=9999, method = "bray")
b.permanova






############### Fungal community #################
##Calculate contribution of each factor to compositional varations estimated by beta diversity (PCoA)
order.plantspecies <- c("Aralia cordata", "Angelica gigas", "Astragalus propinquus" , 
                        "Cynanchum auriculatum", "Cynanchum wilfordii", "Glycyrrhiza uralensis", 
                        "Panax ginseng", "Peucedanum japonicum" )
#control
sample_data(fun.clean.ss)$PlantSpecies

f.meta <- sample_data(fun.clean.log)
f.meta <- data.frame(f.meta)




f.meta$PlantSpecies <- as.character(f.meta$PlantSpecies)



f.meta$PlantSpecies <- factor(f.meta$PlantSpecies)
# Calculate compositional version of the data ####### from core otu analysis script I have to get fun.clean.ss.rel data
# (relative abundances)
fun.clean.ss.rel <- microbiome::transform(fun.CW.CA.ss, "compositional")
t(otu_table(fun.clean.ss.rel))


###

sample_data(fun.clean.log) <- sample_data(f.meta)
sample_data(fun.CW.CA.ss) <- sample_data(f.meta)
sample_data(fun.clean.ss.rel) <-sample_data(f.meta)

fun.clean.log
fun.CW.CA.ss
fun.clean.ss.rel

bray1.fun <-  ordinate(fun.clean.log, 'PCoA', 'bray', ~ PlantSpecies)
sample_data(fun.clean.log) <- sample_data(f.meta)


write.csv(bray1.bac$vectors, "Supplementary Fig. A_PCoA.csv")
write.csv(f.meta, "Source data for PCoA and CAP_bac.csv")
write.csv(f.meta, "Source data for PCoA and CAP_fun.csv") ###f.meta dont forget


##Unconstrained PCoA
#Control vs RA
plot_ordination(fun.CW.CA.ss, bray1.fun, type = "samples", color='PlantSpecies', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c(
    "Cynanchumauriculatum" = "navyblue",
    "Cynanchumwilfordii" = "deeppink"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())

### CAP (Canonical analysis of principals)
##Relative abundance
f.cap.PlantSpecies <- ordinate(fun.clean.ss.rel, "CAP", "bray", ~ PlantSpecies)

perm_anova.ord <- anova.cca(f.cap.PlantSpecies)
perm_anova.ord2 <- permutest(f.cap.PlantSpecies)

## Plotting
# PlantSpecies
plot.f.cap.type     <- plot_ordination(fun.clean.ss.rel, f.cap.PlantSpecies, type = "samples", color='PlantSpecies', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c(
    "Cynanchumauriculatum" = "navyblue",
    "Cynanchumwilfordii" = "deeppink"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type


##Normalized abundance and log transformed
f.cap.PlantSpecies <- ordinate(fun.clean.log, "CAP", "bray", ~ PlantSpecies)

f.cap.PlantSpecies$CCA$Xbar

perm_anova.ord <- anova.cca(f.cap.PlantSpecies)
perm_anova.ord2 <- permutest(f.cap.PlantSpecies)

## Plotting
# PlantSpecies
plot.f.cap.type     <- plot_ordination(fun.clean.log, f.cap.PlantSpecies, type = "samples", color='PlantSpecies', axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 3)+ scale_color_manual(values=c("Cynanchumauriculatum" = "navyblue",
                                                    "Cynanchumwilfordii" = "deeppink"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type





