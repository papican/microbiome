##Calculate contribution of each factor to compositional varations estimated by beta diversity (PCoA)
#bac
b.map <- sample_data(bac.clean.log)
b.map <- data.frame(b.map)
#fun
f.map <- sample_data(fun.clean.log)
f.map <- data.frame(f.map)


##Unconstrained PCoA


bray1.bac <-  ordinate(bac.clean.log, 'PCoA', 'bray')
bray1.fun <-  ordinate(fun.clean.log, 'PCoA', 'bray')


##Bacteria/Sample
plot_ordination(bac.clean.log, bray1.bac, type = "samples", color = "new" , shape = "PlantSpecies", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4) + scale_shape_manual(values=c(18, 17, 16, 15))+ scale_color_manual(values=c("Araliaceae" = "olivedrab4", "Apiaceae" = "deeppink2"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
##Fungi/Sample
plot_ordination(fun.clean.log, bray1.fun, type = "samples", color = "new" , shape = "PlantSpecies", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_shape_manual(values=c(18, 17, 16, 15))+ scale_color_manual(values=c("Araliaceae" = "olivedrab4", "Apiaceae" = "deeppink2"))+
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
bac.clean.ss.rel <- microbiome::transform(bac.AR.AP.ss, "compositional")
fun.clean.ss.rel <- microbiome::transform(fun.AR.AP.ss, "compositional")

f.cap.duration <- ordinate(fun.clean.ss.rel, "CAP", "bray", ~ new)
b.cap.duration <- ordinate(bac.clean.ss.rel, "CAP", "bray", ~ new)

perm_anova.b.ord <- anova.cca(b.cap.duration)
perm_anova.ord2 <- permutest(b.cap.duration)
write.xlsx(perm_anova.b.ord, "AR-AP bac Cap-Rel.xlsx")
perm_anova.f.ord <- anova.cca(f.cap.duration)
perm_anova.ord2 <- permutest(f.cap.duration)
perm_anova.f.ord
write.xlsx(perm_anova.b.ord, "AR-AP fun Cap-Rel.xlsx")


plot.b.cap.type     <- plot_ordination(bac.clean.ss.rel, b.cap.duration, type = "samples", color = "new" , shape = "PlantSpecies", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_shape_manual(values=c(18, 17, 16, 15))+ scale_color_manual(values=c("Araliaceae" = "olivedrab4", "Apiaceae" = "deeppink2"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type

plot.f.cap.type     <- plot_ordination(fun.clean.ss.rel, f.cap.duration, type = "samples", color = "new" , shape = "PlantSpecies", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_shape_manual(values=c(18, 17, 16, 15))+ scale_color_manual(values=c("Araliaceae" = "olivedrab4", "Apiaceae" = "deeppink2"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type





##loglu olan
b.cap.diagnosis <- ordinate(bac.clean.log, "CAP", "bray", ~ new)
f.cap.diagnosis <- ordinate(fun.clean.log, "CAP", "bray", ~ new)


perm_anova.b.ord <- anova.cca(b.cap.diagnosis)
write.xlsx(perm_anova.b.ord, "AR-AP bac Cap-log.xlsx")
perm_anova.b.ord
perm_anova.f.ord <- anova.cca(f.cap.diagnosis)
write.xlsx(perm_anova.b.ord, "AR-AP fun Cap-log.xlsx")
perm_anova.f.ord

##Bacteria/Sample
plot.b.cap.type     <- plot_ordination(bac.clean.log, b.cap.diagnosis, type = "samples", color = "new" , shape = "PlantSpecies", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_shape_manual(values=c(18, 17, 16, 15))+ scale_color_manual(values=c("Araliaceae" = "olivedrab4", "Apiaceae" = "deeppink2"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type
##Fungi/Sample
plot.f.cap.type     <- plot_ordination(fun.clean.log, f.cap.diagnosis, type = "samples", color = "new" , shape = "PlantSpecies", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_shape_manual(values=c(18, 17, 16, 15))+ scale_color_manual(values=c("Araliaceae" = "olivedrab4", "Apiaceae" = "deeppink2"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type

##Normalized abundance
b.cap.diagnosis <- ordinate(bac.clean.nolog, "CAP", "bray", ~ new)
f.cap.diagnosis <- ordinate(fun.clean.nolog, "CAP", "bray", ~ new)

perm_anova.b.ord <- anova.cca(b.cap.diagnosis)
write.xlsx(perm_anova.b.ord, "AR-AP bac nolog.xlsx")

perm_anova.f.ord <- anova.cca(f.cap.diagnosis)
write.xlsx(perm_anova.f.ord, "AR-AP fun nolog.xlsx")



plot.b.cap.type     <- plot_ordination(bac.clean.nolog, b.cap.diagnosis, type = "samples", color = "new" , shape = "PlantSpecies", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_shape_manual(values=c(18, 17, 16, 15))+ scale_color_manual(values=c("Araliaceae" = "olivedrab4", "Apiaceae" = "deeppink2"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.b.cap.type

plot.f.cap.type     <- plot_ordination(fun.clean.nolog, f.cap.diagnosis, type = "samples", color = "new" , shape = "PlantSpecies", axes = c(1,2))+
  theme(plot.title = element_text(size = 20,hjust = 0.5, face="bold")) +
  geom_point(size = 4)+ scale_shape_manual(values=c(18, 17, 16, 15))+ scale_color_manual(values=c("Araliaceae" = "olivedrab4", "Apiaceae" = "deeppink2"))+
  theme(aspect.ratio=1)+theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=15, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  geom_hline(yintercept=0, linetype='dashed', color='black', size = 0.75)+geom_vline(xintercept=0, linetype='dashed', color='black', size = 0.75)+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(), plot.background=element_blank())
plot.f.cap.type

#PERMANOVA

#Bacteria
b.meta <-sample_data(bac.AR.AP.ss)
b.meta <- data.frame(b.meta)
f.meta <- sample_data(fun.AR.AP.ss)
f.meta <- data.frame(f.meta)
#unconstrained
bac.clean.log
b.otu <- otu_table(bac.clean.log)
b.permanova <- adonis2(formula = t(b.otu) ~ (new + PlantSpecies), data = b.meta, permutations=9999, method = "bray")
write.xlsx(b.permanova, "AR-AP Permanova log Bac PS.xlsx")


#CAP/rel
bac.clean.ss.rel
b.otu <- otu_table(bac.clean.ss.rel)
b.permanova <- adonis2(formula = t(b.otu) ~ (new + PlantSpecies), data = b.meta, permutations=9999, method = "bray")
write.xlsx(b.permanova, "AR-AP Permanova Cap-Rel Bac PS.xlsx")


#Nolog
bac.clean.nolog
b.otu <- otu_table(bac.clean.nolog)
b.permanova <- adonis2(formula = t(b.otu) ~ (new + PlantSpecies), data = b.meta, permutations=9999, method = "bray")
write.xlsx(b.permanova, "AR-AP Permanova Nolog Bac PS.xlsx")


##Fungi

#unconstrained
fun.clean.log
b.otu <- otu_table(fun.clean.log)
b.permanova <- adonis2(formula = t(b.otu) ~ (new + PlantSpecies), data = f.meta, permutations=9999, method = "bray")
write.xlsx(b.permanova, "AR-AP Permanova log Fun PS.xlsx")


#CAP/rel
fun.clean.ss.rel
b.otu <- otu_table(fun.clean.ss.rel)
b.permanova <- adonis2(formula = t(b.otu) ~ (new + PlantSpecies), data = f.meta, permutations=9999, method = "bray")
write.xlsx(b.permanova, "AR-AP Permanova Cap-Rel Fun PS.xlsx")

#Nolog
fun.clean.nolog
b.otu <- otu_table(fun.clean.nolog)
b.permanova <- adonis2(formula = t(b.otu) ~ (new + PlantSpecies), data = f.meta, permutations=9999, method = "bray")
write.xlsx(b.permanova, "AR-AP Permanova Nolog Fun PS.xlsx")








