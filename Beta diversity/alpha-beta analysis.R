#beta ve alpha modification
theme_set(theme_bw())
sample_data(bac.clean.ss)$PlantSpecies <- gsub(" ", "", sample_data(bac.clean.ss)$PlantSpecies)
#alpha
# prune OTUs that are not present in at least one sample
GP <- prune_taxa(taxa_sums(fun.clean.ss) > 0, fun.clean.ss)
# Define a human-associated versus non-human categorical variable:
sample_data(GP)$PlantSpecies <- get_variable(GP, "PlantSpecies") %in% c("Panaxginseng", "Araliacordata", "Cynanchumwilfordii", "Cynanchumauriculatum", "Angelicagigas", "Peucedanumjaponicum","Astragaluspropinquus","Glycyrrhizauralensis")
# Add new human variable to sample data:
sample_data(GP)$PlantSpecies <- factor(PlantSpecies)

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
p <- plot_richness(GP, x="PlantSpecies" ,color = "PlantSpecies", measures=alpha_meas)
p<- p+geom_boxplot(data=p$data, aes(x=PlantSpecies, y=value, color=PlantSpecies), alpha=0.1)
p
#fungi
GP <- prune_taxa(taxa_sums(bac.CW.CA.ss) > 0, bac.CW.CA.ss)
# Define a human-associated versus non-human categorical variable:
PlantSpecies <- get_variable(GP, "PlantSpecies") %in% c("Cynanchumwilfordii", "Cynanchumauriculatum")
# Add new human variable to sample data:
sample_data(GP)$PlantSpecies <- factor(sample_data(GP)$PlantSpecies)

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
p <- plot_richness(GP, "PlantSpecies" ,color = "PlantSpecies", measures=alpha_meas)+ geom_boxplot(data=p$data, aes(x=PlantSpecies, y=value, color=PlantSpecies), alpha=0.1)
p
#beta
#CLR transform
ps_clr <- microbiome::transform(bac.CW.CA.ss, "clr")
#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")

#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")


#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)   

sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig)) 

#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(bac.clean.ss, ord_clr, type="SampleID", color="PlantSpecies") + 
  geom_point(size = 3, ) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = PlantSpecies), linetype = 2)+  scale_color_manual(values = c("Cynanchumwilfordii"= "darkorange", "Cynanchumauriculatum"= "darkgreen"))+
theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) 

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
#ADONIS test
vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$PlantSpecies)

#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$PlantSpecies)
dispr

p <- plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", col= c("darkorange", "darkgreen"), lwd=2 )+
points(p, "sites", pch=18, col="darkorange", bg="darkgreen")+
text(p, "centroids", col="black", cex=0.9)
p

boxplot(dispr, main = "", xlab = "")
permutest(dispr)

#fungi
#beta
#CLR transform
ps_clr <- microbiome::transform(fun.clean.ss, "clr")
#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")

#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

#Examine eigenvalues and % prop. variance explained
head(ord_clr$CA$eig)   

sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig)) 

#Scale axes and plot ordination
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(fun.clean.ss, ord_clr, type="SampleID", color="PlantSpecies") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group = PlantSpecies), linetype = 2)

#Generate distance matrix
clr_dist_matrix <- phyloseq::distance(ps_clr, method = "euclidean") 
#ADONIS test
vegan::adonis2(clr_dist_matrix ~ phyloseq::sample_data(ps_clr)$PlantSpecies)

#Dispersion test and plot
dispr <- vegan::betadisper(clr_dist_matrix, phyloseq::sample_data(ps_clr)$PlantSpecies)
dispr

plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")


boxplot(dispr, main = "", xlab = "")
permutest(dispr)

#subsample reads
#Subsample reads
ps_rare <- phyloseq::rarefy_even_depth(fun.clean.ss, rngseed = 123, replace = FALSE)

#Generate distances
ord_unifrac <- ordinate(ps_rare, method = "PCoA", distance = "wunifrac") 
ord_unifrac_un <- ordinate(ps_rare, method = "PCoA", distance = "unifrac")   
#Plot ordinations
a <- plot_ordination(ps_rare, ord_unifrac, color = "PlantSpecies") + geom_point(size = 2)
b <- plot_ordination(ps_rare, ord_unifrac_un, color = "PlantSpecies") + geom_point(size = 2)
cowplot::plot_grid(a, b, nrow = 1, ncol = 2, scale = .9, labels = c("Weighted", "Unweighted"))


###CW-CA

GP <- prune_taxa(taxa_sums(fun.CW.CA.ss) > 0, fun.CW.CA.ss)
# Define a human-associated versus non-human categorical variable:
sample_data(GP)$PlantSpecies <- get_variable(GP, "PlantSpecies") %in% c("Cynanchumwilfordii", "Cynanchumauriculatum")
# Add new human variable to sample data:
sample_data(GP)$PlantSpecies <- factor(sample_data(GP)$PlantSpecies)
GP.chl <- subset_taxa(bac.CW.CA.ss, Genus=="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium")
ntaxa(GP.chl)
plot_tree(GP.chl, color="PlantSpecies", shape="Family", label.tips="Genus", size="Abundance")



alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
p <- plot_richness(GP, x="PlantSpecies" , color = "PlantSpecies" , measures=alpha_meas) + geom_boxplot(data=p$data, aes(x=PlantSpecies, y=value, color=PlantSpecies), alpha=0.1)+
  scale_color_manual(values = c("Cynanchumwilfordii"= "indianred", "Cynanchumauriculatum"= "navyblue"))+
  geom_boxplot(lwd = 0.75 )+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 10,hjust = 0.5, face='bold')) 
p
