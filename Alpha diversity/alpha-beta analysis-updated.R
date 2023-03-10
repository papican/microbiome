#beta ve alpha modification
theme_set(theme_bw())
sample_data(fun.clean.ss)$PlantSpecies <- gsub(" ", "", sample_data(fun.clean.ss)$PlantSpecies)
#alpha
# prune OTUs that are not present in at least one sample
GP <- prune_taxa(taxa_sums(bac.clean.ss) > 0, bac.clean.ss)
# Define a human-associated versus non-human categorical variable:
PlantSpecies <- get_variable(GP, "PlantSpecies") %in% c("Panaxginseng", "Araliacordata", "Cynanchumwilfordii", "Cynanchumauriculatum", "Angelicagigas", "Peucedanumjaponicum","Astragaluspropinquus","Glycyrrhizauralensis")
# Add new human variable to sample data:
sample_data(GP)$PlantSpecies <- factor(PlantSpecies)

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
p <- plot_richness(GP, "PlantSpecies", "PlantSpecies" , measures=alpha_meas)+ geom_boxplot(data=p$data, aes(x=PlantSpecies, y=value, color=PlantSpecies), alpha=0.1)
p
#fungi
GP <- prune_taxa(taxa_sums(fun.clean.ss) > 0, fun.clean.ss)
# Define a human-associated versus non-human categorical variable:
PlantSpecies <- get_variable(GP, "PlantSpecies") %in% c("Panaxginseng", "Araliacordata", "Cynanchumwilfordii", "Cynanchumauriculatum", "Angelicagigas", "Peucedanumjaponicum","Astragaluspropinquus","Glycyrrhizauralensis")
# Add new human variable to sample data:
sample_data(GP)$PlantSpecies <- factor(sample_data(GP)$PlantSpecies)

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
p <- plot_richness(GP, "PlantSpecies", "PlantSpecies" , measures=alpha_meas)+ geom_boxplot(data=p$data, aes(x=PlantSpecies, y=value, color=PlantSpecies), alpha=0.1)
p
#beta
#CLR transform
ps_clr <- microbiome::transform(bac.clean.ss, "clr")
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




