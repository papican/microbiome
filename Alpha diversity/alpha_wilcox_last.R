##starting files
bac.rarefy <- rarefy_even_depth(bac.AR.AP.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
bac.rarefy #3029 taxa #217OTUs were removed because they are no longer present in any sample after random subsampling

bac.rarefy.f <- prune_taxa(taxa_sums(bac.rarefy) > 0, bac.rarefy)

fun.rarefy <- rarefy_even_depth(fun.AR.AP.ss, rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
fun.rarefy #3029 taxa #217OTUs were removed because they are no longer present in any sample after random subsampling

fun.rarefy.f <- prune_taxa(taxa_sums(fun.rarefy) > 0, fun.rarefy)

##metadata
b.meta <- sample_data(bac.AR.AP.ss)
b.meta <- data.frame(b.meta)

b.meta$new <- as.factor(b.meta$new)
b.meta

f.meta <- sample_data(fun.AR.AP.ss)
f.meta <- data.frame(f.meta)

f.meta$new <- as.factor(f.meta$new)
f.meta
########## Single files bac ###############
tab_all.b <- microbiome::alpha(bac.rarefy.f, index = "all")
write.table(tab_all.b, "Alpha_diversity_bac_Araliaceae_Apiaceae.txt", sep = "\t", row.names = TRUE,  quote = TRUE, na = "NA")
tab_all.f <- microbiome::alpha(fun.rarefy.f, index = "all")
write.table(tab_all.f, "Alpha_diversity_fun_Araliaceae_Apiaceae.txt", sep = "\t", row.names = TRUE,  quote = TRUE, na = "NA")

ps2.meta <- b.meta

## using Observed OTUs
ps2.meta$Observed <- tab_all.b$observed 
ps2.meta


##

##Richness (Observed OTUs)
max.diversity <- aggregate(ps2.meta$Observed, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")



# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$Observed
y <- subset(ps2.meta, new=='Apiaceae')$Observed
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.000931

p <- ggplot(data = ps2.meta, aes(x=new, y=Observed, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Observed OTU \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p



##Shannon
ps2.meta$Shannon <- tab_all.b$diversity_shannon 
ps2.meta$Shannon


max.diversity <- aggregate(ps2.meta$Shannon, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")



# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$Shannon
y <- subset(ps2.meta, new=='Apiaceae')$Shannon
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.002953

p <- ggplot(data = ps2.meta, aes(x=new, y=Shannon, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Shannon \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p



##Simpson evenness
ps2.meta$simp <- tab_all.b$evenness_simpson 
ps2.meta$simp

max.diversity <- aggregate(ps2.meta$simp, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")



# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$simp
y <- subset(ps2.meta, new=='Apiaceae')$simp
wilcox.test(x, y, conf.int = TRUE) #0.0003108


p <- ggplot(data = ps2.meta, aes(x=new, y=simp, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Simpson \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p

####invsimp bac
## using  OTUs
ps2.meta$diversity_inverse_simpson <- tab_all.b$diversity_inverse_simpson 
ps2.meta


##

##inverted simp
max.diversity <- aggregate(ps2.meta$diversity_inverse_simpson, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")



# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$diversity_inverse_simpson
y <- subset(ps2.meta, new=='Apiaceae')$diversity_inverse_simpson
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.08298

p <- ggplot(data = ps2.meta, aes(x=new, y=diversity_inverse_simpson, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Inv. Simpson\n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p

####chao bac
## using  OTUs
ps2.meta$chao1 <- tab_all.b$chao1 
ps2.meta


##

##chao1

max.diversity <- aggregate(ps2.meta$chao1, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")

# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$chao1
y <- subset(ps2.meta, new=='Apiaceae')$chao1
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.0001554

p <- ggplot(data = ps2.meta, aes(x=new, y=chao1, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Chao1 \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p
########## CW-CA fun ###############

ps2.meta <- f.meta

## using Observed OTUs
ps2.meta$Observed <- tab_all.f$observed 
ps2.meta


#new
##Richness (Observed OTUs)
max.diversity <- aggregate(ps2.meta$Observed, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")



# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$Observed
y <- subset(ps2.meta, new=='Apiaceae')$Observed
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.958

p <- ggplot(data = ps2.meta, aes(x=new, y=Observed, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Observed OTU \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 10,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p



##Shannon
ps2.meta$Shannon <- tab_all.f$diversity_shannon 
ps2.meta$Shannon


max.diversity <- aggregate(ps2.meta$Shannon, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")



# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$Shannon
y <- subset(ps2.meta, new=='Apiaceae')$Shannon
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.01476

p <- ggplot(data = ps2.meta, aes(x=new, y=Shannon, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Shannon \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p



##Simpson evenness
ps2.meta$simp <- tab_all.f$evenness_simpson 
ps2.meta$simp

max.diversity <- aggregate(ps2.meta$simp, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")



# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$simp
y <- subset(ps2.meta, new=='Apiaceae')$simp
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.01476


p <- ggplot(data = ps2.meta, aes(x=new, y=simp, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Simpson \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p


##chao fun
##chao1
ps2.meta$chao1 <- tab_all.f$chao1 
ps2.meta$chao1

max.diversity <- aggregate(ps2.meta$chao1, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")



# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$chao1
y <- subset(ps2.meta, new=='Apiaceae')$chao1
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.01476


p <- ggplot(data = ps2.meta, aes(x=new, y=chao1, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Chao1 \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p


##inverted simp fun
##inverted simp
ps2.meta$diversity_inverse_simpson <- tab_all.f$diversity_inverse_simpson
ps2.meta$diversity_inverse_simpson

max.diversity <- aggregate(ps2.meta$diversity_inverse_simpson, by = list(ps2.meta$new), max)

colnames(max.diversity) <- c("Group", "MaxDiversity")



# wilcoxon test
x <- subset(ps2.meta, new=='Araliaceae')$diversity_inverse_simpson
y <- subset(ps2.meta, new=='Apiaceae')$diversity_inverse_simpson
wilcox.test(x, y, conf.int = TRUE) #p-value = 0.01476


p <- ggplot(data = ps2.meta, aes(x=new, y=diversity_inverse_simpson, fill = new)) + geom_boxplot(width = 0.8) +
  theme_bw() + scale_fill_manual(values=c("Araliaceae" = "indianred", "Apiaceae" = "seagreen"))+
  geom_point(position='jitter',shape=1, alpha=.5)+theme(aspect.ratio=2)+
  xlab(paste('P = ', wilcox.test(unlist(x), unlist(y), conf.int = TRUE)$p.value))+
  ylab("Inv. Simpson \n") + theme(plot.title = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")

p
