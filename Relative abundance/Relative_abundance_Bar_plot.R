###### 2. make barplot with phyloseq

order.sample <- c("AC_1_1" , "AC_1_2" , "AC_2_1" , "AC_2_2" , "AC_3_1" , "AC_3_2" , "AC_4_1" , "AC_4_2" , "AG_1_1" , "AG_1_2" , "AG_2_1" , "AG_2_2" ,
                  "AG_3_1" , "AG_3_2" , "AG_4_1" , "AG_4_2" , "AP_1_1" , "AP_1_2" , "AP_2_1" , "AP_2_2" , "AP_3_1" , "AP_3_2" , "AP_4_1" , "AP_4_2" ,
                  "CA_1_1" , "CA_1_2" , "CA_2_1" , "CA_2_2" , "CA_3_1" , "CA_3_2" , "CA_4_1" , "CA_4_2" , "CW_1_1" , "CW_1_2" , "CW_2_1" , "CW_2_2" ,
                  "CW_3_1" , "CW_3_2" , "CW_4_1" , "CW_4_2" , "GU_1_1" , "GU_1_2" , "GU_2_1" , "GU_2_2" , "GU_3_1" , "GU_3_2" , "GU_4_1" , "GU_4_2" ,
                  "PG_1_1" , "PG_1_2" , "PG_2_1" , "PG_2_2" , "PG_3_1" , "PG_3_2" , "PG_4_1" , "PG_4_2" , "PJ1_1_1" , "PJ1_1_2" , "PJ1_2_1" , "PJ1_2_2",
                  "PJ1_3_1" , "PJ1_3_2" , "PJ1_4_1" , "PJ1_4_2")
## Phylum level
df.phylum <- bac.clean.ss %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE) %>% 
  psmelt()



df.phylum$Phylum <- as.character(df.phylum$Phylum)
df.phylum$Phylum2 <- df.phylum$Phylum
df.phylum$Phylum2[which(df.phylum$Class=="Alphaproteobacteria")] <- "Alphaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Gammaproteobacteria")] <- "Gammaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Deltaproteobacteria")] <- "Deltaproteobacteria"

unique(df.phylum$Phylum2)
head(df.phylum)

library(forcats) 
df.phylum %<>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))

df.phylum.2 <- df.phylum %>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))

levels(df.phylum$Phylum2)
levels(df.phylum.2$Phylum2) = c(levels(df.phylum.2$Phylum2), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum.2 %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance <5,]$Phylum2 <- 'Low abundance'
unique(df.phylum$Phylum2)

ord <- df.phylum.rel %>% group_by(Phylum2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Phylum2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac =="Low abundance")]
vec.Low <- c("Low abundance")
vec.reorder <- append(vec.Low, vec.order)


df.phylum.rel$Phylum2 <- factor(df.phylum.rel$Phylum2, levels = vec.reorder) 
table(df.phylum.rel$Phylum2)

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=Sample, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Actinobacteriota"= "darkolivegreen", "Bacteroidota"= "darkolivegreen3", "Proteobacteria"= "darkorange",   
                                "Firmicutes"= "darkorchid",    "Verrucomicrobiota"="indianred", "Myxococcota"="blue" , 
                                 "Low abundance"= "gray", "unidentified"="black"  )) +
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Phylum Community Composition by Plant Species \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 2,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.4)

df.phylum.rel.p1


##Class level
df.class <- bac.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.class %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))


levels(df.class$Class)
levels(df.class$Class) = c(levels(df.class$Class), 'Low abundance')

# we need to group by samples
df.class.rel <- df.class %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.class.rel[df.class.rel$RelAbundance < 5,]$Class <- 'Low abundance'
unique(df.class$Class)

ord <- df.class.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Class
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac =="Low abundance")]
vec.Low <- c("Low abundance")
vec.reorder <- append(vec.Low, vec.order)


df.class.rel$Class <- factor(df.class.rel$Class, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.class.rel$Class)
## plot relative abundance
#bar plot
df.class.rel.p1 <- ggplot(df.class.rel, aes(x = Sample, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Alphaproteobacteria"= "darkslategray2", "Gammaproteobacteria"= "coral", "Bacteroidia"= "darkorange",        
                     "Actinobacteria"= "darkolivegreen",   "Clostridia"="violet", "Chlamydiae"="lightgreen",  
 "Bacilli"= "darkorchid","Polyangia"= "indianred","Low abundance"= "gray", "unidentified"="black")) +
  
  xlab('') +
  ylab("Relative abundance (‰)") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 2,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+theme(aspect.ratio=0.4)

df.class.rel.p1


###Order
df.order <- bac.clean.ss %>%
  tax_glom(taxrank = "Order", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.order %<>% mutate(Order = fct_explicit_na(Order, na_level = "unidentified"))


levels(df.order$Order)
levels(df.order$Order) = c(levels(df.order$Order), 'Low abundance')

# we need to group by samples
df.order.rel <- df.order %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.order.rel[df.order.rel$RelAbundance <5,]$Order <- 'Low abundance'
unique(df.order$Order)

ord <- df.order.rel %>% group_by(Order) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Order
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder <- append(vec.Low, vec.order)


df.order.rel$Order <- factor(df.order.rel$Order, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
table((df.order.rel$Order))
#bar plot
df.order.rel.p1 <- ggplot(df.order.rel, aes(x=Sample, y = RelAbundance, fill = Order)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Pseudomonadales"= "darkolivegreen3" ,                  
                               "Sphingomonadales"= "darkolivegreen" ,                  
                                "Rhizobiales"= "darkorange"         ,               
                               "Burkholderiales"= "darksalmon"      ,             
                               "Caulobacterales"= "darkorchid"      ,              
                                     
                               "Pseudonocardiales"= "darkseagreen4" ,                
                                              
                               "Xanthomonadales"= "darkslategray"   ,                 
                                "Cytophagales"= "darkseagreen1"     ,                  
                               "Flavobacteriales"= "deeppink4"      ,            
                               "Micromonosporales"= "indianred" ,                
                             "Haliangiales"= "cornsilk4"          ,            
                               "Enterobacterales"= "dodgerblue4"    ,              
                               "Streptomycetales"= "blueviolet"     ,              
                               "Steroidobacterales"= "lightpink4"   ,             
                               
                               
                               "Frankiales"= "palevioletred"               ,  
                               "Chlamydiales"= "aquamarine"                     , 
                               "Micrococcales"= "coral"                     , 
                               "Clostridiales"= "lightcoral"           ,     
                               "Bacteroidales"= "gold4"                    , 
                               "Streptosporangiales"= "burlywood3"          ,       
                               "Azospirillales"= "firebrick1"               ,     
                               "Propionibacteriales"= "coral4"              ,  
                               "Paenibacillales"= "lightsteelblue2"                ,           
                               
"Corynebacteriales"= "paleturquoise2",                 
            
           
"Bacillales"= "darkseagreen"        ,                
"Low abundance"= "gray", "unidentified"="black")) +

  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Order Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 0.4)

df.order.rel.p1


##Family
df.family <- bac.clean.ss %>%
  tax_glom(taxrank = "Family", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.family %<>% mutate(Family = fct_explicit_na(Family, na_level = "unidentified"))


levels(df.family$Family)
levels(df.family$Family) = c(levels(df.family$Family), 'Low abundance')

# we need to group by samples
df.family.rel <- df.family %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.family.rel[df.family.rel$RelAbundance < 5,]$Family <- 'Low abundance'
unique(df.family$Family)

ord <- df.family.rel %>% group_by(Family) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Family
vec.charac<-as.character(vec)
vec.family <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder.fam <- append(vec.Low, vec.family)


df.family.rel$Family <- factor(df.family.rel$Family, levels = vec.reorder.fam) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.family.rel$Family)
## plot relative abundance
#bar plot
df.family.rel.p1 <- ggplot(df.family.rel, aes(x=Sample, y = RelAbundance, fill = Family)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = 
                      c(
"Pseudomonadaceae"= "darkolivegreen",
"Sphingomonadaceae"= "darkolivegreen3"  ,  
"Rhizobiaceae"= "darkorange",
"Caulobacteraceae"= "darkorchid",
"Pseudonocardiaceae"= "darksalmon",    
"Comamonadaceae"= "darkslategray"     ,  
"Microscillaceae"= "darkslategray2",
"Xanthobacteraceae"= "deeppink4"   ,
"Oxalobacteraceae"= "dodgerblue4", 
"Flavobacteriaceae"= "darkseagreen1"  ,  
"Xanthomonadaceae"= "coral4", 
"Micromonosporaceae"= "indianred"    ,
"Rhodanobacteraceae"= "gold4",
"Haliangiaceae"= "aquamarine",
"Streptomycetaceae"= "yellow",     
"Methylophilaceae"= "lightpink4",
"Steroidobacteraceae"= "lightcoral"  ,
"Alcaligenaceae"= "lightpink3"        ,
"Weeksellaceae"= "burlywood3",
"cvE6"= "coral"                 ,
"Kaistiaceae"= "firebrick1",
"Geodermatophilaceae"= "darkseagreen"  ,   
"Hyphomonadaceae"= "lightsteelblue2",
"Corynebacteriaceae "= "lightskyblue3"   ,  
"Nocardiaceae"= "orangered",
"Inquilinaceae"= "plum2"       ,
"Prevotellaceae"= "plum4",
"Devosiaceae"= "palevioletred"        ,  
"Clostridiaceae"= "peachpuff3",
"Propionibacteriaceae"= "pink"  ,
"Beijerinckiaceae"= "palevioletred4",
"Streptosporangiaceae"= "yellow4"            ,  
"Microbacteriaceae"= "steelblue4",
"Burkholderiaceae"= "turquoise"    ,      
"Paenibacillaceae"= "wheat1",
"Nocardioidaceae"= "violet"           ,
"Mycobacteriaceae"= "thistle",
"Dysgonomonadaceae" ="darkred",
"Intrasporangiaceae"= "rosybrown4",        
"Bacillaceae"= "powderblue",
"Low abundance"= "gray", 
"unidentified"="black"  )) +
  
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Family Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 6,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 0.4)

df.family.rel.p1


  ### Genus
df.genus <- bac.clean.ss %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.genus %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))

df.genus$Phylum<-as.character(df.genus$Phylum)
df.genus$Class<-as.character(df.genus$Class)
df.genus$Order<-as.character(df.genus$Order)
df.genus$Family<-as.character(df.genus$Family)
df.genus$Genus<-as.character(df.genus$Genus)

df.genus$Phylum[is.na(df.genus$Phylum)] <- "unidentified"
df.genus$Class[is.na(df.genus$Class)] <- "unidentified"
df.genus$Order[is.na(df.genus$Order)] <- "unidentified"
df.genus$Family[is.na(df.genus$Family)] <- "unidentified"
df.genus$Genus[is.na(df.genus$Genus)] <- "unidentified"

#df.genus$Genus2 <- ifelse(df.genus$Genus == "unidentified",ifelse(df.genus$Family=="unidentified",paste0(df.genus$Order,'_',"Unidentified genus"),paste0(df.genus$Family,'_',"Unidentified genus")),paste0(df.genus$Genus))
#df.genus$Genus2[grep("Prevotella ", df.genus$Genus2)] <- "Prevotella"
#df.genus$Genus2[grep("Ruminococcus ", df.genus$Genus2)] <- "Ruminococcus"
#df.genus$Genus2[grep("Coprococcus ", df.genus$Genus2)] <- "Coprococcus"
#df.genus$Genus2[grep("Ruminiclostridium ", df.genus$Genus2)] <- "Ruminiclostridium"
#df.genus$Genus2[grep("Tyzzerella ", df.genus$Genus2)] <- "Tyzzerella"


levels(df.genus$Genus) = c(levels(df.genus$Genus), 'Low abundance')

# we need to group by samples
df.genus.rel <- df.genus %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.genus.rel[df.genus.rel$RelAbundance <5,]$Genus <- 'Low abundance'
unique(df.genus$Genus)

ord <- df.genus.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.genus <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder <- append(vec.Low, vec.genus)


df.genus.rel$Genus <- factor(df.genus.rel$Genus, levels = vec.reorder) 



## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.genus.rel$Genus)
## plot relative abundance
#bar plot
df.genus.rel.p1 <- ggplot(df.genus.rel, aes(x=Sample, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values =  c(
"Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"= "darkolivegreen3" , 
"Mycobacterium"= "darkorange",                                   
"Amycolatopsis"= "darkred"       ,                              
"Cutibacterium"= "darkorchid"     ,                                
"Mesorhizobium"= "darksalmon"      ,                               
"Paenibacillus"= "darkseagreen4"   ,                                
"Pseudomonas"= "darkslategray"         ,                                
"Caulobacter"= "darkslategray2"     ,                             
"Nocardioides"= "deeppink4"  ,                                 
"Sphingomonas"= "dodgerblue4"   ,                                   
"Flavobacterium"= "darkseagreen1",                                    
"Chryseobacterium"= "coral4"    ,                              
"Sphingobium"= "indianred"      ,                                 
"Pseudoxanthomonas"= "blueviolet",                                  
"Burkholderia-Caballeronia-Paraburkholderia"= "aquamarine"           ,                                  
"Nonomuraea"= "lemonchiffon4" ,                                     
"cvE6"= "coral"                 ,                             
"Steroidobacter"= "lightcoral"  ,                                   
"Aeromicrobium"= "lightpink3"    ,                                   
"Bacillus"= "burlywood3"       ,                                  
"Ohtaekwangia"= "lightpink4"    ,                                  
"Clostridium_sensu_stricto_1"= "firebrick1"  ,                                  
"Lysobacter"= "darkseagreen" ,                                     
"Corynebacterium"= "lightsteelblue2"     ,                                        
"Pseudonocardia"= "lightskyblue3"     ,                                     
"Cloacibacterium"= "orangered"           ,                                 
"Methylobacterium-Methylorubrum"= "plum2"           ,                             
"Proteiniphilum"= "plum4"       ,                             
"Kaistia"= "palevioletred" ,                                     
"Rhodococcus"= "peachpuff3"         ,                                  
"Haliangium"= "pink"             ,                             
"Inquilinus"= "steelblue4"      ,                                  
"Bosea"= "orchid4"                ,                             

"unidentified"= "black",  "Low abundance"= "gray" )) +
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.4)

df.genus.rel.p1


###Species için deneme
df.species <- bac.clean.ss %>%
  tax_glom(taxrank = "Species", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.species %<>% mutate(Species = fct_explicit_na(Species, na_level = "unidentified"))


levels(df.species$Species)
levels(df.species$Species) = c(levels(df.species$Species), 'Low abundance')

# we need to group by samples
df.species.rel <- df.species %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.species.rel[df.species.rel$RelAbundance < 5,]$Species <- 'Low abundance'
unique(df.species$Species)

ord <- df.species.rel %>% group_by(Species) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Species
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac =="Low abundance")]
vec.Low <- c("Low abundance")
vec.reorder <- append(vec.Low, vec.order)


df.species.rel$Species <- factor(df.species.rel$Species, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.species.rel$Species)
## plot relative abundance
#bar plot
df.species.rel.p1 <- ggplot(df.species.rel, aes(x = Sample, y = RelAbundance, fill = Species)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("uncultured_Rhizobiales"="indianred", "uncultured_Flavobacteriia"="plum", "Gammaproteobacteria_bacterium"="darkorchid",   "unidentified"= "black",  "Low abundance"= "gray")) +
  
  xlab('') +
  ylab("Relative abundance (‰)") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+theme(aspect.ratio=0.4)

df.species.rel.p1

#groupped by PlantSpecies
## Phylum level
df.phylum <- bac.clean.ss %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE) %>% 
  psmelt()



df.phylum$Phylum <- as.character(df.phylum$Phylum)
df.phylum$Phylum2 <- df.phylum$Phylum
df.phylum$Phylum2[which(df.phylum$Class=="Alphaproteobacteria")] <- "Alphaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Gammaproteobacteria")] <- "Gammaproteobacteria"
df.phylum$Phylum2[which(df.phylum$Class=="Deltaproteobacteria")] <- "Deltaproteobacteria"

unique(df.phylum$Phylum2)
head(df.phylum)

library(forcats) 
df.phylum %<>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))

df.phylum.2 <- df.phylum %>% mutate(Phylum2 = fct_explicit_na(Phylum2, na_level = "unidentified"))

levels(df.phylum$Phylum2)
levels(df.phylum.2$Phylum2) = c(levels(df.phylum.2$Phylum2), 'Low abundance')

# we need to group by samples
df.phylum.rel <- df.phylum.2 %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.rel[df.phylum.rel$RelAbundance <5,]$Phylum2 <- 'Low abundance'
unique(df.phylum$Phylum2)

ord <- df.phylum.rel %>% group_by(Phylum2) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Phylum2
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac =="Low abundance")]
vec.Low <- c("Low abundance")
vec.reorder <- append(vec.Low, vec.order)


df.phylum.rel$Phylum2 <- factor(df.phylum.rel$Phylum2, levels = vec.reorder) 
table(df.phylum.rel$Phylum2)

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
#bar plot
df.phylum.rel.p1 <- ggplot(df.phylum.rel, aes(x=PlantSpecies, y = RelAbundance, fill = Phylum2)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Actinobacteriota"= "darkolivegreen", "Bacteroidota"= "darkolivegreen3", "Proteobacteria"= "darkorange",   
                               "Firmicutes"= "darkorchid",     
                               "Low abundance"= "gray", "unidentified"="black"  )) +
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Phylum Community Composition by Plant Species \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 2,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=1.8)

df.phylum.rel.p1


##Class level
df.class <- bac.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.class %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))


levels(df.class$Class)
levels(df.class$Class) = c(levels(df.class$Class), 'Low abundance')

# we need to group by samples
df.class.rel <- df.class %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.class.rel[df.class.rel$RelAbundance < 5,]$Class <- 'Low abundance'
unique(df.class$Class)

ord <- df.class.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Class
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac =="Low abundance")]
vec.Low <- c("Low abundance")
vec.reorder <- append(vec.Low, vec.order)


df.class.rel$Class <- factor(df.class.rel$Class, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.class.rel$Class)
## plot relative abundance
#bar plot
df.class.rel.p1 <- ggplot(df.class.rel, aes(x = PlantSpecies, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Alphaproteobacteria"= "darkslategray2", "Gammaproteobacteria"= "indianred", "Bacteroidia"= "darkorange",        
                               "Actinobacteria"= "darkolivegreen",    
                               "Bacilli"= "darkorchid","Low abundance"= "gray", "unidentified"="black")) +
  
  xlab('') +
  ylab("Relative abundance (‰)") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 2,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+theme(aspect.ratio=1.8)

df.class.rel.p1


###Order
df.order <- bac.clean.ss %>%
  tax_glom(taxrank = "Order", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.order %<>% mutate(Order = fct_explicit_na(Order, na_level = "unidentified"))


levels(df.order$Order)
levels(df.order$Order) = c(levels(df.order$Order), 'Low abundance')

# we need to group by samples
df.order.rel <- df.order %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.order.rel[df.order.rel$RelAbundance <5,]$Order <- 'Low abundance'
unique(df.order$Order)

ord <- df.order.rel %>% group_by(Order) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Order
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder <- append(vec.Low, vec.order)


df.order.rel$Order <- factor(df.order.rel$Order, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
table((df.order.rel$Order))
#bar plot
df.order.rel.p1 <- ggplot(df.order.rel, aes(x=PlantSpecies, y = RelAbundance, fill = Order)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Pseudomonadales"= "darkolivegreen3" ,                  
                              "Sphingomonadales"= "darkolivegreen" ,                  
                              "Rhizobiales"= "darkorange"         ,               
                               "Burkholderiales"= "darksalmon"      ,             
                               "Caulobacterales"= "darkorchid"      ,              
                               
                               "Pseudonocardiales"= "darkseagreen4" ,                
                               
                               "Xanthomonadales"= "darkslategray"   ,                 
                                               
                               "Flavobacteriales"= "deeppink4"      ,            
                               "Micromonosporales"= "indianred" ,                
                                          
                               "Enterobacterales"= "dodgerblue4"    ,              
                              "Streptomycetales"= "blueviolet"     ,              
                           
                               "Micrococcales"= "coral"                     , 
                             
                               "Streptosporangiales"= "burlywood3"          ,       
                                    
                               "Propionibacteriales"= "coral4"              ,  
                               "Paenibacillales"= "lightsteelblue2"                ,           
                               
                               "Corynebacteriales"= "paleturquoise2",                 
                               
                               
                               "Bacillales"= "darkseagreen"        ,                
                               "Low abundance"= "gray", "unidentified"="black")) +
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Order Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 1.8)

df.order.rel.p1


##Family
df.family <- bac.clean.ss %>%
  tax_glom(taxrank = "Family", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.family %<>% mutate(Family = fct_explicit_na(Family, na_level = "unidentified"))


levels(df.family$Family)
levels(df.family$Family) = c(levels(df.family$Family), 'Low abundance')

# we need to group by samples
df.family.rel <- df.family %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.family.rel[df.family.rel$RelAbundance < 5,]$Family <- 'Low abundance'
unique(df.family$Family)

ord <- df.family.rel %>% group_by(Family) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Family
vec.charac<-as.character(vec)
vec.family <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder.fam <- append(vec.Low, vec.family)


df.family.rel$Family <- factor(df.family.rel$Family, levels = vec.reorder.fam) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.family.rel$Family)
## plot relative abundance
#bar plot
df.family.rel.p1 <- ggplot(df.family.rel, aes(x=PlantSpecies, y = RelAbundance, fill = Family)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = 
                      c(
                        "Pseudomonadaceae"= "darkolivegreen",
                       "Sphingomonadaceae"= "darkolivegreen3"  ,  
                       "Rhizobiaceae"= "darkorange",
                       
                       "Pseudonocardiaceae"= "darksalmon",    
                       "Comamonadaceae"= "darkslategray"     ,  
                     "Oxalobacteraceae"= "dodgerblue4", 
                     "Flavobacteriaceae"= "darkseagreen1"  ,  
                     "Xanthomonadaceae"= "coral4", 
                     "Micromonosporaceae"= "indianred"    ,
                     "Rhodanobacteraceae"= "gold4",
                      
                "Streptomycetaceae"= "yellow",     
                       
                "Weeksellaceae"= "burlywood3",
                        
                "Propionibacteriaceae"= "pink"  ,
                      
                "Streptosporangiaceae"= "yellow4"            ,  
                "Microbacteriaceae"= "steelblue4",
                         
                "Paenibacillaceae"= "wheat1",
                "Nocardioidaceae"= "violet"           ,
                 "Mycobacteriaceae"= "thistle",
                      
                     "Low abundance"= "gray", 
                        "unidentified"="black"  )) +
  
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Family Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 6,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 1.8)

df.family.rel.p1


### Genus
df.genus <- bac.clean.ss %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.genus %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))

df.genus$Phylum<-as.character(df.genus$Phylum)
df.genus$Class<-as.character(df.genus$Class)
df.genus$Order<-as.character(df.genus$Order)
df.genus$Family<-as.character(df.genus$Family)
df.genus$Genus<-as.character(df.genus$Genus)

df.genus$Phylum[is.na(df.genus$Phylum)] <- "unidentified"
df.genus$Class[is.na(df.genus$Class)] <- "unidentified"
df.genus$Order[is.na(df.genus$Order)] <- "unidentified"
df.genus$Family[is.na(df.genus$Family)] <- "unidentified"
df.genus$Genus[is.na(df.genus$Genus)] <- "unidentified"

#df.genus$Genus2 <- ifelse(df.genus$Genus == "unidentified",ifelse(df.genus$Family=="unidentified",paste0(df.genus$Order,'_',"Unidentified genus"),paste0(df.genus$Family,'_',"Unidentified genus")),paste0(df.genus$Genus))
#df.genus$Genus2[grep("Prevotella ", df.genus$Genus2)] <- "Prevotella"
#df.genus$Genus2[grep("Ruminococcus ", df.genus$Genus2)] <- "Ruminococcus"
#df.genus$Genus2[grep("Coprococcus ", df.genus$Genus2)] <- "Coprococcus"
#df.genus$Genus2[grep("Ruminiclostridium ", df.genus$Genus2)] <- "Ruminiclostridium"
#df.genus$Genus2[grep("Tyzzerella ", df.genus$Genus2)] <- "Tyzzerella"


levels(df.genus$Genus) = c(levels(df.genus$Genus), 'Low abundance')

# we need to group by samples
df.genus.rel <- df.genus %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.genus.rel[df.genus.rel$RelAbundance <5,]$Genus <- 'Low abundance'
unique(df.genus$Genus)

ord <- df.genus.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Genus
vec.charac<-as.character(vec)
vec.genus <- vec.charac[-which(vec.charac %in% c("Low abundance", "unidentified"))]
vec.Low <- c("Low abundance","unidentified")
vec.reorder <- append(vec.Low, vec.genus)


df.genus.rel$Genus <- factor(df.genus.rel$Genus, levels = vec.reorder) 



## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.genus.rel$Genus)
## plot relative abundance
#bar plot
df.genus.rel.p1 <- ggplot(df.genus.rel, aes(x=PlantSpecies, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values =  c(
    "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"= "darkolivegreen3" , 
    "Mycobacterium"= "darkorange",                                   
    "Amycolatopsis"= "darkred"       ,                              
    "Cutibacterium"= "darkorchid"     ,                                
    "Mesorhizobium"= "darksalmon"      ,                               
    "Paenibacillus"= "darkseagreen4"   ,                                
    "Pseudomonas"= "darkslategray"         ,                                
                              
    "Nocardioides"= "deeppink4"  ,                                 
    "Sphingomonas"= "dodgerblue4"   ,                                   
    "Flavobacterium"= "darkseagreen1",                                    
    "Chryseobacterium"= "coral4"    ,                              
                                   
    "Pseudoxanthomonas"= "blueviolet",                                  
                                  
    "Nonomuraea"= "lemonchiffon4" ,                                     
                                 
                                            
    
    "unidentified"= "black",  "Low abundance"= "gray" )) +
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=1.8)

df.genus.rel.p1


###Species için deneme
df.species <- bac.clean.ss %>%
  tax_glom(taxrank = "Species", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.species %<>% mutate(Species = fct_explicit_na(Species, na_level = "unidentified"))


levels(df.species$Species)
levels(df.species$Species) = c(levels(df.species$Species), 'Low abundance')

# we need to group by samples
df.species.rel <- df.species %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.species.rel[df.species.rel$RelAbundance < 5,]$Species <- 'Low abundance'
unique(df.species$Species)

ord <- df.species.rel %>% group_by(Species) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Species
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac =="Low abundance")]
vec.Low <- c("Low abundance")
vec.reorder <- append(vec.Low, vec.order)


df.species.rel$Species <- factor(df.species.rel$Species, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.species.rel$Species)
###Low abundance  unidentified 
######24952           264 
###NO PLOTTTTT


#Excel
#phylum
df.phylumS.rel.tab <- df.phylum.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.phylumP.rel.tab <- df.phylum.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.phylumS.rel.tab2 <- df.phylum.rel %>% group_by(Sample,PlantSpecies,Phylum2) %>% summarise(sumRA = sum(RelAbundance))
df.phylumP.rel.tab2 <- df.phylum.rel %>% group_by(PlantSpecies, Phylum2) %>% summarise(sumRA = sum(RelAbundance))
#Class
df.classS.rel.tab <- df.class.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.classS.rel.tab2 <- df.class.rel %>% group_by(Sample,PlantSpecies, Class) %>% summarise(sumRA = sum(RelAbundance))
df.classP.rel.tab <- df.class.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.classP.rel.tab2 <- df.class.rel %>% group_by(PlantSpecies, Class) %>% summarise(sumRA = sum(RelAbundance))
#Order
df.orderS.rel.tab <- df.order.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.orderS.rel.tab2 <- df.order.rel %>% group_by(Sample,PlantSpecies, Order) %>% summarise(sumRA = sum(RelAbundance))
df.orderP.rel.tab <- df.order.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.orderP.rel.tab2 <- df.order.rel %>% group_by(PlantSpecies, Order) %>% summarise(sumRA = sum(RelAbundance))
#Family
df.familyS.rel.tab <- df.family.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.familyS.rel.tab2 <- df.family.rel %>% group_by(Sample,PlantSpecies, Family) %>% summarise(sumRA = sum(RelAbundance))
df.familyP.rel.tab <- df.family.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.familyP.rel.tab2 <- df.family.rel %>% group_by(PlantSpecies, Family) %>% summarise(sumRA = sum(RelAbundance))
#Genus
df.genusS.rel.tab <- df.genus.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.genusS.rel.tab2 <- df.genus.rel %>% group_by(Sample,PlantSpecies, Genus) %>% summarise(sumRA = sum(RelAbundance))
df.genusP.rel.tab <- df.genus.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.genusP.rel.tab2 <- df.genus.rel %>% group_by(PlantSpecies, Genus) %>% summarise(sumRA = sum(RelAbundance))

#phylum
write.csv(df.phylumS.rel.tab2, "Bacterial RA table Phylum Rel Abun Sample.csv")
write.csv(df.phylumS.rel.tab, "Bacterial RA table Phylum Rel Abun Sample2.csv")
write.csv(df.phylumP.rel.tab2, "Bacterial RA table Phylum Rel Abun PlantSpecies.csv")
write.csv(df.phylumP.rel.tab, "Bacterial RA table Phylum Rel Abun PlantSpecies2.csv")

#Order

write.csv(df.orderS.rel.tab2, "Bacterial RA table Order Rel Abun Sample.csv")
write.csv(df.orderS.rel.tab, "Bacterial RA table Order Rel Abun Sample2.csv")
write.csv(df.orderP.rel.tab2, "Bacterial RA table Order Rel Abun PlantSpecies.csv")
write.csv(df.orderP.rel.tab, "Bacterial RA table Order Rel Abun PlantSpecies2.csv")

#Class

write.csv(df.classS.rel.tab2, "Bacterial RA table Class Rel Abun Sample.csv")
write.csv(df.classS.rel.tab, "Bacterial RA table Class Rel Abun Sample2.csv")
write.csv(df.classP.rel.tab2, "Bacterial RA table Class Rel Abun PlantSpecies.csv")
write.csv(df.classP.rel.tab, "Bacterial RA table Class Rel Abun PlantSpecies2.csv")

#Family
write.csv(df.familyS.rel.tab2, "Bacterial RA table Family Rel Abun Sample.csv")
write.csv(df.familyS.rel.tab, "Bacterial RA table Family Rel Abun Sample2.csv")
write.csv(df.familyP.rel.tab2, "Bacterial RA table Family Rel Abun PlantSpecies.csv")
write.csv(df.familyP.rel.tab, "Bacterial RA table Family Rel Abun PlantSpecies2.csv")

#Genus
write.csv(df.genusS.rel.tab2, "Bacterial RA table Genus Rel Abun Sample.csv")
write.csv(df.genusS.rel.tab, "Bacterial RA table Genus Rel Abun Sample2.csv")
write.csv(df.genusP.rel.tab2, "Bacterial RA table Genus Rel Abun PlantSpecies.csv")
write.csv(df.genusP.rel.tab, "Bacterial RA table Genus Rel Abun PlantSpecies2.csv")





################ Fungal community #######################grouped by Sampless
## Phylum level
##Grouped by diagnosis
df.phylum.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.phylum.fun %<>% mutate(Phylum = fct_explicit_na(Phylum, na_level = "unidentified"))
levels(df.phylum.fun$Phylum) = c(levels(df.phylum.fun$Phylum), 'Low abundance')

# we need to group by samples
df.phylum.fun.rel <- df.phylum.fun %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.fun.rel[df.phylum.fun.rel$RelAbundance < 5,]$Phylum <- 'Low abundance'
unique(df.phylum.fun$Phylum)

ord.f <- df.phylum.fun.rel %>% group_by(Phylum) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Phylum
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.phylum.fun.rel$Phylum <- factor(df.phylum.fun.rel$Phylum, levels = vec.reorder.f) 


## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
table(df.phylum.fun.rel$Phylum)
#bar plot
df.phylum.fun.rel.p1 <- ggplot(df.phylum.fun.rel, aes(x=Sample, y = RelAbundance, fill = Phylum)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Basidiomycota"="pink3","Ascomycota" ="skyblue2", 
                               
                               
                               "Low abundance" = "light grey", "unidentified" ="black")) +
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 2,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.4)
df.phylum.fun.rel.p1



###class
##Grouped by diagnosis
df.class.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.class.fun %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
levels(df.class.fun$Class) = c(levels(df.class.fun$Class), 'Low abundance')

# we need to group by samples
df.class.fun.rel <- df.class.fun %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.class.fun.rel[df.class.fun.rel$RelAbundance < 5,]$Class <- 'Low abundance'
unique(df.class.fun$Class)

ord.f <- df.class.fun.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Class
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.class.fun.rel$Class <- factor(df.class.fun.rel$Class, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.class.fun.rel$Class)
## plot relative abundance
#bar plot
df.class.fun.rel.p1 <- ggplot(df.class.fun.rel, aes(x=Sample, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Agaricomycetes"= "blueviolet",
                               "Pezizomycetes"= "gold4",
                               "Ustilaginomycetes"= "brown",        
                               "Tremellomycetes"= "darkolivegreen",
                               "Malasseziomycetes"="violet", 
                               "Eurotiomycetes"= "lightsteelblue2"     ,                                        
                               "Leotiomycetes"= "darkorange"     ,                                     
                               "Sordariomycetes"= "dodgerblue4"           ,                                 
                               "Dothideomycetes"= "darkslategray2"           ,                             
                               "Low abundance"= "gray", "unidentified"="black")) +
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 2,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.4)
df.class.fun.rel.p1



##Order
##Grouped by diagnosis
df.order.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Order", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.order.fun %<>% mutate(Order = fct_explicit_na(Order, na_level = "unidentified"))
levels(df.order.fun$Order) = c(levels(df.order.fun$Order), 'Low abundance')

# we need to group by samples
df.order.fun.rel <- df.order.fun %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.order.fun.rel[df.order.fun.rel$RelAbundance <5,]$Order <- 'Low abundance'
unique(df.order.fun$Order)

ord.f <- df.order.fun.rel %>% group_by(Order) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Order
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.order.fun.rel$Order <- factor(df.order.fun.rel$Order, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.order.fun.rel$Order)
## plot relative abundance
#bar plot
df.order.fun.rel.p1 <- ggplot(df.order.fun.rel, aes(x=Sample, y = RelAbundance, fill = Order)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Pleosporales"= "darkolivegreen3" ,                  
                               "Helotiales"= "darkolivegreen" ,                  
                               "Hypocreales"= "darkorange"         ,               
                               "Glomerellales"= "darksalmon"      ,             
                               "Chaetothyriales"= "darkorchid"      ,              
                               "Pezizales"= "darkred"      ,          
                               "Sebacinales"= "darkseagreen4" ,                
                               "Eurotiales"= "darkslategray2"  ,                 
                               "Sordariales"= "darkslategray"   ,                 
                               "Thelebolales"= "darkseagreen1"     ,                  
                               "Ustilaginales"= "deeppink4"      ,            
                               "Cystofilobasidiales"= "indianred" ,                
                               "Malasseziales"= "rosybrown4"          ,            
                               "Capnodiales"= "dodgerblue4"    ,              
                               "Myrmecridiales"= "blueviolet"     ,              
                               "Microascales"= "lightpink4"   ,             
                               "Agaricales"= "palevioletred"               ,  
                               "Polyporales"= "aquamarine"                     , 
                               "Auriculariales"= "coral"                     , 
                               "Cantharellales"= "lightcoral"           ,     
                               "Trichosphaeriales"= "gold4"                    , 
                               "Chaetosphaeriales"= "burlywood3"          ,       
                               "Branch06"= "firebrick1"               ,     
                               "Tremellodendropsidales"= "coral4"              ,  
                               "Venturiales"= "lightsteelblue2"                ,           
                               "Filobasidiales"= "paleturquoise2",                 
                               "Trichosporonales"= "sienna1"          ,            
                               "Low abundance"= "gray", "unidentified"="black")) +
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Order Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 4,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.4)
df.order.fun.rel.p1


##Family
df.family.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Family", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.family.fun %<>% mutate(Family = fct_explicit_na(Family, na_level = "unidentified"))
levels(df.family.fun$Family) = c(levels(df.family.fun$Family), 'Low abundance')

# we need to group by samples
df.family.fun.rel <- df.family.fun %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.family.fun.rel[df.family.fun.rel$RelAbundance <5,]$Family <- 'Low abundance'
unique(df.family.fun$Family)

ord.f <- df.family.fun.rel %>% group_by(Family) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Family
vec.f.charac<-as.character(vec.f)
vec.family.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.family.f)


df.family.fun.rel$Family <- factor(df.family.fun.rel$Family, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.family.fun.rel$Family)
## plot relative abundance
#bar plot
df.family.fun.rel.p1 <- ggplot(df.family.fun.rel, aes(x=Sample, y = RelAbundance, fill = Family)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Phaeosphaeriaceae"= "darkolivegreen",
                               "Nectriaceae"= "darkolivegreen3"  ,  
                               "Helotiaceae"= "darkorange",
                               "Rutstroemiaceae"= "darkred"   ,
                               "Plectosphaerellaceae"= "darkorchid",
                               "Herpotrichiellaceae"= "darksalmon",    
                               "Serendipitaceae"= "darkseagreen4",
                               "Trichocomaceae"= "darkslategray"     ,  
                               "Pyronemataceae"= "darkslategray2",
                               "Pleosporaceae"= "deeppink4"   ,
                               "Didymellaceae"= "dodgerblue4", 
                               "Pseudeurotiaceae"= "darkseagreen1"  ,  
                               "Morosphaeriaceae"= "coral4", 
                               "Helotiales_fam_Incertae_sedis"= "indianred"    ,
                               "Chaetomiaceae"= "gold4",
                               "Bionectriaceae"= "blueviolet"      ,   
                               "Pezizaceae"= "aquamarine",
                               "Ustilaginaceae"= "yellow",     
                               "Mrakiaceae"= "lightpink4",
                               "Malasseziaceae"= "lightcoral"  ,
                               "Cladosporiaceae"= "lightpink3"        ,
                               "Myrmecridiaceae"= "burlywood3",
                               "Sporormiaceae"= "coral"                 ,
                               "Microascaceae"= "firebrick1",
                               "Psathyrellaceae"= "darkseagreen"  ,   
                               "Hypocreales_fam_Incertae_sedis"= "lightsteelblue2",
                               "Ceratobasidiaceae "= "lightskyblue3"   ,  
                               "Hypocreaceae"= "orangered",
                               "Auriculariales_fam_Incertae_sedis"= "plum2"       ,
                               "Didymosphaeriaceae"= "plum4",
                               "Melanommataceae"= "palevioletred"        ,  
                               "Phanerochaetaceae"= "peachpuff3",
                               "Aspergillaceae"= "pink"  ,
                               "Leptosphaeriaceae"= "palevioletred4",
                               "Trichosphaeriaceae"= "orange"            ,  
                               "Polyporaceae"= "steelblue4",
                               "Sordariales_fam_Incertae_sedis"= "turquoise"    ,      
                               "Irpicaceae"= "wheat1",
                               "Sympoventuriaceae"= "violet"           ,
                               "Filobasidiaceae"= "thistle",
                              "Exidiaceae"= "magenta",
                              "Tricholomataceae"= "maroon4"       , 
                              "Trichosporonaceae"= "rosybrown", 
                              "Low abundance"= "gray", 
                              "unidentified"="black"  )) +
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Family Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.3)
df.family.fun.rel.p1

##Genus
df.genus.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.genus.fun %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))
df.genus.fun$Phylum<-as.character(df.genus.fun$Phylum)
df.genus.fun$Class<-as.character(df.genus.fun$Class)
df.genus.fun$Order<-as.character(df.genus.fun$Order)
df.genus.fun$Family<-as.character(df.genus.fun$Family)
df.genus.fun$Genus<-as.character(df.genus.fun$Genus)

df.genus.fun$Phylum[is.na(df.genus.fun$Phylum)] <- "unidentified"
df.genus.fun$Class[is.na(df.genus.fun$Class)] <- "unidentified"
df.genus.fun$Order[is.na(df.genus.fun$Order)] <- "unidentified"
df.genus.fun$Family[is.na(df.genus.fun$Family)] <- "unidentified"
df.genus.fun$Genus[is.na(df.genus.fun$Genus)] <- "unidentified"

#df.genus.fun$Genus2 <- ifelse(df.genus.fun$Genus == "unidentified",ifelse(df.genus.fun$Family=="unidentified",paste0(df.genus.fun$Order,'_',"Unidentified genus"),paste0(df.genus.fun$Family,'_',"Unidentified genus")),paste0(df.genus.fun$Genus))


levels(df.genus.fun$Genus) = c(levels(df.genus.fun$Genus), 'Low abundance')

# we need to group by samples
df.genus.fun.rel <- df.genus.fun %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.genus.fun.rel[df.genus.fun.rel$RelAbundance <5,]$Genus <- 'Low abundance'
unique(df.genus.fun$Genus)

ord.f <- df.genus.fun.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Genus
vec.f.charac<-as.character(vec.f)
vec.genus.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified_Unidentified genus"))]
vec.Low.f <- c("Low abundance","unidentified_Unidentified genus")
vec.reorder.f <- append(vec.Low.f, vec.genus.f)


df.genus.fun.rel$Genus <- factor(df.genus.fun.rel$Genus, levels = vec.reorder.f) 





## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.genus.fun.rel$Genus)
## plot relative abundance
#bar plot
df.genus.fun.rel.p1 <- ggplot(df.genus.fun.rel, aes(x=Sample, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c(
"Paraphoma"= "darkolivegreen",
"Fusarium"= "darkolivegreen3"  ,  
"Tetracladium"= "darkorange",
"Lambertella"= "darkred"   ,
"Plectosphaerella"= "darkorchid",
"Exophiala"= "darksalmon",    
"Talaromyces"= "darkseagreen4",
"Lasiobolidium"= "darkslategray"       ,
"Alternaria"= "darkslategray2",
"Neonectria"= "deeppink4"   ,
"Phaeosphaeria"= "dodgerblue4",
"Acrocalymma"= "darkseagreen1"  ,  
"Pseudogymnoascus"= "coral4",
"Cadophora"= "cornsilk4"    ,
"Clonostachys"= "gold4",
"Moesziomyces"= "blueviolet"        , 
"Malassezia"= "aquamarine",
"Peziza"= "lemonchiffon4",     
"Zopfiella"= "lightpink4",
"Cladosporium"= "lightcoral"  ,
"Myrmecridium"= "snow1", 
"Tausonia"= "lightpink3"        ,
"Preussia"= "burlywood3", 
"Mrakia"= "coral"                 ,
"Cephalotrichum"= "firebrick1",
"Verticillium"= "darkseagreen"    , 
"Psathyrella"= "lightsteelblue2",
"Trichocladium"= "lightskyblue3"     ,
"Trichoderma"= "orangered",
"Thanatephorus"= "plum2"       ,
"Oliveonia"= "plum4",
"Toxicocladosporium"= "palevioletred"        ,  
"Epicoccum"= "peachpuff3",
"Setophoma"= "pink"  ,
"Acremonium"= "palevioletred4",
"Didymosphaeria"= "yellow4"            ,  
"Alpinaria"= "steelblue4",
"Phanerochaete"= "turquoise"   ,       
"Sarocladium"= "wheat1", 
"Setophaeosphaeria"= "violet"           ,
"Leptosphaeria"= "violetred2", 
"Nigrospora"= "yellowgreen" ,     
"Conlarium"= "thistle", 
"Aspergillus"= "thistle4" , 
"Hexagonia"= "magenta", 
"Botryotrichum"= "maroon4"        ,
"Mariannaea"= "rosybrown",
"Scolecobasidium"= "orchid2"   ,    
"Ilyonectria"= "rosybrown2",
"unidentified" ="black",
"Low abundance"="gray"
))+
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 6,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=0.4)
df.genus.fun.rel.p1



################ Fungal community #######################grouped by PlantSpecies
## Phylum level
##Grouped by diagnosis
df.phylum.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Phylum", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.phylum.fun %<>% mutate(Phylum = fct_explicit_na(Phylum, na_level = "unidentified"))
levels(df.phylum.fun$Phylum) = c(levels(df.phylum.fun$Phylum), 'Low abundance')

# we need to group by samples
df.phylum.fun.rel <- df.phylum.fun %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.phylum.fun.rel[df.phylum.fun.rel$RelAbundance < 5,]$Phylum <- 'Low abundance'
unique(df.phylum.fun$Phylum)

ord.f <- df.phylum.fun.rel %>% group_by(Phylum) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Phylum
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.phylum.fun.rel$Phylum <- factor(df.phylum.fun.rel$Phylum, levels = vec.reorder.f) 



## relative abundance with less than 0.5% in sample was labeled 'low abundance'

## plot relative abundance
table(df.phylum.fun.rel$Phylum)
#bar plot
df.phylum.fun.rel.p1 <- ggplot(df.phylum.fun.rel, aes(x=PlantSpecies, y = RelAbundance, fill = Phylum)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Basidiomycota"="pink3","Ascomycota" ="skyblue2", 
                              
                               
                               "Low abundance" = "gray", "unidentified" ="black")) +
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Phylum Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 2,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=1.8)
df.phylum.fun.rel.p1



###class
##Grouped by diagnosis
df.class.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Class", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.class.fun %<>% mutate(Class = fct_explicit_na(Class, na_level = "unidentified"))
levels(df.class.fun$Class) = c(levels(df.class.fun$Class), 'Low abundance')

# we need to group by samples
df.class.fun.rel <- df.class.fun %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.class.fun.rel[df.class.fun.rel$RelAbundance < 5,]$Class <- 'Low abundance'
unique(df.class.fun$Class)

ord.f <- df.class.fun.rel %>% group_by(Class) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Class
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.class.fun.rel$Class <- factor(df.class.fun.rel$Class, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.class.fun.rel$Class)
## plot relative abundance
#bar plot
df.class.fun.rel.p1 <- ggplot(df.class.fun.rel, aes(x=PlantSpecies, y = RelAbundance, fill = Class)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Agaricomycetes"= "blueviolet",
                               "Pezizomycetes"= "gold4",
                                    
                               "Tremellomycetes"= "darkolivegreen",
                               
                                 
                               
                               
                               
                               "Eurotiomycetes"= "lightsteelblue2"     ,                                        
                               "Leotiomycetes"= "darkorange"     ,                                     
                               "Sordariomycetes"= "dodgerblue4"           ,                                 
                               "Dothideomycetes"= "darkslategray2"           ,                             
                                 
                               
                               "Low abundance"= "gray", "unidentified"="black")) +
  
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=1.8)
df.class.fun.rel.p1



##Order
##Grouped by diagnosis
df.order.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Order", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.order.fun %<>% mutate(Order = fct_explicit_na(Order, na_level = "unidentified"))
levels(df.order.fun$Order) = c(levels(df.order.fun$Order), 'Low abundance')

# we need to group by samples
df.order.fun.rel <- df.order.fun %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.order.fun.rel[df.order.fun.rel$RelAbundance <5,]$Order <- 'Low abundance'
unique(df.order.fun$Order)

ord.f <- df.order.fun.rel %>% group_by(Order) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Order
vec.f.charac<-as.character(vec.f)
vec.order.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.order.f)


df.order.fun.rel$Order <- factor(df.order.fun.rel$Order, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.order.fun.rel$Order)
## plot relative abundance
#bar plot
df.order.fun.rel.p1 <- ggplot(df.order.fun.rel, aes(x=PlantSpecies, y = RelAbundance, fill = Order)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Pleosporales"= "darkolivegreen3" ,                  
                               "Helotiales"= "darkolivegreen" ,                  
                               "Hypocreales"= "darkorange"         ,               
                               "Glomerellales"= "darksalmon"      ,             
                               "Chaetothyriales"= "darkorchid"      ,              
                               "Pezizales"= "darkred"      ,          
                               "Sebacinales"= "darkseagreen4" ,                
                               "Eurotiales"= "darkslategray2"  ,                 
                               "Sordariales"= "darkslategray"   ,                 
                               "Thelebolales"= "darkseagreen1"     ,                  
                               "Cystofilobasidiales"= "indianred" ,                
                               "Agaricales"= "palevioletred"               ,  
                               "Cantharellales"= "lightcoral"           ,     
                               "Low abundance"= "gray", "unidentified"="black")) +
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Order Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=1.8)
df.order.fun.rel.p1


##Family
df.family.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Family", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.family.fun %<>% mutate(Family = fct_explicit_na(Family, na_level = "unidentified"))
levels(df.family.fun$Family) = c(levels(df.family.fun$Family), 'Low abundance')

# we need to group by samples
df.family.fun.rel <- df.family.fun %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.family.fun.rel[df.family.fun.rel$RelAbundance <5,]$Family <- 'Low abundance'
unique(df.family.fun$Family)

ord.f <- df.family.fun.rel %>% group_by(Family) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Family
vec.f.charac<-as.character(vec.f)
vec.family.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified"))]
vec.Low.f <- c("Low abundance","unidentified")
vec.reorder.f <- append(vec.Low.f, vec.family.f)


df.family.fun.rel$Family <- factor(df.family.fun.rel$Family, levels = vec.reorder.f) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.family.fun.rel$Family)
## plot relative abundance
#bar plot
df.family.fun.rel.p1 <- ggplot(df.family.fun.rel, aes(x=PlantSpecies, y = RelAbundance, fill = Family)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Phaeosphaeriaceae"= "darkolivegreen",
                              "Nectriaceae"= "darkolivegreen3"  ,  
                              "Helotiaceae"= "darkorange",
                              "Rutstroemiaceae"= "darkred"   ,
                              "Plectosphaerellaceae"= "darkorchid",
                              "Herpotrichiellaceae"= "darksalmon",    
                              "Serendipitaceae"= "darkseagreen4",
                              "Trichocomaceae"= "darkslategray"     ,  
                              "Pyronemataceae"= "darkslategray2",
                              "Pleosporaceae"= "deeppink4"   ,
                              "Didymellaceae"= "dodgerblue4", 
                              "Pseudeurotiaceae"= "darkseagreen1"  ,  
                              "Morosphaeriaceae"= "coral4", 
                              "Helotiales_fam_Incertae_sedis"= "indianred"    ,
                              "Chaetomiaceae"= "gold4",
                              "Bionectriaceae"= "blueviolet"      ,   
                              "Pezizaceae"= "aquamarine",
                              "Mrakiaceae"= "lightpink4",
                              "Sporormiaceae"= "coral"                 ,
                              "Psathyrellaceae"= "darkseagreen"  ,   
                              "Ceratobasidiaceae"= "lightskyblue3"   ,  
                              "Low abundance"= "gray", 
                               "unidentified"="black"  )) +
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Family Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 8,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=1.8)
df.family.fun.rel.p1

##Genus
df.genus.fun <- fun.clean.ss %>%
  tax_glom(taxrank = "Genus", NArm = FALSE) %>% 
  psmelt()


library(forcats) 
df.genus.fun %<>% mutate(Genus = fct_explicit_na(Genus, na_level = "unidentified"))
df.genus.fun$Phylum<-as.character(df.genus.fun$Phylum)
df.genus.fun$Class<-as.character(df.genus.fun$Class)
df.genus.fun$Order<-as.character(df.genus.fun$Order)
df.genus.fun$Family<-as.character(df.genus.fun$Family)
df.genus.fun$Genus<-as.character(df.genus.fun$Genus)

df.genus.fun$Phylum[is.na(df.genus.fun$Phylum)] <- "unidentified"
df.genus.fun$Class[is.na(df.genus.fun$Class)] <- "unidentified"
df.genus.fun$Order[is.na(df.genus.fun$Order)] <- "unidentified"
df.genus.fun$Family[is.na(df.genus.fun$Family)] <- "unidentified"
df.genus.fun$Genus[is.na(df.genus.fun$Genus)] <- "unidentified"

#df.genus.fun$Genus2 <- ifelse(df.genus.fun$Genus == "unidentified",ifelse(df.genus.fun$Family=="unidentified",paste0(df.genus.fun$Order,'_',"Unidentified genus"),paste0(df.genus.fun$Family,'_',"Unidentified genus")),paste0(df.genus.fun$Genus))


levels(df.genus.fun$Genus) = c(levels(df.genus.fun$Genus), 'Low abundance')

# we need to group by samples
df.genus.fun.rel <- df.genus.fun %>%  
  group_by(PlantSpecies) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.genus.fun.rel[df.genus.fun.rel$RelAbundance <5,]$Genus <- 'Low abundance'
unique(df.genus.fun$Genus)

ord.f <- df.genus.fun.rel %>% group_by(Genus) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec.f <- ord.f$Genus
vec.f.charac<-as.character(vec.f)
vec.genus.f <- vec.f.charac[-which(vec.f.charac %in% c("Low abundance","unidentified_Unidentified genus"))]
vec.Low.f <- c("Low abundance","unidentified_Unidentified genus")
vec.reorder.f <- append(vec.Low.f, vec.genus.f)


df.genus.fun.rel$Genus <- factor(df.genus.fun.rel$Genus, levels = vec.reorder.f) 





## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.genus.fun.rel$Genus)
## plot relative abundance
#bar plot
df.genus.fun.rel.p1 <- ggplot(df.genus.fun.rel, aes(x=PlantSpecies, y = RelAbundance, fill = Genus)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = c("Paraphoma"= "darkolivegreen",
                               "Fusarium"= "darkolivegreen3"  ,  
                               "Tetracladium"= "darkorange",
                               "Lambertella"= "darkred"   ,
                               "Plectosphaerella"= "darkorchid",
                               "Exophiala"= "darksalmon",    
                               "Talaromyces"= "darkseagreen4",
                               "Lasiobolidium"= "darkslategray"       ,
                               "Alternaria"= "darkslategray2",
                               "Neonectria"= "deeppink4"   ,
                               "Phaeosphaeria"= "dodgerblue4",
                               "Acrocalymma"= "darkseagreen1"  ,  
                               "Pseudogymnoascus"= "coral4",
                               "Cadophora"= "indianred"    ,
                               "Clonostachys"= "gold4",
                             "Peziza"= "red",     
                             "Zopfiella"= "lightpink4",
                             "Tausonia"= "lightpink3"        ,
                             "Preussia"= "burlywood3", 
                             "Mrakia"= "yellow"                 ,
                             "Psathyrella"= "lightsteelblue2",
                             "Thanatephorus"= "plum2"       ,
                             "unidentified" ="black",
                             "Low abundance"="gray"
  ))+
  xlab('')+
  ylab("Relative abundance (‰)") +
  #ggtitle("Genus Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 8,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio=1.8)
df.genus.fun.rel.p1
#Excel
#phylum
df.phylumS.rel.tab <- df.phylum.fun.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.phylumS.rel.tab2 <- df.phylum.fun.rel %>% group_by(Sample,PlantSpecies, Phylum) %>% summarise(sumRA = sum(RelAbundance))
df.phylumP.rel.tab <- df.phylum.fun.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.phylumP.rel.tab2 <- df.phylum.fun.rel %>% group_by(PlantSpecies, Phylum) %>% summarise(sumRA = sum(RelAbundance))
#Class
df.classS.rel.tab <- df.class.fun.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.classS.rel.tab2 <- df.class.fun.rel %>% group_by(Sample,PlantSpecies, Class) %>% summarise(sumRA = sum(RelAbundance))
df.classP.rel.tab <- df.class.fun.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.classP.rel.tab2 <- df.class.fun.rel %>% group_by(PlantSpecies, Class) %>% summarise(sumRA = sum(RelAbundance))
#Order
df.orderS.rel.tab <- df.order.fun.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.orderS.rel.tab2 <- df.order.fun.rel %>% group_by(Sample,PlantSpecies, Order) %>% summarise(sumRA = sum(RelAbundance))
df.orderP.rel.tab <- df.order.fun.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.orderP.rel.tab2 <- df.order.fun.rel %>% group_by(PlantSpecies, Order) %>% summarise(sumRA = sum(RelAbundance))
#Family
df.familyS.rel.tab <- df.family.fun.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.familyS.rel.tab2 <- df.family.fun.rel %>% group_by(Sample, PlantSpecies,Family) %>% summarise(sumRA = sum(RelAbundance))
df.familyP.rel.tab <- df.family.fun.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.familyP.rel.tab2 <- df.family.fun.rel %>% group_by(PlantSpecies, Family) %>% summarise(sumRA = sum(RelAbundance))
#Genus
df.genusS.rel.tab <- df.genus.fun.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.genusS.rel.tab2 <- df.genus.fun.rel %>% group_by(Sample,PlantSpecies, Genus) %>% summarise(sumRA = sum(RelAbundance))
df.genusP.rel.tab <- df.genus.fun.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.genusP.rel.tab2 <- df.genus.fun.rel %>% group_by(PlantSpecies, Genus) %>% summarise(sumRA = sum(RelAbundance))

#phylum
write.csv(df.phylumS.rel.tab2, "Fungi RA table Phylum Rel Abun Sample.csv")
write.csv(df.phylumS.rel.tab, "Fungi RA table Phylum Rel Abun Sample2.csv")
write.csv(df.phylumP.rel.tab2, "Fungi RA table Phylum Rel Abun PlantSpecies.csv")
write.csv(df.phylumP.rel.tab, "Fungi RA table Phylum Rel Abun PlantSpecies2.csv")

#Order

write.csv(df.orderS.rel.tab2, "Fungi RA table Order Rel Abun Sample.csv")
write.csv(df.orderS.rel.tab, "Fungi RA table Order Rel Abun Sample2.csv")
write.csv(df.orderP.rel.tab2, "Fungi RA table Order Rel Abun PlantSpecies.csv")
write.csv(df.orderP.rel.tab, "Fungi RA table Order Rel Abun PlantSpecies2.csv")

#Class

write.csv(df.classS.rel.tab2, "Fungi RA table Class Rel Abun Sample.csv")
write.csv(df.classS.rel.tab, "Fungi RA table Class Rel Abun Sample2.csv")
write.csv(df.classP.rel.tab2, "Fungi RA table Class Rel Abun PlantSpecies.csv")
write.csv(df.classP.rel.tab, "Fungi RA table Class Rel Abun PlantSpecies2.csv")

#Family
write.csv(df.familyS.rel.tab2, "Fungi RA table Family Rel Abun Sample.csv")
write.csv(df.familyS.rel.tab, "Fungi RA table Family Rel Abun Sample2.csv")
write.csv(df.familyP.rel.tab2, "Fungi RA table Family Rel Abun PlantSpecies.csv")
write.csv(df.familyP.rel.tab, "Fungi RA table Family Rel Abun PlantSpecies2.csv")

#Genus
write.csv(df.genusS.rel.tab2, "Fungi RA table Genus Rel Abun Sample.csv")
write.csv(df.genusS.rel.tab, "Fungi RA table Genus Rel Abun Sample2.csv")
write.csv(df.genusP.rel.tab2, "Fungi RA table Genus Rel Abun PlantSpecies.csv")
write.csv(df.genusP.rel.tab, "Fungi RA table Genus Rel Abun PlantSpecies2.csv")



#Species için deneme
df.species <- bac.clean.ss %>%
  tax_glom(taxrank = "Species", NArm = FALSE) %>% 
  psmelt()

library(forcats) 
df.species %<>% mutate(Species = fct_explicit_na(Species, na_level = "unidentified"))


levels(df.species$Species)
levels(df.species$Species) = c(levels(df.species$Species), 'Low abundance')

# we need to group by samples
df.species.rel <- df.species %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.species.rel[df.species.rel$RelAbundance < 5,]$Species <- 'Low abundance'
unique(df.species$Species)

ord <- df.species.rel %>% group_by(Species) %>% summarise(Abundance = sum(Abundance)) %>% arrange(Abundance)
vec <- ord$Species
vec.charac<-as.character(vec)
vec.order <- vec.charac[-which(vec.charac =="Low abundance")]
vec.Low <- c("Low abundance")
vec.reorder <- append(vec.Low, vec.order)


df.species.rel$Species <- factor(df.species.rel$Species, levels = vec.reorder) 

## relative abundance with less than 0.5% in sample was labeled 'low abundance'
table(df.species.rel$Species)
## plot relative abundance
#bar plot
df.species.rel.p1 <- ggplot(df.species.rel, aes(x = Sample, y = RelAbundance, fill = Species)) + 
  geom_bar(stat="identity", width = 0.8, position = 'stack') +
  #scale_fill_discrete() +
  scale_fill_manual(values = my_color_Class) +
  
  xlab('') +
  ylab("Relative abundance (‰)") +
  #ggtitle("Class Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 3,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(breaks=seq(0,1000,100))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+theme(aspect.ratio=0.4)

df.species.rel.p1


#Genus
df.speciesS.rel.tab <- df.species.rel %>% group_by(Sample) %>% summarise(sumRA = sum(RelAbundance))
df.speciesS.rel.tab2 <- df.species.rel %>% group_by(Sample,PlantSpecies, Species) %>% summarise(sumRA = sum(RelAbundance))
df.speciesP.rel.tab <- df.species.rel %>% group_by(PlantSpecies) %>% summarise(sumRA = sum(RelAbundance))
df.speciesP.rel.tab2 <- df.species.rel %>% group_by(PlantSpecies, Species) %>% summarise(sumRA = sum(RelAbundance))
#Genus
write.csv(df.speciesS.rel.tab2, "Bacteria RA table Species Rel Abun Sample.csv")
write.csv(df.speciesS.rel.tab, "Bacteria RA table Species Rel Abun Sample2.csv")
write.csv(df.speciesP.rel.tab2, "Bacteria RA table Species Rel Abun PlantSpecies.csv")
write.csv(df.speciesP.rel.tab, "Bacteria RA table Species Rel Abun PlantSpecies2.csv")
