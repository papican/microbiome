# Prepearing files for qiime2 environment
  #Biom file import
  qiime tools import \
   --input-path asv_table_final.biom \
   --type 'FeatureTable[Frequency]' \
   --input-format BIOMV210Format \
   --output-path table.qza

  #tree file import
  qiime tools import \
  --input-path tree.nwk \
  --output-path tree.qza \
  --type 'Phylogeny[Unrooted]'

  #fasta file import
  qiime tools import \
  --input-path dna-sequences.fna \
  --output-path sequences.qza \
  --type 'FeatureData[Sequence]'


  #mapped taxa file
  qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-nb-classifier.qza \
  --i-reads sequences.qza \
  --o-classification taxonomy.qza



#core taxa with qiime2
qiime feature-table core-features \
--i-table table.qza \
--p-min-fraction 0.4 \
--p-max-fraction 1 \
--p-steps 12 \
--o-visualization cor_viz.qzv



#Creating phyloseq object from biom files artifacts
library(dplyr)
library(forcats)
library(metagenomeSeq)
library(vegan)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
#library(seqtime)
library(agricolae)
library(RColorBrewer)
library(xlsx)
library(magrittr)
library(indicspecies)
library(Hmisc)
library(igraph)
library(qgraph)
library(randomForest)
#library(multifunc)
library(OTUtable)
library(FSA)
library(rcompanion)



phylum_colors <- c(
  "gray",'black', "#DA5724", "#5F7FC7","#508578", "#CD9BCD", "orange",
  "#5F7FC7","#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#D1A33D", "#8A7C64", "#599861","#5E738F"
)



my_color_collection <- c(
  "#CBD588", "#5F7FC7", "orange", "#AD6F3B", "#673770",
  "#D14285", "#652926", "#C84248", "#8569D5", "#5E738F",
  "#D1A33D", "#8A7C64", "#599861","#616163", "#FFCDB2",
  "#6D9F71", "#242F40",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")

my_color_Class <- c(
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "orange", "#5F7FC7", "#CBD588", "#AD6F3B", "#673770",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',

  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")



my_color_OTU <- c(
  "gray",'black', "#5F7FC7", "orange",  "#AD6F3B",
  "#673770","#D14285", "#652926", "#C84248",  "#8569D5",
  "#5E738F","#D1A33D", "#8A7C64", "#599861","#616163",
  "#FFCDB2", "#242F40", "#6D9F71",  "#CCA43B", "#F92A82",
  "#ED7B84", "#7EB77F", "#DEC4A1", "#E5D1D0", '#0E8482',
  '#C9DAEA', '#337357', '#95C623', '#E55812', '#04471C',
  '#F2D7EE', '#D3BCC0', '#A5668B', '#69306D', 'navy',
  '#1A535C', '#4ECDC4', 'orange', '#FF6B6B', "orchid1",
  'cyan2', '#FFF275', 'springgreen', '#FF3C38', '#A23E48',
  '#000000', '#CF5C36', '#EEE5E9', '#7C7C7C', '#EFC88B',

  '#2E5266', '#6E8898', '#9FB1BC', '#D3D0CB', '#E2C044',
  '#5BC0EB', '#FDE74C', '#9BC53D', '#E55934', '#FA7921',
  "#CD9BCD", "#508578", "#CBD588","#CBD588", "#5F7FC7",
  "orange",   "#AD6F3B", "#673770","#D14285", "#652926",
  "#C84248",  "#8569D5", "#5E738F","#D1A33D", "#8A7C64",
  "#599861","#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", "#DEC4A1",
  "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', '#95C623',
  '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', '#A5668B',
  '#69306D', '#0E103D', '#1A535C', '#4ECDC4', '#F7FFF7',
  '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange")



















# Set plotting theme
theme_set(theme_bw())


### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
bac_phylo=import_biom("asv_table_final.biom")

### merge with metadata
# Import sample metadata
## maybe Gyeryueng data didn't work because it wasn't transformed to json file

## in metadata erase # (This step is essential)
map <- read.table(file = 'sample_metadata.tsv', sep = '\t', header = TRUE)
map <- sample_data(map)


head(map)
dim(map)
summary(map)
str(map)

summary(map)
colnames(map)
rownames(map)
nrow(map)

# Assign rownames to be Sample ID's
map$SampleID
rownames(map) <- map$SampleID
rownames(map)
dim(map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
phy_tree = read_tree("tree.nwk")
phy <- merge_phyloseq(bac_phylo, map, phy_tree)

class(phy)
phy


## changing rank names
colnames(tax_table(phy)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phy

#### get unassigned vectors
#### get CP and MT phyloseq obj and vectors
# (1) CP
phy.cp <- subset_taxa(phy, Order == "o__Chloroplast")
phy.cp <- subset_taxa(phy, Family == "f__Chloroplast")
phy.cp <- subset_taxa(phy, Genus == "g__Chloroplast")

vec.cp <- rownames(otu_table(phy.cp))
length(rownames(otu_table(phy.cp))) ## 10 otus of CP
vec.cp


# (2) MT
phy.mt <- subset_taxa(phy, Family == "f__Mitochondria")

vec.mt <- rownames(otu_table(phy.mt))
tax_table(phy.mt)
length(rownames(otu_table(phy.mt))) ## -> 110 



# (3) Unassigned
#unique(tax_table(phy)[,'Kingdom']) ## only bacteria, then no need to exclude
#phy.un <- subset_taxa(phy, Kingdom == "Unassigned")
#vec.un <- rownames(otu_table(phy.un))
#tax_table(phy.un)
#length(rownames(otu_table(phy.un)))


# (4) Archea
#unique(tax_table(phy)[,'Kingdom']) ## only bacteria, then no need to exclude
#phy.ar <- subset_taxa(phy, Kingdom == "d__Archaea")
#vec.ar <- rownames(otu_table(phy.ar))
#tax_table(phy.ar)
#length(rownames(otu_table(phy.ar)))


### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### clean taxa
phy.clean <- pop_taxa(phy, c(vec.cp, vec.mt))

phy  ## 3367
phy.clean ## 3247


#### We will also remove the "d__" patterns for cleaner labels
tax_table(phy.clean)[,colnames(tax_table(phy.clean))] <- gsub(tax_table(phy.clean)[,colnames(tax_table(phy.clean))],pattern="[a-z]__",replacement="")



### filter otu with total count of 20? (in all samples)
### later we can implement

phy.clean.otu <- otu_table(phy.clean)
head(phy.clean.otu)
df.clean.otu <- data.frame(phy.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')





## Remove reads with over 400 bp
library(seqinr)

bac.seq <- read.fasta(file = "dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")
## max length and min length
getLength(bac.seq)
min(getLength(bac.seq)) #278
max(getLength(bac.seq)) #456
bac.seq[which(getLength(bac.seq)>390)]

otu_over_400bp <- attr(bac.seq[which(getLength(bac.seq)>400)], "names")

otu_over_400bp
phy.clean
phy.clean <- pop_taxa(phy.clean,otu_over_400bp)
phy.clean  ## 3246
sum(otu_table(phy.clean)) # 3315727

summarize_phyloseq(phy.clean)
phy.clean



#last time clean bacteria from archea

bac.clean.ss <- subset_taxa(phy.clean, Kingdom != "Archaea")
bac.clean.ss <- phyloseq::filter_taxa(bac.clean.ss, function(x) sum(x) != 0, TRUE)

# (4)#unidentified issue
#tax.bac<-tax_table(bac.clean.ss)
#tax.bac <-data.frame(tax.bac)
#unidentified.bac.phylum<-rownames(tax.bac)[is.na(tax.bac$Phylum)]

#unidentified.bac.phylum.seq<-bac.seq[which(attr(bac.seq, "names")%in% unidentified.bac.phylum)]
#write.fasta(unidentified.bac.phylum.seq, names(unidentified.bac.phylum.seq), file.out = "unidentified.bac.phylum.seq.fasta")



#bac.clean.ss <- pop_taxa(bac.clean.ss, remove.bacterial.otu)
#summarize_phyloseq(bac.clean.ss)
#bac.clean.ss #3246!!!!!!!!



## Designating OTU id
bac.list <- bac.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

#########################Fungal community
###Fungal community
## Fungal community

### setting input and output path
### We can then load the biom file with phyloseq function import_biom. We extract the OTU table with OTU abundances and the taxonomy table from the resulting phyloseq object.
fun_phylo=import_biom("asv_table_final.biom")


### merge with metadata
# Import sample metadata
## maybe Gyeryueng data didn't work because it wasn't transformed to json file

## in metadata erase # (This step is essential)
f.map <- read.table(file = 'sample_metadata.tsv', sep = '\t', header = TRUE)
f.map <- sample_data(f.map)

head(f.map)
dim(f.map)
summary(f.map)
str(f.map)

summary(f.map)
colnames(f.map)
rownames(f.map)
nrow(f.map)

# Assign rownames to be Sample ID's
f.map$SampleID
rownames(f.map) <- f.map$SampleID
rownames(f.map)
dim(f.map)
# Merge biom data object with sample metadata + tree data(this is rooted!)
fun_tree = read_tree("tree.nwk")
fun <- merge_phyloseq(fun_phylo, f.map, fun_tree)

class(fun)
fun   ## 956 otus 

## changing rank names
colnames(tax_table(fun)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

fun


#### get unassigned vectors

# (3) Unassigned
#unique(tax_table(fun)[,'Kingdom']) ##ONLY FUNGİ SO NO NEED TO EXCLUDE "k__Chromista", "k__Plantae", "Unassigned"
#tax_table(fun)[,'Kingdom'
#fun.un <- subset_taxa(fun, Kingdom %in% c("Unassigned","k__Chromista","k__Plantae"))
#vec.un <- rownames(otu_table(fun.un))
#tax_table(fun.un)
#length(rownames(otu_table(fun.un))) ##  43

### exclude those vectors
fun  

sample_names(fun)
sample_variables(fun)

### pop taxa application
## get rid of CP and MT otus
## pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### Clean it up!!!
#fun.clean <- pop_taxa(fun, c(vec.un))

fun.clean <- fun #135 samples 1598 OTUs
fun.clean
#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(fun.two)[,colnames(tax_table(fun.two))] <- gsub(tax_table(fun.two)[,colnames(tax_table(fun.two))],pattern="[A-Z]_[0-9]__",replacement="")
# fun.test4 <- fun.two %>% psmelt() 
# fun.test4


tax_table(fun.clean)[,colnames(tax_table(fun.clean))] <- gsub(tax_table(fun.clean)[,colnames(tax_table(fun.clean))],pattern="[a-z]__",replacement="")

#' sample_data(fun.clean)$SampleID <- factor(sample_data(fun.clean)$SampleID, levels =target_PAB)

tax_table(fun.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

fun.clean #956 taxa
str(fun.clean)
otu_table(fun.clean)

# ## fix it in fun.clean object!!! pop_taxa does the work
# fun.clean <- pop_taxa(fun.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(fun.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement 
fun.clean.otu <- otu_table(fun.clean)
head(fun.clean.otu)
df.clean.otu <- data.frame(fun.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


# # how about we delete major taxa in the negative sequence?
# negative <- df.clean.otu[,c("OTU",'F204','F205','F206')]
# negative
# negative$total <- apply(negative[,-1], 1, sum)
# negative_0 <- subset(negative, negative$total > 0) 
# neg.otu <- negative_0$OTU  ### 58 otu to be eliminated
# neg.otu
# length(neg.otu)
# tax_table(fun.clean)[neg.otu]
# library(xlsx)
# write.xlsx(tax_table(fun.clean)[neg.otu], 'fungal neg.otu_taxonomy_old DB.xlsx')
# 

sum(otu_table(fun)) ##3624315
sum(otu_table(fun.clean)) ## 3624315



##### get rid of otu of less than 100 reads
library(seqinr)
fun.seq <- read.fasta(file = "dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")
min(getLength(fun.seq)) #163 
max(getLength(fun.seq)) #387
otu_less_than_100bp <- attr(fun.seq[which(getLength(fun.seq)<100)], "names")

fun.clean.ss <- pop_taxa(fun.clean,otu_less_than_100bp)
fun.clean.ss ## 956## same because no seq less then 100bp
sum(otu_table(fun.clean.ss)) ## 3624315

## Fungal OTUs which are 'unidentifed' at the phylum level
#tax.fun<-tax_table(fun.clean.ss)
#tax.fun <-data.frame(tax.fun)
#unidentified.phylum<-rownames(tax.fun)[is.na(tax.fun$Phylum)]

#unidentified.phylum.seq<-fun.seq[which(attr(fun.seq, "names")%in% unidentified.phylum)]
#write.fasta(unidentified.phylum.seq, names(unidentified.phylum.seq), file.out = "unidentified.phylum.seq.fasta")

## Remove unidentified OTUs based on BLAST result, this section did not performed on Bacterial data, the sequences retrieved from unidentified fasta file that you send
remove.fungal.otu <-c ("058e45742c60a6818c2ab64b22c0c2cb",
                       "9be48d5296dbf4505f48bf3c520aeb5d",
                       "d623a39f97eb1f0f59f4216d92f5ea3d",
                       "09777868382bc4303d7524923d7af7b5",
                       "7d66db88a8752c04cd7e4309dae54d9b",
                       "bc249cfe2032df93d377234010b2347b",
                       "4704e4262a3720a8edb161dbb0afc11e",
                       "8fcf67e5d0f7587dc1ca05520bb4aeb5",
                       "c3c940c1d572143ebcfc0bec62217d53",
                       "98607934c6a92804511c567b9c20171b",
                       "ac0c515d8aeaed2233c406b97283559f",
                       "02cdb14e84cc50b11919e8616fdb03be",
                       "0d2b3f7ab38e883849cb84328058b244",
                       "f92a1fa9e07e7f21d81fd582b0fcdefb",
                       "5f5a28d945d0e6d34be017bdc113f464",
                       "7bfddf6dfdb5c4e0736f595755370a10",
                       "e368180b4ee847da71af1f1898dc6f78",
                       "db9e249b13246882f64f7b46a26f4d9b",
                       "36a54f37e92e9f8f3576c9b7d11bb112",
                       "303ad96a9ed59c99c4edaf1a8ddb325c",
                       "e92972f3a614ade16b4f55d353a67efe",
                       "6be93ca2bb84cc629e883f8dffd9ddd2",
                       "10bc93f77666d669dc37205f348eef5c",
                       "6acac390d4cbec364aecf6020245351d",
                       "f5c4af27ecf964418e8fb3228eacf9f3",
                       "7c82f58ac77a94dd8c32c17e3c8f3767",
                       "fe0e2b0e722dcde915646dbc635bc9a7",
                       "f4b51a2413c58d968c2c251476a88c07",
                       "9ff4d482dc75346e3888f2fb5a39d960",
                       "2c82c9def899e2f4fa9bdd447ca7f54f",
                       "0893e66cb445ff3bb03018d748044f4d",
                       "f246be2ef83132cc2848ab75b9c6e40c",
                       "21d4f1ee7197327393c7d110b5ea2ada",
                       "5ba4848b3f56e9ccb129060eec7e0339",
                       "2bb0130d0900f34b1de5f6ca265374eb",
                       "defe625ebde27e72b22d6fd9de55ac9d",
                       "bdab322ca62cfc43368cf5d357409aa2",
                       "8f5810f52138b9b1ded3c23608b47b04",
                       "77267a929db9e20f484816f654bb07c3",
                       "f08d429f583ab7630079f55616c7c3df",
                       "8067591fce42915f400c20ed6c90f55d",
                       "79d704f5da691f29bd25d2ab5ed16167",
                       "5cf5aaf869e31ecb1f716c92b73cd1cb",
                       "ce79fe2eb8822b4e6197e8209643f6af",
                       "94f96e35ba44d89ca1107f1e8cdf84b3",
                       "21979cd2b409319bf1e0c585878fac3a",
                       "047d152baacc959647754e0e80495177",
                       "e6d074f1ca53ce4cfd5c082a4fca8424",
                       "451768ca78fe73eb525a535437d23f1f",
                       "95cd5a97764dfe098c07fc374c2c6f3d",
                       "07946677abe536056f150ae6fafb2150",
                       "993b983e836f16289f4031e744113ef8",
                       "ca00ab029d53c553292413db8b16c35b",
                       "753b2c5485f3b1492f1ef3a3ae861916",
                       "148f4022c0e2767313ecc5ec30880dfa",
                       "8509dc40bc44a5c90223a5187e037a80",
                       "14b11f5f48e8c708d538aaa52e1fe87c",
                       "86bf1c55ef35800a60e053295a1cb5d1",
                       "4c659b3237e1ff7bdfd12a51f261a947",
                       "92418d6d7ef2d661053266038f87897b",
                       "e30f33b1cbc7e363adca5cec1e37ee09",
                       "f220a224508e503e0d4a39d71647dbdb",
                       "511532e39d7b6ca3e3ae8d1babc8b525",
                       "687053a30fa278069a5bdde66ba12bba",
                       "b9b2879f0e4aeddc5960d84b1c859feb",
                       "aa65b913cecb070ae2c86b92176b3631",
                       "e77587d48323886118083590e7dedb72",
                       "5b395f70d83d534530da509d08438532",
                       "dee9a56b833d461761539cafadce181c",
                       "181dda9ab7da808c4ffc676ba2431b71",
                       "c2c87af4d81faffcc286286dbe2a78c0",
                       "cbc7727faa6d9a9757d3238e9c05500e",
                       "9a13ef020ca2982a2997c54cfebb921a",
                       "35f3a4fca8492cc7839aa3aa2c92f479",
                       "01115335e9e04af711b22bbfd73ca31b",
                       "3b09509132a4f02716ce620cb63edfd3",
                       "ab36e1a754a9679ea84e2151036acbca",
                       "c7f8eca5ec46b8bbd15acd982b9b6b31",
                       "fe4c69025d1452a73d0cbd1d32d224df",
                       "2b32db85afb33d0724a8e40b0d7a4235",
                       "98d82e174693fd390ab3f59e6b4d574b",
                       "0d4bef21e4cc8c73d25aaf2ba35bcf63",
                       "b073e8024a4336921a33a2840ec15123",
                       "8d58617cca69d167e4f0c96a66f0c033",
                       "691c2cd60790223cfcb20b91a847c8f8",
                       "b7ee50478a7ed85b69b273c36ac466e1",
                       "1e5e2da93e1e9739768bb46842a7a197",
                       "e3987284444d0574c46b78d63c6a686f",
                       "cd5cf5f04c626628f79c8d3d50e62c94",
                       "5be2b83894a05f6b026ae7cfb5dc673c",
                       "de42f070dff4a4c0dac4ddcd18dfe920",
                       "d05f82b8eb06615daf94ac1068bcee96",
                       "5fdf7cdafc703696088b97f889db30e4",
                       "ec6fd5e007ab2d0534b2f142d580f3eb",
                       "03f0504dee0935a6bbb14bc4ee30e885",
                       "6adcb4226b8ec9c3aac9b84c16bb97bc",
                       "c351b1aa7b2c59667c08db2b878c474c",
                       "cc4e7ceb802643d9ea51030dbcc38ddc",
                       "c913977bd4c0bdbfa1b6c18db2e30299",
                       "3ebdba339289b44624ccd4a6585aeed8",
                       "e32baf7d65717d2f29043cf93b37d105",
                       "5623b47c4c3b87b37b8aeaf52c14ac75",
                       "d452246bda8c617341edaa86dc72ae6b",
                       "60e668d58485847250c2da4702d42379",
                       "6d9f9a143f79f154fd67c0b0e9082e1d",
                       "04ea526654d02532ca1a8247ccbfa993",
                       "537cac483fe5bd3fae7f80b56179c88e",
                       "4fdf84935b2f44dcea536edef87ef369",
                       "7bececaca2d7ed0bafd396ef30ff245e",
                       "3e0f1af5a7fb346918ebdc2c09e80191",
                       "188891a204844d61b175ba4fc814f449",
                       "c36a2fcdbe2e592178af284ddf9b33fb",
                       "da19a81e6ae2577bc8d2a35ce2ff1a9e",
                       "3e0f1af5a7fb346918ebdc2c09e80191",
                       "188891a204844d61b175ba4fc814f449",
                       "c36a2fcdbe2e592178af284ddf9b33fb"
)
                        

fun.clean.ss <- pop_taxa(fun.clean.ss, remove.fungal.otu)
summarize_phyloseq(fun.clean.ss)
fun.clean.ss #874 taxa left after this. 843 taxa left after the feedback dated 02.09.2022

## Designating OTU id
library(phyloseq
        )
bac.list <- bac.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
fun.list <- fun.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))




bac.list$number <- paste0('B',1:dim(bac.list)[1])
bac.list

bac.list$OTU_id <- ifelse(is.na(bac.list$Genus),ifelse(is.na(bac.list$Family),paste0(bac.list$number,'_o_',bac.list$Order),paste0(bac.list$number,'_f_',bac.list$Family)),paste0(bac.list$number,'_',bac.list$Genus))
bac.list$OTU_id

fun.list$number <- paste0('F',1:dim(fun.list)[1])
fun.list

fun.list$OTU_id <- ifelse(is.na(fun.list$Genus),ifelse(is.na(fun.list$Family),paste0(fun.list$number,'_o_',fun.list$Order),paste0(fun.list$number,'_f_',fun.list$Family)),paste0(fun.list$number,'_',fun.list$Genus))
fun.list$OTU_id


bac.list

fun.list

otu.list <- rbind(bac.list,fun.list)
dim(otu.list)
write.csv(otu.list,'otu_id.csv')

write.csv(bac.list, 'bac_otu_list.csv')
write.csv(fun.list, 'fun_otu_list.csv')
OTU_id.list <- rbind(bac.list[c('OTU','OTU_id')],fun.list[c('OTU','OTU_id')])
OTU_id.list$OTU_id

#seperated files
#CW-CA
sample_data(bac.clean.ss)$PlantSpecies <- gsub(" ", "", sample_data(bac.clean.ss)$PlantSpecies)

bac.CW.CA.ss <- subset_samples(bac.clean.ss, PlantSpecies %in% c('Cynanchumwilfordii','Cynanchumauriculatum'))

bac.list <- bac.CW.CA.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

bac.list$number <- paste0('B',1:dim(bac.list)[1])
bac.list

bac.list$OTU_id <- ifelse(is.na(bac.list$Genus),ifelse(is.na(bac.list$Family),paste0(bac.list$number,'_o_',bac.list$Order),paste0(bac.list$number,'_f_',bac.list$Family)),paste0(bac.list$number,'_',bac.list$Genus))
bac.list$OTU_id

