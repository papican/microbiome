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
#remove.bacterial.otu <-c ("8170867ec0b84b57a10a4c40134177c6",
                       "58a15608ac27fe347c20b8575a7a98db",
                       "9abe16354bae1697b20bf6b8b41853f3",
                       "108f40371081a148d16b9b10331ef6a5",
                       "894fb605a9863f550d4a8257bb4b8f85",
                       "8e1dcba358617dc87ba3c52ab1b9083f",
                       "05ab884db890b8c86bcebb56119f32ed",
                       "d750ca49b5b9bfd089fccd8fb5276832",
                       "9c99e394a810585245e3ccb32b2a60c6",
                       "24fa36f95795be6a7a4cacb1b19d9225",
                       "f3afdf8dd75c7e949735241f6d31c1e6",
                       "49232c23e27e46b6ce604f29f75e1da5",
                       "479d653598095238acfd9b677630b2cc",
                       "67e527038cb3ef53d7d1006a63ed9a01",
                       "6fb9b381bd7ef6e66d89e0102ba304e2",
                       "7c9cc30edf93be3d852506770821e29d",
                       "8cf8cfb5f3f8deb095f0fabdc714fb80",
                       "ae1d4b213eefc4c1ca7833ca1d9c63f9",
                       "e6572e70487dde3651fa06553a00b38f",
                       "8136f17c5fe82a1f42a480206f8722a1",
                       "8f58c2a233948004f0537ad6bea6453f",
                       "eebd3c3e373ac778fe644681d316dcbe",
                       "6b9a8700d8f2ed7461b4973139646c29",
                       "f0e33662fdae349689d8f88d5f5f7c40",
                       "c176f6158b62160064477f3c57cace47",
                       "31b98b44271348ed218ded5d8cbbcaf3",
                       "347016993f7528ea0beb87f0aca7a356",
                       "72192680d2296ba564eb35839a05454e",
                       "2e9922672d16b9fe025b6f2a38fffded",
                       "8e74ed7c2840d0dbce33afa2e69d4ce8",
                       "d3fd6410af36c0b6fc81f5de3ef693aa",
                       "ab98212d6425a9ce7daa69febcadb1f7",
                       "8ed5b000abc6045f7c17ef303c8ea0c5",
                       "905576b9d12d00671814f872e339cded",
                       "7c3a42faf90bfdd268fc7ee6f74e5923",
                       "9e6131b1048cbd963bc80f6c64727f31",
                       "88988cd8c01d2e1ceb8d18c2202944eb",
                       "d71acf7d1720df4bb4dd6b4aa9fd2f71",
                       "b068e1e736dc4b547bca78e4e57eb3fb",
                       "fbf884d5f172b17d9ec9203f4e106a2a",
                       "3f62d5fd38864ba490d123068f36f945",
                       "adc9f2efda1680061da9e43cfca72100",
                       "425cfa3f7684fc63e1edbb0dd7827229",
                       "aa5d8adc08016b3ffa885f97546e535c",
                       "d9633d6e8a9ef8fa3c3a5752757fb517",
                       "868c2e669e41dd494a5fb0e4f1a15ca8",
                       "dcc323752ebda91f52c005254f1f8bd3",
                       "0eb3f44b97f9d2acf2811b0216f4c077",
                       "54f2723cd9dac0b7a81a1c2e853e4631",
                       "d18111ff74102d57c84c5d65befe7135",
                       "d8a61c9e613b1900869d2c3203331158",
                       "26ca1e6bd899f8a3cadd83163e2e73b7",
                       "438203050d1568fdcb88c844123b64a9",
                       "4f9e6b905dba28082ccf6390fd179949",
                       "93235defc85aa63a5a98485618b82d4d",
                       "e1f47ddc3d9a519f956a43018478892c",
                       "f06f8216b07dc463f7a63f8df84f17af",
                       "9231207243d233cc02af76a1099e33bb",
                       "5c7664893f9c241397cf7529aac0dbfc",
                       "42642e3d263a909b7a93625a21e5831b",
                       "126fc06ff24d5faeaebcfa4d8cc68431",
                       "630854e031cc15f800f2373cd3e1eea4",
                       "79493292015e4e1e880d91bc0617b700",
                       "339542f3d6e7701067dc8653d8952280",
                       "9b1e2d02955de5678ae968d50c2641f4",
                       "2dc354aec224c8288201ca4259c5bbbf",
                       "85b6e6f8d0c577f41a5d8250b0e6be7b",
                       "1089bc1e429cb4d867d114c91a0d1510",
                       "555670f4a6185340054f7518fea546a3",
                       "1e94f0f32eb871af1d68ada9f3579d1b",
                       "bc6885cf0cec168c6819df3a61864b3d",
                       "99335c91eb6c8cef52252c3b7d02822d",
                       "c3bb05b3bcbe5f502e5b5d81054ca05b",
                       "452619c4cd4ef958c7f8b6858e323f98",
                       "82474995fcecba1f0102b309b4fc200c",
                       "ed665db12bfe17c856007627d0a6723d",
                       "c14da2ec74127055c8e7e47a2d5f931d",
                       "15040f745eb7d7d784c66d0c6a1f8aef",
                       "0a92b248966c863456a2dbc447f1d4f3",
                       "f746d9fe81d80f86e9c90860767a8176",
                       "61dbe07f5bc05e2482db01b8ef31d8f1",
                       "bfbbca9a9f578a3f3f321cd14c76bb50",
                       "81a0834989fa0cef2b056b7ee764457c",
                       "95125b6da93839ba0dfa68420a309af9",
                       "415b9bdb80001f206740d1d280e1ef35",
                       "17b60c856c7c53a85b27c2216cc6654c",
                       "777f56768eb1ced08accc4f477f2eb6a",
                       "91dcbb7790515d396f6e927d17dd78f0",
                       "c77e430cc028132a910072b550732b41",
                       "ea7fd439bf3b051678ed30e67ef754d2",
                       "d076363e5a2f953d339a6f6141d9e35a",
                       "5be8d468f66e1e0ac0a9eb5f13a047f1",
                       "9cba0a5ae153c5066f530e2d7d7505cf",
                       "859637e91a450ec0d1afef9d114a610b",
                       "64189151dc20d7a88f3fc2a02af70d9e",
                       "75a14122d5c99a42c188a5fa83125048",
                       "532080936345b2fc84d1166fab9a0a7e",
                       "192f873cb2ac5b8b0ae1e7c3cd5c3649",
                       "69d00252284f63355fb0ab01839201f4",
                       "d33fb86a3c47bd1fc45a59c7d15332bb",
                       "636cdf6d90eae767d5a17baf2c3b4f10",
                       "c3edc79c4dd0de1ea83291be2889294e",
                       "17e7f2fe268a0500abe51bf818043063",
                       "5161dfc6049c4e010f3d8f00314f131d",
                       "1d5f3282b9fb9f78554c32fa1d0da7ed",
                       "cbe170f406701182617ef4417a9c9e57",
                       "859bc2fb3555a4d2c51d9ea50d128a86",
                       "c59ad3a29ef4989a41d78aa0c0582c63",
                       "5a2423686ab5ec8ccdc8e689cd38c3bb",
                       "4e9becf073027e46b4b950e438cd9ade",
                       "6170fa8781e2e1eb40fe5fcb1153e4b2",
                       "cdb81b45d07b6be488086946bde75626",
                       "8c971491314bd1bbac71c568d18764c8",
                       "bcb58f308b703c2cd0a539e11a69c6cd",
                       "62eea8e4973348ada6ed719cf6f42eea",
                       "08081220367d39c41c432ae8e1749764",
                       "8a3f12d4960f0283fb58033aff43fbf8",
                       "6e6aff216cb09f39f29b377bca29d051",
                       "f4d8a1bb729d7b528f7c1794c2346a82",
                       "84e590192f62ac27fa87a66a47a8c45a",
                       "e68ca9509639ab6bd9ad6868d07fd28f",
                       "7d976b59a01ead16d633b7afd3fb3d35",
                       "42c7da5d1d6f00b80690bfeb3b84137f",
                       "a3c71d9c597b0cd5c7040ef481d45bba",
                       "a4eb124ec4b49b0bc6201f91eb54c133",
                       "36ee116bc572a746958a2a170657873f",
                       "fd85e2cd6b166f29853ca16993de00fd",
                       "608f776fe20cbc7c5fe634fec7b384f4",
                       "75a8c33cff292313363c80137c25fe98",
                       "70056810b010e5785e5cf9a3d09931fd",
                       "f1793be216aa5e911b24eb6ba27b42df",
                       "98c81e8ceababc27f1ca86c7a43db221",
                       "d6b824508c0c787416faa7f75dfeec6c",
                       "609c435df2b6d13fdf74492302880000",
                       "84ca3c84ea7312d7f3fcf5b91afe956f",
                       "40471b3821b868ac00b6e4c4f16ac53e",
                       "4ba0350b746181fce61958727ffe0748",
                       "3e7635b56c46d7826029789b3f01d772",
                       "16a4d0d0be8ad7f9b1530bac66ba5f6e",
                       "8822390f60e14fbe327d2524a1323323",
                       "9b90655aecb332986ab5a161a9a1927f",
                       "4ef86403815d192d12ccc4dc1581e6c0",
                       "3aa70a7b9a796bebe853b4a56c0eaa43",
                       "5ff650b383aa6af0d7f4888d335dffea",
                       "82aadee94872aa7cb6c13e03965acfa1",
                       "610f10c12dbc55317a098e3af57f8bcc",
                       "b2de2d5cc5a4f01c41572b3075a9945e",
                       "07fab17bcc3d2fd7b1a205904e7ca6a6",
                       "b23c42af7174273c406579b907faa290",
                       "15fb05aceaeb9aaa0abae9981394c03d",
                       "16d93c9b7052ac18d6b3cd56a5c7fd4b",
                       "849021ccc47cb1309a3e24b8b961355d",
                       "50a83355f57f902104023b43dfeb038f",
                       "c3fcc4ba8bd2ee99f37846f54d500a5c",
                       "2ef97cadad700e4485b48ff0c0d9456f",
                       "7b1c0033a0cbe762de0017b3af0a3d1d",
                       "06de47dcfbf8c29b945903b9a72d74b2",
                       "b5f8299918ec2396b95bcbad173c5d34",
                       "4f2674c458b6e63086e03205575c6e57",
                       "0692c4768e9fdb885ec0f5675cc30d5a",
                       "5cbbcb3ba1eda5837b24d2cbd6b20f9a",
                       "9c88d6b2306b376fc057064e8976f8ab",
                       "dba45b80bbc9c6af15b78e393f99a862",
                       "a8ac5957dfe680a2fd28a380257f387f",
                       "7042563ca256dd7a440a93480fce183d",
                       "615a1cb7e8a94a5c29593e29f33b5410",
                       "109bfdb7a38d31bed2027b41853f759d",
                       "090d38f1ce7e3d5ddad2eca60f3c67b8",
                       "0aebbc1f9416e3821161424f74b8df6f",
                       "a5326ec03fe203d8c1d22172540df399",
                       "617ffffe49f80741554d81f1f48b850f",
                       "0e57bc265f888c369d84f9c586614d16",
                       "b04c606f42eccad86ff3fc6de23516d5",
                       "c3fadb181c1621f07d588927ac3d112c",
                       "69c63693df05f83e6c463f805f4c8343",
                       "bc2ee788052b5fdfb1508bb9bfbc585c",
                       "ff0b772ddad6c01f956076d6b3580c65",
                       "450ffd4e4a86a43069e70ae8269cfe45",
                       "55fa6f8c879d5cd0fb6c6db1ce3d7549",
                       "2a727685aeb5f76ade35279b17eae9ee",
                       "386f649a482fff8ca221e31a1871ef66",
                       "c4ad77bf34381ff428ab7b37e7cbfbc9",
                       "45910cbaebbc11c832a3d75ed2dc7eae",
                       "b968bb78046899f8871a49796968d8dc",
                       "4cc9d2ddcdc65b6cef77a43accc35871",
                       "e8213417cf045dc2ebdafb68315fa348",
                       "764a322e74cb614189bfaedf165e1992",
                       "229011a4db87cb7453cd5989fb1337d0",
                       "57fb4b6fb2c80a651016a751ccea2a0c",
                       "29b65148ef231cba5ec3e2cade7d1f91",
                       "602d294a97f2fe3fc9f885b0091d19eb",
                       "c884b06617f74fdcf585c1a902089b15",
                       "9a977462b027bce7606d1d5d06a38ebf",
                       "9ae34714e5230bfd26c7000136a8ba29",
                       "1812af4dbad33a8b67bed66725262bb6",
                       "6ce20f0eebe1ef6443ed1ca05d91f1ff",
                       "bc2885aaaf3bdd28c31f27cc05c7a4db",
                       "daebcc65b1ff091e12d74f22bda290d7",
                       "d7ca4745b500e47931b435ca89b97112",
                       "0e911357b454b96676904abe9e0d6bae",
                       "0e8b5097cec2354b95605488aed1ab91",
                       "619967c69e3f52acb9df670ba7214ce3",
                       "46809ea68169156b812fd24a52b67843",
                       "e1410b32bee770dfb4f8122546cc3880",
                       "3f70f82c023cbe586c1877089c6ad48b",
                       "2477c856225bb17c588ad804fe1a38a3",
                       "0884208d9e91c097d6f457a1d37c2c78",
                       "d732313abd629bf5e122290b12e588a5",
                       "92661763e92696b76eb6751a5e8d7c04",
                       "259365a636305e699b9137488256a767",
                       "2f72c2a56e136ad16e0616326a1be4ca",
                       "44aea0f356e4e0e205af89f3596ef9b7",
                       "ed79612d7964f87ec71a750e92493a4b",
                       "54130e3e2654ee431b9133a90f985012",
                       "956126e26a123c4119274c986619135c",
                       "1fe06ecadf05711c46f089a102bb7f67",
                       "c900c53726929a6f4c88829012c312c6",
                       "d67d9861df2ba69cabc3142dfc3afe38",
                       "6edb774bcd601723d04baf47fa6f00b5",
                       "f146bcf9eb19bc96441c64f1dfe28190",
                       "b9a052aaca8a99a8cb54d210102709fa",
                       "b3d5d3ca13995075c24923445d3b6264",
                       "b053cd4bc21114844ffd3b2d78f6e628",
                       "bf23ad79f678037da6a61bca18ea73fe",
                       "74ddaaac89e220c9a470662861059cac",
                       "5d60f8053034dba291000b888fb0042c",
                       "28a09e67e5a9197cb2b23344e5390e0c",
                       "386459840c5e3cb34ca61ac530b62864",
                       "d47c48ddf2e145603e78193412444f5d",
                       "0f6966d83aa1028c4b00316c616d8e45",
                       "52735e81128256dd66cc26125a22244b",
                       "31c005c497dd79bdbe3accf889b46c0b",
                       "3d8b371b9832fbc0be246f459a0fe593" )


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
                       "0d4bef21e4cc8c73d25aaf2ba35bcf63"
)
                        

fun.clean.ss <- pop_taxa(fun.clean.ss, remove.fungal.otu)
summarize_phyloseq(fun.clean.ss)
fun.clean.ss #874 taxa left after this. 

## Designating OTU id
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

write.csv(bac.list, 'bac_otu_list.xlsx')
write.csv(bac.list, 'fun_otu_list.xlsx')
OTU_id.list <- rbind(bac.list[c('OTU','OTU_id')],fun.list[c('OTU','OTU_id')])
OTU_id.list$OTU_id
