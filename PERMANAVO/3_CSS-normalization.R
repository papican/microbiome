### Normalization

## Let's do normalization with CSS
## phyloseq to metagenomeSeq

#### Bacteria ###


bac.clean.ss

## Remove residual taxa that do not have any sequences
#Bacteria
#sum(taxa_sums(bac.clean.ss) == 0) #result =0
#taxa_sums(bac.clean.ss)

## CODE for CSS normalization using preloaded data
sort(sample_sums(bac.clean.ss))

bac.clean.filt <- bac.clean.ss

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization


p <- cumNormStatFast(met.bac.clean) #Calculates the percentile for which to sum counts up to and scale by.
p #0.5
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog <- bac.clean.filt
otu_table(bac.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log <- bac.clean.filt
otu_table(bac.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)


### Fungi ###
fun.clean.filt <- fun.clean.ss

# (filt.sample <- sample_sums(fun.clean.filt) > 1000)
sum(sample_sums(fun.clean.filt) <= 1000)  ## 0 sample discarded
#fun.clean.over1000 <- prune_samples(filt.sample, fun.clean.filt)
#fun.clean.over1000  ## 

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization


p <- cumNormStatFast(met.fun.clean)
p #0.5
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog <- fun.clean.filt
otu_table(fun.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log <- fun.clean.filt
otu_table(fun.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

fun.clean.log #874 taxa 64 samples
bac.clean.log #3246 taxa and 64 samples




#### CW-CA trial bac###


bac.CW.CA.ss

## Remove residual taxa that do not have any sequences
#Bacteria
#sum(taxa_sums(bac.clean.ss) == 0) #result =0
#taxa_sums(bac.clean.ss)

## CODE for CSS normalization using preloaded data
sort(sample_sums(bac.CW.CA.ss))

bac.clean.filt <- bac.CW.CA.ss

# met.phy.clean <- phyloseq_to_metagenomeSeq(phy.clean.prune)
met.bac.clean <- phyloseq_to_metagenomeSeq(bac.clean.filt)

# normalization


p <- cumNormStatFast(met.bac.clean) #Calculates the percentile for which to sum counts up to and scale by.
p #0.5
met.bac.norm <- cumNorm(met.bac.clean, p =p)

# returns normalized factors for each sample
normFactors(met.bac.norm)
sort(normFactors(met.bac.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.bac.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.bac.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
bac.clean.nolog <- bac.clean.filt
otu_table(bac.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

bac.clean.log <- bac.clean.filt
otu_table(bac.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)



### CW-CA Fungi ###
fun.clean.filt <- fun.CW.CA.ss

# (filt.sample <- sample_sums(fun.clean.filt) > 1000)
sum(sample_sums(fun.clean.filt) <= 1000)  ## 0 sample discarded
#fun.clean.over1000 <- prune_samples(filt.sample, fun.clean.filt)
#fun.clean.over1000  ## 

met.fun.clean <- phyloseq_to_metagenomeSeq(fun.clean.filt)

# normalization


p <- cumNormStatFast(met.fun.clean)
p #0.5
met.fun.norm <- cumNorm(met.fun.clean, p =p)

# returns normalized factors for each sample
normFactors(met.fun.norm)
sort(normFactors(met.fun.norm))

# To export normalized count matrix
met.CSS.log <- MRcounts(met.fun.norm, norm = T, log = T)
met.CSS.nolog <- MRcounts(met.fun.norm, norm = T, log = F)
## this was log2!!

## back to the phyloseq file
fun.clean.nolog <- fun.clean.filt
otu_table(fun.clean.nolog) <- otu_table(met.CSS.nolog, taxa_are_rows = TRUE)

fun.clean.log <- fun.clean.filt
otu_table(fun.clean.log) <- otu_table(met.CSS.log, taxa_are_rows = TRUE)

fun.clean.log #874 taxa 64 samples
