# check the data 
##Relative abundance


############
##Read OTU tables and convert to phyloseq object
bac.clean.ss


##phyloseq files
bac.clean.ss

fun.clean.ss


# Calculate compositional version of the data
# (relative abundances)
bac.clean.ss.rel <- microbiome::transform(bac.clean.ss, "compositional")
t(otu_table(bac.clean.ss.rel))

fun.clean.ss.rel <- microbiome::transform(fun.clean.ss, "compositional")
t(otu_table(fun.clean.ss.rel))

# Filter the data to include only healthy subjects
#pseq.1 <- subset_samples(phy.clean.l, ibd_subtype == "HC" & timepoint == "1") 
#print(pseq.1)
# keep only taxa with positive sums
#pseq.2 <- prune_taxa(taxa_sums(pseq.1) > 0, pseq.1)

#print(pseq.2)

#Relative population frequencies; at 1% compositional abundance threshold:
# 
# head(prevalence(bac.clean.ss.rel, detection = 0.01, sort = TRUE))
# 
# head(prevalence(fun.clean.ss.rel, detection = 0.01, sort = TRUE))
# 

#We can see that only OTU ids are listed with no taxonomic information. Absolute population frequencies (sample count):


#Core microbiota analysis
#If you only need the names of the core taxa, do as follows. This returns the taxa that exceed the given prevalence and detection thresholds.


core.bac.50 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 50/100)
core.bac.50
[1] "408af2b2f704cbffb6bf8244cbbffc02" "4d49755c33e4090e62cf211f0f771f71"
[3] "12f536207ed71bd9bd113ce252a826f5" "05ab884db890b8c86bcebb56119f32ed"
[5] "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6"
[7] "8e1dcba358617dc87ba3c52ab1b9083f" "444a211bd5cb8bcbf6bd29f57f552286"
[9] "76e17fad933f1d388073d173734d6bd9" "9304494ad5b82f0ad6c585d39f883f82"
[11] "8f3dec6ee90c1c720a3434ce393a3ac7" "7b616034f3ca5cd6fd664d04efa99d8e"
[13] "34d9c708d20438326cb8ad7bb8574b25" "d5a97f89689286ae2e98fc103a7b21e0"
[15] "b1854c63cb99e8c9cf57ea46a0dd9bc7" "34f407e3d79a2b08a3c9ac4c3275f192"
[17] "84f7c22212eaaa840a909947bd53b5a0" "2457d3a95c5d465e4a5c8247fa7b6976"
[19] "c0c52d61eab3d65a001a1a11aacbb1f3" "d3a42a4240197f4e44457d9fb21469d9"
[21] "d4337fb69d0c52d5ffada5b4fcb6d6b7" "7081f2bd5e858632f3c21b7f70101af6"
[23] "f2a15a4db58c4d18503fe81dc8bc9449" "cf15d93871a5b0ccc29932f59aee0941"
[25] "66012bb025b201b70d944c0ae82fa0d1" "96d0f816130691843827c7bfc9d13bc8"

core.bac.60 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 60/100)
core.bac.60  
[1] "4d49755c33e4090e62cf211f0f771f71" "12f536207ed71bd9bd113ce252a826f5"
[3] "05ab884db890b8c86bcebb56119f32ed" "58a15608ac27fe347c20b8575a7a98db"
[5] "8170867ec0b84b57a10a4c40134177c6" "444a211bd5cb8bcbf6bd29f57f552286"
[7] "76e17fad933f1d388073d173734d6bd9" "9304494ad5b82f0ad6c585d39f883f82"
[9] "7b616034f3ca5cd6fd664d04efa99d8e" "d5a97f89689286ae2e98fc103a7b21e0"
[11] "34f407e3d79a2b08a3c9ac4c3275f192" "84f7c22212eaaa840a909947bd53b5a0"
[13] "c0c52d61eab3d65a001a1a11aacbb1f3" "d3a42a4240197f4e44457d9fb21469d9"
[15] "7081f2bd5e858632f3c21b7f70101af6" "f2a15a4db58c4d18503fe81dc8bc9449"
[17] "cf15d93871a5b0ccc29932f59aee0941" "96d0f816130691843827c7bfc9d13bc8"
core.bac.70 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 70/100)
core.bac.70
[1] "4d49755c33e4090e62cf211f0f771f71" "12f536207ed71bd9bd113ce252a826f5"
[3] "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6"
[5] "76e17fad933f1d388073d173734d6bd9" "9304494ad5b82f0ad6c585d39f883f82"
[7] "7b616034f3ca5cd6fd664d04efa99d8e" "34f407e3d79a2b08a3c9ac4c3275f192"
[9] "84f7c22212eaaa840a909947bd53b5a0" "c0c52d61eab3d65a001a1a11aacbb1f3"
[11] "d3a42a4240197f4e44457d9fb21469d9" "7081f2bd5e858632f3c21b7f70101af6"
[13] "cf15d93871a5b0ccc29932f59aee0941" "96d0f816130691843827c7bfc9d13bc8"
core.bac.75 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 75/100)
core.bac.75
[1] "4d49755c33e4090e62cf211f0f771f71" "12f536207ed71bd9bd113ce252a826f5" "58a15608ac27fe347c20b8575a7a98db"
[4] "8170867ec0b84b57a10a4c40134177c6" "9304494ad5b82f0ad6c585d39f883f82" "34f407e3d79a2b08a3c9ac4c3275f192"
[7] "d3a42a4240197f4e44457d9fb21469d9" "96d0f816130691843827c7bfc9d13bc8"

core.bac.80 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 80/100)
core.bac.80
[1] "4d49755c33e4090e62cf211f0f771f71" "58a15608ac27fe347c20b8575a7a98db"
[3] "8170867ec0b84b57a10a4c40134177c6" "9304494ad5b82f0ad6c585d39f883f82"
[5] "34f407e3d79a2b08a3c9ac4c3275f192" "96d0f816130691843827c7bfc9d13bc8"
core.bac.85 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 85/100)
core.bac.85
[1] "4d49755c33e4090e62cf211f0f771f71" "58a15608ac27fe347c20b8575a7a98db"
[3] "8170867ec0b84b57a10a4c40134177c6" "9304494ad5b82f0ad6c585d39f883f82"
core.bac.90 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 90/100)
core.bac.90
[1] "8170867ec0b84b57a10a4c40134177c6" "9304494ad5b82f0ad6c585d39f883f82"
core.bac.95 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 95/100)
core.bac.95
character(0)


tax_core_bac <- subset(bac.list, OTU%in% core.bac.70)
tax_core_bac <- subset(bac.list, OTU%in% core.bac.75)
tax_core_bac <- subset(bac.list, OTU%in% core.bac.80)
tax_core_bac <- subset(bac.list, OTU%in% core.bac.85)
tax_core_bac <- subset(bac.list, OTU%in% core.bac.90)


write.csv(tax_core_bac,"bac Core OTUs-70% prevalence.csv")
write.csv(tax_core_bac,"bac Core OTUs-75% prevalence.csv")
write.csv(tax_core_bac,"bac Core OTUs-80% prevalence.csv")
write.csv(tax_core_bac,"bac Core OTUs-85% prevalence.csv")
write.csv(tax_core_bac,"bac Core OTUs-90% prevalence.csv")

tax_core_bac <- subset(bac.list, OTU%in% core.bac.70)




tax_core_ra_fun <- subset(fun.list, OTU%in% core.fun.70.ra)
tax_core_fun <- subset(fun.list, OTU%in% core.fun.70)
write.table(tax_core_fun, "fungal 70% core.tsv", sep='\t', quote=F)

otu_fun.rel<-otu_table(fun.clean.ss.rel)
otu_fun.rel.core <- subset(otu_fun.rel, rownames(otu_fun.rel)%in%core.fun.70)
write.table(otu_fun.rel.core, "fungal abundance table of 70% core-RA.tsv", sep='\t', quote=F)

otu_fun.normal<-otu_table(fun.clean.nolog)
otu_fun.normal.core <- subset(otu_fun.normal, rownames(otu_fun.normal)%in%core.fun.70.ra)
write.table(otu_fun.normal.core, "fungal normalized abundance table of 70% core-RA.tsv", sep='\t', quote=F)

## Core's correlation with factors
### Abundance table
b.meta <- sample_data(bac.clean.nolog)
b.meta <- data.frame(b.meta)
otu_bac.rel.core.t <- t(otu_bac.rel.core)
b.meta.abun <- merge(b.meta, otu_bac.rel.core.t, by = "row.names")

otu_fun.rel.core.t <- t(otu_fun.rel.core)
f.meta.abun <- merge(f.meta, otu_fun.rel.core.t, by = "row.names")

## Correlation between OTU abundance and each quantitative factors
b.meta.abun$PlantSpecies <- as.factor(b.meta.abun$PlantSpecies)
cor_5.otu <- Hmisc::rcorr(b.meta.abun$PlantSpecies, b.meta.abun$`4d49755c33e4090e62cf211f0f771f71`, type="spearman")
cor_5.otu$r
x         y
x 1.0000000 0.3900736
y 0.3900736 1.0000000
cor_5.otu$P
x           y
x          NA 0.001440919
y 0.001440919          NA
cor_5.otu <- Hmisc::rcorr(b.meta.abun$PlantSpecies, b.meta.abun$`58a15608ac27fe347c20b8575a7a98db`, type="spearman")
cor_5.otu$r
x         y
x 1.0000000 0.1553775
y 0.1553775 1.0000000
cor_5.otu$P
x         y
x        NA 0.2202052
y 0.2202052        NA

cor_5.otu <- Hmisc::rcorr(b.meta.abun$PlantSpecies, b.meta.abun$`8170867ec0b84b57a10a4c40134177c6`, type="spearman")
cor_5.otu$r
x         y
x 1.0000000 0.1539712
y 0.1539712 1.0000000

cor_5.otu$P
x         y
x        NA 0.2244607
y 0.2244607        NA

cor_5.otu <- Hmisc::rcorr(b.meta.abun$PlantSpecies, b.meta.abun$`9304494ad5b82f0ad6c585d39f883f82`, type="spearman")
cor_5.otu$r
x          y
x  1.0000000 -0.6065845
y -0.6065845  1.0000000

cor_5.otu$P
x           y
x          NA 1.08025e-07
y 1.08025e-07          NA




############CA-CW trial
##Read OTU tables and convert to phyloseq object
bac.CW.CA.ss


##phyloseq files
bac.CW.CA.ss

fun.CW.CA.ss


# Calculate compositional version of the data
# (relative abundances)
bac.clean.ss.rel <- microbiome::transform(bac.CW.CA.ss, "compositional")
t(otu_table(bac.clean.ss.rel))

fun.clean.ss.rel <- microbiome::transform(fun.CW.CA.ss, "compositional")
t(otu_table(fun.clean.ss.rel))




core.bac.50 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 50/100)
core.bac.50
[1] "021c5b891400b6d9b57112db2890f06c" "4d49755c33e4090e62cf211f0f771f71" "12f536207ed71bd9bd113ce252a826f5"
[4] "6b9a8700d8f2ed7461b4973139646c29" "8f58c2a233948004f0537ad6bea6453f" "108f40371081a148d16b9b10331ef6a5"
[7] "894fb605a9863f550d4a8257bb4b8f85" "f3afdf8dd75c7e949735241f6d31c1e6" "49232c23e27e46b6ce604f29f75e1da5"
[10] "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6" "9abe16354bae1697b20bf6b8b41853f3"
[13] "8e1dcba358617dc87ba3c52ab1b9083f" "76e17fad933f1d388073d173734d6bd9" "9304494ad5b82f0ad6c585d39f883f82"
[16] "34d9c708d20438326cb8ad7bb8574b25" "684b79b9e58e119cc5584d975e9bcf78" "d5a97f89689286ae2e98fc103a7b21e0"
[19] "34f407e3d79a2b08a3c9ac4c3275f192" "84f7c22212eaaa840a909947bd53b5a0" "2457d3a95c5d465e4a5c8247fa7b6976"
[22] "c0c52d61eab3d65a001a1a11aacbb1f3" "77e5ce2cfb4e91a0056d507afdbbf68d" "7081f2bd5e858632f3c21b7f70101af6"
[25] "ee902da0182b78aea19a6cb6c851a7af" "4cbbab0891db355f05a10bdaf05fff34" "cf15d93871a5b0ccc29932f59aee0941"
[28] "a97e5e031cb783f618adfaa58aba963d" "96d0f816130691843827c7bfc9d13bc8"

core.bac.60 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 60/100)
core.bac.60  
[1] "4d49755c33e4090e62cf211f0f771f71" "12f536207ed71bd9bd113ce252a826f5" "8f58c2a233948004f0537ad6bea6453f"
[4] "108f40371081a148d16b9b10331ef6a5" "894fb605a9863f550d4a8257bb4b8f85" "f3afdf8dd75c7e949735241f6d31c1e6"
[7] "49232c23e27e46b6ce604f29f75e1da5" "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6"
[10] "8e1dcba358617dc87ba3c52ab1b9083f" "76e17fad933f1d388073d173734d6bd9" "9304494ad5b82f0ad6c585d39f883f82"
[13] "684b79b9e58e119cc5584d975e9bcf78" "d5a97f89689286ae2e98fc103a7b21e0" "34f407e3d79a2b08a3c9ac4c3275f192"
[16] "84f7c22212eaaa840a909947bd53b5a0" "2457d3a95c5d465e4a5c8247fa7b6976" "77e5ce2cfb4e91a0056d507afdbbf68d"
[19] "7081f2bd5e858632f3c21b7f70101af6" "ee902da0182b78aea19a6cb6c851a7af" "4cbbab0891db355f05a10bdaf05fff34"
[22] "a97e5e031cb783f618adfaa58aba963d" "96d0f816130691843827c7bfc9d13bc8"
core.bac.70 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 70/100)
core.bac.70
[1] "4d49755c33e4090e62cf211f0f771f71" "108f40371081a148d16b9b10331ef6a5" "894fb605a9863f550d4a8257bb4b8f85"
[4] "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6" "9304494ad5b82f0ad6c585d39f883f82"
[7] "684b79b9e58e119cc5584d975e9bcf78" "d5a97f89689286ae2e98fc103a7b21e0" "34f407e3d79a2b08a3c9ac4c3275f192"
[10] "84f7c22212eaaa840a909947bd53b5a0" "7081f2bd5e858632f3c21b7f70101af6" "a97e5e031cb783f618adfaa58aba963d"
[13] "96d0f816130691843827c7bfc9d13bc8"
core.bac.75 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 75/100)
core.bac.75
[1] "4d49755c33e4090e62cf211f0f771f71" "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6"
[4] "9304494ad5b82f0ad6c585d39f883f82" "84f7c22212eaaa840a909947bd53b5a0" "7081f2bd5e858632f3c21b7f70101af6"
[7] "a97e5e031cb783f618adfaa58aba963d"

core.bac.80 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 80/100)
core.bac.80
[1] "4d49755c33e4090e62cf211f0f771f71" "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6"
[4] "9304494ad5b82f0ad6c585d39f883f82" "84f7c22212eaaa840a909947bd53b5a0" "7081f2bd5e858632f3c21b7f70101af6"
[7] "a97e5e031cb783f618adfaa58aba963d"
core.bac.85 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 85/100)
core.bac.85
[1] "4d49755c33e4090e62cf211f0f771f71" "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6"
[4] "9304494ad5b82f0ad6c585d39f883f82" "7081f2bd5e858632f3c21b7f70101af6"
core.bac.90 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 90/100)
core.bac.90
[1] "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6" "9304494ad5b82f0ad6c585d39f883f82"
core.bac.95 <- core_members(bac.clean.ss.rel, detection = 0, prevalence = 95/100)
core.bac.95
[1] "58a15608ac27fe347c20b8575a7a98db" "8170867ec0b84b57a10a4c40134177c6" "9304494ad5b82f0ad6c585d39f883f82"
character(0)


tax1_core_bac <- subset(bac.list, OTU%in% core.bac.70)
tax2_core_bac <- subset(bac.list, OTU%in% core.bac.75)
tax3_core_bac <- subset(bac.list, OTU%in% core.bac.80)
tax4_core_bac <- subset(bac.list, OTU%in% core.bac.85)
tax5_core_bac <- subset(bac.list, OTU%in% core.bac.90)
tax6_core_bac <- subset(bac.list, OTU%in% core.bac.95)

write.csv(tax1_core_bac,"bac CW-CA Core OTUs-70% prevalence.csv")
write.csv(tax2_core_bac,"bac CW-CA Core OTUs-75% prevalence.csv")
write.csv(tax3_core_bac,"bac CW-CA Core OTUs-80% prevalence.csv")
write.csv(tax4_core_bac,"bac CW-CA Core OTUs-85% prevalence.csv")
write.csv(tax5_core_bac,"bac CW-CA Core OTUs-90% prevalence.csv")
write.csv(tax6_core_bac,"bac CW-CA Core OTUs-90% prevalence.csv")


#FUNGi cw-ca
core.fun.50 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 50/100)
core.fun.50
[1] "251c3304e559314cbe1ddd1196b66fe3" "fbc512a33651bfc19ff67d9a19a26be2" "7902a69089f69b1838947eac2af665ce"
[4] "722bb1dea751bf6f51a025823c04d3c7" "b7d56027f65de48259ecde66b95a2247" "a362cbc76a0a9b2bd9b241c0212ffc5e"
[7] "15de6498f020a65a0c78712f16910dd5" "080d4dea4c8475ae76ffdee07dc2c28c" "fd8adf73b2a4df1a36a0cb3f59668edb"
[10] "3ba279c8f85f27396ed9451d5421468b" "727d121a2c994b24135dad1dd97db133" "90165c66ecbed09057411d35ee7d616e"
[13] "0c4765daaaf31a064c8780ceb6575eaf" "12643d885423aaf82af592b55c660f19" "40747166cfcfc3d95b8946d2f46b6b99"
[16] "6d34111b2197daa594aaf8d86b5df65c" "e1bc523f62796664a63dac4dc7a51495" "2852c9b3a8a7199cc0d30d8356543b41"
[19] "3e8cf843ae1e99d8a13f65659dcc103d" "bc99afc5cb72e66fe6d101b4aa059f01" "23128e0c6dfe2c3229f1fed66a838989"
[22] "e0f04ceb7064e1ab035e0f31ef74d1f7" "5b54c5c00abaaffd4fc3de04c0c3fe99" "d324a225f173b5e8e905208d3f5b360d"

core.fun.60 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 60/100)
core.fun.60  
[1] "fbc512a33651bfc19ff67d9a19a26be2" "7902a69089f69b1838947eac2af665ce" "722bb1dea751bf6f51a025823c04d3c7"
[4] "b7d56027f65de48259ecde66b95a2247" "a362cbc76a0a9b2bd9b241c0212ffc5e" "15de6498f020a65a0c78712f16910dd5"
[7] "080d4dea4c8475ae76ffdee07dc2c28c" "fd8adf73b2a4df1a36a0cb3f59668edb" "727d121a2c994b24135dad1dd97db133"
[10] "90165c66ecbed09057411d35ee7d616e" "0c4765daaaf31a064c8780ceb6575eaf" "12643d885423aaf82af592b55c660f19"
[13] "40747166cfcfc3d95b8946d2f46b6b99" "6d34111b2197daa594aaf8d86b5df65c" "e1bc523f62796664a63dac4dc7a51495"
[16] "2852c9b3a8a7199cc0d30d8356543b41" "3e8cf843ae1e99d8a13f65659dcc103d" "bc99afc5cb72e66fe6d101b4aa059f01"
[19] "e0f04ceb7064e1ab035e0f31ef74d1f7" "5b54c5c00abaaffd4fc3de04c0c3fe99" "d324a225f173b5e8e905208d3f5b360d"
core.fun.70 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 70/100)
core.fun.70
[1] "fbc512a33651bfc19ff67d9a19a26be2" "722bb1dea751bf6f51a025823c04d3c7" "b7d56027f65de48259ecde66b95a2247"
[4] "a362cbc76a0a9b2bd9b241c0212ffc5e" "15de6498f020a65a0c78712f16910dd5" "080d4dea4c8475ae76ffdee07dc2c28c"
[7] "fd8adf73b2a4df1a36a0cb3f59668edb" "727d121a2c994b24135dad1dd97db133" "90165c66ecbed09057411d35ee7d616e"
[10] "0c4765daaaf31a064c8780ceb6575eaf" "40747166cfcfc3d95b8946d2f46b6b99" "6d34111b2197daa594aaf8d86b5df65c"
[13] "2852c9b3a8a7199cc0d30d8356543b41" "3e8cf843ae1e99d8a13f65659dcc103d" "d324a225f173b5e8e905208d3f5b360d"
core.fun.75 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 75/100)
core.fun.75
[1] "fbc512a33651bfc19ff67d9a19a26be2" "b7d56027f65de48259ecde66b95a2247" "a362cbc76a0a9b2bd9b241c0212ffc5e"
[4] "15de6498f020a65a0c78712f16910dd5" "080d4dea4c8475ae76ffdee07dc2c28c" "90165c66ecbed09057411d35ee7d616e"
[7] "0c4765daaaf31a064c8780ceb6575eaf" "40747166cfcfc3d95b8946d2f46b6b99" "6d34111b2197daa594aaf8d86b5df65c"
[10] "3e8cf843ae1e99d8a13f65659dcc103d" "d324a225f173b5e8e905208d3f5b360d"

core.fun.80 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 80/100)
core.fun.80
[1] "fbc512a33651bfc19ff67d9a19a26be2" "b7d56027f65de48259ecde66b95a2247" "a362cbc76a0a9b2bd9b241c0212ffc5e"
[4] "15de6498f020a65a0c78712f16910dd5" "080d4dea4c8475ae76ffdee07dc2c28c" "90165c66ecbed09057411d35ee7d616e"
[7] "0c4765daaaf31a064c8780ceb6575eaf" "40747166cfcfc3d95b8946d2f46b6b99" "6d34111b2197daa594aaf8d86b5df65c"
[10] "3e8cf843ae1e99d8a13f65659dcc103d" "d324a225f173b5e8e905208d3f5b360d"
core.fun.85 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 85/100)
core.fun.85
[1] "fbc512a33651bfc19ff67d9a19a26be2" "b7d56027f65de48259ecde66b95a2247" "15de6498f020a65a0c78712f16910dd5"
[4] "080d4dea4c8475ae76ffdee07dc2c28c" "90165c66ecbed09057411d35ee7d616e" "40747166cfcfc3d95b8946d2f46b6b99"
[7] "3e8cf843ae1e99d8a13f65659dcc103d" "d324a225f173b5e8e905208d3f5b360d"
core.fun.90 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 90/100)
core.fun.90
[1] "fbc512a33651bfc19ff67d9a19a26be2" "b7d56027f65de48259ecde66b95a2247" "15de6498f020a65a0c78712f16910dd5"
[4] "3e8cf843ae1e99d8a13f65659dcc103d" "d324a225f173b5e8e905208d3f5b360d"
core.fun.95 <- core_members(fun.clean.ss.rel, detection = 0, prevalence = 95/100)
core.fun.95
[1] "b7d56027f65de48259ecde66b95a2247" "15de6498f020a65a0c78712f16910dd5" "3e8cf843ae1e99d8a13f65659dcc103d"
[4] "d324a225f173b5e8e905208d3f5b360d"


tax1_core_fun <- subset(fun.list, OTU%in% core.fun.70)
tax2_core_fun <- subset(fun.list, OTU%in% core.fun.75)
tax3_core_fun <- subset(fun.list, OTU%in% core.fun.80)
tax4_core_fun <- subset(fun.list, OTU%in% core.fun.85)
tax5_core_fun <- subset(fun.list, OTU%in% core.fun.90)
tax6_core_fun <- subset(fun.list, OTU%in% core.fun.95)



write.csv(tax1_core_fun,"fun CW-CA Core OTUs-70% prevalence.csv")
write.csv(tax2_core_fun,"fun CW-CA Core OTUs-75% prevalence.csv")
write.csv(tax3_core_fun,"fun CW-CA Core OTUs-80% prevalence.csv")
write.csv(tax4_core_fun,"fun CW-CA Core OTUs-85% prevalence.csv")
write.csv(tax5_core_fun,"fun CW-CA Core OTUs-90% prevalence.csv")
write.csv(tax6_core_fun,"fun CW-CA Core OTUs-95% prevalence.csv")

