###PERMANOVA
## ALL
bac.clean.log
b.otu <- otu_table(bac.clean.log)
b.meta <- sample_data(bac.clean.log)
b.meta <- data.frame(b.meta)

b.meta$PlantSpecies <- as.factor(b.meta$PlantSpecies)
b.permanova <- adonis2(formula = t(b.otu) ~ (PlantSpecies), data = b.meta, permutations=9999, method = "bray")
b.permanova


#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999

#adonis2(formula = t(b.otu) ~ (PlantSpecies), data = b.meta, permutations = 9999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#PlantSpecies  7   8.4414 0.41499 5.6751  1e-04 ***
#Residual     56  11.8996 0.58501                  
#Total        63  20.3410 1.00000                  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

fun.clean.log
f.otu <- otu_table(fun.clean.log)
f.meta <- sample_data(fun.clean.log)
f.meta <- data.frame(f.meta)


f.meta$PlantSpecies <- as.factor(f.meta$PlantSpecies)
f.permanova <- adonis2(formula = t(f.otu) ~ (PlantSpecies), data = f.meta, permutations=9999, method = "bray")
f.permanova

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999

#adonis2(formula = t(f.otu) ~ (PlantSpecies), data = f.meta, permutations = 9999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#PlantSpecies  7   9.1559 0.53786 9.3109  1e-04 ***
#  Residual     56   7.8668 0.46214                  
#Total        63  17.0227 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#cw-ca
bac.clean.log
b.otu <- otu_table(bac.clean.log)
b.meta <- sample_data(bac.clean.log)
b.meta <- data.frame(b.meta)

b.meta$PlantSpecies <- as.factor(b.meta$PlantSpecies)
b.meta$SampleID <- as.factor(b.meta$SampleID)

b.permanova <- adonis2(formula = t(b.otu) ~ (PlantSpecies), data = b.meta, permutations=9999, method = "bray")
b.permanova

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999

#adonis2(formula = t(b.otu) ~ (PlantSpecies), data = b.meta, permutations = 9999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#PlantSpecies  1   1.3568 0.34712 7.4434  1e-04 ***
#  Residual     14   2.5520 0.65288                  
#Total        15   3.9088 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
b.permanova <- adonis2(formula = t(b.otu) ~ (SampleID), data = b.meta, permutations=9999, method = "bray")
b.permanova

#adonis2(formula = t(b.otu) ~ (SampleID), data = b.meta, permutations = 9999, method = "bray")
#Df SumOfSqs R2 F Pr(>F)
#Model    15   3.9088  1         
#Residual  0   0.0000  0         
#Total    15   3.9088  1 


b.permanova <- adonis2(formula = t(b.otu) ~ (SampleID+PlantSpecies), data = b.meta, permutations=9999, method = "bray")
b.permanova
#adonis2(formula = t(b.otu) ~ (SampleID + PlantSpecies), data = b.meta, permutations = 9999, method = "bray")
#Df SumOfSqs R2 F Pr(>F)
#Model    15   3.9088  1         
#Residual  0   0.0000  0         
#Total    15   3.9088  1 

fun.clean.log
f.otu <- otu_table(fun.clean.log)
f.meta <- sample_data(fun.clean.log)
f.meta <- data.frame(f.meta)


f.meta$PlantSpecies <- as.factor(f.meta$PlantSpecies)
f.meta$SampleID <- as.factor(f.meta$SampleID)


f.permanova <- adonis2(formula = t(f.otu) ~ (PlantSpecies), data = f.meta, permutations=9999, method = "bray")
f.permanova
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 9999

#adonis2(formula = t(f.otu) ~ (PlantSpecies), data = f.meta, permutations = 9999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)    
#PlantSpecies  1  0.87378 0.31145 6.3324  3e-04 ***
#  Residual     14  1.93180 0.68855                  
#Total        15  2.80558 1.00000                  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

f.permanova <- adonis2(formula = t(f.otu) ~ (SampleID), data = f.meta, permutations=9999, method = "bray")
f.permanova
#No residual component

#adonis2(formula = t(f.otu) ~ (SampleID), data = f.meta, permutations = 9999, method = "bray")
#Df SumOfSqs R2 F Pr(>F)
#Model    15   2.8056  1         
#Residual  0   0.0000  0         
#Total    15   2.8056  1  
f.permanova <- adonis2(formula = t(f.otu) ~ (PlantSpecies+SampleID), data = f.meta, permutations=9999, method = "bray")
f.permanova
#No residual component

#adonis2(formula = t(f.otu) ~ (PlantSpecies + SampleID), data = f.meta, permutations = 9999, method = "bray")
#Df SumOfSqs R2 F Pr(>F)
#Model    15   2.8056  1         
#Residual  0   0.0000  0         
#Total    15   2.8056  1  