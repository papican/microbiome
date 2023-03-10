knitr::opts_chunk$set(echo = TRUE)


library(purrr)
library(dplyr)
library(tidyr)
library(seqinr)
#devtools::install_github('mhahsler/rBLAST')
library(rBLAST)
#after diff abundant otu analysis
rightside=fun.list
leftside=mydata
leftside_tibbled <- tibble::rownames_to_column(leftside, "ASV")

rightside_tibbled <- tibble::rownames_to_column(rightside, "ASV")

joined=left_join(leftside_tibbled,rightside_tibbled,by='OTU',copy=FALSE)
joined
write.csv(joined,'CW-CA_diff_fungi_OTU.csv',row.names=F)


#getting the sequences of diff otus


# You should install the "seqinr" package first to run the code described below.
bac.seq <- read.fasta(file = "dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")



PG_AC_fun <- read.csv(file = "AR-AP_Bacteria_diff_OTU.csv")

enrichedASVs_CA <- PG_AC_fun$OTU[which(PG_AC_fun$threshold == "Apiaceae")]
difBacSeq.CA<-bac.seq[which(attr(bac.seq, "names") %in% enrichedASVs_CA)]
write.fasta(difBacSeq.CA,bac.list$OTU_id[which(bac.list$OTU %in%attr(difBacSeq.CA, "names"))], "Enriched_ASVs_bacterial_sequences_Apiaceae.fasta")

##getting rid of whitespaces in fasta fles


##### Construct isolate sequence-based curated DB
CuDB<-makeblastdb("Apiaceae_Bacteria.txt", dbtype = "nucl")
bl.b <- blast(db="Apiaceae_Bacteria.txt")
bl.b
CuDB<-makeblastdb("Apiaceae_Fungi.txt", dbtype = "nucl")
bl.f <- blast(db="Apiaceae_Fungi.txt")
bl.f


#### Load differentially abundant bacterial ASV sequences
bac.asv.seq<-readDNAStringSet("Enriched_ASVs_bacterial_sequences_Apiaceae.fasta")
fun.asv.seq<-readDNAStringSet("Enriched_ASVs_fungal_sequences_Apiaceae.fasta")


### BLAST
##bacterial
predic.blast.CW <- predict(bl.b, bac.asv.seq)
write.csv(predic.blast.CW,"diff_asv_blast_result_bacterial_Apiaceae.csv")

##Fungal
predic.blast.CW.fun <- predict(bl.f, fun.asv.seq )
fungi_diff_filtered <- predic.blast.CW.fun %>% filter(Alignment.Length > 100)
write.csv(fungi_diff_filtered,"diff_asv_blast_result_fungi_Apiaceae.csv")










table(predic.blast.CW$SubjectID)

df.blast.p1 <- ggplot(predic.blast.CW, aes(SubjectID, fill=SubjectID)) + 
  geom_bar(stat="count", width = 0.5) +
  #scale_fill_discrete() +
  scale_fill_manual(values = PG_AC_color) +
  
  
  xlab('')+
  ylab("Number of Blast Hits") +
  #ggtitle("Family Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
   scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, by = 40))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 0.2)

df.blast.p1

table(fungi_diff_filtered$SubjectID)

df.blast.p1 <- ggplot(predic.blast.CW.fun, aes(SubjectID, fill=SubjectID)) + 
  geom_bar(stat="count", width = 0.5) +
  #scale_fill_discrete() +
  scale_fill_manual(values = PG_AC_color) +
  xlab('')+
  ylab("Number of Blast Hits") +
  #ggtitle("Family Community Composition by Sample \n") +
  ## adjust positions
  guides(fill = guide_legend(nrow = 5,reverse = T))+
  theme(legend.position="bottom") +
  theme(plot.title = element_text(size = 20,hjust = 0.5, face='bold')) + 
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 13,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.4,size=10, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10))+
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.background=element_blank(),panel.border=element_blank(), plot.background=element_blank())+
  theme(aspect.ratio = 0.2)

df.blast.p1
#colors
##CW-CA
CW_CA_color <- c("#5F7FC7", "orange",  "#AD6F3B",
"#673770","#D14285", "#652926", "#C84248",  "#8569D5",
"#5E738F","#D1A33D", "#8A7C64", "#599861","#616163",
"#FFCDB2", "#242F40","indianred")
PG_AC_color <- c("#6D9F71",  "#CCA43B", "#F92A82",
  "#ED7B84", "#7EB77F", "#DEC4A1", "#E5D1D0", '#0E8482',
  '#C9DAEA', '#337357', '#95C623', '#E55812', '#04471C',
  '#F2D7EE', '#D3BCC0', '#A5668B', '#69306D', 'navy',
  '#1A535C', '#4ECDC4', 'orange', '#FF6B6B', "orchid1",
  'cyan2',"#599861","#616163",
  "#FFCDB2", "#242F40","indianred")


