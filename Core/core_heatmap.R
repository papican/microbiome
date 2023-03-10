


### core microbiome heatmap

#bacteria
# get_df_rel <- function(phy.clean.ss.5, b.dom){
df.otu <- bac.clean.ss %>% psmelt()
df.otu$Sample <- factor(df.otu$Sample, levels = rev(b.meta$SampleID))
head(df.otu)

# we need to group by samples
df.otu.rel <- df.otu %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.selected.rel <- df.otu.rel %>% filter(OTU %in% core.bac.90)
excel.bac <- df.selected.rel %>%   group_by(OTU, RelAbundance) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

# install.packages('xlsx')
library(xlsx)
#write.xlsx(as.data.frame(excel.bac) , 'core_bacteria.xlsx')

b.add <- merge(excel.bac, bac.list, by ="OTU")
b.order<- b.add %>% dplyr::select(OTU,total.x, OTU_id)%>% group_by(OTU_id)%>% summarise(Totals=sum(total.x)) %>% arrange(desc(Totals))
b.order.vec <- b.order$OTU_id

b.add <- as.data.frame(b.add)
df.rel <- df.selected.rel %>% dplyr::select(OTU,Sample,PlantSpecies,RelAbundance,Abundance) %>% arrange(desc(Abundance))
df.rel <- merge(df.rel,bac.list,by='OTU')
head(df.rel)

df.rel$OTU_id <- factor(df.rel$OTU_id, levels=b.order.vec) 
head(df.rel)
str(df.rel)
# df.rel$id_plus <- paste0(df.rel$id,' ',df.rel$Genus)
df.rel

colnames(df.rel)
#   return(df.selected.rel)
# }



plot_core_heatmap <- function(df_bwil, colfunc){
  textcol <- "black"
  
  ggplot(df_bwil, aes(x=Sample, y=OTU_id, fill=RelAbundance)) +
    geom_tile(colour="black",size=.4) +
    scale_y_discrete(expand=c(0,0), position = "right")+
    scale_x_discrete(expand=c(0,0), position = "bottom")+
    theme_grey(base_size=8)+
    labs(x="",y="",title="Fungi Core OTUs")+
    #scale_fill_gradient(low="blue", high='red', na.value="grey90") +
    scale_fill_gradientn(colours =colfunc(10), na.value="grey90")+
    # theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.3,size=7)) +
    guides(fill = guide_colorbar( barwidth = 0.5, barheight = 10, title="RA", frame.linetype = 1, frame.colour = "#000000"))+
    coord_fixed()+ #set base size for all font elements
    theme_grey(base_size=15)+
    #theme options
    theme(axis.text=element_text(face="bold"),
          #set thickness of axis ticks
          axis.ticks=element_line(size=1),
          #remove plot background
          plot.background=element_blank(),
          #remove plot border
          panel.border=element_blank()) +
    #theme options
    theme(legend.margin = grid::unit(0,"cm"),
          #change legend text properties
          legend.text=element_text(colour=textcol,size=7,face="bold"),
          #change legend key height
          legend.key.height=grid::unit(0.8,"cm"),
          #set a slim legend
          legend.key.width=grid::unit(0.2,"cm"),
          #set x axis text size and colour
          axis.text.x=element_text(colour=textcol, angle = 90,hjust=0.9,vjust=0.4, face="bold"),
          #set y axis text colour and adjust vertical justification
          axis.text.y=element_text(vjust = 0.2,colour=textcol, face="bold"),
          #change axis ticks thickness
          axis.ticks=element_line(size=0.4),
          #change title font, size, colour and justification
          plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
          #remove plot background
          plot.background=element_blank(),
          #remove plot border
          panel.border=element_blank())
}
colfunc.bac <- colorRampPalette(c("#FFFFFF","#0a6600"))
plot_core_heatmap(df.rel, colfunc.bac)


b.meta <- sample_data(bac.clean.ss)
b.meta <- data.frame(b.meta)


## (2) Fungi prevalence



# get_df_rel <- function(phy.clean.ss.5, b.dom){

df.otu.fun <- fun.its1.clean.ss %>% psmelt()

lb <- c("F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12",
        "F13","F14","F15","F16","F17","F18","F19","F20","F21","F22","F23",
        "F24","F25","F26","F27","F28","F29","F30","F36","F37","F38","F39",
        "F40","F41","F42","F43", "F44","F45","F46","F47","F48","F49","F50",
        "F51","F52","F53","F54","F55","F56","F57","F58","F59","F60","F61",
        "F62","F63","F64","F65","F71","F72","F73","F74","F75","F76","F77",
        "F78","F79","F81","F82","F83","F84","F87","F88","F89","F90","F91",
        "F93","F94","F96","F97","F98","F99","F100","F101","F102","F103","F104",
        "F105","F106","F107","F108","F109","F110","F111","F112","F113",
        "F114","F115","F116","F117","F118","F119","F120","F121","F122","F123",
        "F124","F125","F126","F128","F129","F130","F131","F132","F133","F134",
        "F135","F136","F137","F138","F139","F140","F141","F142","F143",
        "F144","F145")


df.otu.fun$Sample <- factor(df.otu.fun$Sample, levels = rev(lb))

# we need to group by samples
df.otu.rel.fun <- df.otu.fun %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.selected.rel.fun <- df.otu.rel.fun %>% filter(OTU %in% core.fun.its1.70)
excel.fun <- df.selected.rel.fun %>% dplyr::select(OTU,Sample,Phylum, Class, Order,Family,Genus,Species,Abundance, RelAbundance) %>%  group_by(OTU, RelAbundance) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

# install.packages('xlsx')
library(xlsx)
#write.xlsx(as.data.frame(excel.bac) , 'core_bacteria.xlsx')

f.add <- merge(excel.fun, fun.list.its1.OTU_id, by ="OTU")
f.order<- f.add %>% dplyr::select(OTU,total, OTU_id)%>% group_by(OTU_id)%>% summarise(Totals=sum(total)) %>% arrange(desc(Totals))
f.order.vec <- f.order$OTU_id

df.rel.fun <- df.selected.rel.fun %>% dplyr::select(OTU,Sample,RelAbundance,Abundance) %>% arrange(desc(Abundance))
df.rel.fun <- merge(df.rel.fun,fun.list.its1.OTU_id,by='OTU')
head(df.rel)

df.rel.fun$OTU_id <- factor(df.rel.fun$OTU_id, levels=f.order.vec) 
head(df.rel.fun)
str(df.rel.fun)

colfunc.fun <- colorRampPalette(c("#FFFFFF","#336633"))

plot_core_heatmap(df.rel.fun, colfunc.fun)

####CW-CA trial
b.meta <- sample_data(fun.CW.CA.ss)
b.meta <- data.frame(f.meta)
# get_df_rel <- function(phy.clean.ss.5, b.dom){
df.otu <- fun.CW.CA.ss %>% psmelt()
df.otu$Sample <- factor(df.otu$Sample, levels = rev(b.meta$SampleID))
head(df.otu)

# we need to group by samples
df.otu.rel <- df.otu %>%  
  group_by(Sample) %>%                         # Filter out at absolute read of 20       
  mutate(RelAbundance = Abundance*1000/sum(Abundance))  # Transform to rel. abundance

df.selected.rel <- df.otu.rel %>% filter(OTU %in% core.fun.95)
excel.bac <- df.selected.rel %>%   group_by(OTU, RelAbundance) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

# install.packages('xlsx')
library(xlsx)
#write.xlsx(as.data.frame(excel.bac) , 'core_bacteria.xlsx')

b.add <- merge(excel.bac, fun.list, by ="OTU")
b.order<- b.add %>% dplyr::select(OTU,total.x, OTU_id)%>% group_by(OTU_id)%>% summarise(Totals=sum(total.x)) %>% arrange(desc(Totals))
b.order.vec <- b.order$OTU_id

b.add <- as.data.frame(b.add)
df.rel <- df.selected.rel %>% dplyr::select(OTU,Sample,PlantSpecies,RelAbundance,Abundance) %>% arrange(desc(Abundance))
df.rel <- merge(df.rel,fun.list,by='OTU')
head(df.rel)

df.rel$OTU_id <- factor(df.rel$OTU_id, levels=b.order.vec) 
head(df.rel)
str(df.rel)
# df.rel$id_plus <- paste0(df.rel$id,' ',df.rel$Genus)
df.rel

colnames(df.rel)
#   return(df.selected.rel)
# }
plot_core_heatmap(df.rel, colfunc.bac)
