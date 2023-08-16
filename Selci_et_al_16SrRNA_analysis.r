library(tidyverse)
library(ggrepel)
library(phyloseq)
library(vegan)
library(viridis) #Viridis color library
library(microbiome) # data analysis and visualisation
library(plotly)
library("IRdisplay")
library(ggpmisc) #to use stat_poly_eq
library(ggpubr)
library("gridExtra")
library("factoextra")
library(rioja) # plotting poackages for tabular bubbleplots
library(reshape2)
library(RColorBrewer) # nice color options
library(randomcoloR) #generate random colors
library(ggthemes)
library(Biostrings)

# set size of the plot
options(repr.plot.width=14, repr.plot.height=12)

theme_glab <- function(base_size = 35,
                    base_family = "",
                    base_line_size = base_size / 180,
                    base_rect_size = base_size / 180) {
   
    font <- "Helvetica" #assign font family up front
   
    theme_bw(base_size = base_size,
                base_family = base_family,
                base_line_size = base_line_size) %+replace%
    theme(
        legend.background =  element_blank(),
        legend.title =       element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65),
                                         hjust = 0),
        legend.text =        element_text(color = rgb(100, 100, 100, maxColorValue = 255),
                                          size = rel(0.65)),
        legend.key.size =    unit(0.8, "lines"),
     
      plot.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        hjust = 0),
       
      axis.title = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
      axis.text = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.65)),
       
      plot.caption = element_text(
        color = rgb(100, 100, 100, maxColorValue = 255),
        size = rel(0.7),
        hjust = 1),
       
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      panel.border = element_rect(fill = NA, colour = rgb(100, 100, 100, maxColorValue = 255)),

     
      complete = TRUE
    )
}

path_df<-read.csv("../tables/Table1_sites_R.csv", header=TRUE, sep=",")

path_df

#Temp/pH samples distribution
ggplot(path_df, aes(x=temp,y=ph)) + 
geom_point(size=8,shape=21,aes(fill=type),stroke=.3) + 
#geom_text(aes(label=id)) +
scale_fill_manual(values = c("#2a7fffff","#ff5555ff","#66CD00"),labels=c("N","R","R/N")) +
labs(title="", x="Temperature (°C)", y="pH") + 
scale_x_continuous(breaks = c(10,20,30,40,50,60,70,80,90,100)) +
scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)) +
theme_glab()  + 
theme(aspect.ratio = 1/2)
#ggsave("../figures/figure1_ph_temp.svg",width=14,height=12)

#Phyloseq object

#otu_table 
otu_table_don<-read.csv("../16S_tables/otu_table.csv",header=T,sep=",",row.names=1)
otu_table_cr17<-otu_table(as.matrix(otu_table_don),taxa_are_rows = TRUE)

#tax_table
taxa_table_don<-read.csv("../16S_tables/taxonomy_table.csv",header=T,sep=",",row.names=1)
taxa_cr17<-tax_table(as.matrix(taxa_table_don))

#tree
bac_tree <- read_tree("../16S_tables/bac_tree.tre")

#refseq
bac_refseq <- readDNAStringSet("../16S_tables/bac_repset.fasta")

sample_namecr17<-read.csv("../16S_tables/cr17_sample.csv",header=T,sep=",",row.names=1)

#ASV with less than 5 reads were filtered
bac_data<-phyloseq(sample_data(sample_namecr17),
                   otu_table(otu_table_cr17),
                   tax_table(taxa_cr17), 
                   phy_tree(bac_tree),
                   refseq(bac_refseq))

bac_data #already normalized by the library median size 

## Subset and aggreate the dataset based on taxonomic level to ease computational work downstream
# Check the number of aggregated taxa
length(get_taxa_unique(bac_data, taxonomic.rank = "Genus"))
length(get_taxa_unique(bac_data, taxonomic.rank = "Family"))
length(get_taxa_unique(bac_data, taxonomic.rank = "Order"))
length(get_taxa_unique(bac_data, taxonomic.rank = "Class"))
length(get_taxa_unique(bac_data, taxonomic.rank = "Phyla"))

# Transform Bacteria abundance to relative abundances for plotting and some stats
bac_ra = transform_sample_counts(bac_data, function(x){x / sum(x)})

## Agglomerate at a specific taxonomic level at the Genus level
bac_ra_genus = tax_glom(bac_ra, "Genus", NArm = TRUE)
bac_ra_family = tax_glom(bac_ra, "Family", NArm = TRUE)
bac_ra_order = tax_glom(bac_ra, "Order", NArm = TRUE)
bac_ra_class = tax_glom(bac_ra, "Class", NArm = TRUE)
bac_ra_phyla = tax_glom(bac_ra, "Phyla", NArm = TRUE)

# Combined abundance plot of potential metabolic genera combined at the genus level
# Subset specific taxa involved in selected metabolisms for downstream analysis
bac_patho <- subset_taxa(bac_ra_genus, 
Genus == "Acinetobacter" |
Genus == "Pseudomonas" |
Genus == "Abiotrophia" |
Genus == "Achromobacter" |
Genus == "Acinetobacter" |
Genus == "Actinobacillus" |
Genus == "Arcanobacterium" |
Genus == "Arcobacter" |
Genus == "Babesia" |
Genus == "Bacillus" |
Genus == "Bartonella" |
Genus == "Bordetella" |
Genus == "Borrelia" |
Genus == "Brodetella" |
Genus == "Brucella" |
Genus == "Burkholderia" |
Genus == "Campylobacter" |
Genus == "Capnocytophaga" |
Genus == "Chlamydia" |
Genus == "Clostridium" |
Genus == "Comamonas" |
Genus == "Corynebacterium" |
Genus == "Coxiella" |
Genus == "Cronobacter" |
Genus == "Deinococcus" |
Genus == "Dermatophilus" |
Genus == "Ehrlichia" |
Genus == "Enterococcus" |
Genus == "Erysipelothrix" |
Genus == "Escherichia" |
Genus == "Escherichia/Shigella" |
Genus == "Flavobacterium" |
Genus == "Francisella" |
Genus == "Gardnerella" |
Genus == "Granulicatella" |
Genus == "Haemophilus" |
Genus == "Hafnia" |
Genus == "Halomonas" |
Genus == "Helicobacter" |
Genus == "Klebsiella" |
Genus == "Kocuria" |
Genus == "Lawsonia" |
Genus == "Legionella" |
Genus == "Leptospira" |
Genus == "Listeria" |
Genus == "Merkel_cell" |
Genus == "Micrococcus" |
Genus == "Morganella" |
Genus == "Mycobacterium" |
Genus == "Mycoplasma" |
Genus == "Neisseria" |
Genus == "Nocardia" |
Genus == "Pasteurella" |
Genus == "Photobacterium" |
Genus == "Plesiomonas" |
Genus == "Propionibacterium" |
Genus == "Proteus" |
Genus == "Providencia" |
Genus == "Pseudomonas" |
Genus == "Rhodococcus" |
Genus == "Rickettsiae" |
Genus == "Roseomonas" |
Genus == "Rothia" |
Genus == "Salmonella" |
Genus == "Serratia" |
Genus == "Shewanella" |
Genus == "Shigella" |
Genus == "Sphaerophorus" |
Genus == "Staphylococcus" |
Genus == "Stenotrophomonas" |
Genus == "Streptococcus" |
Genus == "Treponema" |
Genus == "Vibrio" |
Genus == "Yersinia"
                      ) # Subsetting only the putative pathogens

sample_data(bac_patho)

#excluding poas: PL and PG
bac_patho<-subset_samples(bac_patho,  ! name %in% "PG" )
bac_patho<-subset_samples(bac_patho,  ! name %in% "PL" )

bac_ra_genus
bac_patho


# Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <-subset_samples(bac_patho, type %in% "S")
summary_otu <- as.matrix((otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
write.csv(summary_comb, "../16S_tables/bac_patho_ra_S.csv")

save.image()

samples_order<-list("S1",
"S2",
"S3",
"S4",
"S5",
"S6",
"S7",
"S8",
"S9",
"S10",
"S11",
"S12",
"S13",
"S14",
"S15",
"S16",
"S17",
"S18",
"S19",
"S20",
"S21",
"S22")

subset_samples(bac_patho, type %in% "F")
subset_samples(bac_patho, type %in% "S")

# Putative Pathogens
path_f_s<-plot_bar(bac_patho, fill="Genus", x="id.1", facet_grid = ~type, title="") +
theme_glab() +
labs(x="",y="Relative Abundance (%)") +
theme(legend.position = "none") +
theme(plot.title = element_text(size=8),legend.position = "bottom",axis.text.x = element_text(angle = 90, 
vjust = 0.6, hjust=1.2),axis.text=element_text(size=12),
legend.key.size = unit(6, 'mm'),legend.text = element_text(size=12))

path_f_s$data$id.1 <- factor(path_f_s$data$id.1, levels = samples_order)

print(path_f_s)
#ggsave("../figures/barplot_F_S.svg", width=14, height=12)

# Combined abundance plot of potential metabolic genera combined at the genus level
# Subset specific taxa involved in selected metabolisms for downstream analysis
bac_patho2 <- subset_taxa(bac_data, 
Genus == "Acinetobacter" |
Genus == "Pseudomonas" |
Genus == "Abiotrophia" |
Genus == "Achromobacter" |
Genus == "Acinetobacter" |
Genus == "Actinobacillus" |
Genus == "Arcanobacterium" |
Genus == "Arcobacter" |
Genus == "Babesia" |
Genus == "Bacillus" |
Genus == "Bartonella" |
Genus == "Bordetella" |
Genus == "Borrelia" |
Genus == "Brodetella" |
Genus == "Brucella" |
Genus == "Burkholderia" |
Genus == "Campylobacter" |
Genus == "Capnocytophaga" |
Genus == "Chlamydia" |
Genus == "Clostridium" |
Genus == "Comamonas" |
Genus == "Corynebacterium" |
Genus == "Coxiella" |
Genus == "Cronobacter" |
Genus == "Deinococcus" |
Genus == "Dermatophilus" |
Genus == "Ehrlichia" |
Genus == "Enterococcus" |
Genus == "Erysipelothrix" |
Genus == "Escherichia" |
Genus == "Escherichia/Shigella" |
Genus == "Flavobacterium" |
Genus == "Francisella" |
Genus == "Gardnerella" |
Genus == "Granulicatella" |
Genus == "Haemophilus" |
Genus == "Hafnia" |
Genus == "Halomonas" |
Genus == "Helicobacter" |
Genus == "Klebsiella" |
Genus == "Kocuria" |
Genus == "Lawsonia" |
Genus == "Legionella" |
Genus == "Leptospira" |
Genus == "Listeria" |
Genus == "Merkel_cell" |
Genus == "Micrococcus" |
Genus == "Morganella" |
Genus == "Mycobacterium" |
Genus == "Mycoplasma" |
Genus == "Neisseria" |
Genus == "Nocardia" |
Genus == "Pasteurella" |
Genus == "Photobacterium" |
Genus == "Plesiomonas" |
Genus == "Propionibacterium" |
Genus == "Proteus" |
Genus == "Providencia" |
Genus == "Pseudomonas" |
Genus == "Rhodococcus" |
Genus == "Rickettsiae" |
Genus == "Roseomonas" |
Genus == "Rothia" |
Genus == "Salmonella" |
Genus == "Serratia" |
Genus == "Shewanella" |
Genus == "Shigella" |
Genus == "Sphaerophorus" |
Genus == "Staphylococcus" |
Genus == "Stenotrophomonas" |
Genus == "Streptococcus" |
Genus == "Treponema" |
Genus == "Vibrio" |
Genus == "Yersinia"
                      ) # Subsetting only the putative pathogens

save.image()

# Putative Pathogens - Fluids
path_gen_ra1<-plot_bar(subset_samples(bac_patho, type %in% "F"), fill="Genus", x="id.1") +
    theme_glab() +
    scale_y_continuous(labels=c(0,20,40,60,80)) +
    labs(x="",y="Relative Abundance") +
    theme(plot.title = element_text(size=8),legend.position = "right",axis.text.x = element_text(angle = 90, 
    vjust = 0.6, hjust=1.2),axis.text=element_text(size=12),
    legend.key.size = unit(6, 'mm'),legend.text = element_text(size=12)) +
    theme(legend.position = "bottom")

path_gen_ra1$data$id.1 <- factor(path_gen_ra1$data$id.1, levels = samples_order)

path_gen_ra1

# Transform Bacteria abundance to relative abundances for plotting and some stats
bac_ra_genus2 = transform_sample_counts(bac_patho, function(x){x / sum(x)})

# Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <-subset_samples(bac_ra_genus2, use %in% "N")
summary_otu <- as.matrix((otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
write.csv(summary_comb, "../16S_tables/bac_patho_ra_N.csv")

# Putative Pathogens
path_gen_ra2<-plot_bar(subset_samples(bac_ra_genus2, type %in% "F"), fill="Genus", x="id.1") +
    theme_glab() +
    scale_y_continuous(labels=c(0,25,50,75,100)) +
    labs(x="", y="Rleative Abundance (%)") +
    theme(plot.title = element_text(size=8),legend.position = "right",axis.text.x = element_text(angle = 90, 
    vjust = 0.6, hjust=1.2),axis.text=element_text(size=12),
    legend.key.size = unit(6, 'mm'),legend.text = element_text(size=12)) +
    theme(legend.position = "bottom")

path_gen_ra2$data$id.1 <- factor(path_gen_ra2$data$id.1, levels = samples_order)

path_gen_ra2

ggarrange(path_gen_ra1,path_gen_ra2, nrow=2, ncol=1, widths = c(2,2),
  heights = c(1,2), common.legend = TRUE)

#ggsave("../figures/barplot_order.svg", width=14, height=14)

# Saving a .csv with combined taxonomy and ASVs abundances for the selected taxonomic level
summary_object <-bac_ra_genus2
summary_otu <- as.matrix((otu_table(summary_object)))
summary_tax <- as.matrix(tax_table(summary_object))
summary_comb <- cbind(summary_otu, summary_tax)
write.csv(summary_comb, "../16S_tables/bac_patho_ra.csv")

bac_data_F<-subset_samples(bac_data, type %in% "F")

# Putative Pathogens - Sediments
path_gen_ra_sed<-plot_bar(subset_samples(bac_patho, type %in% "S"), fill="Genus", x="id.1") +
    theme_glab() +
    scale_y_continuous(labels=c(0,0.1,0.2,0.3,0.4)) +
    labs(x="",y="Relative Abundance (%)") +
    theme(plot.title = element_text(size=8),legend.position = "right",axis.text.x = element_text(angle = 90, 
    vjust = 0.6, hjust=1.2),axis.text=element_text(size=12),
    legend.key.size = unit(6, 'mm'),legend.text = element_text(size=12)) +
    theme(legend.position = "bottom")

path_gen_ra_sed$data$id.1 <- factor(path_gen_ra_sed$data$id.1, levels = samples_order)

path_gen_ra_sed
#ggsave("../figures/barplot_order_sediment.svg", width=14, height=12)

# Antibiograms results:

eucast_df<-read.csv("../antibiograms_stat/dataset/eucast.csv", header=T, sep=",")
clsi_df<-read.csv("../antibiograms_stat/dataset/clsi.csv", header=T, sep=",")
generic_df<-read.csv("../antibiograms_stat/dataset/generic.csv", header=T, sep=",")

eucast_df$antibiotics<-factor(eucast_df$antibiotics,levels=c("M","G","Cp","Cf","MDR"))
colnames(eucast_df)<- c("Antibiotic","S5","S4","S2","S3","PAO1")
eucast_df

#eucast_S5
s5_eucast<-ggplot(eucast_df, aes(x=Antibiotic,y=S5)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S5",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s5_eucast

#eucast_S4
s4_eucast<-ggplot(eucast_df, aes(x=Antibiotic,y=S4)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S4",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s4_eucast

#eucast_S2
s2_eucast<-ggplot(eucast_df, aes(x=Antibiotic,y=S2)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S2",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s2_eucast

#eucast_S3
s3_eucast<-ggplot(eucast_df, aes(x=Antibiotic,y=S3)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S3",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s3_eucast

#eucast_pao1
pao1_eucast<-ggplot(eucast_df, aes(x=Antibiotic,y=PAO1)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="PAO1",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

pao1_eucast

ggarrange(s5_eucast,s4_eucast,s2_eucast,s3_eucast,pao1_eucast,nrow=1,ncol=5, 
          common.legend=TRUE, legend="right") +
theme(aspect.ratio = 0.4) 
ggsave("../antibiograms_stat/eucast_MS.svg", width = 20, height = 10)

clsi_df$antibiotics<-factor(clsi_df$antibiotics,levels=c("M","G","Cp","Cf","MDR"))
colnames(clsi_df)<- c("Antibiotic","S5","S4","S2","S3","PAO1")
clsi_df

#clsi_S5
s5_clsi<-ggplot(clsi_df, aes(x=Antibiotic,y=S5)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S5",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s5_clsi

#clsi_S4
s4_clsi<-ggplot(clsi_df, aes(x=Antibiotic,y=S4)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S4",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s4_clsi

#clsi_S2
s2_clsi<-ggplot(clsi_df, aes(x=Antibiotic,y=S2)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S2",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s2_clsi

#clsi_S3
s3_clsi<-ggplot(clsi_df, aes(x=Antibiotic,y=S3)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S3",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s3_clsi

#clsi_pao1
pao1_clsi<-ggplot(clsi_df, aes(x=Antibiotic,y=PAO1)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="PAO1",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

pao1_clsi

ggarrange(s5_clsi,s4_clsi,s2_clsi,s3_clsi,pao1_clsi,nrow=1,ncol=5, 
          common.legend=TRUE, legend="right") +
theme(aspect.ratio = 0.4) 
#ggsave("../antibiograms_stat/clsi_MS.svg", width = 20, height = 10)

colnames(generic_df)<- c("Antibiotic","S5","S4","S2","S3","PAO1")
generic_df

#generic_S5
s5_generic<-ggplot(generic_df, aes(x=Antibiotic,y=S5)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S5",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s5_generic

#generic_S4
s4_generic<-ggplot(generic_df, aes(x=Antibiotic,y=S4)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S4",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s4_generic

#generic_S2
s2_generic<-ggplot(generic_df, aes(x=Antibiotic,y=S2)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S2",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s2_generic

#generic_S3
s3_generic<-ggplot(generic_df, aes(x=Antibiotic,y=S3)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="S3",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

s3_generic

#generic_pao1
pao1_generic<-ggplot(generic_df, aes(x=Antibiotic,y=PAO1)) +
geom_bar(stat="identity",aes(fill=Antibiotic), color="black") + 
labs(y="Antibiotic Resistant Strains (%)") + 
scale_fill_viridis_d() + 
coord_cartesian(ylim = c(0, 100)) + 
labs(title="",subtitle="PAO1",x="") + 
theme(axis.text.x=element_text(face="bold", size="15")) +
theme_glab(base_size = 25)

pao1_generic

ggarrange(s5_generic,s4_generic,s2_generic,s3_generic,pao1_generic,nrow=1,ncol=5, 
          common.legend=TRUE, legend="right") +
theme(aspect.ratio = 0.4)
#ggsave("../antibiograms_stat/generic_MS.svg", width = 20, height = 10)

library("gplots")
library("corrplot")
library(RColorBrewer)

eucast_matrix<-read.csv("../antibiograms_stat/dataset/eucast_antibiogram_matrix.csv",head=TRUE,sep=",",row.names=1)
eucast_matrix<-as.matrix(eucast_matrix)
head(eucast_matrix)

#svg("../figures/svg/eucast_corr.svg",width=8,height=16)
corrplot(eucast_matrix,method="square",tl.col="black", tl.srt=90,rect.lwd = .1,addgrid.col = "white",
        col=c("white", "black"))
#dev.off()

clsi_matrix<-read.csv("../antibiograms_stat/dataset/clsi_antibiogram_matrix.csv",head=TRUE,sep=",",row.names=1)
clsi_matrix<-as.matrix(clsi_matrix)
head(clsi_matrix)

#svg("../figures/svg/clsi_corr.svg",width=8,height=16)
corrplot(clsi_matrix,method="square",
                    tl.col="black", 
                    tl.srt=90,rect.lwd = .1,addgrid.col = "white",
         
        col=c("white", "black"))
#dev.off()


#Metagenomic data: Kaiju results - Genus level

kaiju_cr17<-read.csv("../kaiju/kaiju_env_data_cr17.csv", header=TRUE, sep=",", row.names=1)
head(kaiju_cr17)

cr17_genus<-readRDS("../../../cr17_cr18_metagenomes/kaiju/cr17_reads_based/kaiju_summaries/rds_otu_table/cr17_genus_otu_table.rds")

tax_tab_genus<-data.frame(row.names(t(cr17_genus)))
tax_tab_genus<-tax_table(tax_tab_genus)
colnames(tax_tab_genus)<-"Genus"
row.names(tax_tab_genus)<-tax_tab_genus[,1]

cr17_prok_genus<-phyloseq(sample_data(kaiju_cr17),
                   otu_table(cr17_genus,taxa_are_rows = F),
                   tax_table(tax_tab_genus))
cr17_prok_genus

# Combined abundance plot of potential metabolic genera combined at the genus level
# Subset specific taxa involved in selected metabolisms for downstream analysis
bac_patho_kaiju <- subset_taxa(cr17_prok_genus, 
Genus == "Acinetobacter" |
Genus == "Pseudomonas" |
Genus == "Abiotrophia" |
Genus == "Achromobacter" |
Genus == "Acinetobacter" |
Genus == "Actinobacillus" |
Genus == "Arcanobacterium" |
Genus == "Arcobacter" |
Genus == "Babesia" |
Genus == "Bacillus" |
Genus == "Bartonella" |
Genus == "Bordetella" |
Genus == "Borrelia" |
Genus == "Brodetella" |
Genus == "Brucella" |
Genus == "Burkholderia" |
Genus == "Campylobacter" |
Genus == "Capnocytophaga" |
Genus == "Chlamydia" |
Genus == "Clostridium" |
Genus == "Comamonas" |
Genus == "Corynebacterium" |
Genus == "Coxiella" |
Genus == "Cronobacter" |
Genus == "Deinococcus" |
Genus == "Dermatophilus" |
Genus == "Ehrlichia" |
Genus == "Enterococcus" |
Genus == "Erysipelothrix" |
Genus == "Escherichia" |
Genus == "Escherichia/Shigella" |
Genus == "Flavobacterium" |
Genus == "Francisella" |
Genus == "Gardnerella" |
Genus == "Granulicatella" |
Genus == "Haemophilus" |
Genus == "Hafnia" |
Genus == "Halomonas" |
Genus == "Helicobacter" |
Genus == "Klebsiella" |
Genus == "Kocuria" |
Genus == "Lawsonia" |
Genus == "Legionella" |
Genus == "Leptospira" |
Genus == "Listeria" |
Genus == "Merkel_cell" |
Genus == "Micrococcus" |
Genus == "Morganella" |
Genus == "Mycobacterium" |
Genus == "Mycoplasma" |
Genus == "Neisseria" |
Genus == "Nocardia" |
Genus == "Pasteurella" |
Genus == "Photobacterium" |
Genus == "Plesiomonas" |
Genus == "Propionibacterium" |
Genus == "Proteus" |
Genus == "Providencia" |
Genus == "Pseudomonas" |
Genus == "Rhodococcus" |
Genus == "Rickettsiae" |
Genus == "Roseomonas" |
Genus == "Rothia" |
Genus == "Salmonella" |
Genus == "Serratia" |
Genus == "Shewanella" |
Genus == "Shigella" |
Genus == "Sphaerophorus" |
Genus == "Staphylococcus" |
Genus == "Stenotrophomonas" |
Genus == "Streptococcus" |
Genus == "Treponema" |
Genus == "Vibrio" |
Genus == "Yersinia"
                      ) # Subsetting only the putative pathogens
bac_patho_kaiju

# Putative Pathogens
kaiju_path_gen<-
    plot_bar(subset_samples(bac_patho_kaiju, type %in% "F"), fill="Genus", x="site") +
    theme_glab() +
#   scale_y_continuous(labels=c(0,20,40,60,80)) +
    labs(x="",y="Relative Abundance (%)") +
    theme(plot.title = element_text(size=8),legend.position = "right",axis.text.x = element_text(angle = 90, 
    vjust = 0.6, hjust=1.2),axis.text=element_text(size=12),
    legend.key.size = unit(6, 'mm'),legend.text = element_text(size=12)) +
    theme(legend.position = "bottom")

kaiju_path_gen$data$site <- factor(kaiju_path_gen$data$site, levels = samples_order)

kaiju_path_gen
ggsave("../kaiju/kaiju_ra_genus.svg", width=14, height=10)

write.csv(otu_table(bac_patho_kaiju), "../kaiju/kaiju__path_summary.csv")

#Metagenomic data: Contigs blast against CARD 
card_samp<-read.csv("../card_blast/SAMPLE_DATA_BLAST_results_CR_ARGs.csv", header=TRUE, sep=",")
rownames(card_samp)<-card_samp[,1]
card_samp<-sample_data(card_samp)
card_samp

card_otu<-as.matrix(read.csv("../card_blast/OTU_TABLE_BLAST_results_CR_ARGs.csv", header=TRUE, sep=",", row.names=1))
card_otu<-otu_table(card_otu, taxa_are_rows=TRUE)
card_otu

card_tax<-as.matrix(read.csv("../card_blast/TAX_TABLE_BLAST_results_CR_ARGs.csv", header=TRUE, sep=","))
rownames(card_tax)<- card_tax[,1]
card_tax<-tax_table(card_tax)
card_tax

card_phy<-phyloseq(card_samp,card_otu,card_tax)
card_phy

# Transform Bacteria abundance to relative abundances for plotting and some stats
card_phy_ra = transform_sample_counts(card_phy, function(x){x / sum(x)})
card_phy_ra

sample_data(card_phy_ra)

sample_order2<-c("S2","S3","S4","S5","S6","S8","S13","S16","S18")

CustomPalette <- palette(brewer.pal(n = 10, name = "RdYlGn")) # Change the palette and number of steps accordingly. For more than 11 color switch to ColorRamp

library("paletteer")

length(table(tax_table(card_phy)[,"ARG"])) # Get the total number of taxa

# ARG absolute abudance
card_arg<-plot_bar(card_phy, fill="ARG",x="id1") +
    theme_glab() +
#   scale_fill_manual(values=palette) +
#   scale_y_continuous(labels=c(0,25,50,75,100)) +
    labs(x="",y="Absolute Abundance") +
    theme(plot.title = element_text(size=8),legend.position = "right",axis.text.x = element_text(angle = 90, 
    vjust = 0.6, hjust=1.2),axis.text=element_text(size=12),
    legend.key.size = unit(6, 'mm'),legend.text = element_text(size=12)) +
    theme(legend.position = "bottom")

card_arg$data$id1 <- factor(card_arg$data$id1, levels = sample_order2)

card_arg

# ARG RA
card_arg_ra<-plot_bar(card_phy_ra, fill="ARG",x="id1") +
    theme_glab() +
    scale_y_continuous(labels=c(0,25,50,75,100)) +
    labs(x="",y="Relative Abundance (%)") +
    theme(plot.title = element_text(size=8),legend.position = "right",axis.text.x = element_text(angle = 90, 
    vjust = 0.6, hjust=1.2),axis.text=element_text(size=12),
    legend.key.size = unit(6, 'mm'),legend.text = element_text(size=12)) +
    theme(legend.position = "bottom")

card_arg_ra$data$id1 <- factor(card_arg_ra$data$id1, levels = sample_order2)

card_arg_ra

ggarrange(card_arg,card_arg_ra, nrow=2, ncol=1, widths = c(2,2),
  heights = c(1,2), common.legend = TRUE)

#ggsave("../figures/card_blast_barplot.svg", width=14, height=14)

save.image()
