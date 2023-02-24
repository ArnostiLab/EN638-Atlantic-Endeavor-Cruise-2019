#Here we have the code that produces bar charts and bubble plots
#of bacterial community composition in amended and unamended mesocosms
#from each station and depth sampled during the 2019 Endeavor cruise (EN638).
#All of the community composition data examined here was compiled by Greta
#Giljan.


#Set the working directory to community_comp_bar_graphs_input

#Load packages
library(vegan)
library(maptools)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(dplyr)
library(plyr)
library(reshape2)
library(grid)
library(readxl)
library(colorspace)


#Load ggpattern and magick packages so that patterns can be added over colors
library(ggpattern)
library(magick)

#Load ggh4x to enable the facet_nested function, which allows graphs to have
#multiple, stacked facets
library(ggh4x)

###################################################################################                                                                                #
#                            Read in data for mesocosms                           #                                                                                #
###################################################################################

InputData_LV_dA_Stagged_5Perc <- read_xlsx("EN638_LV_dA_norm_5Perc_stagged.xlsx")
InputData_LV_dB_Stagged_5Perc <- read_xlsx("EN638_LV_dB_norm_5Perc_stagged.xlsx")
InputData_LV_dC_Stagged_5Perc <- read_xlsx("EN638_LV_dC_norm_5Perc_averraged_stagged.xlsx")

###################################################################################
#                       Create a manual color scale function                      #
###################################################################################

#This assigns a color to each taxa and is used in place of scale_fill_manual; 
#make sure all taxa from all data sets are included here so that this function can 
#be used to create all graphs 

scale_fill_bacteria <- function(...){
  ggplot2:::manual_scale(
    'fill', 
    values = setNames(c("red", "red1", "red2", "red3", "darkorange1", "darkorange2", "darkorange3", "red1", "red2", "red3", "darkorange1",
                        "darkorange2", "darkorange3", "red1", "red2", "red3", "red4", "darkorange1", "darkorange2", "darkorange3",  "red2",
                        "darkorange3", "red3", "red2", "darkorange4", "red1", "red2", "darkorange3", 	"darkorange4", "red1", "red3", "red1",
                        "darkorange1", "darkorange2", "darkorange3", "red1", "red3", "darkorange2", "red1", "darkorange1", "darkorange2",
                        "darkorange1",
                        
                        "cornflowerblue", "dodgerblue1", "dodgerblue3", "blue1", "blue2", "blue3", "dodgerblue1", "dodgerblue2", "dodgerblue3",
                        "blue1", "blue2", "blue3", "dodgerblue1", "dodgerblue2", "dodgerblue3", "dodgerblue4", "dodgerblue1", "dodgerblue2", 
                        "blue1", "blue2", "blue1", "blue4", "dodgerblue2", "dodgerblue3", "dodgerblue1", "blue2", "dodgerblue3", "dodgerblue1", 
                        
                        "yellow", "gold1", "gold2", "gold3", "yellow1", "yellow2", "yellow3", "gold1", "gold2", "gold3", "yellow1", "yellow2",
                        "yellow3", "gold1", "yellow3", "yellow2", "yellow4", "yellow3", "gold2", "yellow3", "gold2", "gold2",
                        
                        "darkblue",
                        
                        "hotpink", "hotpink3", "deeppink1",
                        
                        "yellowgreen", "chartreuse", "chartreuse3",
                        
                        "tan2",
                        
                        "turquoise1",
                        
                        "grey77", "grey47",
                        
                        "blueviolet", "darkorchid1",
                        
                        "khaki", "khaki3", 
                        
                        "chocolate", "chocolate4", "chocolate1", "chocolate3", 
                        
                        "springgreen1", "springgreen3", "springgreen4",
                        
                        "grey20"), 
                      
                      c("Gammaproteobacteria <5%", "Coxiella", "Alteromonas", "Salinimonas", "Colwellia", "Pseudoalteromonas", "Vibrio", "Marine Methylotrophic Group 3",
                        "Methylophaga", "Alcanivorax", "Pseudomonadales KI89A clade", "Marinobacter", "Amphritea", "Marinobacterium", "Oceanospirillum", "Oleiphilus",
                        "OM43 clade", "SAR86 clade", "Spongiibacter", "Zhongshania", "Thalassotalea", "Colwelliaceae uncultured", "Pseudohongiella", "Moritella",
                        "Psychrobium", "P13-46", "Saccharospirillaceae uncultured", "Sinobacterium", "SUP05 cluster", "Catenococcus", "Candidatus Endoecteinascidia", 
                        "Limnobacter", "Marinomonas", "Acinetobacter", "Pseudomonas", "Oleispira", "Thalassolituus", "KI89A clade", "Nitrincolaceae uncultured",
                        "Halomonas", "Woeseia", "Aliivibrio",
                        
                        "Alphaproteobacteria <5%", "Ascidiaceihabitans", "Planktomarina", "Micavibrionaceae uncultured", "Rhodobacteraceae", "Lacimonas", "Shimia",
                        "Sulfitobacter", "Thalassobius", "Rhodobacteraceae uncultured", "AEGEAN-169 marine group", "Rickettsiales S25-593", "SAR11 clade I", 
                        "SAR11 clade I uncultured", "SAR11 clade Ia", "SAR11 clade Ib", "SAR11 clade II", "SAR11 clade IV", "Amylibacter", "Dinoroseobacter",
                        "SAR116 clade", "Sagittulla", "SAR11 clade IV", "Methylobacterium-Methylorobrum", "Magnetospiraceae uncultured", "SD2D12", "SM2D12", "Methylorubrum",
                        
                        "Bacteroidota <5%", "Saprospiraceae uncultured", "Aurantivirga", "Dokdonia", "Flavicella", "Flavobacterium", "Leeuwenhoekiella",
                        "Lutibacter", "Maribacter", "Polaribacter", "Tenacibaculum", "Flavobacteriaceae uncultured", "Winogradskyella", "Cryomorphaceae uncultured",
                        "Formosa", "NS4 marine group", "NS5 marine group", "NS9 marine group", "Reichenbachiella", "Pseudofulvibacter", "Arenibacter", "Flavirhabdus",
                        
                        "SAR324 clade(Marine group B)",
                        
                        "Actinobacteria <5%", "Candidatus Actinomarina", "Sva0996 marine group",
                        
                        "Cyanobacteria <5%", "Prochlorococcus MIT9313", "Synechococcus CC9902",
                        
                        "Dadabacteriales",
                        
                        "Marinimicrobia (SAR406 clade)",
                        
                        "Chloroflexi <5%", "SAR202 clade",
                        
                        "Firmicutes <5%", "Paenibacillus",
                        
                        "Patescibacteria <5%", "JGI 0000069-P22",
                        
                        "Planctomycetes <5%", "CL500-3", "Phycisphaeraceae uncultured", "Blastopirellula",
                        
                        "Verrucomicrobia <5%", "DEV007", "Roseibacillus",
                        
                        "other bacteria <5%")), 
    ...
  )
}

###################################################################################
#                               Choosing patterns                                 #
###################################################################################

#This code prints a figure displaying possible magick patterns and their names that
#can be used in the following manual pattern scale function.
#From https://coolbutuseless.github.io/package/ggpattern/articles/pattern-magick.html

df1 <- data.frame(
  x    = rep(1:6, 9),
  y    = rep(1:9, each=6),
  name = gridpattern::names_magick,
  stringsAsFactors = FALSE
)


ggplot(df1) + 
  geom_tile_pattern(
    aes(x, y, pattern_type = I(name)),
    pattern       = 'magick',
    pattern_scale = 1.5,
    pattern_fill  = 'black', 
    width         = 0.9, 
    height        = 0.9
  ) + 
  geom_label(aes(x+0.4, y+0.4, label = name), hjust = 1, vjust = 1) + 
  theme_void() + 
  labs(
    title = "All the possible magick pattern names"
  ) +
  coord_fixed(1)

###################################################################################
#                     Create a manual pattern scale function                      #
###################################################################################

#This assigns a pattern to each taxa and is used in place of 
#scale_pattern_type_manual; make sure all taxa from all data sets are included here
#so that this function can be used to create all graphs. The value 'gray100' leaves
#an assigned taxa blank.

scale_pattern_type_bacteria <- function(...){
  ggpattern:::manual_scale(
    'pattern_type',
    values = setNames(c('gray100', 'gray100', 'vertical', 'gray100', 'bricks', 'verticalsaw', 'crosshatch', 'verticalrightshingle', 
                        'gray100', 'horizontal3', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100',
                        'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'hs_diagcross',
                        'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 
                        'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100',
                        'gray100', 'gray100', 'gray100',
                        
                        'horizontalsaw', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100',
                        'rightshingle', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 
                        'gray100', 'hexagons', 'circles', 'gray100', 'gray100', 'gray100', 'gray100', 
                        'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100',
                        
                        'verticalbricks', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 
                        'gray100', 'gray100', 'smallfishscales', 'left45', 'right45', 'gray100', 'gray100',
                        'verticalsaw', 'gray100', 'gray100', 'gray100', 'gray100', 'right30', 'gray100', 'gray100',
                        
                        'gray100',
                        
                        'gray100', 'gray100', 'gray100',
                        
                        'gray100', 'hs_fdiagonal', 'gray100', 
                        
                        'gray100',
                        
                        'gray100',
                        
                        'gray100', 'hs_cross',
                        
                        'gray100', 'gray100',
                        
                        'gray100', 'gray100',
                        
                        'gray100', 'gray100', 'gray100', 'gray100',
                        
                        'gray100', 'vertical2', 'gray100',
                        
                        'gray100'),
                      
                      c("Gammaproteobacteria <5%", "Coxiella", "Alteromonas", "Salinimonas", "Colwellia", "Pseudoalteromonas", "Vibrio", "Marine Methylotrophic Group 3",
                        "Methylophaga", "Alcanivorax", "Pseudomonadales KI89A clade", "Marinobacter", "Amphritea", "Marinobacterium", "Oceanospirillum", "Oleiphilus",
                        "OM43 clade", "SAR86 clade", "Spongiibacter", "Zhongshania", "Thalassotalea", "Colwelliaceae uncultured", "Pseudohongiella", "Moritella",
                        "Psychrobium", "P13-46", "Saccharospirillaceae uncultured", "Sinobacterium", "SUP05 cluster", "Catenococcus", "Candidatus Endoecteinascidia", 
                        "Limnobacter", "Marinomonas", "Acinetobacter", "Pseudomonas", "Oleispira", "Thalassolituus", "KI89A clade", "Nitrincolaceae uncultured",
                        "Halomonas", "Woeseia", "Aliivibrio",
                        
                        "Alphaproteobacteria <5%", "Ascidiaceihabitans", "Planktomarina", "Micavibrionaceae uncultured", "Rhodobacteraceae", "Lacimonas", "Shimia",
                        "Sulfitobacter", "Thalassobius", "Rhodobacteraceae uncultured", "AEGEAN-169 marine group", "Rickettsiales S25-593", "SAR11 clade I", 
                        "SAR11 clade I uncultured", "SAR11 clade Ia", "SAR11 clade Ib", "SAR11 clade II", "SAR11 clade IV", "Amylibacter", "Dinoroseobacter",
                        "SAR116 clade", "Sagittulla", "SAR11 clade IV", "Methylobacterium-Methylorobrum", "Magnetospiraceae uncultured", "SD2D12", "SM2D12", "Methylorubrum",
                        
                        "Bacteroidota <5%", "Saprospiraceae uncultured", "Aurantivirga", "Dokdonia", "Flavicella", "Flavobacterium", "Leeuwenhoekiella",
                        "Lutibacter", "Maribacter", "Polaribacter", "Tenacibaculum", "Flavobacteriaceae uncultured", "Winogradskyella", "Cryomorphaceae uncultured",
                        "Formosa", "NS4 marine group", "NS5 marine group", "NS9 marine group", "Reichenbachiella", "Pseudofulvibacter", "Arenibacter", "Flavirhabdus",
                        
                        "SAR324 clade(Marine group B)",
                        
                        "Actinobacteria <5%", "Candidatus Actinomarina", "Sva0996 marine group",
                        
                        "Cyanobacteria <5%", "Prochlorococcus MIT9313", "Synechococcus CC9902",
                        
                        "Dadabacteriales",
                        
                        "Marinimicrobia (SAR406 clade)",
                        
                        "Chloroflexi <5%", "SAR202 clade",
                        
                        "Firmicutes <5%", "Paenibacillus",
                        
                        "Patescibacteria <5%", "JGI 0000069-P22",
                        
                        "Planctomycetes <5%", "CL500-3", "Phycisphaeraceae uncultured", "Blastopirellula",
                        
                        "Verrucomicrobia <5%", "DEV007", "Roseibacillus",
                        
                        "other bacteria <5%")),
    
    ...
  )
}

###################################################################################
#            Create barchart - Amended dA (DCM) 5% Genus level                    #
###################################################################################

#Factor data and define levels
InputData_LV_dA_Stagged_5Perc$IncubationTime <- factor(InputData_LV_dA_Stagged_5Perc$IncubationTime, levels=c("0", "1", "2", "3", "4"), labels=c("0", "3", "5", "10", "15"))
InputData_LV_dA_Stagged_5Perc$Substrate_f <- factor(InputData_LV_dA_Stagged_5Perc$Substrate, levels=c("Pullulan", "Laminarin", "Xylan", "Fucoidan", "Arabinogalactan", "Chondroitin", "Control"), 
                                                                                             labels=c("Pullulan", "Laminarin", "Xylan", "Fucoidan", "Arabinogalactan", "Chondroitin", "Amended Meso"))
InputData_LV_dA_Stagged_5Perc$TaxonShort_f <- factor(InputData_LV_dA_Stagged_5Perc$TaxonShort, levels=c("Gammaproteobacteria", "Gammaproteobacteria <5%",
                                                                                                        "Alteromonas", 
                                                                                                        "Colwellia", "Colwelliaceae uncultured", "Moritella", "Pseudoalteromonas", "Vibrio",
                                                                                                        "Marine Methylotrophic Group 3", "Methylophaga", "Alcanivorax",
                                                                                                        "KI89A clade", "Marinobacter", "Acinetobacter", "Amphritea", 
                                                                                                        "Nitrincolaceae uncultured", "Pseudomonas", "Oleispira",
                                                                                                        "OM43 clade",
                                                                                                        "SAR86 clade", 
                                                                                                        
                                                                                                        "Alphaproteobacteria", "Alphaproteobacteria <5%",
                                                                                                        "Micavibrionaceae uncultured", "Amylibacter", "Ascidiaceihabitans", "Lacimonas", 
                                                                                                        "Planktomarina", "Shimia", "Sulfitobacter",  
                                                                                                        "Rhodobacteraceae uncultured", "AEGEAN-169 marine group", 
                                                                                                        "SAR11 clade I uncultured", "SAR11 clade Ia", "SAR11 clade Ib", "SAR11 clade II",
                                                                                                        
                                                                                                        "Bacteroidota", "Bacteroidota <5%",
                                                                                                        "Saprospiraceae uncultured", "Arenibacter",
                                                                                                        "Flavicella", "Flavobacterium", "Formosa", 
                                                                                                        "Maribacter", "NS5 marine group", "Polaribacter", "Pseudofulvibacter",
                                                                                                        "Tenacibaculum", "Flavobacteriaceae uncultured", 
                                                                                                        
                                                                                                        "SAR324 clade", "SAR324 clade(Marine group B)", 
                                                                                                        
                                                                                                        "Actinobacteria", "Actinobacteria <5%",
                                                                                                        "Candidatus Actinomarina",  
                                                                                                        
                                                                                                        "Cyanobacteria", "Cyanobacteria <5%",
                                                                                                        "Prochlorococcus MIT9313", "Synechococcus CC9902", 
                                                                                                        
                                                                                                        "Chloroflexi", "Chloroflexi <5%",
                                                                                                        "SAR202 clade",
                                                                                                        
                                                                                                        "Verrucomicrobia", "Verrucomicrobia <5%",
                                                                                                        "DEV007", "Roseibacillus",
                                                                                                        
                                                                                                        "other bacteria <5%"))


InputData_LV_dA_Stagged_5Perc$Station_f <- factor(InputData_LV_dA_Stagged_5Perc$Station, levels = c("18", "19", "20"))
InputData_LV_dA_Stagged_5Perc$Depth_f <- factor(InputData_LV_dA_Stagged_5Perc$Depth, levels = c("A"), labels=c("DCM"))

#Filter data - keep only Amended Mesos
InputData_LV_dA_Meso_Filt <- InputData_LV_dA_Stagged_5Perc %>% 
  filter(Substrate_f %in% c("Amended Meso"))

#Relabel Mesocosm numbers to Amend1, Amend2, Amend3, and Unamend
InputData_LV_dA_Meso_Filt$Mesocosm_f <- factor(InputData_LV_dA_Meso_Filt$Mesocosm, levels=c("1", "2", "3", "4"), labels=c("Amend1", "Amend2", "Amend3", "Unamend"))

#Creating the bar chart with the colors and patterns overlapping
Barchart_LV_Meso_stagged_dA_Genus_5Perc <- ggplot(InputData_LV_dA_Meso_Filt, aes(x=IncubationTime, y=ReadAbundance, fill=TaxonShort_f)) +
  
  geom_bar_pattern(position = "stack",
                   stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_spacing = 0.05,
                   pattern = 'magick',
                   aes(pattern_type = TaxonShort_f)) +
  scale_fill_bacteria(guide = "none") +
  scale_pattern_type_bacteria(guide = "none") +
  scale_y_continuous(limits=c(0,101)) +
  facet_grid(vars(Station_f), vars(Mesocosm_f)) + 
  theme_test(base_size = 16) +
  labs(x = "Incubation time [d]", y = "Relative read abundance [%]")

print(Barchart_LV_Meso_stagged_dA_Genus_5Perc)

#Saving as an editable pdf
ggsave(filename="EN638_LV_Mesos_stagged_dA_Genus_5Perc_with_patterns.pdf", plot=Barchart_LV_Meso_stagged_dA_Genus_5Perc, width=12, height=10, device=cairo_pdf)

#For use when you want a legend:
#scale_fill_bacteria(guide = guide_legend(name = "taxa")) +
#scale_pattern_type_bacteria(guide = guide_legend(override.aes = list(fill = "white")))

#Note that the control here is actually samples directly from the 20 L amended mesocosm, not from a smaller incubation run in parallel with the 
#FLA-PS incubations that contains no substrate. 

#The controls (which were relabeled to Amended Mesos) have timepoints of 0, 3, 5, 10, 15 (no 15 d sample for stn 18 dA). FLA-PS incubations,
#which are fluorescently-labeled polysaccharide incubations created with water sampled from Mesocosm Amend1, have timepoints of 
#This means that the timepoints for the Control here are also different than those for the
#2, 5, 7, 12, 17

###################################################################################
#       Create barchart - Amended dB (Oxygen minimum) 5% Genus level              #
###################################################################################
#Factor data and define levels
InputData_LV_dB_Stagged_5Perc$IncubationTime <- factor(InputData_LV_dB_Stagged_5Perc$IncubationTime, levels=c("0", "1", "2", "3", "4"), labels=c("0", "3", "5", "10", "15"))
InputData_LV_dB_Stagged_5Perc$TaxonShort_f <- factor(InputData_LV_dB_Stagged_5Perc$TaxonShort, levels=c("Gammaproteobacteria <5%", "Coxiella",
                                                                                                        "Alteromonas", 
                                                                                                        "Colwellia", "Colwelliaceae uncultured", "Moritella", "Pseudoalteromonas", "Vibrio",
                                                                                                        "Alcanivorax",
                                                                                                        "KI89A clade", "Marinobacter",  
                                                                                                        "Nitrincolaceae uncultured", "SUP05 cluster", "Oleispira",
                                                                                                        "SAR86 clade", "Marinomonas", "Halomonas", "Woeseia",
                                                                                                        
                                                                                                        "Alphaproteobacteria <5%",
                                                                                                        "Amylibacter", 
                                                                                                        "Shimia",  
                                                                                                        "Rhodobacteraceae uncultured", 
                                                                                                        "SAR11 clade I uncultured", "SAR11 clade Ia", "SAR11 clade Ib", "SAR11 clade II",
                                                                                                        "SAR11 clade IV", "Magnetospiraceae uncultured",
                                                                                                        
                                                                                                        "Bacteroidota <5%", "Aurantivirga",
                                                                                                        "NS5 marine group", "Polaribacter", "Pseudofulvibacter",
                                                                                                        "Tenacibaculum", "Flavobacteriaceae uncultured", "Flavirhabdus",
                                                                                                        
                                                                                                        "SAR324 clade", "SAR324 clade(Marine group B)", 
                                                                                                        
                                                                                                        "Actinobacteria <5%", "Sva0996 marine group",
                                                                                                        
                                                                                                        "Marinimicrobia (SAR406 clade)",
                                                                                                        
                                                                                                        "Chloroflexi <5%",
                                                                                                        "SAR202 clade",
                                                                                                        
                                                                                                        "Firmicutes <5%", "Paenibacillus",
                                                                                                        
                                                                                                        "Patescibacteria <5%", "JGI 0000069-P22",
                                                                                                        
                                                                                                        "Planctomycetes <5%", "Phycisphaeraceae uncultured", "Blastopirellula",
                                                                                                       
                                                                                                        "Verrucomicrobia <5%",
                                                                                                        "DEV007", "Roseibacillus",
                                                                                                        
                                                                                                        "other bacteria <5%"))


InputData_LV_dB_Stagged_5Perc$Station_f <- factor(InputData_LV_dB_Stagged_5Perc$Station, levels = c("18", "19", "20"))
InputData_LV_dB_Stagged_5Perc$Depth_f <- factor(InputData_LV_dB_Stagged_5Perc$Depth, levels = c("B"), labels=c("O2 Min"))

#Using dplyr to relabel Bacteroidetes (it should be Bacteroidota) and SAR324 (it should be SAR324 clade)
#This way the taxonomic names used is consistent across all data sets, which is important when we create the legend and bubble plot later in the code
InputData_LV_dB_Stagged_5Perc <- InputData_LV_dB_Stagged_5Perc %>%
  mutate(TaxonShort_f = recode(TaxonShort_f, "Bacteroidetes" = 'Bacteroidota', "Bacteroidetes <5%" = 'Bacteroidota <5%'))

InputData_LV_dB_Stagged_5Perc <- InputData_LV_dB_Stagged_5Perc %>%
  mutate(Phylum = recode(Phylum, "Bacteroidetes" = 'Bacteroidota'))

InputData_LV_dB_Stagged_5Perc <- InputData_LV_dB_Stagged_5Perc %>%
  mutate(Phylum = recode(Phylum, "SAR324" = 'SAR324 clade'))


#Relabel Mesocosm numbers to Amend1, Amend2, Amend3, and Unamend
InputData_LV_dB_Stagged_5Perc$Mesocosm_f <- factor(InputData_LV_dB_Stagged_5Perc$Mesocosm, levels=c("1", "2", "3", "4"), labels=c("Amend1", "Amend2", "Amend3", "Unamend"))

#Creating the bar chart with the colors and patterns overlapping
Barchart_LV_Meso_stagged_dB_Genus_5Perc <- ggplot(InputData_LV_dB_Stagged_5Perc, aes(x=IncubationTime, y=ReadAbundance, fill=TaxonShort_f)) +
  
  geom_bar_pattern(position = "stack",
                   stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_spacing = 0.05,
                   pattern = 'magick',
                   aes(pattern_type = TaxonShort_f)) +
  scale_fill_bacteria(guide = "none") +
  scale_pattern_type_bacteria(guide = "none") +
  scale_y_continuous(limits=c(0,101)) +
  facet_grid(vars(Station_f), vars(Mesocosm_f)) + 
  theme_test(base_size = 16) +
  labs(x = "Incubation time [d]", y = "Relative read abundance [%]")

print(Barchart_LV_Meso_stagged_dB_Genus_5Perc)

#Saving as an editable pdf
ggsave(filename="EN638_LV_Mesos_stagged_dB_Genus_5Perc_with_patterns.pdf", plot=Barchart_LV_Meso_stagged_dB_Genus_5Perc, width=12, height=10, device=cairo_pdf)

#For use when you want a legend:
#scale_fill_bacteria(guide = guide_legend(name = "taxa")) +
#scale_pattern_type_bacteria(guide = guide_legend(override.aes = list(fill = "white")))

#Note that the control here is actually samples directly from the 20 L amended mesocosm, not from a smaller incubation run in parallel with the 
#FLA-PS incubations that contains no substrate. 

#The controls (which were relabeled to Amended Mesos) have timepoints of 0, 3, 5, 10, 15 (no 15 d sample for stn 18 dA). FLA-PS incubations,
#which are fluorescently-labeled polysaccharide incubations created with water sampled from Mesocosm Amend1, have timepoints of 
#This means that the timepoints for the Control here are also different than those for the
#2, 5, 7, 12, 17

#No FLA-PS incubations for this depth, only mesocosms

###################################################################################
#           Create barchart - Amended dC (Bottom) 5% Genus level                  #
###################################################################################

#Factor data and define levels
InputData_LV_dC_Stagged_5Perc$IncubationTime <- factor(InputData_LV_dC_Stagged_5Perc$IncubationTime, levels=c("0", "1", "2", "3", "4"), labels=c("0", "3", "5", "10", "15"))
InputData_LV_dC_Stagged_5Perc$Substrate_f <- factor(InputData_LV_dC_Stagged_5Perc$Substrate, levels=c("Pullulan", "Laminarin", "Xylan", "Fucoidan", "Arabinogalactan", "Chondroitin", "Control"),
                                                                                             labels=c("Pullulan", "Laminarin", "Xylan", "Fucoidan", "Arabinogalactan", "Chondroitin", "Amended Meso"))
InputData_LV_dC_Stagged_5Perc$TaxonShort_f <- factor(InputData_LV_dC_Stagged_5Perc$TaxonShort, levels=c("Gammaproteobacteria", "Gammaproteobacteria <5%",
                                                                                                        "Alteromonas", 
                                                                                                        "Colwellia", "Colwelliaceae uncultured", "Moritella", "Pseudoalteromonas", 
                                                                                                        "Psychrobium", "Aliivibrio", "Vibrio", "Methylophaga",
                                                                                                        "Alcanivorax", "Marinobacter", "Marinomonas",  
                                                                                                        "Nitrincolaceae uncultured", "Oleiphilus", "Pseudomonas", "Oleispira",
                                                                                                        
                                                                                                        "Alphaproteobacteria", "Alphaproteobacteria <5%",
                                                                                                        "Sulfitobacter", "SM2D12", 
                                                                                                        "SAR11 clade I uncultured", "SAR11 clade Ia", "SAR11 clade Ib", "SAR11 clade II", 
                                                                                                        
                                                                                                        "Bacteroidetes", "Bacteroidetes <5%",
                                                                                                        "Reichenbachiella", "Arenibacter", "Flavirhabdus",
                                                                                                        "Polaribacter", "Pseudofulvibacter",
                                                                                                        "Tenacibaculum", "Flavobacteriaceae uncultured", 
                                                                                                        
                                                                                                        "SAR324 clade", "SAR324 clade(Marine group B)",
                                                                                                        
                                                                                                        "Marinimicrobia", "Marinimicrobia (SAR406 clade)", 
                                                                                                        
                                                                                                        "Chloroflexi", "Chloroflexi <5%",
                                                                                                        "SAR202 clade",
                                                                                                        
                                                                                                        "Verrucomicrobia", "Verrucomicrobia <5%",
                                                                                                        "Roseibacillus", 
                                                                                                        
                                                                                                        "other bacteria <5%"))

InputData_LV_dC_Stagged_5Perc$Station_f <- factor(InputData_LV_dC_Stagged_5Perc$Station, levels = c("18", "19", "20"))
InputData_LV_dC_Stagged_5Perc$Depth_f <- factor(InputData_LV_dC_Stagged_5Perc$Depth, levels = c("C"), labels=c("Bottom"))

#Using dplyr to relabel Bacteroidetes (it should be Bacteroidota), Marinimonas (it should be Marinimicrobia), and SAR324 (it should be SAR324 clade)
#This way the taxonomic names used is consistent across all data sets, which is important when we create the legend and bubble plot later in the code
InputData_LV_dC_Stagged_5Perc <- InputData_LV_dC_Stagged_5Perc %>%
  mutate(TaxonShort_f = recode(TaxonShort_f, "Bacteroidetes" = 'Bacteroidota', "Bacteroidetes <5%" = 'Bacteroidota <5%'))

InputData_LV_dC_Stagged_5Perc <- InputData_LV_dC_Stagged_5Perc %>%
  mutate(Phylum = recode(Phylum, "Bacteroidetes" = 'Bacteroidota'))

InputData_LV_dC_Stagged_5Perc <- InputData_LV_dC_Stagged_5Perc %>%
  mutate(Phylum = recode(Phylum, "Marinimonas" = 'Marinimicrobia'))

InputData_LV_dC_Stagged_5Perc <- InputData_LV_dC_Stagged_5Perc %>%
  mutate(Phylum = recode(Phylum, "SAR324" = 'SAR324 clade'))

#Filter data - keep only Amended Mesos
InputData_LV_dC_Meso_Filt <- InputData_LV_dC_Stagged_5Perc %>% 
  filter(Substrate_f %in% c("Amended Meso"))

#Relabel Mesocosm numbers to Amend1, Amend2, Amend3, and Unamend
InputData_LV_dC_Meso_Filt$Mesocosm_f <- factor(InputData_LV_dC_Meso_Filt$Mesocosm, levels=c("1", "2", "3", "4"), labels=c("Amend1", "Amend2", "Amend3", "Unamend"))

#Creating the bar chart with the colors and patterns overlapping
Barchart_LV_Meso_stagged_dC_Genus_5Perc <- ggplot(InputData_LV_dC_Meso_Filt, aes(x=IncubationTime, y=ReadAbundance, fill=TaxonShort_f)) +
  
  geom_bar_pattern(position = "stack",
                   stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_spacing = 0.05,
                   pattern = 'magick',
                   aes(pattern_type = TaxonShort_f)) +
  scale_fill_bacteria(guide = "none") +
  scale_pattern_type_bacteria(guide = "none") +
  scale_y_continuous(limits=c(0,101)) +
  facet_grid(vars(Station_f), vars(Mesocosm_f)) + 
  theme_test(base_size = 16) +
  labs(x = "Incubation time [d]", y = "Relative read abundance [%]")

print(Barchart_LV_Meso_stagged_dC_Genus_5Perc)

#Saving as an editable pdf
ggsave(filename="EN638_LV_Mesos_stagged_dC_Genus_5Perc_with_patterns.pdf", plot=Barchart_LV_Meso_stagged_dC_Genus_5Perc, width=12, height=10, device=cairo_pdf)


#Note that the control here is actually samples directly from the 20 L amended mesocosm, not from a smaller incubation run in parallel with the 
#FLA-PS incubations that contains no substrate. This means that the timepoints for the Control here are also different than those for the
#FLA-PS incubations - they should be 0, 3, 5, 10, 15 


###################################################################################
#    Create bubble plot - Unamended mesocosms and averaged amended mesocosms      #
###################################################################################

#The above sections of code should be run before this section to ensure the data is
#correctly labeled. Unlike the above sections, this code will create a bubble plot
#rather than a barchart.

#Remove extra columns not included in all data sets 
Amend_dA_filt <- subset(InputData_LV_dA_Meso_Filt, select=-c(Substrate, Substrate_f))
Amend_dC_filt <- subset(InputData_LV_dC_Meso_Filt, select=-c(Substrate, Substrate_f))

#Combining the subsets above
Total_community_comp <- rbind(Amend_dA_filt, InputData_LV_dB_Stagged_5Perc, Amend_dC_filt)

#Factor and define levels for Mesocosm_f
Total_community_comp$Mesocosm_f <- factor(Total_community_comp$Mesocosm_f, levels = c("Amend1", "Amend2", "Amend3", "Unamend"))

#Creating a Mesocosm_type column that lists whether the mesocosm was amended or unamended
Total_community_comp_mean <- Total_community_comp %>% mutate(Mesocosm_type =
                                 case_when(Mesocosm_f == "Amend1" ~ "Amended", 
                                           Mesocosm_f == "Amend2" ~ "Amended",
                                           Mesocosm_f == "Amend3" ~ "Amended",
                                           Mesocosm_f == "Unamend" ~ "Unamended")
)


#Create an average across the three amended mesocosms from each station/depth using dplyr
Total_community_comp_mean <- Total_community_comp_mean %>% group_by(Station_f, Depth_f, Mesocosm_type, IncubationTime, Phylum, TaxonShort_f) %>% 
  summarise_at("ReadAbundance", mean)

#Remove TaxonShort_f with Abundance lower than 5%
Total_community_comp_mean_filt <- subset(Total_community_comp_mean, ReadAbundance>=5)

#Factor and define Mesocosm_type
Total_community_comp_mean_filt$Mesocosm_type <- factor(Total_community_comp_mean_filt$Mesocosm_type, levels=c("Amended", "Unamended"))

#Choosing a color scheme
#colours = c("#F4DE52FF", "#2BA872", "#421C57")
#lam_color = c("#FFC425", "grey60", "grey60")
Meso_colors = c("grey60", "#FFC425")

#Make the bubble plot - first, group the families by Phylum
xx = ggplot(Total_community_comp_mean_filt, aes(x=IncubationTime, y = TaxonShort_f)) + 
  geom_point(aes(size = ReadAbundance, fill = Depth_f), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(5, 100), range = c(1,4), breaks = c(5,25,50,75)) + 
  labs(x= "Time (days)", y = "", size = "Relative Abundance", fill = "")  + 
  theme_classic() +
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 5, angle = 0, vjust = 0.3, hjust = 0.5), 
        axis.text.y = element_text(colour = "black", size = 5), 
        axis.title.x = element_text(colour = "black", size = 7),
        axis.line = element_line(colour = "black", size = 0.2),
        axis.ticks = element_line(size = 0.1),
        legend.text = element_text(size = 6, face ="bold", colour ="black"), 
        legend.title = element_text(size = 6, face = "bold"), 
        legend.box = "horizontal",
        strip.text = element_text(size=6),
        panel.background = element_blank(), panel.border = element_rect(colour = NA, fill = NA, size = 0.3), 
        panel.grid.major.y = element_line(color="gray70", size=0.1),
        panel.spacing = unit(0, "mm"), legend.position = "right") +  
  scale_fill_discrete_qualitative(palette = "Harmonic", guide = "none") +
  scale_y_discrete(limits = rev(levels(Total_community_comp$ReadAbundance)))+
  facet_nested(Phylum ~ Mesocosm_type + Depth_f + Station_f, scales="free_y", space="free_y", switch="y")


#Phylum MUST be included in facet_nested before the ~
#Otherwise the code combines lower taxonomic levels with the same names even if they're in different phylums
#For example, the SAR202 (Class) ambiguous taxa (family) and the Ambiguous taxa (class) ambiguous taxa (family) end up combined


#Changing legend placement, line width, etc
xy <- xx + 
  theme(strip.placement = "outside",                                # Place facet labels outside y axis labels.
        strip.background = element_rect(linewidth=0.1, fill = "white"),  # Make facet label background white.
        legend.position = "bottom", legend.justification = "right") # Move legend position                 

print(xy)

#Export
ggsave(filename="Bubbleplot_abundance_5percent.pdf", plot=xy, width=10, height=7, device=cairo_pdf)


#Fixed orientation of Phylum names, added shading, and changed the time points to the correct number of days in Illustrator.

###################################################################################
#       Other options for graphing overlapping colors and patterns                #
###################################################################################

fill <-                                                             c("grey99", "red",
                                                                      "red1", "red2", "red3",
                                                                      "darkorange1", "darkorange2", "darkorange3", 
                                                                      "red1", "red2", "red3", 
                                                                      "darkorange1", "darkorange2", "darkorange3", 
                                                                      "red1", "red2", "red3",
                                                                      "darkorange1", "darkorange2", "darkorange3",
                                                                      
                                                                      "grey99", "cornflowerblue",
                                                                      "blue1", "blue2", "blue3",
                                                                      "dodgerblue1", "dodgerblue2", "dodgerblue3", 
                                                                      "blue1", "blue2", "blue3",
                                                                      "dodgerblue1", "dodgerblue2", "dodgerblue3", "dodgerblue4", 
                                                                      
                                                                      "grey99", "yellow", 
                                                                      "gold1","gold2","gold3", 
                                                                      "yellow1", "yellow2", "yellow3",
                                                                      "gold1","gold2","gold3", 
                                                                      "yellow1", "yellow2",
                                                                      
                                                                      "grey99", "darkblue",
                                                                      
                                                                      "grey99", "hotpink", 
                                                                      "hotpink3",
                                                                      
                                                                      "grey99", "yellowgreen",
                                                                      "chartreuse", "chartreuse3",
                                                                      
                                                                      "grey99", "blueviolet",
                                                                      "darkorchid1", 
                                                                      
                                                                      "grey99", "Springgreen1",
                                                                      "springgreen3",
                                                                      "springgreen4",
                                                                      
                                                                      "grey20")

pattpatt <- c('vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw', 'vertical', 'bricks', 'hs_diagcross', 'verticalsaw', 
              'crosshatch', 'horizontal3', 'horizontalsaw')


#Another way to create a graph with overlapping colors and patterns. This particular code is great for when you don't have an
#unreasonable number of categories, since it requires manual input of colors and patterns rather than using a function.
#Unlike the code we used above, because the colors and patterns are not separate functions, this code generates a graph with 
#a single legend, rather than 2 legends as the graphs above do, with one for the colors and one for the patterns.

Barchart_LV_FLAPS_stagged_dA_Genus_5Perc <- ggplot(subset(InputData_LV_dA_Stagged_5Perc, Mesocosm %in% c("1")), aes(x=IncubationTime, y=ReadAbundance, fill=TaxonShort_f)) +
  
  geom_bar_pattern(position = "stack",
                   stat = "identity",
                   color = "black", 
                   pattern_fill = "black",
                   pattern_spacing = 0.05,
                   pattern = 'magick',
                   aes(pattern_type = TaxonShort_f)) +
  scale_fill_manual(values = fill) +
  scale_pattern_type_manual(values = pattpatt) +
  scale_y_continuous(limits=c(0,101)) +
  facet_grid(Station_f~Substrate_f) + 
  theme_test(base_size = 16) +
  labs(x = "Incubation time [d]", y = "Relative read abundance [%]") +
  guides()


print(Barchart_LV_FLAPS_stagged_dA_Genus_5Perc)


