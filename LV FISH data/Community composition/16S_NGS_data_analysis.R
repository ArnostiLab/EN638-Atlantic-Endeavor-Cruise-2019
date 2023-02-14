#Here we have the code that produces bar charts and bubble plots
#of bacterial community composition in amended and unamended mesocosms
#from each station and depth sampled during the 2019 Endeavor cruise (EN638)

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


#Load ggpattern and magick packages so that patterns can be added over colors
library(ggpattern)
library(magick)

#Load ggh4x to enable the facet_nested function, which allows graphs to have
#multiple, stacked facets
library(ggh4x)

###################################################################################                                                                                #
#                        Read in data for unamended mesocosms                     #                                                                                #
###################################################################################

InputData_dB_Stagged_5Perc <- read_xlsx("EN638_norm_dB_5Perc_averraged_stagged.xlsx")
InputData_dD_Stagged_5Perc <- read_xlsx("EN638_norm_dD_5Perc_averraged_stagged.xlsx")

###################################################################################                                                                                
#                         Read in data for amended mesocosms                      #                                                                                
###################################################################################

InputData_LV_dA_Stagged_5Perc <- read_xlsx("EN638_LV_dA_norm_5Perc_stagged.xlsx")
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
                        "blue1", "blue2", "blue1", "blue4", "dodgerblue2", "dodgerblue3", "dodgerblue1", "blue2", "dodgerblue1", 
                        
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
                        "SAR116 clade", "Sagittulla", "SAR11 clade IV", "Methylobacterium-Methylorobrum", "Magnetospiraceae uncultured", "SD2D12", "Methylorubrum",
                        
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
                        'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100', 'gray100',
                        
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
                        "SAR116 clade", "Sagittulla", "SAR11 clade IV", "Methylobacterium-Methylorobrum", "Magnetospiraceae uncultured", "SD2D12", "Methylorubrum",
                        
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
InputData_LV_dA_Stagged_5Perc$IncubationTime <- factor(InputData_LV_dA_Stagged_5Perc$IncubationTime, levels=c("0", "1", "2", "3", "4"), labels=c("2", "5", "7", "12", "17"))
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

#Creating the bar chart with the colors and patterns overlapping
Barchart_LV_FLAPS_stagged_dA_Genus_5Perc <- ggplot(subset(InputData_LV_dA_Stagged_5Perc, Mesocosm %in% c("1")), aes(x=IncubationTime, y=ReadAbundance, fill=TaxonShort_f)) +
  
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
  facet_grid(vars(Station_f), vars(Substrate_f)) + 
  theme_test(base_size = 16) +
  labs(x = "Incubation time [d]", y = "Relative read abundance [%]")

print(Barchart_LV_FLAPS_stagged_dA_Genus_5Perc)

#Saving as an editable pdf
ggsave(filename="EN638_LV_FLAPS_stagged_dA_Genus_5Perc_with_patterns_legend.pdf", plot=Barchart_LV_FLAPS_stagged_dA_Genus_5Perc, width=12, height=10, device=cairo_pdf)

#For use when you want a legend:
#scale_fill_bacteria(guide = guide_legend(name = "taxa")) +
#scale_pattern_type_bacteria(guide = guide_legend(override.aes = list(fill = "white")))

#Note that the control here is actually samples directly from the 20 L amended mesocosm, not from a smaller incubation run in parallel with the 
#FLA-PS incubations that contains no substrate. This means that the timepoints for the Control here are also different than those for the
#FLA-PS incubations - they should be 0, 3, 5, 10, 15 (no 15 d sample for stn 18 dA)

###################################################################################
#           Create barchart - Amended dC (Bottom) 5% Genus level                  #
###################################################################################

#Factor data and define levels
InputData_LV_dC_Stagged_5Perc$IncubationTime <- factor(InputData_LV_dC_Stagged_5Perc$IncubationTime, levels=c("0", "1", "2", "3", "4"), labels=c("2", "5", "7", "12", "17"))
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
#This way the taxonomic names used is consistent across all datasets, which is important when we create the legend and bubble plot later in the code
InputData_LV_dC_Stagged_5Perc <- InputData_LV_dC_Stagged_5Perc %>%
  mutate(TaxonShort_f = recode(TaxonShort_f, "Bacteroidetes" = 'Bacteroidota', "Bacteroidetes <5%" = 'Bacteroidota <5%'))

InputData_LV_dC_Stagged_5Perc <- InputData_LV_dC_Stagged_5Perc %>%
  mutate(Phylum = recode(Phylum, "Bacteroidetes" = 'Bacteroidota'))

InputData_LV_dC_Stagged_5Perc <- InputData_LV_dC_Stagged_5Perc %>%
  mutate(Phylum = recode(Phylum, "Marinimonas" = 'Marinimicrobia'))

InputData_LV_dC_Stagged_5Perc <- InputData_LV_dC_Stagged_5Perc %>%
  mutate(Phylum = recode(Phylum, "SAR324" = 'SAR324 clade'))

#Creating the bar chart with the colors and patterns overlapping
Barchart_LV_FLAPS_stagged_dC_Genus_5Perc <- ggplot(subset(InputData_LV_dC_Stagged_5Perc, Mesocosm %in% c("1")), aes(x=IncubationTime, y=ReadAbundance, fill=TaxonShort_f)) +
  
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
  facet_grid(vars(Station_f), vars(Substrate_f)) + 
  theme_test(base_size = 16) +
  labs(x = "Incubation time [d]", y = "Relative read abundance [%]")

print(Barchart_LV_FLAPS_stagged_dC_Genus_5Perc)

#Saving as an editable pdf
ggsave(filename="EN638_LV_FLAPS_stagged_dC_Genus_5Perc_with_patterns.pdf", plot=Barchart_LV_FLAPS_stagged_dC_Genus_5Perc, width=12, height=10, device=cairo_pdf)


#Note that the control here is actually samples directly from the 20 L amended mesocosm, not from a smaller incubation run in parallel with the 
#FLA-PS incubations that contains no substrate. This means that the timepoints for the Control here are also different than those for the
#FLA-PS incubations - they should be 0, 3, 5, 10, 15 


###################################################################################
#            Create barchart - Unamended dB (DCM) 5% Genus level                  #
###################################################################################

#Factor data and define levels
InputData_dB_Stagged_5Perc$IncubationTime <- factor(InputData_dB_Stagged_5Perc$IncubationTime, levels=c("0", "1", "3", "10"), labels=c("0", "1", "3", "10"))
InputData_dB_Stagged_5Perc$Substrate_f <- factor(InputData_dB_Stagged_5Perc$Substrate, levels=c("Pullulan", "Laminarin", "Xylan", "Fucoidan", "Arabinogalactan", "Chondroitin", "Control"))
InputData_dB_Stagged_5Perc$TaxonShort_f <- factor(InputData_dB_Stagged_5Perc$TaxonShort, levels=c("Gammaproteobacteria <5%",
                                                                                                    "Alteromonas", 
                                                                                                    "Colwellia", "Colwelliaceae uncultured", "Thalassotalea", "Pseudoalteromonas", "Vibrio",
                                                                                                    "Marine Methylotrophic Group 3", "Methylophaga", "Pseudohongiella",
                                                                                                    "SAR86 clade", 
                                                                                                    
                                                                                                    "Alphaproteobacteria <5%",
                                                                                                    "Amylibacter", "Dinoroseobacter", "Lacimonas", 
                                                                                                    "Shimia", 
                                                                                                    "Rhodobacteraceae uncultured", "AEGEAN-169 marine group", 
                                                                                                    "SAR11 clade I uncultured", "SAR11 clade Ia", "SAR11 clade Ib", "SAR11 clade II",
                                                                                                    
                                                                                                    "Bacteroidota <5%",
                                                                                                    "Cryomorphaceae uncultured",
                                                                                                    "Flavicella", "Formosa", 
                                                                                                    "Maribacter", "NS4 marine group", "NS5 marine group", "Polaribacter",
                                                                                                    "Tenacibaculum", "Flavobacteriaceae uncultured", "NS9 marine group",
                                                                                                    
                                                                                                    "SAR324 clade(Marine group B)", 
                                                                                                    
                                                                                                    "Actinobacteria <5%",
                                                                                                    "Candidatus Actinomarina",  
                                                                                                    
                                                                                                    "Cyanobacteria <5%",
                                                                                                    "Prochlorococcus MIT9313", "Synechococcus CC9902", 
                                                                                                    
                                                                                                    "Marinimicrobia (SAR406 clade)", 
                                                                                                    
                                                                                                    "Chloroflexi <5%",
                                                                                                    "SAR202 clade",
                                                                                                    
                                                                                                    "Planctomycetes <5%", 
                                                                                                    "CL500-3",
                                                                                                    
                                                                                                    "Verrucomicrobia <5%",
                                                                                                    "DEV007", 
                                                                                                    
                                                                                                    "other bacteria <5%"))

InputData_dB_Stagged_5Perc$Station_f <- factor(InputData_dB_Stagged_5Perc$Station, levels = c("18", "19", "20"))
InputData_dB_Stagged_5Perc$Depth_f <- factor(InputData_dB_Stagged_5Perc$Depth, levels = c("B"), labels=c("DCM"))

#Using dplyr to relabel Bacteroitoda (it should be Bacteroidota)
InputData_dB_Stagged_5Perc <- InputData_dB_Stagged_5Perc %>%
  mutate(Phylum = recode(Phylum, "Bacteroitoda" = 'Bacteroidota'))


#Creating the bar chart with the colors and patterns overlapping
Barchart_Unamend_FLAPS_stagged_dB_Genus_5Perc <- ggplot(InputData_dB_Stagged_5Perc, aes(x=IncubationTime, y=RelativeAbundance, fill = TaxonShort_f)) + 
  
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
  facet_grid(vars(Station_f), vars(Substrate_f)) + 
  theme_test(base_size = 16) +
  labs(x = "Incubation time [d]", y = "Relative read abundance [%]")

print(Barchart_Unamend_FLAPS_stagged_dB_Genus_5Perc)

#Saving as an editable pdf
ggsave(filename="EN638_Unamended_FLAPS_stagged_DCM_Genus_5Perc_with_patterns.pdf", plot=Barchart_Unamend_FLAPS_stagged_dB_Genus_5Perc, width=12, height=10, device=cairo_pdf)


#Note that the control here is from an incubation run in parallel with the FLA-PS 
#incubations that contains no substrate, so the time points are the same as the FLA-PS time points

###################################################################################
#           Create barchart - Unamended dD (Bottom) 5% Genus level                #
###################################################################################

#Factor data and define levels
InputData_dD_Stagged_5Perc$IncubationTime <- factor(InputData_dD_Stagged_5Perc$IncubationTime, levels=c("0", "1", "3", "10"), labels=c("0", "1", "3", "10"))
InputData_dD_Stagged_5Perc$Substrate_f <- factor(InputData_dD_Stagged_5Perc$Substrate, levels=c("Pullulan", "Laminarin", "Xylan", "Fucoidan", "Arabinogalactan", "Chondroitin", "Control"))
InputData_dD_Stagged_5Perc$TaxonShort_f <- factor(InputData_dD_Stagged_5Perc$TaxonShort, levels=c("Gammaproteobacteria <5%",
                                                                                                  "Limnobacter", "Alteromonas", "Salinimonas",
                                                                                                  "Colwellia", "Colwelliaceae uncultured", "Pseudoalteromonas", "Vibrio",
                                                                                                  "Marine Methylotrophic Group 3", "Alcanivorax",
                                                                                                  "Marinomonas", "Acinetobacter", "Pseudomonas", 
                                                                                                  "Oleispira", "Saccharospirillaceae uncultured", "Thalassolituus", 
                                                                                                  "SAR86 clade", "Sinobacterium",
                                                                                                  
                                                                                                  "Alphaproteobacteria <5%",
                                                                                                  "Shimia", "Sulfitobacter", "Methylobacterium-Methylorubrum", 
                                                                                                  "AEGEAN-169 marine group", 
                                                                                                  "SAR11 clade I uncultured", "SAR11 clade Ia", "SAR11 clade Ib", "SAR11 clade II",
                                                                                                  
                                                                                                  "Bacteroidota <5%",
                                                                                                  "Cryomorphaceae uncultured",
                                                                                                  
                                                                                                  "SAR324 clade(Marine group B)",
                                                                                                  
                                                                                                  "Marinimicrobia (SAR406 clade)",
                                                                                                  
                                                                                                  "Chloroflexi <5%",
                                                                                                  "SAR202 clade",
                                                                                                  
                                                                                                  "Verrucomicrobia <5%",
                                                                                                  "DEV007", 
                                                                                                  
                                                                                                  "other bacteria <5%"))

InputData_dD_Stagged_5Perc$Station_f <- factor(InputData_dD_Stagged_5Perc$Station, levels = c("18", "19", "20"))
InputData_dD_Stagged_5Perc$Depth_f <- factor(InputData_dD_Stagged_5Perc$Depth, levels = c("D"), labels=c("Bottom"))

#Creating the bar chart with the colors and patterns overlapping
Barchart_Unamend_FLAPS_stagged_dD_Genus_5Perc <- ggplot(InputData_dD_Stagged_5Perc, aes(x=IncubationTime, y=RelativeAbundance, fill=TaxonShort_f)) +
  
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
  facet_grid(vars(Station_f), vars(Substrate_f)) + 
  theme_test(base_size = 16) +
  labs(x = "Incubation time [d]", y = "Relative read abundance [%]")

print(Barchart_Unamend_FLAPS_stagged_dD_Genus_5Perc)

#Saving as an editable pdf
ggsave(filename="EN638_Unamended_FLAPS_stagged_Bottom_Genus_5Perc_with_patterns.pdf", plot=Barchart_Unamend_FLAPS_stagged_dD_Genus_5Perc, width=12, height=10, device=cairo_pdf)

#Note that the control here is from an incubation run in parallel with the FLA-PS 
#incubations that contains no substrate, so the time points are the same as the FLA-PS time points

###################################################################################
#    Create bubble plot - Laminarin, amended mesocosms, and controls only         #
###################################################################################

#The above sections of code should be run before this section to ensure the data is
#correctly labeled. Unlike the above sections, this code will create a bubble plot
#rather than a barchart.

#Subsetting the data to include only laminarin and control incubations
Amend_dA_filt <- InputData_LV_dA_Stagged_5Perc %>% 
  filter(Substrate_f %in% c("Laminarin", "Amended Meso"))

Amend_dC_filt <- InputData_LV_dC_Stagged_5Perc %>% 
  filter(Substrate_f %in% c("Laminarin", "Amended Meso"))

Unamend_dB_filt <- InputData_dB_Stagged_5Perc %>% 
  filter(Substrate_f %in% c("Laminarin", "Control"))

Unamend_dD_filt <- InputData_dD_Stagged_5Perc %>% 
  filter(Substrate_f %in% c("Laminarin", "Control"))

#For the amended mesocosms, we only want to include amended mesocosm 1; in the above
#code, the other mesocosms are filtered out during the graphing step (i.e., Mesocosm %in% c("1"))
#Note that we don't have to do this for the unamended mesocosms because there was only one
#Unamended mesocosm from each station/depth, while there were 3 amended mesocosm for each 
#station/depth, with the FLA-PS incubations done using water from amended mesocosm 1
Amend_dA_filt <- Amend_dA_filt %>% 
  filter(Mesocosm %in% c("1"))

Amend_dC_filt <- Amend_dC_filt %>% 
  filter(Mesocosm %in% c("1"))

#To limit the amount of timepoints displayed in this graph, we'll remove the 17 d timepoint from the amended
#mesocosms
Amend_dA_filt <- Amend_dA_filt %>% 
  filter(IncubationTime %in% c("2", "5", "7", "12"))

Amend_dC_filt <- Amend_dC_filt %>% 
  filter(IncubationTime %in% c("2", "5", "7", "12"))

#Adding a Treatment column to each data subset so that amended and unamended can be 
#distinguished from one another
Amend_dA_filt$Treatment <- "Amended"
Amend_dC_filt$Treatment <- "Amended"
Unamend_dB_filt$Treatment <- "Unamended"
Unamend_dD_filt$Treatment <- "Unamended"

#Change column names ReadAbundance (Amended) and RelativeAbundance (Unamended) to Abundance
#so that the columns in different data sets match and they can be combined into one data set
Amend_dA_filt <- dplyr::rename(Amend_dA_filt, Abundance = ReadAbundance)

Amend_dC_filt <- dplyr::rename(Amend_dC_filt, Abundance = ReadAbundance)

Unamend_dB_filt <- dplyr::rename(Unamend_dB_filt, Abundance = RelativeAbundance)

Unamend_dD_filt <- dplyr::rename(Unamend_dD_filt, Abundance = RelativeAbundance)

#Remove extra columns not included in all data sets 
Amend_dA_filt <- subset(Amend_dA_filt, select=-c(BC, Mesocosm))
Amend_dC_filt <- subset(Amend_dC_filt, select=-c(BC, Mesocosm))

#Changing the IncubationTime from days to timepoints to minimize the number of blank
#days on the bubbleplot
Amend_dA_filt$IncubationTime_f <- factor(Amend_dA_filt$IncubationTime, levels = c("2", "5", "7", "12"), labels=c("0", "1", "2", "3"))
Amend_dC_filt$IncubationTime_f <- factor(Amend_dC_filt$IncubationTime, levels = c("2", "5", "7", "12"), labels=c("0", "1", "2", "3"))
Unamend_dB_filt$IncubationTime_f <- factor(Unamend_dB_filt$IncubationTime, levels = c("0", "1", "3", "10"), labels=c("0", "1", "2", "3"))
Unamend_dD_filt$IncubationTime_f <- factor(Unamend_dD_filt$IncubationTime, levels = c("0", "1", "3", "10"), labels=c("0", "1", "2", "3"))

#Combining the subsets above
Total_community_comp <- rbind(Amend_dA_filt, Amend_dC_filt, Unamend_dB_filt, Unamend_dD_filt)

#Factor and define levels for Treatment and IncubationTime
Total_community_comp$Treatment_f <- factor(Total_community_comp$Treatment, levels = c("Amended", "Unamended"))
Total_community_comp$Substrate_f <- factor(Total_community_comp$Substrate_f, levels = c("Laminarin", "Amended Meso", "Control"), labels = c("Lam", "Meso", "Cont"))

#Remove TaxonShort_f with Abundance lower than 5%
Total_community_comp <- subset(Total_community_comp, Abundance>=5)

#Choosing a color scheme
colours = c("#F4DE52FF", "#2BA872", "#421C57")
lam_color = c("#FFC425", "grey60", "grey60")


#Make the bubble plot - first, group the families by Phylum
xx = ggplot(Total_community_comp, aes(x=IncubationTime_f, y = TaxonShort_f)) + 
  geom_point(aes(size = Abundance, fill = Substrate_f), alpha = 1, shape = 21) + 
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
  scale_fill_manual(values = lam_color, guide = FALSE) + 
  scale_y_discrete(limits = rev(levels(Total_community_comp$Abundance)))+
  facet_nested(Phylum ~ Treatment_f + Depth_f + Station_f + Substrate_f, scales="free_y", space="free_y", switch="y")


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
ggsave(filename="Bubbleplot_abundance_5percent.pdf", plot=xy, width=9, height=7, device=cairo_pdf)


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


