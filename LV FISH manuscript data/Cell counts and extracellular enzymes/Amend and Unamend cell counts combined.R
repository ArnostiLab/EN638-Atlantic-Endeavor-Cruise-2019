
#Load packages
library("plyr")
library("ggplot2")
library("scales")
library("dplyr")
library("tidyr")
library("reshape2")
library("ggh4x")
library("patchwork")
library("tidyverse")
library("ggbreak")
#devtools::install_github("teunbrand/ggh4x")

###################################################################################                                                                                
#                     Read in the amended cell count data                         #                                                                                
###################################################################################
Cell <- read.csv("EN638_LV_substrate_cellcounts_all_forgraphs.csv")

#Calculating percent of selfish bacteria
avg <- Cell %>%
  group_by(stn, depth_name, substrate, time_d, replicate, cells_ml, substratecells_ml) %>%
  summarise(selfish_percent = (substratecells_ml/cells_ml)*100)

#Calculating total averages and standard deviations from replicate samples for total cell counts and selfish counts
selfish_Amend <- avg %>%
  group_by(stn, depth_name, substrate, time_d) %>%
  summarise(mean.cells=mean(cells_ml), sd.cells=sd(cells_ml), 
            mean.selfish=mean(substratecells_ml), sd.selfish=sd(substratecells_ml),
            mean.percent=mean(selfish_percent), sd.percent=sd(selfish_percent))

#Ordering Factors
selfish_Amend$stn <- factor(selfish_Amend$stn, levels = c("stn18","stn19", "stn20"), labels=c("Station 18", "Station 19", "Station 20"))
selfish_Amend$depth_name <- factor(selfish_Amend$depth_name, levels = c("DCM", "bottom"))
selfish_Amend$substrate <- factor(selfish_Amend$substrate, levels = c("pul", "lam", "xyl", "fuc", "ara", "chn"))

#Ordering factors (time_day) and renaming them to be the time since water collection so that they are on the same time frame as the Unamended
selfish_Amend$time_d <- factor(selfish_Amend$time_d, levels = c("0", "3", "5", "10", "15", "30"), 
                     labels = c("2", "5", "7", "12", "17", "32"))

#Adding a column noting that these are amended cell counts
selfish_Amend$Treatment <- "Amended"

###################################################################################                                                                                
#                    Read in the unamended cell count data                        #                                                                                
###################################################################################
selfish_Unamend <- read.csv("EN638_LV_Unamend_cellcounts.csv")

#Ordering Factors
selfish_Unamend$stn <- factor(selfish_Unamend$stn, levels = c("stn18","stn19", "stn20"), labels=c("Station 18", "Station 19", "Station 20"))
selfish_Unamend$depth_name <- factor(selfish_Unamend$depth_name, levels = c("DCM", "bottom"))
selfish_Unamend$substrate <- factor(selfish_Unamend$substrate, levels = c("pul", "lam", "xyl", "fuc", "ara", "chn", "control"))

#Ordering factors (time_day) and changing them to the correct times
selfish_Unamend$time_d <- factor(selfish_Unamend$time_d, levels = c("0.2", "1", "2", "3", "5", "7", "10", "15", "30"), 
                        labels = c("0", "1", "1", "3", "5", "7", "10", "15", "30"))

#Adding a column noting that these are unamended cell counts
selfish_Unamend$Treatment <- "Unamended"

###################################################################################                                                                                
#      Combining amended and unamended cell counts into a single dataframe        #                                                                                
###################################################################################
selfish_all <- bind_rows(selfish_Amend, selfish_Unamend)

#Remove 15, 30, and 32 day timepoints since there's no data
selfish_all <- filter(selfish_all, time_d != "15")
selfish_all <- filter(selfish_all, time_d != "30")
selfish_all <- filter(selfish_all, time_d != "32")

#Ordering factors
selfish_all$time_d <- factor(selfish_all$time_d, levels = c("0", "1", "2", "3", "5", "7", "10", "12", "17"))
selfish_all$depth_name <- factor(selfish_all$depth_name, levels = c("DCM", "bottom"), labels=c("DCM", "Bottom"))

###################################################################################                                                                                
#   Read in the extracellular enzymatic activity data from amended mesocosms      #                                                                                
###################################################################################
LV <- read.csv("FlaRates_EN638_LV_Amend_Means_Corrected_Final.csv")

#Creating columns from row names in the LV dataframe
LV_edit <- LV %>%
  separate(X, c("Station", "Depth", "LV", "SpecificTreatment", "Substrate", "Timepoint"), "-")

#Deleting columns that we won't be using in the LV dataframe
Amend_EEA <- subset(LV_edit, select=-c(X.1, LV))

#Filtering the dataset; we only need depths dA and dC. Even though amendments were
#done at another depth (dB) selfish incubations were done only with dA and dC water.
Amend_EEA <- Amend_EEA %>% 
  select(Station, Depth, SpecificTreatment, Substrate, Timepoint, rate.x.nM.hr, rate.1.nM.hr, rate.2.nM.hr, rate.3.nM.hr, mean.rate.nM.hr, sd.rate.nM.hr,
         kcrate.x.nM.hr, kcrate.1.nM.hr, kcrate.2.nM.hr, kcrate.3.nM.hr, mean.kcrate.nM.hr, sd.kcrate.nM.hr) %>%
  filter(Depth %in% c("dA", "dC"))

#Here, we only use the data from mesocosm A1, because the selfish incubations were 
#only done in mesocosm A1. So now we'll filter out mesocosms A2 and A3 as well.
Amend_EEA <- Amend_EEA %>% 
  select(Station, Depth, SpecificTreatment, Substrate, Timepoint, mean.kcrate.nM.hr, sd.kcrate.nM.hr) %>%
  filter(!SpecificTreatment %in% c("A2", "A3"))

#Ordering Factors
#Note: the t0 timepoint in the amended mesocosms was 2 days after HMW-OM addition to 
#amended mesocosms, so t0 is 2 d.
Amend_EEA$Station <- factor(Amend_EEA$Station, levels = c("stn18", "stn19", "stn20"), labels=c("Station 18", "Station 19", "Station 20"))
Amend_EEA$Depth <- factor(Amend_EEA$Depth, levels = c("dA",  "dC"), labels=c("DCM", "Bottom"))
Amend_EEA$Substrate <- factor(Amend_EEA$Substrate, levels = c("pul", "lam", "xyl", "fuc", "ara", "chn"))
Amend_EEA$Timepoint <- factor(Amend_EEA$Timepoint, levels = c("t0", "t1", "t2", "t3", "t4", "t5"), labels=c("2", "5", "7", "12", "17", "31"))

###################################################################################                                                                                
#   Read in the extracellular enzymatic activity data from unamended mesocosms    #                                                                                
###################################################################################
#FLA samples from the unamended mesocosms were run with the bulk data; they're pulled 
#from the bulk dataframe below.
Bulk <- read.csv("FlaRatesFINAL_Corrected_Sept2021.csv")

#Creating columns from row names
Bulk_edit <- Bulk %>%
  separate(X, c("Station", "Depth", "SpecificTreatment", "Substrate", "Timepoint"), "-")

#Deleting columns that we won't be using in the bulk dataframe
Bulk_edit <- subset(Bulk_edit, select=-c(X.1))

#Filter bulk so that only dA and dC unamended mesocosms remain; these have 7 timepoints
Bulk_filt <- Bulk_edit %>% 
  select(Station, Depth, SpecificTreatment, Substrate, Timepoint, rate.x.nM.hr, rate.1.nM.hr, rate.2.nM.hr, rate.3.nM.hr, mean.rate.nM.hr, sd.rate.nM.hr,
         kcrate.x.nM.hr, kcrate.1.nM.hr, kcrate.2.nM.hr, kcrate.3.nM.hr, mean.kcrate.nM.hr, sd.kcrate.nM.hr) %>%
  filter(Depth %in% c("d2", "d7"))

#Filter bulk to get stn 18 d6 (there was no d7 at station 18)
stn18d6 <- Bulk_edit %>% 
  select(Station, Depth, SpecificTreatment, Substrate, Timepoint, rate.x.nM.hr, rate.1.nM.hr, rate.2.nM.hr, rate.3.nM.hr, mean.rate.nM.hr, sd.rate.nM.hr,
         kcrate.x.nM.hr, kcrate.1.nM.hr, kcrate.2.nM.hr, kcrate.3.nM.hr, mean.kcrate.nM.hr, sd.kcrate.nM.hr) %>%
  filter(Station == "stn18" & Depth == "d6")

#Combine these two dataframes
Unamend_EEA <- rbind(Bulk_filt, stn18d6)

#Remove station 17
Unamend_EEA <- filter(Unamend_EEA, Station != "stn17")

#In the unamend dataframe, change bulk in the SpecificTreatment column to U and change d2, d6, and d7 to dA and dC
Unamend_EEA$SpecificTreatment <- factor(Unamend_EEA$SpecificTreatment, levels = c("bulk"), labels=c("U"))
Unamend_EEA$Depth <- factor(Unamend_EEA$Depth, levels = c("d2", "d6", "d7"), labels=c("dA", "dC", "dC"))

#Ordering Factors
#Note: the t0 timepoint in the amended mesocosms was 2 days after HMW-OM addition to amended mesocosms
#dB unamended mesocosms were started with FLA immediately, not 2 days after HMW-OM addition like the amended mesocosms
#MPI did dA unamended and dC unamended mesocosms for FLA incubations with selfish bacteria, so these have some extra timepoints
#.2 is 5 hours
Unamend_EEA$Station <- factor(Unamend_EEA$Station, levels = c("stn18", "stn19", "stn20"), labels=c("Station 18", "Station 19", "Station 20"))
Unamend_EEA$Depth <- factor(Unamend_EEA$Depth, levels = c("dA", "dC"), labels=c("DCM", "Bottom"))
Unamend_EEA$Substrate <- factor(Unamend_EEA$Substrate, levels = c("pul", "lam", "xyl", "fuc", "ara", "chn"))
Unamend_EEA$Timepoint <- factor(Unamend_EEA$Timepoint, levels = c("t0", "t1", "t2", "t3", "t4", "t5", "t6", "t7"), labels=c("0", "0.2", "1", "3", "7", "10", "15", "30"))

###################################################################################                                                                                
#                                       Figure 1a                                 #                                                                                
###################################################################################
#Graphing cell counts, selfish uptake, and extracellular enzymatic activities in
#unamended DCM mesocosms.

#Subset the cell count data by depth so that they can be graphed individually
Selfish_DCM <- subset(selfish_all, depth_name=="DCM")
Selfish_Bottom <- subset(selfish_all, depth_name=="Bottom")

#Further subset DCM into unamended mesocosms only 
Selfish_DCM_Unamend <- subset(Selfish_DCM, Treatment=="Unamended")

#Remove the 5 d timepoint for which there is no data
Selfish_DCM_Unamend <- subset(Selfish_DCM_Unamend, time_d!="5")


#Add a (blank) 15 and 30 day timepoint so that the cell count data has the same number
#of timepoints on the x-axis as the extracellular enzymatic activity data.
#To do this, we have to first create a function so that we don't get the 
#"Can't add rows to grouped data frames" error.
add_summary_rows <- function(.data, ...) {
  group_modify(.data, function(x, y) bind_rows(x, summarise(x, ...)))
}

Selfish_DCM_Unamend <- Selfish_DCM_Unamend %>% 
  add_summary_rows(
    time_d = "15"
  )

Selfish_DCM_Unamend <- Selfish_DCM_Unamend %>% 
  add_summary_rows(
    time_d = "30"
  )

#Order factors again
Selfish_DCM_Unamend$time_d <- factor(Selfish_DCM_Unamend$time_d, levels = c("0", "1", "3", "7", "10", "15", "30"), 
                     labels = c("0", "1", "3", "7", "10", "15", "30"))


#Plot cell count line graph
Unamend_DCM_cell_line <-ggplot(Selfish_DCM_Unamend, aes(x=time_d, y=mean.cells, group=substrate))+geom_line(aes(color=substrate))+geom_point(aes(color=substrate))+
  scale_y_continuous(breaks=c(0, 1000000, 2000000, 3000000, 4000000, 5000000), limits=c(0, 5000000), labels = scientific)+
  facet_grid(cols = vars(stn))+labs(x="Time (days)", y=("Total cellular abundance (cells/mL)"))+
  scale_color_manual(values=c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF", "control" = "grey60"))+
  geom_errorbar(aes(ymin=mean.cells-sd.cells, ymax=mean.cells+sd.cells), width = 0.1, size = 0.5)+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_blank(), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text = element_text(size=10, color="black"))

print(Unamend_DCM_cell_line)


#DCM selfish percentage bar graph
#Remove xylan
Selfish_DCM_Unamend_percent <- Selfish_DCM_Unamend %>% filter(substrate != "xyl")

#Plot %selfish uptake
Unamend_DCM_selfish_percent <-ggplot(Selfish_DCM_Unamend_percent,aes(x=time_d, y=mean.percent, ymin=mean.percent-sd.percent, ymax=mean.percent+sd.percent, fill=substrate))+
  geom_bar(position=position_dodge(), stat="identity", width=.9, colour = "black", size=0.2)+
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30), limits=c(0, 30))+
  facet_grid(cols = vars(stn))+labs(x="Time (days)", y=("Percent of selfish bacteria"))+
  scale_fill_manual(values = c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF"))+
  geom_errorbar(position=position_dodge(width = 0.9), size = 0.5) +
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_blank(), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text = element_blank())

print(Unamend_DCM_selfish_percent)


#DCM extracellular enzymatic activities bar graph
#Filter the data first; remove bottom water.
Unamend_DCM_EEA <- Unamend_EEA %>% 
  select(Station, Depth, SpecificTreatment, Substrate, Timepoint, mean.kcrate.nM.hr, sd.kcrate.nM.hr) %>%
  filter(Depth %in% c("DCM"))

#Order factors
Unamend_DCM_EEA$Timepoint <- factor(Unamend_DCM_EEA$Timepoint, levels = c("0", "0.2", "1", "3", "7", "10", "15", "30"), labels=c("0", "0", "1", "3", "7", "10", "15", "30"))
Unamend_DCM_EEA$Substrate <- factor(Unamend_DCM_EEA$Substrate, levels = c("pul", "lam", "xyl", "fuc", "ara", "chn"))

#Plot EEA
Unamend_DCM_EEA_bar <-ggplot(Unamend_DCM_EEA, aes(x=Timepoint, y=mean.kcrate.nM.hr, ymin=mean.kcrate.nM.hr-sd.kcrate.nM.hr, ymax=mean.kcrate.nM.hr+sd.kcrate.nM.hr, fill=Substrate))+
  geom_bar(position=position_dodge(), stat="identity", width=1, colour = "black", size=0.2) +
  facet_grid(cols = vars(Station))+
  labs(x="Time (days)", y=substitute(paste("Hydrolysis rate (nmol ",L^-1, h^-1,")"))) +
  scale_fill_manual(values = c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF")) +
  geom_errorbar(position=position_dodge(1), size = 0.5) +
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_text(size=9, colour="black", angle=0, hjust=0.5, vjust=1), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_text(size=10), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text.x = element_blank())


print(Unamend_DCM_EEA_bar)

#Combine these graphs using patchwork
Fig_1a <- Unamend_DCM_cell_line / Unamend_DCM_selfish_percent / Unamend_DCM_EEA_bar 

print(Fig_1a)

ggsave(filename="EN638_LV_Fig_1a.pdf", plot=Fig_1a, width=8, height=8, device=cairo_pdf)

###################################################################################                                                                                
#                                       Figure 2a                                 #                                                                                
###################################################################################
#Graphing cell counts, selfish uptake, and extracellular enzymatic activities in
#unamended Bottom water mesocosms.

#The cell count data have already been subsetted by depth (in the section for 
#figure 1a) so that they can be graphed individually; here, we'll further subset 
#the bottom water cell count data so that it only includes unamended mesocosms 
Selfish_Bottom_Unamend <- subset(Selfish_Bottom, Treatment=="Unamended")

#Remove the 5 d timepoint for which there is no data
Selfish_Bottom_Unamend <- subset(Selfish_Bottom_Unamend, time_d!="5")


#Now we'll use the function that we created in section 1a to add a 15 and 30 d timepoint
#to the dataframe so that it will match that of the enzymatic activity dataset
Selfish_Bottom_Unamend <- Selfish_Bottom_Unamend %>% 
  add_summary_rows(
    time_d = "15"
  )

Selfish_Bottom_Unamend <- Selfish_Bottom_Unamend %>% 
  add_summary_rows(
    time_d = "30"
  )

#Order factors again
Selfish_Bottom_Unamend$time_d <- factor(Selfish_Bottom_Unamend$time_d, levels = c("0", "1", "3", "7", "10", "15", "30"), 
                                     labels = c("0", "1", "3", "7", "10", "15", "30"))


#Plot cell count line graph
Unamend_Bottom_cell_line <-ggplot(Selfish_Bottom_Unamend, aes(x=time_d, y=mean.cells, group=substrate))+geom_line(aes(color=substrate))+geom_point(aes(color=substrate))+
  scale_y_continuous(breaks=c(0, 25000, 50000, 75000, 100000), limits=c(0, 100000), labels = scientific)+
  facet_grid(cols = vars(stn))+labs(x="Time (days)", y=("Total cellular abundance (cells/mL)"))+
  scale_color_manual(values=c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF", "control" = "grey60"))+
  geom_errorbar(aes(ymin=mean.cells-sd.cells, ymax=mean.cells+sd.cells), width = 0.1, size = 0.5)+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_blank(), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text = element_text(size=10, color="black"))

print(Unamend_Bottom_cell_line)


#Bottom water selfish percentage bar graph
#Remove xylan
Selfish_Bottom_Unamend_percent <- Selfish_Bottom_Unamend %>% filter(substrate != "xyl")

#Plot %selfish uptake
Unamend_Bottom_selfish_percent <-ggplot(Selfish_Bottom_Unamend_percent,aes(x=time_d, y=mean.percent, ymin=mean.percent-sd.percent, ymax=mean.percent+sd.percent, fill=substrate))+
  geom_bar(position=position_dodge(), stat="identity", width=.9, colour = "black", size=0.2)+
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30), limits=c(0, 30))+
  facet_grid(cols = vars(stn))+labs(x="Time (days)", y=("Percent of selfish bacteria"))+
  scale_fill_manual(values = c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF"))+
  geom_errorbar(position=position_dodge(width = 0.9), size = 0.5) +
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_blank(), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text = element_blank())

print(Unamend_Bottom_selfish_percent)


#Bottom water extracellular enzymatic activities bar graph
#Filter the data first; remove DCM water.
Unamend_Bottom_EEA <- Unamend_EEA %>% 
  select(Station, Depth, SpecificTreatment, Substrate, Timepoint, mean.kcrate.nM.hr, sd.kcrate.nM.hr) %>%
  filter(Depth %in% c("Bottom"))

#Order factors
Unamend_Bottom_EEA$Timepoint <- factor(Unamend_Bottom_EEA$Timepoint, levels = c("0", "0.2", "1", "3", "7", "10", "15", "30"), labels=c("0", "0", "1", "3", "7", "10", "15", "30"))
Unamend_Bottom_EEA$Substrate <- factor(Unamend_Bottom_EEA$Substrate, levels = c("pul", "lam", "xyl", "fuc", "ara", "chn"))

#Plot EEA
Unamend_Bottom_EEA_bar <-ggplot(Unamend_Bottom_EEA, aes(x=Timepoint, y=mean.kcrate.nM.hr, ymin=mean.kcrate.nM.hr-sd.kcrate.nM.hr, ymax=mean.kcrate.nM.hr+sd.kcrate.nM.hr, fill=Substrate))+
  geom_bar(position=position_dodge(), stat="identity", width=1, colour = "black", size=0.2) +
  facet_grid(cols = vars(Station))+
  labs(x="Time (days)", y=substitute(paste("Hydrolysis rate (nmol ",L^-1, h^-1,")"))) +
  scale_fill_manual(values = c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF")) +
  geom_errorbar(position=position_dodge(1), size = 0.5) +
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_text(size=9, colour="black", angle=0, hjust=0.5, vjust=1), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_text(size=10), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text.x = element_blank())


print(Unamend_Bottom_EEA_bar)

#Combine these graphs using patchwork
Fig_2a <- Unamend_Bottom_cell_line / Unamend_Bottom_selfish_percent / Unamend_Bottom_EEA_bar 

print(Fig_2a)

ggsave(filename="EN638_LV_Fig_2a.pdf", plot=Fig_2a, width=8, height=8, device=cairo_pdf)

###################################################################################                                                                                
#                                       Figure 1b                                 #                                                                                
###################################################################################
#Graphing cell counts, selfish uptake, and extracellular enzymatic activities in
#amended DCM mesocosms.


#Subset the cell count data by depth so that they can be graphed individually
Selfish_DCM <- subset(selfish_all, depth_name=="DCM")
Selfish_Bottom <- subset(selfish_all, depth_name=="Bottom")

#Further subset DCM into amended mesocosms only 
Selfish_DCM_Amend <- subset(Selfish_DCM, Treatment=="Amended")

#Add a (blank) 31 day timepoint so that the cell count data has the same number
#of timepoints on the x-axis as the extracellular enzymatic activity data.
#To do this, we have to first create a function so that we don't get the 
#"Can't add rows to grouped data frames" error.
add_summary_rows <- function(.data, ...) {
  group_modify(.data, function(x, y) bind_rows(x, summarise(x, ...)))
}

Selfish_DCM_Amend <- Selfish_DCM_Amend %>% 
  add_summary_rows(
    time_d = "31"
  )

#Order factors again
Selfish_DCM_Amend$time_d <- factor(Selfish_DCM_Amend $time_d, levels = c("2", "5", "7", "12", "17", "31"))


#Plot cell count line graph
DCM_cell_line <-ggplot(Selfish_DCM_Amend, aes(x=time_d, y=mean.cells, group=substrate))+geom_line(aes(color=substrate))+geom_point(aes(color=substrate))+
  scale_y_continuous(breaks=c(0, 1000000 ,2000000, 3000000, 4000000, 5000000), limits=c(0, 5000000), labels = scientific)+
  facet_grid(cols = vars(stn))+labs(x="Time (days)", y=("Total cellular abundance (cells/mL)"))+
  scale_color_manual(values=c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF", "control" = "grey60"))+
  geom_errorbar(aes(ymin=mean.cells-sd.cells, ymax=mean.cells+sd.cells), width = 0.1, size=0.5)+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_blank(), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text = element_text(size=10, color="black"))

print(DCM_cell_line)


#DCM selfish percentage bar graph
#Remove xylan
DCM_selfish <- Selfish_DCM_Amend %>% filter(substrate != "xyl")

#Plot %selfish uptake
DCM_selfish_percent <-ggplot(DCM_selfish, aes(x=time_d, y=mean.percent, ymin=mean.percent-sd.percent, ymax=mean.percent+sd.percent, fill=substrate))+
  geom_bar(position=position_dodge(), stat="identity", width=.9, colour = "black", size=0.2)+
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30), limits=c(-2, 30))+
  facet_grid(cols = vars(stn))+
  labs(x="Time (days)", y=("Percent of selfish bacteria"))+
  scale_fill_manual(values = c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF"))+
  geom_errorbar(position=position_dodge(width = 0.9), size = 0.5) +
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_blank(), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text = element_blank())

print(DCM_selfish_percent)


#DCM extracellular enzymatic activities bar graph
#Filter the data first; remove bottom water.
DCM_EEA <- Amend_EEA %>% 
  select(Station, Depth, SpecificTreatment, Substrate, Timepoint, mean.kcrate.nM.hr, sd.kcrate.nM.hr) %>%
  filter(Depth %in% c("DCM"))

#Plot EEA
DCM_EEA_bar <-ggplot(DCM_EEA, aes(x=Timepoint, y=mean.kcrate.nM.hr, ymin=mean.kcrate.nM.hr-sd.kcrate.nM.hr, ymax=mean.kcrate.nM.hr+sd.kcrate.nM.hr, fill=Substrate))+
  geom_bar(position=position_dodge(), stat="identity", width=1, colour = "black", size=0.2)+
  facet_grid(cols = vars(Station))+
  labs(x="Time (days)", y=substitute(paste("Hydrolysis rate (nmol ",L^-1, h^-1,")")))+
  scale_fill_manual(values = c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF"))+
  geom_errorbar(position=position_dodge(1), size = 0.5)+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_text(size=9, colour="black", angle=0, hjust=0.5, vjust=1), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_text(size=10), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text.x = element_blank())


print(DCM_EEA_bar)

#Combine these graphs using patchwork
Fig_1b <- DCM_cell_line / DCM_selfish_percent / DCM_EEA_bar 
Fig_1b + plot_layout(guides = 'collect')

print(Fig_1b)

ggsave(filename="EN638_LV_Fig_1b.pdf", plot=Fig_1b, width=8, height=8, device=cairo_pdf)


###################################################################################                                                                                
#                                       Figure 2b                                 #                                                                                
###################################################################################
#Graphing cell counts, selfish uptake, and extracellular enzymatic activities in
#amended Bottom water mesocosms.

#The cell count data have already been subsetted by depth (in the section for 
#figure 1b) so that they can be graphed individually; here, we'll further subset 
#the bottom water cell count data so that it only includes amended mesocosms 
Selfish_Bottom_Amend <- subset(Selfish_Bottom, Treatment=="Amended")

#Now we'll use the function that we created in section 1b to add a 31 d timepoint
#to the dataframe so that it will match that of the enzymatic activity dataset
Selfish_Bottom_Amend <- Selfish_Bottom_Amend %>% 
  add_summary_rows(
    time_d = "31"
  )

#Order factors again
Selfish_Bottom_Amend$time_d <- factor(Selfish_Bottom_Amend $time_d, levels = c("2", "5", "7", "12", "17", "31"))


#Plot cell count line graph
Bottom_cell_line <-ggplot(Selfish_Bottom_Amend, aes(x=time_d, y=mean.cells, group=substrate))+geom_line(aes(color=substrate))+geom_point(aes(color=substrate))+
  scale_y_continuous(breaks=c(0, 1000000, 2000000, 3000000, 4000000, 5000000), limits=c(0, 5000000), labels = scientific)+
  facet_grid(cols = vars(stn))+labs(x="Time (days)", y=("Total cellular abundance (cells/mL)"))+
  scale_color_manual(values=c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF", "control" = "grey60"))+
  geom_errorbar(aes(ymin=mean.cells-sd.cells, ymax=mean.cells+sd.cells), width = 0.1, size=0.5)+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_blank(), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text = element_text(size=10, color="black"))

print(Bottom_cell_line)


#Bottom selfish percentage bar graph
#Remove xylan
Bottom_selfish <- Selfish_Bottom_Amend %>% filter(substrate != "xyl")

#Plot %selfish uptake
Bottom_selfish_percent <-ggplot(Bottom_selfish,aes(x=time_d, y=mean.percent, ymin=mean.percent-sd.percent, ymax=mean.percent+sd.percent, fill=substrate))+
  geom_bar(position=position_dodge(), stat="identity", width=.9, colour = "black", size=0.2)+
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30), limits=c(-2, 30))+
  facet_grid(cols = vars(stn))+
  labs(x="Time (days)", y=("Percent of selfish bacteria"))+
  scale_fill_manual(values = c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF"))+
  geom_errorbar(position=position_dodge(width = 0.9), size = 0.5)+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_blank(), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text.x = element_blank())

print(Bottom_selfish_percent)


#Bottom water extracellular enzymatic activities bar graph
#Filter the data first; remove DCM.
Bottom_EEA <- Amend_EEA %>% 
  select(Station, Depth, SpecificTreatment, Substrate, Timepoint, mean.kcrate.nM.hr, sd.kcrate.nM.hr) %>%
  filter(Depth %in% c("Bottom"))

#Plot EEA
Bottom_EEA_bar <-ggplot(Bottom_EEA, aes(x=Timepoint, y=mean.kcrate.nM.hr, ymin=mean.kcrate.nM.hr-sd.kcrate.nM.hr, ymax=mean.kcrate.nM.hr+sd.kcrate.nM.hr, fill=Substrate))+
  geom_bar(position=position_dodge(), stat="identity", width=1, colour = "black", size=0.2)+
  facet_grid(cols = vars(Station))+
  labs(x="Time (days)", y=substitute(paste("Hydrolysis rate (nmol ",L^-1, h^-1,")")))+
  scale_fill_manual(values = c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF"))+
  geom_errorbar(position=position_dodge(1), size = 0.5) +
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_text(size=9, colour="black", angle=0, hjust=0.5, vjust=1), 
                   axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_text(size=10), axis.title.y = element_text(size=10),
                   axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text.x = element_blank())

print(Bottom_EEA_bar)

#Combine these graphs using patchwork
Fig_2b <- Bottom_cell_line / Bottom_selfish_percent / Bottom_EEA_bar 

print(Fig_2b)

ggsave(filename="EN638_LV_Fig_2b.pdf", plot=Fig_2b, width=8, height=8, device=cairo_pdf)


###################################################################################                                                                                
#                                       Figure 3                                  #                                                                                
###################################################################################
#Total selfish cells (line graph); this provides line graphs on selfish bacterial
#abundance, grouped by depth. Both amended (solid line) and unamended (dotted line)
#are included in this plot. 

#Subset by depth so that the y-axis can be adjusted for each
Selfish_DCM <- subset(selfish_all, depth_name=="DCM")
Selfish_Bottom <- subset(selfish_all, depth_name=="Bottom")

#Create the line graph for DCM mesocosms
DCM_selfish_line <-ggplot(Selfish_DCM, aes(x=time_d, y=mean.selfish, group=interaction(substrate, Treatment)))+geom_line(aes(color=substrate, linetype=Treatment), size=1)+geom_point(aes(color=substrate))+
  scale_y_continuous(breaks=c(0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000), limits=c(0, 1000000), labels = scientific)+
  facet_grid(depth_name ~ stn)+labs(x="Time (days) since water collection", y=("Selfish cellular abundance (cells/mL)"))+
  scale_color_manual(values=c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF", "control" = "grey60"))+
  geom_errorbar(aes(ymin=mean.selfish-sd.selfish, ymax=mean.selfish+sd.selfish), width = 0.1)+
  theme_bw()+theme(panel.border = element_rect(fill=NA, colour = "black"), panel.grid.minor = element_blank(), 
                        panel.grid.major = element_blank(), axis.text.x  = element_blank(), 
                        axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                        axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                        axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1), 
                        strip.text = element_text(size=10, color="black"))


print(DCM_selfish_line)

#Create the line graph for bottom water mesocosms; this graph contains a y-axis break
Bottom_selfish_line <-ggplot(Selfish_Bottom, aes(x=time_d, y=mean.selfish, group=interaction(substrate, Treatment)))+geom_line(aes(color=substrate, linetype=Treatment), size=1)+geom_point(aes(color=substrate))+
  scale_y_continuous(breaks=c(0,10000,20000,30000,40000,50000,200000,300000,400000,500000,600000), limits=c(-5, 600000), labels = scientific)+
  scale_y_break(c(50000, 200000), scales = .5)+
  facet_grid(depth_name ~ stn)+labs(x="Time (days) since water collection", y=("Selfish cellular abundance (cells/mL)"))+
  scale_color_manual(values=c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF", "control" = "grey60"))+
  geom_errorbar(aes(ymin=mean.selfish-sd.selfish, ymax=mean.selfish+sd.selfish), width = 0.1)+
  theme_bw()+theme(panel.border = element_rect(fill=NA, colour = "black"), panel.grid.minor = element_blank(), 
                        panel.grid.major = element_blank(), axis.text.x  = element_text(size=9, colour="black", angle=0, hjust=0.5, vjust=1), 
                        axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                        axis.title.x = element_text(size=10), axis.title.y = element_text(size=10),
                        axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1), 
                        strip.text.x = element_blank(), strip.text.y = element_text(size=10, color="black"),
                        axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank())


print(Bottom_selfish_line)

#Combine DCM and bottom line graphs
Combo <- DCM_selfish_line / Bottom_selfish_line + plot_layout(widths = c(1,0.15)) 
Combo + plot_layout(guides = 'collect')

print(Combo)

ggsave(filename="EN638_LV_selfish_counts_line_Fig_3.pdf", plot=Combo, width=12, height=8, device=cairo_pdf)

#Edits (removing double strip.text.y labels, removing second y axis labels, removing lines in the 
#y-axis break, removing double legends) completed in Illustrator

###################################################################################                                                                                
#                                       Figure 4                                  #                                                                                
###################################################################################
#Total selfish cells (line graph) - no laminarin or pullulan. Same line graphs on 
#selfish bacterial abundance as above, but without laminarin or pullulan so that
#selfish uptake of other substrates can be more easily seen.

#Removing laminarin and pullulan
selfish_no_lam_pul_for_line <- subset(selfish_all, substrate != "lam")
selfish_no_lam_pul_for_line <- subset(selfish_no_lam_pul_for_line, substrate != "pul")

#Subset by depth so that the y-axis can be adjusted for each
Selfish_no_lam_DCM <- subset(selfish_no_lam_pul_for_line, depth_name=="DCM")
Selfish_no_lam_Bottom <- subset(selfish_no_lam_pul_for_line, depth_name=="Bottom")

#Create the line graph for DCM mesocosms
DCM_no_lam_selfish_line <-ggplot(Selfish_no_lam_DCM, aes(x=time_d, y=mean.selfish, group=interaction(substrate, Treatment)))+geom_line(aes(color=substrate, linetype=Treatment), size=1)+geom_point(aes(color=substrate))+
  scale_y_continuous(breaks=c(0,10000,20000,30000,40000,50000,60000), limits=c(0, 60000), labels = scientific)+
  facet_grid(depth_name ~ stn)+labs(x="Time (days) since water collection", y=("Selfish cellular abundance (cells/mL)"))+
  scale_color_manual(values=c("xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF", "control" = "grey60"))+
  geom_errorbar(aes(ymin=mean.selfish-sd.selfish, ymax=mean.selfish+sd.selfish), width = 0.1)+
  theme_bw()+theme(panel.border = element_rect(fill=NA, colour = "black"), panel.grid.minor = element_blank(), 
                        panel.grid.major = element_blank(), axis.text.x  = element_blank(), 
                        axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                        axis.title.x = element_blank(), axis.title.y = element_text(size=10),
                        axis.ticks.x = element_blank(), axis.ticks.y = element_line(colour="#000000", size=0.1), 
                        strip.text = element_text(size=10, color="black"))


print(DCM_no_lam_selfish_line)

#Create the line graph for bottom water mesocosms
Bottom_no_lam_selfish_line <-ggplot(Selfish_no_lam_Bottom, aes(x=time_d, y=mean.selfish, group=interaction(substrate, Treatment)))+geom_line(aes(color=substrate, linetype=Treatment), size=1)+geom_point(aes(color=substrate))+
  scale_y_continuous(breaks=c(0,10000,20000), limits=c(0, 20000), labels = scientific)+
  facet_grid(depth_name ~ stn)+labs(x="Time (days) since water collection", y=("Selfish cellular abundance (cells/mL)"))+
  scale_color_manual(values=c("xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF"))+
  geom_errorbar(aes(ymin=mean.selfish-sd.selfish, ymax=mean.selfish+sd.selfish), width = 0.1)+
  theme_bw()+theme(panel.border = element_rect(fill=NA, colour = "black"), panel.grid.minor = element_blank(), 
                        panel.grid.major = element_blank(), axis.text.x  = element_text(size=9, colour="black", angle=0, hjust=0.5, vjust=1), 
                        axis.text.y = element_text(size=9, colour="black"), plot.title = element_text(hjust = 0.5), 
                        axis.title.x = element_text(size=10), axis.title.y = element_text(size=10),
                        axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1), 
                        strip.text.x = element_blank(), strip.text.y = element_text(size=10, color="black"))


print(Bottom_no_lam_selfish_line)

#Combine DCM and bottom line graphs
No_lam_Combo <- DCM_no_lam_selfish_line / Bottom_no_lam_selfish_line 
print(No_lam_Combo)

ggsave(filename="EN638_LV_selfish_counts_no_lam_pul_line_Fig_4.pdf", plot=No_lam_Combo, width=12, height=8, device=cairo_pdf)


###################################################################################                                                                                
#                                       Figure S2                                 #                                                                                
###################################################################################
#Graphing percent selfish uptake (stagged bars)

#Remove xylan
selfish_all_for_bar <- selfish_all %>% filter(substrate != "xyl")

#Plot 
All_selfish_percent_stagged <-ggplot(selfish_all_for_bar,aes(x=time_d, y=mean.percent, fill=substrate))+
  geom_bar(position="stack", stat="identity", width=1, colour = "black", size=0.2)+
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30), limits=c(-2, 30))+
  facet_grid(depth_name + Treatment ~ stn)+
  labs(x="Time (days) since water collection", y=("Percent of selfish bacteria"))+
  scale_fill_manual(values = c("pul" = "#0000FF","lam" = "#FFC425","xyl" = "#E90000","fuc" = "#007E00","ara" = "#333333","chn" = "#66CCFF"))+
  theme_bw()+theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                   axis.text.x  = element_text(size=11, colour="black", angle=0, hjust=0.5, vjust=1), 
                   axis.text.y = element_text(size=11, colour="black"), plot.title = element_text(hjust = 0.5), 
                   axis.title.x = element_text(size=12), axis.title.y = element_text(size=12),
                   axis.ticks.x = element_line(colour="#000000", size=0.1), axis.ticks.y = element_line(colour="#000000", size=0.1),
                   strip.text = element_text(size=12, color="black"))

print(All_selfish_percent_stagged)

ggsave(filename="EN638_LV_selfish_counts_stacked_percent_Fig_S2.pdf", plot=All_selfish_percent_stagged, width=12, height=10, device=cairo_pdf)




