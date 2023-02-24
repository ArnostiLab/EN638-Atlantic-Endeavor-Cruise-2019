#Script for 16s rRNA tag analysis. Adapted from Greta Reintjes GitHub package NGSTAG-master. 


#Set the working directory to NMDS_input

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

#################################################################################################################################
#Here starts the part for statistical analysis. Input is the total community composition for each sample and the 
#environmental data that are associated to the samples. 
#The individual script sections from 'Interactive data submission' to 'NMDS' need to be run in order. Then, scripts 
#for NMDS plots from 'NMDS of surface' on can be run individually.
#################################################################################################################################
#########################################################
#Interactive data submission for statistical analysis
#########################################################

#Here, data need to be in .csv format, the separator is a ; and the decimal separator is a ,

#Community data submission
  comm <<- read.csv(file.choose(), header = TRUE, sep=",", row.names = 1, dec = ".") 

#Environmental data submission 
  contex<<- read.csv(file.choose(), header = TRUE, sep=",", row.names = 1, dec = ".")

#Factors data submission 
  biome <<- read.csv(file.choose(), sep="," ,header = TRUE, row.names = 1, dec = ".")

#########################################################
#                  Diversity indices 
#########################################################

x <- comm
y <- contex

#species number 
  sn<-specnumber(x) 						
    png(filename="Richness.png", width=1500, height=1450, res=250)
    par( bty="n")		
    plot(sn, main="Richness", col="black", xlab="", ylab="") # no x and y labels, x and y limits set.
    lines(sn,col="grey20",lwd=3) 		        		# join points with pint line, width 3
    mtext(side=1,line=2.5,cex=1,"Sample") 			# add text under side 1 of plot
    mtext(side=2,line=2.5,cex=1,"Richness") 		# add text under side 2 of plot	
    text(sn, row.names(x), cex=0.6, pos=3, col="grey20") 	# add text to dots, pos is location of text
    dev.off()
#frequencies of species
  snfreq<- specnumber(x, MARGIN=2) 			
    png(file="Genus_Frequency.png", width=1500, height=1450, res=250)
    par( bty="n")					
    plot(snfreq, main="Genus Frequency", col="black", xlab="", ylab="")
    mtext(side=1,line=2.5,cex=1,"Genus Frequencey") 
    mtext(side=2,line=2.5,cex=1,"Sample") 			      
    dev.off()
# shannon-weaver
  h<-diversity(x)						
    png(file="Shannon_Weaver.png", width=1500, height=1450, res=250)	
    par( bty="n")						
    plot(h, main="Shannon-Weaver", col="black", xlab="", ylab="") 
    mtext(side=1,line=2.5,cex=1,"Sample") 			    
    mtext(side=2,line=2.5,cex=1,"Shannon Index") 	
    text(h, row.names(x), cex=0.6, pos=1, col="grey20")
    dev.off()
# simpsons
  s<-diversity(x,index="simpson")
    png(file="Simpson.png", width=1500, height=1450, res=250)
    par( bty="n")	
    plot(s, main="Simpsons", col="black", xlab="", ylab="")
    mtext(side=1,line=2.5,cex=1,"Sample") 			
    mtext(side=2,line=2.5,cex=1,"Simpsons Index") 		
    text(s, row.names(x), cex=0.6, pos=1, col="grey20")
    dev.off()
# inverse simpsons 
  i<-diversity(x, index="invsimpson")
    png(file="Inverse_Simpson.png", width=1500, height=1450, res=250)
    par( bty="n")	
    plot(i, main="inverse Simpsons", col="black", xlab="", ylab="")
    mtext(side=1,line=2.5,cex=1,"Sample") 			
    mtext(side=2,line=2.5,cex=1,"Inverse Simpsons Index") 
    text(i, row.names(x), cex=0.6, pos=1, col="grey20")
    dev.off()
#Evenness
  j<-h/log(specnumber(x))
    png(file="Eveness_bacteria.png", width=1500, height=1450, res=250)
    par( bty="n")	
    plot(j, main="Evenness", col="black" , xlab="", ylab="") # ylim may cause issue change to just plot(j)
    mtext(side=1,line=2.5,cex=1,"Depth") 			
    mtext(side=2,line=2.5,cex=1,"Eveness") 			
    text(j, row.names(x), cex=0.6, pos=1, col="grey20")
    dev.off()
#species richness per factor (province)
    attach(biome)# attach dataframe with factor!!
    png(file="Species_richness_with_depth.png", width=1500, height=1450, res=200)
    boxplot(specnumber(comm) ~ Depth, xlab = "Depth")
    dev.off()
    
    png(file="Species_richness_between_IncubationTime.png", width=1500, height=1450, res=200)
    boxplot(specnumber(comm) ~ IncubationTime, xlab = "Incubation time [h]")
    dev.off()
    
    png(file="Species_richness_between_Substrates.png", width=1500, height=1450, res=200)
    boxplot(specnumber(comm) ~ Substrate, xlab = "Substrate")
    dev.off()

#Results return values for all diversity indecies in diversity_results.csv file
  results<-list(sn=sn, h=h, s=s, i=i, j=j)
  write.table(results, file = "diversity_results.csv", sep = "\t", row.names=TRUE, col.names=TRUE)
  return(results)

#########################################################
#                NGS species occurrence
#########################################################

#How many sites does each species occur in.
  spc.pres<-apply(x>0,2,sum)    # to get number of presences for each species. (1,0 presence and absence)
    png(filename="Species_Occurance.png", width=1500, height=1450, res=250)
    plot(sort(spc.pres), main="Species Occurance", xlab='Cumulative Number of Genera',ylab='Number of Sites', col="black", pch=16)
    dev.off()
    png(filename="Species_Occurance(log).png", width=1500, height=1450, res=250)
    plot(sort(spc.pres),log='x',main="Species Occurance(log)", xlab='Cumulative Number of Genera (log)',ylab='Number of Sites', col="black", pch=16) 					#if skewed put y axis on log scale
    dev.off()
#Abundant species (enter number of site (e.g. 35 or more))  
  spc.pres29<-	spc.pres[spc.pres>=29]
#Rare species (that occur 1 or less sites)
  spc.pres1<-	spc.pres[spc.pres<=1]    
  
#mean abundance of each species  
  tmp <- 	apply(x,2,sum)	
  spc.mean <- tmp/spc.pres     
    png(filename="Mean Coverage When a Genera Occurs.png", width=1500, height=1450, res=250)
    plot(sort(spc.mean),main="Mean Coverage When a Genera Occurs", xlab="Cumulative Number of Sites",ylab="Mean Abundance", col="black", pch=16)
    dev.off()
#mean abundance of each species against number of sites it occur in
    png(filename="Species Abundance (high).png", width=1500, height=1450, res=250)
    plot(spc.pres,spc.mean, main="Cumulative Count of Genera Against Mean Abundance", xlab="Cumulative Count of Sites", ylab="Mean Abundance", col="black", cex=1, pch=16)	
    dev.off()
  
#optional identify which are the dominant organism using listpts function
#need to click near points of interest and use ESC to stop identification of points
    plot(spc.pres,spc.mean, main="Cumulative Count of Genera Against Mean Abundance", xlab="Cumulative Count of Sites", ylab="Mean Abundance", col="grey60", cex=1)
    listpts <- identify(spc.pres,spc.mean,names(x), cex=0.3) 	# list identify abundance species
    dev.copy(png, filename="Species Abundance with Names (high).png", width=1500, height=1450, res=250)
    dev.off()
  
    plot(spc.pres,spc.mean, main="Cumulative Count of Genera Against Mean Abundance", xlab="Cumulative Count of Sites", ylab="Mean Abundance", col="grey60", cex=1)
    listpts2 <- identify(spc.pres,spc.mean,names(x), cex=0.3) 	# list identify abundance species
    dev.copy(png, filename="Species Abundance with Names(low).png", width=1500, height=1450, res=250)
    dev.off()  		
  
#Calculate species area relationship
  plt.pres<-apply(x>0,1,sum) # to calculate the number of species in each site (for 1 means for column)
  plt.sum <- apply(x,1,sum) 	# to calculate the total number of species on each site
    png(filename="De Relationship between number of Genera per site and total area.png", width=1500, height=1450, res=250)
    plot(plt.pres,plt.sum, main="Number of Genera  Against Abundance of Genera ", xlab="Number of Genera per Site", ylab="Abundance of Genera per Site", col="black")  	# see the relationship between number of species/plot and total number of species 
    res=lm(plt.sum~plt.pres)
    res
    abline(res) #draw line for relationship between species/area
    dev.off()
#calculate top 10 species from occurance  
  ctmp<-apply(x,2, sum)
  ctmp<-sort(ctmp, decreasing=TRUE)
  top10<-ctmp[1:10]
  
results<-c(spc.pres29=spc.pres29, spc.pres1=spc.pres1, listpts=listpts, top10=top10)
  write.table(results, file="species_occurance_results.csv", sep = "\t", row.names=TRUE, col.names=TRUE)
  return(results)

#########################################################
#                  Cluster analysis
#########################################################

  sample = row.names(x)
  
# bray curtis tree calculation
  bray_cluster<-vegdist(x, "bray", na.rm=TRUE)  # calculate Bray-Curtis distance among samples
  treesin=hclust(bray_cluster, method="single") 
  treeave=hclust(bray_cluster, method="average")  # cluster communities using average-linkage algorithm
  treecom=hclust(bray_cluster, method="complete")
  
    png(file="Tree_braycurtis.png", width=1500, height=1450, res=200)
    par(mfrow=c(3,1))
    plot(treesin, hang=-1, main="single")
    plot(treeave, hang=-1, main="average")  # cluster communities using average-linkage algorithm
    plot(treecom, hang=-1, main="complete")
    dev.off()
  
# Jaccard tree calculation
  jaccard_cluster<-vegdist(x, "jaccard", na.rm=TRUE)  # calculate Jaccard distance among samples
  treesin=hclust(jaccard_cluster, method="single")
  treeave=hclust(jaccard_cluster, method="average")
  treecom=hclust(jaccard_cluster, method="complete")
  
    png(file="Tree_jaccard.png", width=1500, height=1450, res=200)
    par(mfrow=c(3,1))
    plot(treesin, hang=-1, main="single")
    plot(treeave, hang=-1, main="average")
    plot(treecom, hang=-1, main="complete")
    dev.off()

    # bootstrap tree   
  t_vector=t(x)     # transpose data set
  library(pvclust)
  
  tree=pvclust(t_vector, method.hclust="ward", nboot=1000, method.dist="euclidean")
    png(file="Tree_bootstrapped_all_mesos.png", width=1500, height=1450, res=200)
    par(mfrow=c(1,1))
    plot(tree, hang=-1, main="boot") 								
    pvrect(tree, alpha=.95, border=2, lwd=3)		# bootstrap							
    dev.off()
  
#########################################################
#                   NMDS
#########################################################
  
#Calculate Bray-Curtis distances among samples based on species P/A/abundance 
  bray_all <- vegdist(x, method="bray")  
  mds_bray<-metaMDS(x, k=2, distance="bray") # robustness
    png(file="stressplot.png", width=1500, height=1450, res=200)   # Assess goodness of ordination fit (stress plot)
    stressplot(mds_bray, bray_all, pch=".", p.col="black", l.col="black") 
    dev.off()
#Calculate Jaccard distances among samples based on species P/A/abundance   
  jaccard_all <- vegdist(x, method="jaccard")  
  mds_jaccard<-metaMDS(x, K=2, distance="jaccard")
    png(file="stressplot_jaccard.png", width=1500, height=1450, res=200)
    stressplot(mds_jaccard, jaccard_all, pch="." , p.col="black", l.col="black")  
  dev.off()
  
  y<-contex

  rank<-rankindex(y, x, c("euc","man","bray","jac","kul")) 
  ef<-envfit(mds_bray, y, permu=999)

  
  #########################################################
  #                 NMDS of Main Mesocosms
  #########################################################
  #Plot NMDS basic plot 
  
  png(file="NMDS_EN638_main_mesos_Basic_vegan_Bray.png", width=1500, height=1450, res=200)
  par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
  plot(mds_bray, display = "sites",   type = "p", xlim=c(-1,1), ylim=c(-1, 1))       # plot site scores as text
  dev.off()  
  
  #Plot NMDS showing specific factors (Station with polygon)  
  
  y<-contex
  
  rank<-rankindex(y, x, c("euc","man","bray","jac","kul")) 
  ef<-envfit(mds_bray, y, permu=999)
  
  png(file="NMDS_EN638_main_mesos_Station_Poly_Bray.png", width=1500, height=1500, res=200)
  
  setEPS()
  postscript("NMDS_EN638_main_mesos_Station_Poly_Bray.eps")
  
  par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
  mds.fig <- ordiplot(mds_bray, type = "none", xlim=c(-2.0,2.0), ylim=c(-2, 2), cex.axis=1)   # don't plot anything yet
  points(mds.fig, "sites", pch = 19, cex = 2, col = "red", select = Station == 
           "18")   # plot just the samples, colour by habitat, pch=19 means plot a circle
  points(mds.fig, "sites", pch = 19, cex = 2, col = "green", select = Station == 
           "19")
  points(mds.fig, "sites", pch = 19, cex = 2, col = "blue", select = Station == 
           "20")
  ordihull(mds_bray, Station, draw = "polygon", label = TRUE, cex=1.5)
  #plot(ef,p.max=0.05,cex=0.8, col="blue")
  legend("topright", legend = c("Station 18", "Station 19", "Station 20"), # select cluster number 
         col = c("red", "green", "blue"), bty = "n", cex = 1.5, pch = 19, title ="Station") # add confidence ellipses around habitat types
  
  dev.off()
  
  #Plot NMDS showing specific factors (Amendment_type with polygon)  
  
  png(file="NMDS_EN638_main_mesos_Amendment_type_Poly_Bray.png", width=1500, height=1500, res=200)
  
  setEPS()
  postscript("NMDS_EN638_main_mesos_Amendment_type_Poly_Bray.eps")
  
  par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
  mds.fig <- ordiplot(mds_bray, type = "none", xlim=c(-2.0,2.0), ylim=c(-2, 2), cex.axis=1)   # don't plot anything yet
  points(mds.fig, "sites", pch = 19, cex = 2, col = "gold", select = Amendment_type == 
           "Amended")   # plot just the samples, colour by habitat, pch=19 means plot a circle
  points(mds.fig, "sites", pch = 19, cex = 2, col = "steelblue1", select = Amendment_type == 
           "Unamended")
 
  # points(mds.fig, "sites", pch = 19, cex = 2, col = "black", select = Amendment_type == 
  #          "Control")
  
  legend("topright", legend = c("Amended", "Unamended"), # select cluster number 
         col = c("gold", "steelblue1"), bty = "n", cex = 1.5, pch = 19, title = "Amendment")# add confidence ellipses around habitat types
  ordihull(mds_bray, Amendment_type, draw = "polygon", label = TRUE, cex=1.5) # add confidence ellipses around habitat types
  
  dev.off()
  
  #Plot NMDS showing specific factors (Incubation time with polygon)  
  
  png(file="NMDS_EN638_main_mesos_IncubationTime_Poly_Bray.png", width=1500, height=1500, res=200)
  
  setEPS()
  postscript("NMDS_EN638_main_mesos_IncubationTime_Poly_Bray.eps")
  
  par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
  mds.fig <- ordiplot(mds_bray, type = "none", xlim=c(-2.0,2.0), ylim=c(-2, 2), cex.axis=1)   # don't plot anything yet
  points(mds.fig, "sites", pch = 19, cex = 2, col = "grey10", select = IncubationTime == 
           "0")
  points(mds.fig, "sites", pch = 19, cex = 2, col = "grey20", select = IncubationTime == 
           "3")
  points(mds.fig, "sites", pch = 19, cex = 2, col = "grey40", select = IncubationTime == 
           "5")
  points(mds.fig, "sites", pch = 19, cex = 2, col = "grey60", select = IncubationTime == 
           "10")
  points(mds.fig, "sites", pch = 19, cex = 2, col = "grey80", select = IncubationTime == 
           "15")
  #plot(ef,p.max=0.05,cex=0.8, col="blue")
  legend("topright", legend = c("0", "3", "5","10", "15"), # select cluster number 
         col = c("grey10", "grey20", "grey40", "grey60", "grey80"), bty = "n", cex = 1.5, pch = 19, title = "IncubationTime")
  ordihull(mds_bray, IncubationTime, draw = "polygon", label = TRUE, cex=1) # add confidence ellipses around habitat types
  
  dev.off()
  
  #Plot NMDS showing specific factors (Depth with polygon)  
  
  png(file="NMDS_EN638_main_mesos_Depth_Poly_Bray.png", width=1500, height=1500, res=200)
  
  setEPS()
  postscript("NMDS_EN638_main_mesos_Depth_Poly_Bray.eps")
  
  par(oma=c(1,1,.5,.5), mar=c(1.5,1.5,0,0))
  mds.fig <- ordiplot(mds_bray, type = "none", xlim=c(-2.0,2.0), ylim=c(-2, 2), cex.axis=1)   # don't plot anything yet
  points(mds.fig, "sites", pch = 19, cex = 2, col = "red", select = Depth == 
           "A")
  points(mds.fig, "sites", pch = 19, cex = 2, col = "green", select = Depth == 
           "B")
  points(mds.fig, "sites", pch = 19, cex = 2, col = "blue", select = Depth == 
           "C")
  #plot(ef,p.max=0.05,cex=0.8, col="blue")
  legend("topright", legend = c("A", "B", "C"), # select cluster number 
         col = c("red", "green", "blue"), bty = "n", cex = 1.5, pch = 19, title = "Depth")
  ordihull(mds_bray, Depth, draw = "polygon", label = TRUE, cex=1) # add confidence ellipses around habitat types
  
  dev.off()
  
  ### ANOSIM analysis ####
  ano_Station = anosim(comm, y$Station, distance = "bray", permutations = 9999)
  ano_Station
  
  ano_Amendment_type = anosim(comm, y$Amendment_type, distance = "bray", permutations = 9999)
  ano_Amendment_type
  
  ano_IncTime = anosim(comm, y$IncubationTime, distance = "bray", permutations = 9999)
  ano_IncTime
  
  ano_Depth = anosim(comm, y$Depth, distance = "bray", permutations = 9999)
  ano_Depth
  
  