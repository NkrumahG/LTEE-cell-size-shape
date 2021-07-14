rm(list=ls())
library(tidyverse)
library(cowplot)
library(dplyr)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(grid)
library(datarium)
library(ggpubr)
library(rstatix)
library(car)


setwd("/Users/nkrumahgrant/Desktop/Cell Morphology/Manuscript-Cell-size-full-submission/Clones-mixed-populations/")
#####functions####

#mod.boxplot: constructs boxplots with fifth to ninety-fifth percentiles as whisker points
mod.boxplot <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#####Data import#####

###### Clones: coulter counter data from stationary phase cultures#####
Coult <-  list.files(pattern = "*clones.CSV") #119 coulter files -  missing Ara−2 50K (may need to get one more measurement)

#Data Wrangling
CoultCol.names <- c( "x", "seconds", "fL", "um", "population", "generation", "block")
Coulter.data <- as.data.frame(setNames(do.call(rbind,Map('cbind', lapply(Coult, read.csv, skip = 474), V4=Coult, V5=Coult, V6=Coult)), CoultCol.names)) %>%  subset(fL > 0 & fL < 6)

#Add section header here where I talk about how I created the data frame 

#Extrating generation information from the file name
Coulter.data$generation <- as.factor(gsub("_.*","",Coulter.data$generation))

#Extracting population information from the file name and change data order 
Coulter.data$population <- as.factor(gsub(".*K_|_2018.*","",Coulter.data$population))
Coulter.data$population <- factor(Coulter.data$population, levels = c("606", "M1", "M2", "M3", "M4", "M5", "M6", "607","P1","P2","P3","P4","P5","P6"))

#Assigning a block number based upon date completed
Coulter.data$block <- as.factor(gsub("^.*2018-10-04.*$","1",Coulter.data$block))
Coulter.data$block <- as.factor(gsub("^.*2018-10-30.*$","2",Coulter.data$block))
Coulter.data$block <- as.factor(gsub("^.*2018-10-31.*$","3",Coulter.data$block))
Coulter.data$block <- as.factor(gsub("^.*2018-11-01.*$","4",Coulter.data$block))
Coulter.data$block <- as.factor(gsub("^.*2018-11-08.*$","5",Coulter.data$block))
Coulter.data$block <- as.factor(gsub("^.*2018-11-09.*$","6",Coulter.data$block))

#Create a new column that has the population and not the block 
Coulter.data$PopGen <- as.factor(paste(Coulter.data$population,Coulter.data$generation, sep="_"))

#Setting factor levels so that data is plotted in the correct order 
Coulter.data$generation <- factor(Coulter.data$generation, c("0K", "2K", "10K", "50K"))

#Creating a column with Freq of counts for each population at each generation 
PopGen.Freq <- as.data.frame(table(Coulter.data$PopGen)) 
names(PopGen.Freq)[1] = "PopGen"
Coulter.data <- left_join(Coulter.data, PopGen.Freq)

###### Clones: segger data from stationary phase cultures#####

#Data Wrangling
Micro <- list.files(pattern = "*.txt") #114 microscopy files - Missing: 3 from rep 1, 1 from rep 2 amd 2 from rep 3
Segger.data <- as.data.frame(do.call(rbind,Map('cbind', lapply(Micro, read.csv))))

#List of colnames so that I am able to establish what columns to include in the final dataset
#ExclSegg <- as.list(colnames(Segger.data))
#View(ExclSegg)

#Excluding all unnecessary columns in the Segger data set
Segger.data<- Segger.data[,c(1:3, 5, 33, 41)] #INCLUDE THE COLUMNS WITH LONG AND SHORT AXIS.... 

#Creating new columns with folder information so that I can gsub indentifying information
#for each of the population and replicate
Segger.data$population <- Segger.data$Folder
Segger.data$generation <- Segger.data$Folder

#Total number of cells that went into the analysis?
celltotal <- Segger.data %>%  group_by(population, generation) %>% 
  summarise (n=n())

sum(celltotal$n) #87811 

#Assigning a block number based upon date completed
#Minus populations were completed in the first three sets and plus populations in the latter 3
#Both ancestors were included in each rep
Segger.data$Folder <- as.factor(gsub("^.*20181004.*$","1",Segger.data$Folder))
Segger.data$Folder<-  as.factor(gsub("^.*20181030.*$","2",Segger.data$Folder))
Segger.data$Folder <- as.factor(gsub("^.*20181031.*$","3",Segger.data$Folder))
Segger.data$Folder <- as.factor(gsub("^.*20181101.*$","4",Segger.data$Folder))
Segger.data$Folder <- as.factor(gsub("^.*20181108.*$","5",Segger.data$Folder))
Segger.data$Folder <- as.factor(gsub("^.*20181109.*$","6",Segger.data$Folder))

#Extracting population information 
Segger.data$population<- as.factor(gsub(".*_ecoli_", "",Segger.data$population))
Segger.data$population<- as.factor(gsub("_2K.*|_10K.*|_50K.*", "",Segger.data$population))

#Extracting generation information
Segger.data$generation <- (str_extract(Segger.data$generation, "(?<=_)\\d+K$"))
Segger.data$generation <- as.factor(ifelse(is.na(Segger.data$generation),"0K", Segger.data$generation))

#Setting factor levels so that data is plotted in the correct order 
Segger.data$generation <- factor(Segger.data$generation, c("0K", "2K", "10K", "50K"))

#Section where I abstract the medians of the median

#For loop for calculating the median of the medians for each unique population.
#1) calculate the median for each replicate
#2) Take the median of the medians and return that value 
#3) Linear regression of medians of the segger data against the medians of the coulter data

#1) Calculate the median for each replicate

#Made a "New" Column that is a concatenation of population, generation and block number seprated 
#by an underscore
Segger.data$New <- as.factor(paste(Segger.data$population,Segger.data$generation,Segger.data$Folder, sep="_"))

#Creating a data frame with summary stats for median Segger data 
Seg.Medians <- as.data.frame(Segger.data %>% 
                               group_by(New) %>% 
                               dplyr::summarize(SegMedian = as.numeric(median(MeshVolume, na.rm =T))))

#Create new column in Seg.Medians dataframe with name that removes block factor level 
Seg.Medians$New2 <- (str_extract(Seg.Medians$New, (".*\\d+_\\d+K") ))

#Calculating median of medians for Segger data 
Seg.Med.Medians <- as.data.frame(Seg.Medians %>% 
                                   group_by(New2) %>% 
                                   dplyr::summarize(SegMed.Medians = median(SegMedian, na.rm = T)))

#Calculating average of medians for Segger data 
Seg.Avg.Medians <- as.data.frame(Seg.Medians %>% 
                                   group_by(New2) %>% 
                                   dplyr::summarize(SegMed.Medians = mean(SegMedian, na.rm = T)))

#Made a new column in the Coulter.data frame that is a concatenation of columns in the dataset
Coulter.data$New <- as.factor(paste(Coulter.data$population,Coulter.data$generation,Coulter.data$block, sep="_"))

#Creating a data frame with summary stats for median Coulter data 
#What is the median cell volume for each replicate?
Coul.Medians <- as.data.frame(Coulter.data %>% 
                                group_by(New) %>% 
                                dplyr::summarize(CoulMedian = as.numeric(median(fL, na.rm =T))))

#Create new column in Coul.Medians dataframe with name that removes block factor level 
Coul.Medians$New2 <- (str_extract(Coul.Medians$New, (".*\\d+_\\d+K") ))

#Now create two dataframes where I:
#1) calculate the median of the median cell volumes
#2) calculate the average of the median cell volumes. 

#Calculating median of medians for Coulter data 
Coul.Med.Medians <- as.data.frame(Coul.Medians %>% 
                                    group_by(New2) %>% 
                                    dplyr::summarize(CoulMed.Medians = median(CoulMedian, na.rm = T)))

#Calculating average of the medians for the coulter data
Coul.Avg.Medians <- as.data.frame(Coul.Medians %>% 
                                    group_by(New2) %>% 
                                    dplyr::summarize(CoulAvg.Medians = mean(CoulMedian, na.rm = T)))

#Coul.Med.Medians$MixedOrClone <- paste("Clone")
#In the following lines of code I change how replicates are named in the Coulter dataframe so that it matches the 
#Segger naming scheme. In this way I can merge the dataframes together, matching by rows, so that I can perform the
#linear regression. 
Coul.Med.Medians$New3 <- paste("ARA", Coul.Med.Medians$New2,sep="_")
Coul.Med.Medians$New3 <- ifelse(Coul.Med.Medians$New3 == "ARA_606_0K", "REL606_0K", Coul.Med.Medians$New3)
Coul.Med.Medians$New3 <- ifelse(Coul.Med.Medians$New3 == "ARA_607_0K", "REL607_0K", Coul.Med.Medians$New3)

#Indexing the dataframe above to contain only the columns of interest and then reordering columns
Coul.Med.Medians <- (Coul.Med.Medians[,2:3])
Coul.Med.Medians <- Coul.Med.Medians[,c(2,1)]

#Merging Median of Median dataframes together. 
#First I will need to give the dfs one column with same name. First two lines. 
#The following line of code merges the data by same column name 
colnames(Coul.Med.Medians)[1] <- "population"
colnames(Seg.Med.Medians)[1] <- "population"

Coul.Seg.MedianDF <- left_join(Coul.Med.Medians, Seg.Med.Medians)

#Create a new column in data frame with generation as factor 
Coul.Seg.MedianDF$Generation <- (str_extract(Coul.Seg.MedianDF$population, "(?<=_)\\d+K$"))

#Setting facctor levels 
Coul.Seg.MedianDF$Generation <-factor(Coul.Seg.MedianDF$Generation ,levels = c("0K", "2K","10K","50K"))

#Correlation plot 

#Generate correlation plot (GGplot)
#Fit a linear model to data
coef(lm(SegMed.Medians ~ CoulMed.Medians , data = Coul.Seg.MedianDF ))

#Plot the median value of medians for each population. Unbiased as opposed to the mean that 
#can skew things with strong outliers.

#size 6.33x5.03
ggplot(Coul.Seg.MedianDF, aes(CoulMed.Medians, SegMed.Medians,color = Generation)) + 
  geom_point(size=2) +
  geom_smooth(method = lm, se = F, colour = "red", lty = 2) +
  labs(y= "Cell volume microscopy (a.u.)", x= "Cell volume Coulter counter (fL)") +
  scale_y_continuous(breaks = c(1000,2000,3000,4000,5000,6000,7000,8000) ,limits = c(1200,8000)) +
  scale_x_continuous(breaks = c(0,0.25,.5,.75,1,1.25,1.5,1.75,2.0), limits = c(0.2,1.75)) +
  scale_color_discrete(name = "Generation", labels = c("0k", "2k", "10k", "50k")) +
  #theme(legend.position = "bottom") +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(
          fill = "white",
          size = 1)) +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = margin(l = 5, r = 5))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = margin(b = 5, t =5))) +
  theme(panel.border = element_blank(), axis.line = element_line())

#Correlation test between median of median cell volumes 
ggqqplot(Coul.Seg.MedianDF$CoulMed.Medians, ylab = "Coulter median of medians cell volume (fL)")
ggqqplot(Coul.Seg.MedianDF$SegMed.Medians, ylab = "Microscopy median of medians cell volume (fL)")

shapiro.test(Coul.Seg.MedianDF$CoulMed.Medians) #W = 0.82914, p-value = 4.343e-05; not normally distributed; can't use Pearsons
shapiro.test(Coul.Seg.MedianDF$SegMed.Medians) #W = 0.72109, p-value = 3.66e-07

#Highly significant difference between the median cell volume of clones and the whole populations 
cor.test(x=Coul.Seg.MedianDF$CoulMed.Medians,y=Coul.Seg.MedianDF$SegMed.Medians, method = "kendall", exact = FALSE) 
#tau = 0.5494663, z = 4.8531, p-value = 1.215e-06

#####Histograms#####

#Coulter data

## using geom density because the transparency and smoothing of the histograms,
##the count argument in the aes is multiplied by .0125 so that the y axis reflects n observations
### in the data set.

#Change how populations are named
Coulter.data$population <- gsub("M", "Ara−",Coulter.data$population)
Coulter.data$population <- gsub("P", "Ara+",Coulter.data$population)
Coulter.data$population <- gsub("606", "REL606",Coulter.data$population)
Coulter.data$population <- gsub("607", "REL607",Coulter.data$population)

#Coulter data distrbutions 
ggplot(Coulter.data, aes(x=fL, color = generation, fill = generation, y= .0125*..count..)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +theme(
          panel.background = element_rect(fill = "white"),
          plot.margin = margin(2, 2, 2, 2, "cm"),
          plot.background = element_rect(
            fill = "white",
            colour = "black",
            size = 1 
          )) +
  geom_density(alpha = 0.5, position = "identity") + 
  xlim(c(.15,2)) +
  scale_y_continuous() + 
  facet_wrap(~population, nrow = 2)  

####### Using Segger.data and Coulter.data files to create summary tables for the quantiles
#I used to construct boxplots. I then take the median of each quantile, for each population and use that
#to construct new boxplots. 

#Coulter data
Coult.quants <- tapply(Coulter.data$fL, Coulter.data$New, mod.boxplot)  #function mod.boxplot calculates the quantiles 
Coult.quants <- as.data.frame(do.call(rbind, Coult.quants)) 
Coult.quants <- rownames_to_column(Coult.quants, var = "rowid")
colnames(Coult.quants) <- c("Population","Fifth","Twenty_fifth","Median","Seventy_fifth","Ninety_fifth")
Coult.quants$Population <- as.factor((str_extract(Coult.quants$Population, (".*\\d+_\\d+K"))))
Coult.quants$Generation <-  (str_extract(Coult.quants$Population, "(?<=_)\\d+K$")) 

#View(Coult.quants)
#Calculate MEDIANS of each quantile for each population sample and then construct a dataframe with these values
#This is to keep things consistent with the later work. And to keep things consistent with median of median 
#analyses. Didn't change assigned variables so that down stream code will work with the variable names.  
#Note made: 05/09/2020

#Dataframe w/ medians of each quantile
Coult.quants.means <- Coult.quants %>% 
  group_by(Population) %>% 
  summarise(Fifth=median(Fifth), Twenty_fifth=median(Twenty_fifth), Median=median(Median), 
            Seventy_fifth=median(Seventy_fifth), Ninety_fifth=median(Ninety_fifth))

#Here I caluclute the grand means of the mean of the median cell volumes
#In simpler terms, I calculate the means of the median cell volume for each
#population. I then take those values as my estimate and then caluclate the grand mean
#of those estimates at each time point.

#Summary table 1 in supplement. Values are written in text so may be unnecessary. 

#05/09/2020 instaed of taking the mean of the medians I will be looking at the median 
#of the medians. This differs from v01 of this code. where I calculated the mean of the medians. 

Coult.quants.means.tab <- full_join(Coult.quants %>% filter(Generation == "0K") %>%
  group_by(Generation) %>%
  summarise(mean.Medians=mean(Median), S.D.=sd(Median), n=n()) %>% 
  mutate(CIlow = mean.Medians - abs(qt(0.05/2, n-1)*(S.D./sqrt(n)))) %>%
  mutate(CIHigh = mean.Medians + abs(qt(0.05/2, n-1)*(S.D./sqrt(n)))) %>%
  droplevels(),
  
Coult.quants.means.tab <- as.data.frame(Coult.quants %>% #[-c(36:38),] %>% #excluding Ara−3
  group_by(Population, Generation) %>% 
  filter(Generation != "0K") %>% 
  summarise(median.median = median(Median)) %>% #here take the median of the medians to keep consistent 05/09/20
  group_by(Generation) %>% 
  summarise(mean.Medians=mean(median.median), S.D.=sd(median.median), n=n()) %>% 
  mutate(CIlow = mean.Medians - abs(qt(0.05/2, n-1)*(S.D./sqrt(n)))) %>% 
  mutate(CIHigh = mean.Medians + abs(qt(0.05/2, n-1)*(S.D./sqrt(n))))) %>% 
  droplevels())

Coult.quants.means.tab$Generation <- factor(Coult.quants.means.tab$Generation, c("0K", "2K", "10K", "50K"))

#Grand mean of median cell volume for clones 
grand.mean.vol.clones <- as.data.frame(full_join(Coult.quants %>% filter(Generation == "0K") %>% 
   group_by(Population, Generation) %>% 
   transmute(median.median = Median) %>%  #here I call the medians "median of median" for the ancestor. Treating the 12 estimates as I did the other populations
  droplevels(), 
Coult.quants[-c(36:38),] %>% #excluding Ara−3
 group_by(Population, Generation) %>% 
  filter(Generation != "0K") %>% 
  summarise(median.median = median(Median)) %>% #The means are caluclated within the ggplot using the stat_summary function. 
  droplevels()))

#Order factor levels 
grand.mean.vol.clones$Generation <- factor(grand.mean.vol.clones$Generation, c("0K", "2K", "10K", "50K"))

#Comparisons interested in for t.tests
#gm.vol.clones_comparisons <- list(c("0K", "2K"), c("2K","10K"), c("10K", "50K"))
#plot size 6.33x5.03
grand.mean.vol.clones %>% 
  ggplot(aes(y = median.median, x = Generation)) +
  labs(x= "Generation", y= "Cell volume of clones (fL)") +
  scale_y_continuous(breaks = seq(.2,1.4,.2)) +
  scale_x_discrete(labels = c("0k", "2k", "10k", "50k")) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "point", position="identity",size=2) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.2) +
  geom_bracket(xmin = "0K", xmax= "2K", label = "0.0067", y.position = 0.65)+
  geom_bracket(xmin = "2K", xmax= "10K", label = "0.1725", y.position = 0.75)+
  geom_bracket(xmin = "10K", xmax= "50K", label = "0.0052", y.position = 1.0)+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black",size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) 
 
#####Statistical analyses comparing differences in the grand mean of median cell volumes of clones#####
grand.mean.vol.clones <- full_join(Coult.quants %>% filter(Generation == "0K") %>% 
                                     group_by(Population, Generation) %>% 
                                     transmute(median.median = Median) %>% 
                                     droplevels(), 
                                   
                                   Coult.quants %>% #[-c(36:38),] %>% #exclude Ara−3
                                     group_by(Population, Generation) %>% 
                                     filter(Generation != "0K") %>%
                                     summarise(median.median = median(Median)) %>%
                                     droplevels())

#Normality test
with(grand.mean.vol.clones, shapiro.test(median.median[Generation == "0K"])) #p-value = 0.0003829
with(grand.mean.vol.clones, shapiro.test(median.median[Generation == "2K"])) #p-value = 0.03757
with(grand.mean.vol.clones, shapiro.test(median.median[Generation == "10K"])) #p-value = 0.1456
with(grand.mean.vol.clones, shapiro.test(median.median[Generation == "50K"])) #p-value = 0.1741; Ara−3 excluded

#filter dataset to include only the data for the comparisons to be made with the t.test. 

#Bartlett's test for the homogenity of variance 
bartlett.test(median.median ~ Generation, data = grand.mean.vol.clones) #df = 3, p-value = 0.05803; variances homogenous

#3 one-tailed tests: I am concerned with whether preceding generations have greater cell
#volumes than time before.

gm.vol.clones2K <- grand.mean.vol.clones %>% filter(Generation == "0K" | Generation == "2K") %>% 
                   droplevels()

gm.vol.clones10K <- grand.mean.vol.clones %>% filter(Generation == "2K" | Generation == "10K") %>% 
                    droplevels()

gm.vol.clones10K$Generation <- factor(gm.vol.clones10K$Generation, levels = c("2K", "10K")) 

gm.vol.clones50K <- grand.mean.vol.clones %>% filter(Generation == "10K" | Generation == "50K") %>% 
                    droplevels()

gm.vol.clones0_10K <- grand.mean.vol.clones %>% filter(Generation == "0K" | Generation == "10K") %>% 
  droplevels()

gm.vol.clones0_50K <- grand.mean.vol.clones %>% filter(Generation == "0K" | Generation == "50K") %>% 
  droplevels()

#View(gm.vol.clones0_50K)

#####Statistical analyses for figure 3#####
#T.tests; all one tailed. Median of medians. 
t.test(median.median ~ Generation, data = gm.vol.clones2K, alternative = "less", var.equal = T, paired = T) #t = -2.9436, df = 11, p-value = 0.00668
t.test(median.median ~ Generation, data = gm.vol.clones10K, alternative = "less", var.equal = T, paired = T) #t = -0.9869, df = 11, p-value = 0.1725
t.test(median.median ~ Generation, data = gm.vol.clones50K, alternative = "less", var.equal = T, paired = T) 
#Add %>% filter(Population != "M3_50K" & Population != "M3_10K") behind gm.vol.clones50K to get p-value when Ara−3 is excluded at 10k and 50k generations. 
#t = -2.79, df = 11, p-value = 0.008794 #including Ara−3
#t = -3.1422, df = 10, p-value = 0.005236  #excluding Ara−3 at 10K and 50K

t.test(median.median ~ Generation, data = gm.vol.clones0_10K, alternative = "less", var.equal = T, paired = T) #t = -3.4871, df = 22, p-value = 0.001044
t.test(median.median ~ Generation, data = gm.vol.clones0_50K, alternative = "less", var.equal = T, paired = T)
#Add %>% filter(Population != "M3_50K" , median.median != "0.647") after "gm.vol.clones0_50K" (line immediately above) to get 
#p-value when Ara−3 is excluded at 50k and the highest value for the ancestor. 
#t = -4.3246, df = 11, p-value = 0.0006027 including Ara−3
#t = -5.5905, df = 10, p-value = 0.0001154 excluding Ara−3 and the greatest ancestor measurement 

#Extract Generation Sampled from 
Coult.quants.means$Generation <-  (str_extract(Coult.quants.means$Population, "(?<=_)\\d+K$")) 

#Extract Population name 
Coult.quants.means$PopGen <- as.factor(gsub("_0K.*|_2K.*|_10K.*|_50K.*", "",Coult.quants.means$Population))

#Setting the order of the data for plotting
Coult.quants.means$Generation <-factor(Coult.quants.means$Generation ,levels = c("0K", "2K","10K","50K"))

#Change column order and rename columns 
Coult.quants.means <- Coult.quants.means[,c(1,8,7,2,3,4,5,6)]
colnames(Coult.quants.means) <- c("PopGen","Population","Generation","Fifth","Twenty_fifth","Median","Seventy_fifth","Ninety_fifth")

#Generate boxplots using the means of the three replicates for each population an at each timepoint
Coult.quants.means$Population <- gsub("M", "Ara−",Coult.quants.means$Population)
Coult.quants.means$Population <- gsub("P", "Ara+",Coult.quants.means$Population)
Coult.quants.means$Population <- gsub("606", "REL606",Coult.quants.means$Population)
Coult.quants.means$Population <- gsub("607", "REL607",Coult.quants.means$Population)

Exp.Coult <- Coult.quants.means[rep(1:nrow(Coult.quants.means[1:2,]),each=6),] 
Coult.quants.means2 <- full_join(Coult.quants.means, Exp.Coult)
Coult.quants.means2[1:12,2] <-  paste(c("Ara−1","Ara−2","Ara−3","Ara−4","Ara−5","Ara−6","Ara+1","Ara+2","Ara+3","Ara+4","Ara+5","Ara+6"))

#[-c(1,2),] This index prevents the ancestors from having their own panel, with a single boxplot in it, 
#when using coult quants means data frame. Used during code development. 

Coult.quants.means <- mutate(Coult.quants.means, IQR = Seventy_fifth-Twenty_fifth)
Coult.quants.means <- mutate(Coult.quants.means, IQR_fold_change = (100*IQR/Twenty_fifth))

#####Boxplots of cell size changes in clones####
#Recall these boxplots are generated using the median of each quantile. 

#Coulter counter
#size 9.5x6
ggplot(Coult.quants.means2, aes(x=Population, fill = Generation)) +
  geom_boxplot(aes(ymin=Fifth, lower=Twenty_fifth, middle=Median, upper=Seventy_fifth, ymax=Ninety_fifth),
               stat = "identity") + 
  labs(y="Cell volume (fL)") +
  facet_wrap(~Population, scales ="free_x", nrow=2) +
  scale_y_continuous(breaks = seq(0.0,9.6,0.8), limits = c(0,4.8)) +
  coord_cartesian(ylim=c(0,3.2)) +
  scale_fill_discrete(labels = c("0k", "2k", "10k", "50k")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(
          fill = "white",
          size = 1)) +
  theme(axis.text.y = element_text(color = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(legend.text = element_text(color = "black", size = 14)) +
  theme(legend.title = element_text( size = 14)) +
  theme(strip.text = element_text( size=14)) 
  

#####New section header here#####

#Segger data distributions
Segger.data$population <- gsub("ARA_M", "Ara−",Segger.data$population)
Segger.data$population <- gsub("ARA_P", "Ara+",Segger.data$population)

ggplot(Segger.data, aes(x=(LongAxis/ShortAxis), color = generation, fill = generation, y=.0125*..count..))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +theme(
          panel.background = element_rect(fill = "white"),
          plot.margin = margin(2, 2, 2, 2, "cm"),
          plot.background = element_rect(
            fill = "white",
            colour = "black",
            size = 1 
          )) +
  geom_density(alpha = 0.5, position = "identity") +
  facet_wrap(~population, nrow = 2)

#Generating dataframes for segger data
Seg.quants <- Segger.data %>% select(population, generation, New, MeshVolume)
Seg.quants <- as.tibble(Seg.quants[complete.cases(Seg.quants), ])
Seg.quants <- (as_tibble(tapply(Seg.quants$MeshVolume, Seg.quants$New, mod.boxplot)))
Seg.quants <- as.data.frame(t(Seg.quants))
Seg.quants <- rownames_to_column(Seg.quants, var = "rowid")
colnames(Seg.quants) <- c("Population","Fifth","Twenty_fifth","Median","Seventy_fifth","Ninety_fifth")
Seg.quants$Population <- as.factor((str_extract(Seg.quants$Population, ("(?<=L|_)\\w+\\d+K"))))
Seg.quants$Generation <-  (str_extract(Seg.quants$Population, "(?<=_)\\d+K$")) 

#Caluclate medians in dataframe
Seg.quants.means <- Seg.quants %>% 
  group_by(Population) %>% 
  summarise(Fifth=median(Fifth), Twenty_fifth=median(Twenty_fifth), Median=median(Median), #Boxplots will now plot the medians
            Seventy_fifth=median(Seventy_fifth), Ninety_fifth=median(Ninety_fifth))

#Identifying information for the data
Seg.quants.means$Generation <-  (str_extract(Seg.quants.means$Population, "(?<=_)\\d+K$")) 

#Renaming column 
Seg.quants.means <- Seg.quants.means %>% rename(PopGen = Population)

#Extract Population name 
Seg.quants.means$Population <- as.factor(gsub("_0K.*|_2K.*|_10K.*|_50K.*", "",Seg.quants.means$PopGen))

#Setting the order of the data for plotting
Seg.quants.means$Generation <-factor(Seg.quants.means$Generation ,levels = c("0K", "2K","10K","50K"))

#Change column order and rename columns 
Seg.quants.means <- Seg.quants.means[,c(1,8,7,2,3,4,5,6)]
colnames(Seg.quants.means) <- c("PopGen","Population","Generation","Fifth","Twenty_fifth","Median","Seventy_fifth","Ninety_fifth")

#Change how populations are named
Seg.quants$Population <- gsub("M", "Ara−",Seg.quants$Population)
Seg.quants$Population <- gsub("P", "Ara+",Seg.quants$Population)
Seg.quants$Population <- gsub("606", "REL606",Seg.quants$Population)
Seg.quants$Population <- gsub("607", "REL607",Seg.quants$Population)

Seg.quants.means$Population <- gsub("M", "Ara−",Seg.quants.means$Population)
Seg.quants.means$Population <- gsub("P", "Ara+",Seg.quants.means$Population)
Seg.quants.means$Population <- gsub("606", "REL606",Seg.quants.means$Population)
Seg.quants.means$Population <- gsub("607", "REL607",Seg.quants.means$Population)

#Made this one so that I can include the ancestor boxplot in each evolved clone panel
#Exp -> I expanded the rows correspondig to the two ancestors 6 times, 1 for each population
#Then I full joined the two data sets
#For visualization purposes!!!!!

Exp.Seg <- Seg.quants.means[rep(1:nrow(Seg.quants.means[1:2,]),each=6),]  
Seg.quants.means2 <- full_join(Seg.quants.means, Exp.Seg)
Seg.quants.means2[1:12,2] <-  paste(c("Ara−1","Ara−2","Ara−3","Ara−4","Ara−5","Ara−6","Ara+1","Ara+2","Ara+3","Ara+4","Ara+5","Ara+6"))

#Boxplots for micoscoopy
#size 9.5x6
Seg.quants.means2 %>%
  ggplot(aes(x=Population, fill = Generation)) +
  geom_boxplot(aes(ymin=Fifth, lower=Twenty_fifth, middle=Median, upper=Seventy_fifth, ymax=Ninety_fifth),
               stat = "identity") + 
  labs(y= "Microscopy volume (a.u.)") +
  scale_y_continuous(breaks = c(1400,2500,3500,4500,5500,6500,7500,8000), limits = c(1400,7500)) +
  scale_fill_discrete(labels = c("0k", "2k", "10k", "50k")) +
  facet_wrap(~Population, scales ="free_x", nrow=2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(panel.background = element_rect(fill = "white"),
          plot.margin = margin(1, 1, 1, 1, "cm"),
          plot.background = element_rect(
            fill = "white",
            size = 1)) +
  theme(axis.text.y = element_text(color = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(legend.text = element_text(color = "black", size = 14)) +
  theme(legend.title = element_text( size = 14)) +
  theme(strip.text = element_text( size=14))

#Minor point in manuscript about Ara−3 being bigger than the rest of the group.
#Coult.quants.means2 %>%  filter(Generation == "50K", Population != "Ara−3") %>%  summarise(mean = mean(Median)) 
#0.766 average excluding Ara−3. Ara−3 has a mean median cell volume of 1.54.

###### Mixed populations: coulter counter data from stationary phase cultures#####

#Ara+4 and REL607 labeled file as 01/09/2020 as this was the replicate I needed these remade for. Data collected on 01/13/2020
#REL607 labeled as 11/21/2019 was collected on 11/22/2019
CoultMP <-  list.files(pattern = "*MP.CSV") 

#load in data from Wiser et al. Science 2013 
FitData <- read.csv("Wiser.etal.Science2013-fit.data.csv")

#skip removes coulter counter summary data from the .csv files 
CoultCol.namesMP <- c( "x", "seconds", "fL", "um", "population", "generation", "block") 
Coulter.dataMP <- as.data.frame(setNames(do.call(rbind,Map('cbind', lapply(CoultMP, read.csv, skip = 274), V4=CoultMP, V5=CoultMP, V6=CoultMP)), CoultCol.namesMP))  %>%  subset(fL > 0 & fL < 6) #see note

#Extrating generation informaiton from the file name
Coulter.dataMP$generation <- as.factor(gsub("_.*","",Coulter.dataMP$generation))

#Extracting population information from the file name and change data order 
Coulter.dataMP$population <- as.factor(gsub(".*K_|_2019.*","",Coulter.dataMP$population))
Coulter.dataMP$population <- as.factor(gsub(".*K_|_2020.*","",Coulter.dataMP$population))
Coulter.dataMP$population <- factor(Coulter.dataMP$population, levels = c("606", "M1", "M2", "M3", "M4", "M5", "M6", "607","P1","P2","P3","P4","P5","P6"))

#Assigning a block number based upon date completed
Coulter.dataMP$block <- as.factor(gsub("^.*2019-11-21.*$","1",Coulter.dataMP$block))
Coulter.dataMP$block <- as.factor(gsub("^.*2019-11-22.*$","2",Coulter.dataMP$block))
Coulter.dataMP$block <- as.factor(gsub("^.*2019-11-26.*$","3",Coulter.dataMP$block))
Coulter.dataMP$block <- as.factor(gsub("^.*2019-11-27.*$","4",Coulter.dataMP$block))
Coulter.dataMP$block <- as.factor(gsub("^.*2020-01-09.*$","5",Coulter.dataMP$block))
Coulter.dataMP$block <- as.factor(gsub("^.*2020-01-10.*$","6",Coulter.dataMP$block))

#Create a new column that has the population and not the block 
Coulter.dataMP$PopGen <- as.factor(paste(Coulter.dataMP$population,Coulter.dataMP$generation, sep="_"))

#Setting factor levels so that data is plotted in the correct order 
Coulter.dataMP$generation <- factor(Coulter.dataMP$generation, c("0K", "2K", "10K", "50K"))

#Creating a column with Freq of counts for each population at each generation.
PopGen.Freq <- as.data.frame(table(Coulter.dataMP$PopGen)) 
names(PopGen.Freq)[1] = "PopGen"
Coulter.dataMP <- left_join(Coulter.dataMP, PopGen.Freq)

#Section where I abstract the medians of the median

#For loop for calculating the median of the medians for each unique population.
#Does what I want but rearranges data
#I have all these replicates of data I want to 
#1) calculate the median for each replicate
#2) Take the median of the medians and return that value 
#3) Linear regression of medians of the segger data against the medians of the coulter data

#1) Calculate the median for each replicate

#Made a new column in the Coulter.dataMP frane that is a concatenation of columns in the dataset
Coulter.dataMP$New <- as.factor(paste(Coulter.dataMP$population,Coulter.dataMP$generation,Coulter.dataMP$block, sep="_"))


#Creating a data frame with summary stats for median Coulter data 
Coul.MediansMP <- as.data.frame(Coulter.dataMP %>% 
                                  group_by(New) %>% 
                                  dplyr::summarize(CoulMedianMP = as.numeric(median(fL, na.rm =T))))
#View(Coul.MediansMP)
#Create new column in Coul.Medians dataframe with name that removes block factor level 
Coul.MediansMP$New2 <- (str_extract(Coul.MediansMP$New, (".*\\d+_\\d+K") ))

#Calculating median of medians for Coulter data 
Coul.Med.MediansMP <- as.data.frame(Coul.MediansMP %>% 
                group_by(New2) %>% 
                 dplyr::summarize(CoulMed.MediansMP = median(CoulMedianMP, na.rm = T)))

#Change how populations are named 
Coul.Med.MediansMP$New3 <- paste("ARA", Coul.Med.MediansMP$New2,sep="_")
Coul.Med.MediansMP$New3 <- ifelse(Coul.Med.MediansMP$New3 == "ARA_606_0K", "REL606_0K", Coul.Med.MediansMP$New3)
Coul.Med.MediansMP$New3 <- ifelse(Coul.Med.MediansMP$New3 == "ARA_607_0K", "REL607_0K", Coul.Med.MediansMP$New3)

#Indexing the dataframe above to contain only the columns of interest and then reordering columns
Coul.Med.MediansMP <- (Coul.Med.MediansMP[,2:3])
Coul.Med.MediansMP <- Coul.Med.MediansMP[,c(2,1)]

colnames(Coul.Med.MediansMP)[1] <- "population"

#####New Section#####
#Here I generate correlation plots for the following. 
#A) cell sizes of mixed populations and the clones... how well does cell size of mixed population predict cell size of clone
#B) cell size of mixed population and MP fitness.. how well does mixed population cell size predict fitness
#C) cell size of clones and MP fitness... how well does clone cell size predict fitness
#D) cell size of clones and fitness of clones using data collected from my aerobic/anaerobic study

#Create dataframe containing columns (3) of measurements for the median of the medians of, 1) mixed populations, 2) clones and
# 3) fitness (Wiser et al. Science 2013)

Test <- as.data.frame(left_join(Coul.Med.MediansMP,Coul.Med.Medians)) 

#changed the name of column from "population" to "New2" in this way I could merge the two dataframes by a common variable
#I will be changing the column name again when comes time to generate actual dataframes 
#Will also change the name of the dataframe from "Test" to something more meaningful.

colnames(FitData)[colnames(FitData)=="population"] <- "New2"

FitData$New3 <- paste("ARA", FitData$New2,sep="_")
FitData$New3 <- ifelse(FitData$New3 == "ARA_606_0K", "REL606_0K", FitData$New3)
FitData$New3 <- ifelse(FitData$New3 == "ARA_607_0K", "REL607_0K", FitData$New3)

#Indexing the dataframe above to contain only the columns of interest and then reordering columns
FitData <- (FitData[,2:4])
FitData <- FitData[,c(3,2,1)]

colnames(FitData)[1] <- "population"

Test <- left_join(Test, FitData)

Test$Generation <-  factor((str_extract(Test$population, "(?<=_)\\d+K$")))

Test$Generation <- factor(Test$Generation, c("0K", "2K", "10K", "50K"))

#####Correlation plot between median cell volumes of clones and mixed populations#####
#A)

#size 6.33x5.03
ggplot(Test, aes(x=CoulMed.MediansMP, y=CoulMed.Medians, colour = factor(Generation))) + 
  geom_smooth(method = lm, se = F, colour = "red", lty = 2) + #geom smooth adds 95% confidence interval 
  geom_point(size=2) +
  scale_y_continuous(breaks = c(0,0.25,.5,.75,1.00,1.25,1.5,1.75) ,limits = c(0.2,1.8)) +
  scale_x_continuous(breaks = c(0,0.25,.5,.75,1.00,1.25,1.5, 1.75), limits = c(0.2,1.8)) +
  scale_color_discrete(name = "Generation", labels = c("0k","2k","10k","50k")) +
  labs(y= "Cell volume of clones (fL)", x= "Cell volume of whole populations (fL)") +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(-.3, 1, 0, 1, "cm"),
        plot.background = element_rect(
          fill = "white",
          size = 1)) +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), axis.line = element_line())

#Correlation test between median of median cell volumes 
ggqqplot(Test$CoulMed.Medians, ylab = "Median of medians cell volume of clones")
ggqqplot(Test$CoulMed.MediansMP, ylab = "Median of medians cell volume of whole populations")

shapiro.test(Test$CoulMed.Medians) #W = 0.82914, p-value = 4.343e-05; not normally distributed; can't use Pearsons
shapiro.test(Test$CoulMed.MediansMP) #W = 0.9736, p-value = 0.4972

#Highly significant difference between the median cell volume of clones and the whole populations 
cor.test(x=Test$CoulMed.MediansMP,y=Test$CoulMed.Medians, method = "kendall", exact = FALSE) #including Ara−3
#tau = 0.4896798, z = 4.3251, p-value = 1.525e-05


Test2 <- Test %>% filter(population != "ARA_M3_50K") %>% 
  droplevels()
cor.test(Test2$CoulMed.MediansMP, y=Test2$CoulMed.Medians, method = "kendall", exact = FALSE) #excluding Ara−3
#tau = 0.479339, z = 4.1725, p-value = 3.012e-05

#####Correlation between mixed population fitness and cell size#####
#size 6.33 x 5.03
ggplot(Test, aes(x=CoulMed.MediansMP, y=fitnessMP, colour = factor(Generation))) + 
  geom_smooth(method = lm, se = F, colour = "red", lty = 2) + #geom smooth adds 95% confidence interval 
  geom_point(size=2) +
  scale_y_continuous(breaks = c(1,1.2,1.4,1.6,1.8,2.0,2.2) ,limits = c(.95,2.2)) +
  scale_x_continuous(breaks = c(0,0.15,.30,.45,.6,.75,.90,1.05,1.20), limits = c(0.2,1.30)) +
  scale_colour_discrete(name = "Generation", labels = c("0k","2k","10k","50k")) +
  labs(x= "Cell volume of whole populations (fL)") +
  ylab(expression(paste("Relative fitness of whole populations" ," (",italic("W"),")"))) +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(-.3, 0, 0, 1, "cm"),
        plot.background = element_rect(
          fill = "white",
          size = 1)) +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), axis.line = element_line())

#Normality assumptions
ggqqplot(Test$fitnessMP, ylab = "Fitness of whole populations")
ggqqplot(Test$CoulMed.MediansMP, ylab = "Median of medians cell volume of whole populations")

shapiro.test(Test$fitnessMP) #W = 0.96067, p-value = 0.2543
shapiro.test(Test$CoulMed.MediansMP) #W = 0.9736, p-value = 0.4972

#Correlation between fitnes and the median cell volume of whole populations 
cor.test(x=Test$CoulMed.MediansMP, y=Test$fitnessMP, method = "kendall", exact = FALSE) 
#tau = 0.6066015, z = 5.0409, p-value = 4.634e-07

#####Correlation plot between clones and fitness of mixed populations#####
#size 6.33 x 5.03
ggplot(Test, aes(x=CoulMed.Medians, y=fitnessMP, colour = factor(Generation))) + 
  geom_smooth(method = lm, se = F, colour = "red", lty = 2) + #geom smooth adds 95% confidence interval 
  geom_point(size=2) +
  scale_colour_discrete(name = "Generation", labels = c("0k","2k","10k","50k")) +
  scale_y_continuous(breaks = c(1,1.2,1.4,1.6,1.8,2.0,2.2) ,limits = c(.95,2.2)) +
  scale_x_continuous(breaks = c(0,.2,.3,.4,.5,.6,.7,.8,.9,1,1.1), limits = c(0.2,1.10)) +
  # geom_label_repel(aes(label = population),
  #                  box.padding   = 0.35, 
  #                  point.padding = 0.5,
  #                  segment.color = 'grey50')+
  labs(x= "Cell volume of clones (fL)") +
  ylab(expression(paste("Relative fitness of whole populations" ," (",italic("W"),")"))) +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(-.3, 1, 0, 1, "cm"),
        plot.background = element_rect(
          fill = "white",
          size = 1)) +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), axis.line = element_line())

#Normality assumptions
ggqqplot(Test$fitnessMP, ylab = "Fitness of whole populations")
ggqqplot(Test$CoulMed.Medians, ylab = "Median of medians cell volume of clones")

shapiro.test(Test$fitnessMP) #W = 0.96067, p-value = 0.2543
shapiro.test(Test$CoulMed.Medians) #W = 0.82914, p-value = 4.343e-05

#Correlation between fitnes of whole populations and the median cell volume of clones 
cor.test(x=Test$CoulMed.Medians, y=Test$fitnessMP, method = "kendall", exact = FALSE) 
#tau= 0.4567352 , z = 3.7955, p-value = 0.0001474

#####Correlation plot between clones and fitness of clones (aero/ana data)#####
#Do not include in the manuscript. Examining for purposes of another manuscript. 
#size = 6.33 x 5.03 
ggplot(Test, aes(x=CoulMed.Medians, y=fitnessClones, colour = factor(Generation))) + 
  geom_smooth(method = lm, se = F, colour = "red", lty = 2) + #geom smooth adds 95% confidence interval 
  scale_color_discrete(name = "Generation", labels = c("0k", "2k", "10k", "50k")) +
  geom_point(size=2) +
  scale_y_continuous(breaks = c(1,1.2,1.4,1.6,1.8,2.0,2.2) ,limits = c(.95,1.8)) +
  scale_x_continuous(breaks = c(0,.2,.3,.4,.5,.6,.7,.8,.9,1,1.1), limits = c(0.2,1.10)) +
  labs(x= "Cell volume of clones (fL)") +
  ylab(expression(paste("Relative fitness of clones" ," (",italic("W"),")"))) +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(-.3, 1, 0, 1, "cm"),
        plot.background = element_rect(
          fill = "white",
          size = 1)) +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), axis.line = element_line())

#Normaility tests 
ggqqplot(Test$CoulMed.Medians, ylab = "Cell volume of clones")
ggqqplot(Test$fitnessClones, ylab = "Relative fitness of clones" )

shapiro.test(Test$CoulMed.Medians) #W = 0.82914, p-value = 4.343e-05; not normally distributed
shapiro.test(Test$fitnessClones) #W = 0.94955, p-value = 0.129

cor.test(x=Test$CoulMed.Medians, y=Test$fitnessClones, method = "kendall", exact = FALSE)
#tau = 0.3658769, z = 2.9908, p-value = 0.002783


#Questions
#A) There is a weak correlation between cell volume of individual clone and mixed population
#B) There is a moderate correlation between the cell volume and fitness when considering mixed populations 
#C) There is a borderline weak correlation (.49) between the cell volume and fitness when considering clones
#D) Negative (exected result). Mike used mixed populations in his Sciecne 2013 paper. 

#####Geom density plots for the mixed population cell size distributions#####
#Change how populations are named
Coulter.dataMP$population <- gsub("M", "Ara−",Coulter.dataMP$population)
Coulter.dataMP$population <- gsub("P", "Ara+",Coulter.dataMP$population)
Coulter.dataMP$population <- gsub("606", "REL606",Coulter.dataMP$population)
Coulter.dataMP$population <- gsub("607", "REL607",Coulter.dataMP$population)

ggplot(Coulter.dataMP, aes(x=fL, color = generation, fill = generation, y= .0125*..count..)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +theme(
          panel.background = element_rect(fill = "white"),
          plot.margin = margin(2, 2, 2, 2, "cm"),
          plot.background = element_rect(
            fill = "white",
            colour = "black",
            size = 1 
          )) +
  geom_density(alpha = 0.5, position = "identity") + 
  xlim(c(.15,2)) +
  scale_y_continuous() + 
  facet_wrap(~population, nrow = 2)  

####### Using Coulter.dataMP files to create summary tables for the quantiles
#I used to construct boxplots. I then take the median of each quantile, for each population and use that
#to construct new boxplots. 

#Coulter data (mixed populations)
Coult.quantsMP <- (as_tibble(tapply(Coulter.dataMP$fL, Coulter.dataMP$New, mod.boxplot)))  #function "mod.boxplot" calculates the quantiles 
Coult.quantsMP <- as.data.frame(t(Coult.quantsMP))
Coult.quantsMP <- rownames_to_column(Coult.quantsMP, var = "rowid")
colnames(Coult.quantsMP) <- c("Population","Fifth","Twenty_fifth","Median","Seventy_fifth","Ninety_fifth")
Coult.quantsMP$Population <- as.factor((str_extract(Coult.quantsMP$Population, (".*\\d+_\\d+K"))))
Coult.quantsMP$Generation <-  (str_extract(Coult.quantsMP$Population, "(?<=_)\\d+K$")) 

#Calculate medians of each quantile for each population sample and then construct a dataframe with these values
Coult.quants.meansMP <- Coult.quantsMP %>% 
  group_by(Population) %>% 
  summarise(Fifth=median(Fifth), Twenty_fifth=median(Twenty_fifth), Median=median(Median), 
            Seventy_fifth=median(Seventy_fifth), Ninety_fifth=median(Ninety_fifth))

#Extract Generation Sampled from 
Coult.quants.meansMP$Generation <-  (str_extract(Coult.quants.meansMP$Population, "(?<=_)\\d+K$")) 

#Extract Population name 
Coult.quants.meansMP$PopGen <- as.factor(gsub("_0K.*|_2K.*|_10K.*|_50K.*", "",Coult.quants.meansMP$Population))

#Setting the order of the data for plotting
Coult.quants.meansMP$Generation <-factor(Coult.quants.meansMP$Generation ,levels = c("0K", "2K","10K","50K"))

#Change column order and rename columns 
Coult.quants.meansMP <- Coult.quants.meansMP[,c(1,8,7,2,3,4,5,6)]
colnames(Coult.quants.meansMP) <- c("PopGen","Population","Generation","Fifth","Twenty_fifth","Median","Seventy_fifth","Ninety_fifth")

#Generate boxplots using the means of the three replicates for each population an at each timepoint
Coult.quants.meansMP$Population <- gsub("M", "Ara−",Coult.quants.meansMP$Population)
Coult.quants.meansMP$Population <- gsub("P", "Ara+",Coult.quants.meansMP$Population)
Coult.quants.meansMP$Population <- gsub("606", "REL606",Coult.quants.meansMP$Population)
Coult.quants.meansMP$Population <- gsub("607", "REL607",Coult.quants.meansMP$Population)

#Made this one so that I can include the ancestor boxplot in each evolved clone panel
#Exp -> I expanded the rows correspondig to the two ancestors 6 times, 1 for each population
#Then I full joined the two data sets
#For figure purposes!!!!!

Exp.CoultMP <- Coult.quants.meansMP[rep(1:nrow(Coult.quants.meansMP[1:2,]),each=6),] 
Coult.quants.means2MP <- full_join(Coult.quants.meansMP, Exp.CoultMP)
Coult.quants.means2MP[1:12,2] <-  paste(c("Ara−1","Ara−2","Ara−3","Ara−4","Ara−5","Ara−6","Ara+1","Ara+2","Ara+3","Ara+4","Ara+5","Ara+6"))

#[-c(1,2),] when using coult quants means I had to include this to exclude the ancestors from having 
#their own panel, with one plot in it. 

Coult.quants.meansMP <- mutate(Coult.quants.meansMP, IQR = Seventy_fifth-Twenty_fifth)
Coult.quants.meansMP <- mutate(Coult.quants.meansMP, IQR_fold_change = (100*IQR/Twenty_fifth))

#Boxplots for coulter counter using means of each quantile to construct them. 
#size9.5 x 6 

ggplot(Coult.quants.means2MP, aes(x=Population, fill = Generation)) +
  geom_boxplot(aes(ymin=Fifth, lower=Twenty_fifth, middle=Median, upper=Seventy_fifth, ymax=Ninety_fifth),
               stat = "identity") + 
  labs(y="Cell volume (fL)") +
  facet_wrap(~Population, scales ="free_x", nrow=2) +
  scale_y_continuous(breaks = seq(0,9.6,.8), limits = c(0,3.3)) +
  scale_fill_discrete(labels = c("0k", "2k", "10k", "50k")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  # axis.title.y = element_text(color="black", size=14, face="bold") + ##This line increases the size y axis title
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.background = element_rect(fill = "white",
          size = 1)) +
  theme(axis.text.y = element_text(color = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(legend.text = element_text(color = "black", size = 14)) +
  theme(legend.title = element_text( size = 14)) +
  theme(strip.text = element_text( size=14))

#Minor points in manuscript Ara+1 and Ara−3
# View(Coult.quants.means2 %>% filter(Generation == "10K", Population == "Ara+1"))
# View(Coult.quants.means2MP %>% filter(Generation == "50K", Population == "Ara+1"))


#####Statistical analyses comparing differences in the grand mean of median cell size differences in MP#####
#Grand mean of median cell volume data using instead of summary 
#In this way I will be able to use the compare means summary thus allowing me to add
#statistical comparisons to the plot. 
#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/

grand.mean.vol.MP <- as.data.frame(full_join(Coult.quantsMP %>% filter(Generation == "0K") %>% 
                                     group_by(Population, Generation) %>% 
                                     transmute(median.median = Median) %>% #here I call the medians "median of median" for the ancestor. Treating the 12 estimates as I did the other populations
                                     droplevels(), 
                                   
                                   Coult.quantsMP %>%  #unlike the clones I include Ara−3 
                                     group_by(Population, Generation) %>% 
                                     filter(Generation != "0K") %>% 
                                     summarise(median.median = median(Median)) %>% 
                                     droplevels()))

#Order factor levels 
grand.mean.vol.MP$Generation <- factor(grand.mean.vol.MP$Generation, c("0K", "2K", "10K", "50K"))

#size 6.33x5.03
grand.mean.vol.MP %>% group_by(Generation) %>% 
  ggplot(aes(y = median.median, x = Generation)) +
  labs(x= "Generation", y= "Cell volume of whole populations (fL)") +
  scale_y_continuous(breaks = seq(.2,1.4,.2)) +
  scale_x_discrete(labels = c("0k", "2k", "10k", "50k")) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "point", position="identity",size=2) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.2) +
  geom_bracket(xmin = "0K", xmax= "2K", label = "<<0.0001", y.position = 0.70)+
  geom_bracket(xmin = "2K", xmax= "10K", label = "0.0002", y.position = 0.95)+
  geom_bracket(xmin = "10K", xmax= "50K", label = "0.0001", y.position = 1.25)+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black",size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#####Statistical analyses comparing differences in the grand mean of median cell volumes of whole populations#####

#Normality test
with(grand.mean.vol.MP, shapiro.test(median.median[Generation == "0K"])) #p-value = 0.3172
with(grand.mean.vol.MP, shapiro.test(median.median[Generation == "2K"])) #p-value = 0.5238
with(grand.mean.vol.MP, shapiro.test(median.median[Generation == "10K"])) #p-value = 0.9623
with(grand.mean.vol.MP, shapiro.test(median.median[Generation == "50K"])) #p-value = 0.6654
#Data is normally distributed. 

#Bartlett's test for the homogenity of variance as a group
bartlett.test(median.median ~ Generation, data = grand.mean.vol.MP) #df = 3, p-value = 9.567e-08; 
#variances not homogenous will adjust for variances with each specific T.test below. 

#3 one-tailed tests: I am concerned with whether preceding generations have greater cell
#volumes than time before.

gm.vol.2KMP <- grand.mean.vol.MP %>% filter(Generation == "0K" | Generation == "2K") %>% 
  droplevels()

gm.vol.10KMP <- grand.mean.vol.MP %>% filter(Generation == "2K" | Generation == "10K") %>% 
  droplevels()

gm.vol.10KMP$Generation <- factor(gm.vol.10KMP$Generation, levels = c("2K", "10K")) 

gm.vol.50KMP <- grand.mean.vol.MP %>% filter(Generation == "10K" | Generation == "50K") %>% 
  droplevels()

gm.vol.0_10KMP <- grand.mean.vol.MP %>% filter(Generation == "0K" | Generation == "10K") %>% 
  droplevels()

gm.vol.0_50KMP <- grand.mean.vol.MP %>% filter(Generation == "0K" | Generation == "50K") %>% 
  droplevels()

#Barlett tests on the variances in each of the comparisons.

bartlett.test(median.median ~ Generation, data = gm.vol.2KMP) #p-value = 0.0001669
bartlett.test(median.median ~ Generation, data = gm.vol.10KMP) #p-value = 0.0784 
bartlett.test(median.median ~ Generation, data = gm.vol.50KMP) #p-value = 0.2992
bartlett.test(median.median ~ Generation, data = gm.vol.0_10KMP) #p-value = 6.163e-07
bartlett.test(median.median ~ Generation, data = gm.vol.0_50KMP) #p-value = 2.089e-08

#3 one-tailed tests: I am concerned with whether preceding generations have greater cell
#volumes than time before.

#####Statisticaly analyses for figure 5#####
#T.tests; all one tailed. 
t.test(median.median ~ Generation, data = gm.vol.2KMP, alternative = "less", var.equal = F, paired = T) #t = -11.612, df = 11, p-value = 8.145e-08
t.test(median.median ~ Generation, data = gm.vol.10KMP, alternative = "less", var.equal = T, paired = T) #t = -4.9766, df = 11, p-value = 0.0002088
t.test(median.median ~ Generation, data = gm.vol.50KMP, alternative = "less", var.equal = T, paired = T) #t = -5.4418, df = 11, p-value = 0.0001017
t.test(median.median ~ Generation, data = gm.vol.0_10KMP, alternative = "less", var.equal = F, paired = T) #t = -10.885, df = 11, p-value = 1.574e-07
t.test(median.median ~ Generation, data = gm.vol.0_50KMP, alternative = "less", var.equal = F, paired = T ) #t = -12.327, df = 11, p-value = 4.412e-08

####Statistical analyses on slope changes using mixed populations####

#Does the change in cell size decellerate?
#calculated the slopes using median of medians for mixed populations

#A: slopes for each individual population.

#Here I format data to calculate slopes.
#using median of medians of mixed populations to make this inference. 

SlopeAnalysis <- Coul.Med.MediansMP

SlopeAnalysis$Generation <- as.factor((str_extract(SlopeAnalysis$population, ("\\d+K"))))
SlopeAnalysis$population <- as.factor(gsub("_0K.*|_2K.*|_10K.*|_50K.*", "",SlopeAnalysis$population))

SlopeAnalysis$Generation <- factor(SlopeAnalysis$Generation,(c("0K","2K","10K","50K"))) #Order the factors here 
levels(SlopeAnalysis$Generation) <- c("0", "2,000", "10,000", "50,000")


#Slope group A and B subtracts grand median cell volumes 2,000 from 0  generations 
#volume 
slopeA <- c("ARA_M1", "ARA_M2", "ARA_M3", "ARA_M4", "ARA_M5", "ARA_M6", "REL606")
slopeB <- c("ARA_P1", "ARA_P2", "ARA_P3", "ARA_P4", "ARA_P5", "ARA_P6", "REL607")

slope.grpA <- SlopeAnalysis %>% filter(Generation == "0" | Generation == "2,000", population %in% slopeA) %>% 
spread(key = "Generation", value = "CoulMed.MediansMP") %>% 
  mutate("Tf-Ti" = as.integer(2000), "Volf-Voli" = (`2,000`- SlopeAnalysis[1,2]), GenComp = "Anc to 2K") %>% 
  select(-c(`0`, `2,000`)) %>%  filter(population != "REL606")

slope.grpB <- SlopeAnalysis %>% filter(Generation == "0" | Generation == "2,000", population %in% slopeB) %>% 
  spread(key = "Generation", value = "CoulMed.MediansMP") %>% 
  mutate("Tf-Ti" = as.integer(2000), "Volf-Voli" = (`2,000`- SlopeAnalysis[2,2]), GenComp = "Anc to 2K") %>% 
  select(-c(`0`, `2,000`)) %>%  filter(population != "REL607")

#Slope group C subtracts grand median cell volumes 2,000 from 10,000 generations
slope.grpC <- SlopeAnalysis %>% filter(Generation == "2,000" | Generation == "10,000") %>% 
  spread(key = "Generation", value = "CoulMed.MediansMP") %>% 
  mutate("Tf-Ti" = as.integer(8000), "Volf-Voli" = (`10,000`-`2,000`), GenComp = "2K to 10K") %>% 
  select(-c(`2,000`, `10,000`))

#Slope group D subtracts grand median cell volumes 2,000 from 10,000 generations
slope.grpD <- SlopeAnalysis %>% filter(Generation == "10,000" | Generation == "50,000") %>% 
  spread(key = "Generation", value = "CoulMed.MediansMP") %>% 
  mutate("Tf-Ti" = as.integer(40000), "Volf-Voli" = (`50,000`-`10,000`), GenComp = "10K to 50K") %>% 
  select(-c(`10,000`, `50,000`))

SlopeDataMP <- as.data.frame(full_join(slope.grpA, slope.grpB) %>% 
  full_join(., slope.grpC) %>% 
  full_join(., slope.grpD) %>% mutate(slope = (`Volf-Voli`/`Tf-Ti`))) 

levels(SlopeDataMP$population) <-c("Ara−1", "Ara−2", "Ara−3", "Ara−4", "Ara−5", "Ara−6","Ara+1","Ara+2","Ara+3","Ara+4","Ara+5","Ara+6", "REL606", "REL607")
SlopeDataMP$GenComp <- factor(SlopeDataMP$GenComp, c("Anc to 2K", "2K to 10K", "10K to 50K"))

SlopeDataMP$slope1000 <- SlopeDataMP$slope * 1000

ggplot(SlopeDataMP, aes(y= slope, x= GenComp, colour = population, group = population)) +
  geom_point() +
  geom_path() +
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2, "cm"),
        plot.background = element_rect(
          fill = "white",
          colour = "black",
          size = 1)) 

#B
#Average the slopes across all populations at each comparison 

#Summary table of slope changes 
avgSlopeDataMP <- SlopeDataMP %>% 
  group_by("Tf-Ti", GenComp) %>% 
  summarise(avg.slope=mean(slope*1000), std.dev.slope= sd(slope*1000), n=n()) %>% 
  mutate(CIlow = avg.slope - abs(qt(0.05/2, n-1)*(std.dev.slope/sqrt(n)))) %>% 
  mutate(CIHigh = avg.slope + abs(qt(0.05/2, n-1)*(std.dev.slope/sqrt(n))))

slope.comparisons <- list(c("Anc to 2K", "2K to 10K"), c("2K to 10K", "10K to 50K"))

slope.change.1 <- SlopeDataMP %>% filter(GenComp == "Anc to 2K" | GenComp == "2K to 10K") %>% 
  droplevels()

slope.change.2 <- SlopeDataMP %>% filter(GenComp == "2K to 10K" | GenComp == "10K to 50K") %>% 
  droplevels()

#Normality test; data is normally distributed 
with(SlopeDataMP, shapiro.test(slope1000[`Tf-Ti` == "2000"])) #p-value = 0.5359
with(SlopeDataMP, shapiro.test(slope1000[`Tf-Ti` == "8000"])) #p-value = 0.2217
with(SlopeDataMP, shapiro.test(slope1000[`Tf-Ti` == "40000"])) #p-value = 0.7573

#Averages
with(SlopeDataMP, mean(slope1000[`Tf-Ti` == "2000"])) #0.1699583
with(SlopeDataMP, mean(slope1000[`Tf-Ti` == "8000"])) #0.02341667
with(SlopeDataMP, mean(slope1000[`Tf-Ti` == "40000"])) #0.006877083

#Bartlett's test homogeniety of variances; variances are NOT homogenous. 

bartlett.test(slope1000 ~ GenComp, data = SlopeDataMP) #p-value = 2.499e-09

#####Wilcoxon tests on slope change differrences#####
wilcox_test(slope1000 ~ GenComp, data = slope.change.1, alternative = "greater") #p= 0.0000182
wilcox_test(slope1000 ~ GenComp, data = slope.change.2, alternative = "greater") #p=0.000824

#size = 4 x 6.45 
SlopeDataMP %>% group_by(`Tf-Ti`) %>% 
  ggplot(aes(y = slope1000, x = GenComp)) +
  labs(x= "Generations compared", y= "Average rate of change in cell volume (fL per 1,000 generations)") +
  scale_y_continuous(breaks = seq(0.00,0.20,.05)) +
  scale_x_discrete(labels = c("Anc to 2k", "2k to 10k", "10k to 50k")) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "point", position="identity",size=2) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.2) +
  geom_bracket(xmin = "Anc to 2K", xmax= "2K to 10K", label = "<<0.0001", y.position = 0.21)+
  geom_bracket(xmin = "2K to 10K", xmax= "10K to 50K", label = "0.0008", y.position = 0.05)+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black",size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#####Mean length-to-width changes among the LTEE clones####
#Plotting data; Lollipop chart 7.50 X 5.00 plot size looks beautiful
#Analysis of length to width ratios
#Create dataframe the calculates the mean Lengths and width 
#Making tables
#Running T.tests on the data

Trail1 <- Segger.data %>% filter(generation == "0K")  %>%
  mutate(l.w. = LongAxis/ShortAxis) %>% 
  group_by(population, generation, Folder) %>%
  summarise(lw.ratio=mean(l.w.),n = n()) %>%
  droplevels()

Trail2 <- Segger.data %>%  
  mutate(LongAxis = LongAxis*0.0664, ShortAxis = ShortAxis*0.0664) %>% 
  mutate(l.w. = LongAxis/ShortAxis) %>%
  group_by(population, generation, Folder) %>%
  summarise(mean.lw = mean(l.w.), mean.l =mean(LongAxis), mean.s =mean(ShortAxis), n = n()) %>% #average length/width ratio of n cells per block 
  group_by(population,generation) %>%
  summarise (lw.ratio = mean(mean.lw),mean.length = mean(mean.l), mean.width = mean(mean.s), n=n()) #average length/width ratio of the three blocks

Trail2$Num.Generation <- as.numeric(Trail2$generation)

#Here I create a dummy column, "Ratio.adj" where I change the value of Ara-3 to "5.0000000" 
#Doing it this way allows me to put a "break" in the x-axis
#The true value is retained in the lw.ratio column.

Trail2 <- Trail2[order(Trail2$generation),]
Trail2$Ratio.adj <- Trail2$lw.ratio
Trail2$length.adj <- Trail2$mean.length
Trail2$PopGen <- as.factor(paste(Trail2$population,Trail2$generation,sep="_"))
Trail2[Trail2$PopGen == "Ara−3_50K", "Ratio.adj"] <- 5.0000000
Trail2[Trail2$PopGen == "Ara−3_50K", "length.adj"] <- 4.5000000

#####plot the data  mean length####
#size 5.90 x 7.00
levels(Trail2$generation) <- c("0", "2,000", "10,000", "50,000")
ggplot(Trail2, aes(x=PopGen, y=length.adj, label =mean.length, color = generation)) +
  geom_point(aes(fct_reorder(PopGen, Num.Generation)), fill="black", size=4.50) +
  scale_x_discrete(labels=(Trail2$population)) +
  scale_color_discrete(labels = c("0k","2k","10k","50k")) +
  geom_segment(aes(y = 2.07, 
                   x = `PopGen`, 
                   yend = length.adj, 
                   xend = `PopGen`)) +
  labs(x ="Population")  +
  ylab(expression(paste("Average cell length" ," (",mu,"m) ")))+
  coord_flip() +
  theme_classic()+
  labs(color = "Generation") +
  theme(axis.text.y = element_text(color = "black", size = 14)) +
  theme(axis.text.x = element_text(color = "black", size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text( size = 14)) +
  scale_y_continuous(breaks = c(1,1.5,2,2.5,3,3.5,4,4.5), labels = c("1","1.5","2","2.5","3","3.5","4","6.5"))

#####plot the data  mean width####
  #size 5.90 x 7.00
  levels(Trail2$generation) <- c("0", "2,000", "10,000", "50,000")
  ggplot(Trail2, aes(x=PopGen, y=mean.width, label = mean.width, color = generation)) +
    geom_point(aes(fct_reorder(PopGen, Num.Generation)), fill="black", size=4.50) +
    scale_x_discrete(labels=(Trail2$population)) +
    scale_color_discrete(labels = c("0k","2k","10k","50k")) +
    geom_segment(aes(y = 0.612, 
                     x = `PopGen`, 
                     yend = mean.width, 
                     xend = `PopGen`)) +
    labs(x ="Population")  +
    ylab(expression(paste("Average cell width" ," (",mu,"m) ")))+
    coord_flip() +
    theme_classic()+
    labs(color = "Generation") +
    theme(axis.text.y = element_text(color = "black", size = 14)) +
    theme(axis.text.x = element_text(color = "black", size = 14)) +
    theme(axis.title.x = element_text(size = 14)) +
    theme(axis.title.y = element_text( size = 14)) +
    theme(legend.text = element_text( size = 14)) +
    theme(legend.title = element_text( size = 14)) +
    scale_y_continuous(breaks = c(0.6, 0.7, 0.8, 0.9))
  
#####plot the data  mean length/width####
#size 5.90 x 7.00
levels(Trail2$generation) <- c("0", "2,000", "10,000", "50,000")
ggplot(Trail2, aes(x=PopGen, y=Ratio.adj, label = lw.ratio, color = generation)) +
  geom_point(aes(fct_reorder(PopGen, Num.Generation)), fill="black", size=4.50) +
  scale_x_discrete(labels=(Trail2$population)) +
  scale_color_discrete(labels = c("0k","2k","10k","50k")) +
  geom_segment(aes(y = 3.34, 
                   x = `PopGen`, 
                   yend = Ratio.adj, 
                   xend = `PopGen`)) +
  labs(x ="Population", y= "Mean (length/width)")  +
  coord_flip() +
  theme_classic()+
  labs(color = "Generation") +
  #theme(legend.position = "none") +
  theme(axis.text.y = element_text(color = "black", size = 14)) +
  theme(axis.text.x = element_text(color = "black", size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text( size = 14)) +
  scale_y_continuous(breaks = c(2,2.5,3,3.5,4,4.5,5), labels = c("2","2.5","3","3.5","4","4.5","8"))

#Calculate avg l/w ratios for each time point to generate table with data for stats

Trail3 <- Trail2 %>% filter(population != "Ara−3" | generation != "50,000") %>%  
  filter(generation != "0") %>% 
  group_by(generation) %>% 
  summarise(grand.mean.lw.ratio = mean(lw.ratio), grand.sd.lw.ratio = sd(lw.ratio), n=n()) %>% 
  mutate(CIlow = grand.mean.lw.ratio - abs(qt(0.05/2, n-1)*(grand.sd.lw.ratio/sqrt(n)))) %>% 
  mutate(CIHigh = grand.mean.lw.ratio + abs(qt(0.05/2, n-1)*(grand.sd.lw.ratio/sqrt(n)))) %>% 
  droplevels()

Trail3.1 <- Trail1 %>% filter(generation  == "0K") %>% group_by(generation) %>% 
  summarise(grand.mean.lw.ratio = mean(lw.ratio), grand.sd.lw.ratio = sd(lw.ratio), n=n()) %>% 
  mutate(CIlow = grand.mean.lw.ratio - abs(qt(0.05/2, n-1)*(grand.sd.lw.ratio/sqrt(n)))) %>% 
  mutate(CIHigh = grand.mean.lw.ratio + abs(qt(0.05/2, n-1)*(grand.sd.lw.ratio/sqrt(n)))) %>% 
  droplevels()

levels(Trail3.1$generation)[levels(Trail3.1$generation)== "0K"] <- "0"
Trail3 <- full_join(Trail3.1, Trail3)
Trail3$generation <- factor(Trail3$generation)


#Setting order levels for plotting
Trail3$generation <-factor(Trail3$generation ,levels = c("0", "2,000","10,000","50,000"))

#Trail 3 reports the grand mean and standard deviation for length-to-width changes
#and their associated 95% confidence intervals.I talk about this in the mansucript, 
#so a table might not be needed. Also, the figure shows this data as well. 

#####Plot mean length to width changes by the populations#####

levels(Trail3$generation) <- c("0k", "2k", "10k", "50k")
ggplot(Trail3, aes(y= grand.mean.lw.ratio, x= generation)) +
  geom_point(size = 2) +
  labs(x= "Generation", y= "Average cell aspect ratio (length/width)") +
  geom_errorbar(aes(ymin=CIlow, ymax=CIHigh), colour="black", size=.3, width=.3, position = "dodge") +
  expand_limits(y=c(2.5,4.0)) +
  geom_bracket(xmin = "0k", xmax= "2k", label = "0.0390", y.position = 3.65) + #See statistical analyses beginning on line 1361 for p-values ("label") 
  geom_bracket(xmin = "2k", xmax= "10k", label = "0.9020", y.position = 3.45) +
  geom_bracket(xmin = "10k", xmax= "50k", label = "0.0245", y.position = 3.85) +
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black",  size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#####Statistical analyses on the mean length to width changes experieneced by the populations#####

#Normality test
with(Trail1, shapiro.test(lw.ratio[generation == "0K"])) #p-value = 0.3491 
with(Trail2, shapiro.test(lw.ratio[generation == "2,000"])) #p-value = 0.5370
with(Trail2, shapiro.test(lw.ratio[generation == "10,000"])) #p-value = 0.6952
with(Trail2, shapiro.test(lw.ratio[generation == "50,000"])) #p-value = 0.001143; Ara−3 is included

Trail2_excl_M3 <- Trail2 %>%  filter(PopGen != "Ara−3_50K")
with(Trail2_excl_M3, shapiro.test(lw.ratio[generation == "50,000"])) #p-value = 0.7093; Ara−3 is excluded

#The results of the Shapiro-Wilks normality test indicates that the 
#mean of the median cell volumes are normally distributed. 

#filter dataframe to abstract the 12 ancestor measurements 
Trail_0K <- as.data.frame(Trail1 %>% filter(generation == "0K") %>% 
                            transmute(lw.ratio = lw.ratio, n = 1, PopGen = factor("Ancestor")) %>% 
                            droplevels())
levels(Trail_0K$generation)[levels(Trail_0K$generation)== "0K"] <- "0" 

Trail_2K <- Trail2 %>%  filter(generation == "2,000") %>% 
  select(population, generation, lw.ratio, n, PopGen) %>% 
  droplevels()

Trail_10K <- Trail2 %>%  filter(generation == "10,000") %>% 
  select(population, generation, lw.ratio, n, PopGen) %>% 
  droplevels()

#Excludes Ara−3 from the dataset.
Trail_50K <- Trail2 %>%  filter(generation == "50,000" , population != "Ara−3") %>% 
  select(population, generation, lw.ratio, n, PopGen) %>% 
  droplevels()

#Make data frames on which I will be doing my t.test comparisons. 
#Ancestor to each generation
T0K_2K <- full_join(Trail_0K,Trail_2K)
T0K_10K <- full_join(Trail_0K,Trail_10K)
T0K_50K <- full_join(Trail_0K,Trail_50K)

#Between generations
T2K_10K <- full_join(Trail_2K,Trail_10K)
T10K_50K <- full_join(Trail_10K,Trail_50K)

#Use Bartlett's test to test homogeneity of variances. This test is valid becuase the 
#Shapiro-Wilks normaility test shows that my data is normally distributed. 

#####Bartlett's test#####
bartlett.test(lw.ratio ~ generation, data = T0K_2K) #p-value = 0.01186  use unequal variance call in t-test
bartlett.test(lw.ratio ~ generation, data = T0K_10K) #p-value = 0.2051
bartlett.test(lw.ratio ~ generation, data = T0K_50K) #p-value = 0.03651 use unequal variance call in t-test 
bartlett.test(lw.ratio ~ generation, data = T2K_10K) #p-value = 0.1798
bartlett.test(lw.ratio ~ generation, data = T10K_50K) #p-value = 0.3722

#H0: There is no difference in the mean cell size difference of populations between the two time points.  
#HA: cell size at later generation in comparison is greater than earlier generation in comparison 
#use alternative hypothesis test = "less,"
#raw p = .05/5 = 0.0100

######Stats for figure 12 manuscript#####
t.test(lw.ratio ~ generation, data = T0K_2K, alternative = "two.sided", var.equal = F) #t = 2.2597, df = 15.15, p-value = 0.03899
t.test(lw.ratio ~ generation, data = T0K_10K, alternative = "two.sided", var.equal = T) #t = 3.3062, df = 22, p-value = 0.003215
t.test(lw.ratio ~ generation, data = T0K_50K, alternative = "two.sided", var.equal = F) #t = -0.067158, df = 14.54, p-value = 0.9474
t.test(lw.ratio ~ generation, data = T2K_10K, alternative = "two.sided", var.equal = T) #t = -0.12452, df = 22, p-value = 0.902
t.test(lw.ratio ~ generation, data = T10K_50K, alternative = "two.sided", var.equal = T) #t = -2.4229, df = 21, p-value = 0.02452

#####Surface-to-volume ratio of evolved lines#####

######Calculate surface to volume ratio using Ojkic equation#####
SAratio <- Segger.data %>% group_by(population, generation, Folder) %>%
  filter(MeshVolume != "NaN") %>% 
  mutate(lw.ratio = LongAxis/ShortAxis, gamma = (lw.ratio*pi) * (lw.ratio*pi/4 - pi/12)^-.6666, S = gamma*(MeshVolume^.666), sav = S/MeshVolume)

Ojkic <- function (x) {
   return(2*pi*(x^.66))
}

SAratio %>%  group_by(population) %>%
  filter(population == "Ara−5") %>%
  ggplot(aes(x = MeshVolume, y = S)) +
  geom_point() +
  facet_wrap(~ generation + population)+
  labs(y = "Surface area (a.u.)", x = "Microscopy volume (a.u.)")+
  stat_function(fun = Ojkic, color = "red", lty = 2) +
  geom_smooth(method= lm, colour = "blue", lty = 2) +
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black",  size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#####surface to volume change for individual clones####
#Ancestor
SAratio_anc <- SAratio %>% filter(generation  == "0K") %>% group_by(New, generation, Folder) %>% 
  summarise(mean.SA.V = mean(sav), n=n()) %>% group_by(generation)
SAratio_anc$New <- as.factor((str_extract(SAratio_anc$New, ("REL\\d+"))))

#Evolved
SAratio.evolved <- SAratio %>% 
  group_by(population, generation, Folder) %>%
  summarise(mean.SAV = mean(sav), n = n(), sd.SAV = sd(sav))

SAratio.evolved.anc <- as.data.frame(SAratio.evolved %>%  filter( generation == "0K") %>% 
  transmute(mean.SAV.evolved = mean.SAV, sd.SAV = sd.SAV, n = n()))

SAratio.evolved.plot <- as.data.frame(SAratio.evolved %>% filter(generation != "0K") %>% group_by(population, generation) %>%  
summarise(mean.SAV.evolved = mean(mean.SAV), sd.SAV = sd(mean.SAV), n=n())) #%>%
#   group_by(generation) %>% 
#   summarise(grand.meanSAV = mean(mean.SAV.evolved), grand.sdSAV = sd(mean.SAV.evolved), n = n()) %>% 
#   mutate(CIlow = grand.meanSAV- abs(qt(0.05/2, n-1)*(grand.sdSAV/sqrt(n)))) %>%
#   mutate(CIHigh = grand.meanSAV + abs(qt(0.05/2, n-1)*(grand.sdSAV/sqrt(n))))

SAratio.evolved.plot <- full_join(SAratio.evolved.anc, SAratio.evolved.plot)

levels(SAratio.evolved.plot$generation) <- c("0k", "2k", "10k", "50k")
SAratio.evolved.plot$PopGen <- as.factor(paste(SAratio.evolved.plot$population,SAratio.evolved.plot$generation, sep="_"))

#Assign an id to each unique clone
#Assigns unique ID based on two columns, in my case population and generation 
SAratio.evolved.plot$popID <- cumsum(!duplicated(SAratio.evolved.plot[1:2]))

#Create column called generation ID
SAratio.evolved.plot$genID <- as.numeric(SAratio.evolved.plot$generation)

#Compute arithmetic averages for the horizontal lines in the plot 
with(SAratio.evolved, mean(mean.SAV[generation == "0K"])) #0.4608907 
with(SAratio.evolved, mean(mean.SAV[generation == "2K"])) #0.4298939
with(SAratio.evolved, mean(mean.SAV[generation == "10K"])) #0.4115077
with(SAratio.evolved, mean(mean.SAV[generation == "50K"])) #0.3902947

####SA/V plot for clones, grouped by generation#####
SAratio.evolved$PopGen <- as.factor(paste(SAratio.evolved$population,SAratio.evolved$generation, sep="_"))

#Assign an id to each unique clone
#Assigns unique ID based on two columns, in my case population and generation 
SAratio.evolved$popID <- cumsum(!duplicated(SAratio.evolved[1:2]))

#Create column called generation ID
SAratio.evolved$genID <- as.numeric(SAratio.evolved$generation)

SAratio.evolved %>%
  ggplot(aes(reorder(PopGen, genID), mean.SAV), color = generation) +
  geom_hline(yintercept = 0.4608907, lty = 2, colour = 2) +
  geom_hline(yintercept = 0.4298939, lty = 2, colour = "3") +
  geom_hline(yintercept = 0.4115077, lty = 2, colour = "cyan3") +
  geom_hline(yintercept = 0.3902947, lty = 2, colour = "mediumorchid3") +
  stat_summary(fun.data = mean_cl_normal, geom = "linerange", size=.85, aes(colour = factor(generation))) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "point", position="identity",size=2.5, aes(color = paste(generation)))+
  labs(y = "Surface-to-volume ratio", x = "Clone") +
  scale_color_discrete(name = "Generation") +
  scale_y_continuous(breaks = c(0,0.30,0.35,0.40,0.45,0.50,0.55) ,limits = c(.30,.55))+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black",  size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(color = "black",angle = 90, size = 14, margin = (margin(l = 5, r=5)), vjust = .5))

#####Grand mean surface-to-volume-change#####
SAratio.evolved.plot %>%  
  group_by(generation) %>% 
  ggplot(aes(x= generation, y = mean.SAV.evolved)) +
  geom_hline(yintercept = 0.375, lty = 2, colour = "white") +
  labs(x= "Generation" , y= "Mean SA/V (a.u.)") +
  scale_y_continuous(breaks = seq(.360,.480,.015)) +
  scale_x_discrete(labels = c("0k", "2k", "10k", "50k")) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "point", position="identity",size=2) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.2, position = "identity") +
  geom_bracket(xmin = "0k", xmax= "2k", label = "0.0001", y.position = 0.475)+
  geom_bracket(xmin = "2k", xmax= "10k", label = "0.0058", y.position = 0.45)+
  geom_bracket(xmin = "10k", xmax= "50k", label = "0.0021", y.position = 0.435)+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black",size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#####Statistical analyses on the surface-to-volume changes experienced by the populations#####

#####Normality tests for surface-to-volume changes#####

#Must use 12 ancestral mean estimates for 0k
with(SAratio.evolved.plot, shapiro.test(mean.SAV.evolved[generation == "0k"])) #p-value = 0.3774

#Other generational time-points
SA.ratio.clone.gen <- SAratio.evolved.plot %>% filter(generation != "0k") %>% 
  droplevels()

with(SA.ratio.clone.gen, shapiro.test(mean.SAV.evolved[generation == "2k"])) #p-value = 0.1101
with(SA.ratio.clone.gen, shapiro.test(mean.SAV.evolved[generation == "10k"])) #p-value = 0.2157
with(SA.ratio.clone.gen, shapiro.test(mean.SAV.evolved[generation == "50k"])) #p-value = p-value = 0.03771; Ara−3 is included

#filter dataframe to abstract the 12 ancestor measurements 
SVR_0k <- as.data.frame(SAratio.evolved.plot %>% filter(generation == "0k") %>% 
                            droplevels())

SVR_2k <- SA.ratio.clone.gen %>%  filter(generation == "2k") %>% 
  droplevels()

SVR_10k <- SA.ratio.clone.gen %>%  filter(generation == "10k") %>% 
  droplevels()

SVR_50k <- SA.ratio.clone.gen %>%  filter(generation == "50k") %>% 
  droplevels()

#Make data frames on which I will be doing my t.test comparisons. 
#Ancestor to each generation
SVR_0k.2k <- full_join(SVR_0k, SVR_2k)

SVR_2k.10k <- full_join(SVR_2k, SVR_10k)
SVR_2k.10k$generation <- factor(SVR_2k.10k$generation, levels = c("2k", "10k")) 

SVR_10k.50k <- full_join(SVR_10k, SVR_50k)

#Between generations
SVR_0k.10k <- full_join(SVR_0k, SVR_10k)
SVR_0k.50k <- full_join(SVR_0k, SVR_50k)

#Use Bartlett's test to test homogeneity of variances.
#####Bartlett's test#####
bartlett.test(mean.SAV.evolved ~ generation, data = SVR_0k.2k) #p-value = 0.07401
bartlett.test(mean.SAV.evolved ~ generation, data = SVR_2k.10k) #p-value = 0.1851
bartlett.test(mean.SAV.evolved ~ generation, data = SVR_10k.50k) #p-value = 0.8754
bartlett.test(mean.SAV.evolved ~ generation, data = SVR_0k.10k) #p-value = 0.003023; use Welch's t-test
bartlett.test(mean.SAV.evolved ~ generation, data = SVR_0k.50k) #p-value = 0.001978; use Welch's t-test

#####t-tests comparing differences in SA/V####
t.test(mean.SAV.evolved ~ generation, data = SVR_0k.2k, paired = T, alternative = "greater") #t = 5.205, df = 11, p-value = 0.0001461
t.test(mean.SAV.evolved ~ generation, data = SVR_2k.10k, paired = T, alternative = "greater") #t = 3.0202, df = 11, p-value = 0.005826
t.test(mean.SAV.evolved ~ generation, data = SVR_10k.50k, paired = T, alternative = "greater") #t = 3.587, df = 11, p-value = 0.002133

####Lollipopchart for SA/V ratios####
SAVlolli <- SAratio.evolved.plot %>%  group_by(population, generation,PopGen) %>% 
summarise(mean.sav = mean(mean.SAV.evolved), genID = mean(genID))

SAVlolli<- SAVlolli[order(SAVlolli$generation),]

#5.9 x 7.00
ggplot(SAVlolli, aes(x=PopGen, y=mean.sav, label = mean.sav, color = generation)) +
  geom_point(aes(fct_reorder(PopGen, genID)), fill="black", size=4.50) +
  scale_x_discrete(labels=(SAVlolli$population)) +
  scale_color_discrete(labels = c("0k","2k","10k","50k")) +
  geom_segment(aes(y = 0.461, 
                   x = `PopGen`, 
                   yend = mean.sav, 
                   xend = `PopGen`)) +
  labs(x ="Population", y= "Mean SA/V (a.u.)")  +
  coord_flip() +
  theme_classic()+
  labs(color = "Generation") +
  #theme(legend.position = "none") +
  theme(axis.text.y = element_text(color = "black", size = 14)) +
  theme(axis.text.x = element_text(color = "black", size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text( size = 14)) +
  #theme(legend.justification=c(-3.0,-1), legend.position=c(0,0)) +
  scale_y_continuous(breaks = c(.30,.32,.34,.36,.38,.40,.42,.44,.46,.48))

#####Isometric analysis###### 

#gamma takes into account aspect ratio
# isometry <- SAratio %>% group_by(population, generation) %>% 
#   summarise(mean.volume = mean(MeshVolume), mean.gamma = mean(gamma))

isometry10k <- SAratio %>% filter(generation == "10K") %>% 
  group_by(population, generation) %>% 
  summarise(mean.volume = mean(MeshVolume), mean.gamma = mean(gamma)) %>% 
  transmute(mean.gamma10k = mean.gamma) %>% 
  droplevels()
  
isometry50k <- SAratio %>% filter(generation == "50K") %>% 
    group_by(population, generation) %>% 
    summarise(mean.volume50k = mean(MeshVolume), mean.gamma50k = mean(gamma)) %>% 
    mutate(S_50k = mean.gamma50k*(mean.volume50k^.66), sav50k = S_50k/mean.volume50k) %>% 
    select(-generation) %>% 
    droplevels()

isometry10k_50k <- full_join(isometry10k,isometry50k) %>%  
  mutate(S_10k = mean.gamma10k*(mean.volume50k^.66), savIso = S_10k/mean.volume50k, "50g10" = savIso < sav50k) 
## savIso reflects the SA/V ratio of a cell with the surface area at 10k generation divided by 
## the volume at 50k generations. This is the hypothetical surface area to volume.
## sav50k is the true sa/v ratio at 50k generations. 

isometry10k_50k_long <- gather(isometry10k_50k, comparison, ratio, savIso, sav50k)
isometry10k_50k_long$comparison <- factor(isometry10k_50k_long$comparison, levels = c("savIso", "sav50k")) 

isometry10k_50k_long  %>%  
  group_by(comparison) %>% 
  ggplot(aes(x= comparison, y = ratio)) +
  labs(x= "" , y= "Mean SA/V (a.u.)") +
  scale_y_continuous(breaks = seq(.300,.410,.010)) +
  scale_x_discrete(labels = c("No shape change", "With shape change")) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "point", position="identity",size=2) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.2, position = "identity") +
  geom_bracket(xmin = "savIso", xmax= "sav50k", label = "0.0144", y.position = 0.40)+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black",size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) 

#Normality tests
with(isometry10k_50k_long, shapiro.test(ratio[comparison == "savIso"])) #p-value = 0.4058
with(isometry10k_50k_long, shapiro.test(ratio[comparison == "sav50k"])) #p-value = 0.08442

#Bartlett's test to test homogeneity of variances.
bartlett.test(ratio ~ comparison, data = isometry10k_50k_long ) #p-value = 0.04739

#t.test 
#Is the mean sa/v ratio lower with no shape change than it is with the shape change?
t.test(ratio ~ comparison, data = isometry10k_50k_long, paired = T, alternative = "l") #t = -2.5132, df = 11, p-value = 0.01441

