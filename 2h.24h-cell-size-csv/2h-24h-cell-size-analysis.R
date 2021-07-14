rm(list=ls())
library(tidyverse)
library(cowplot)
library(dplyr)
library(gridExtra)
library(reshape2)
library(ggrepel)
library(grid)
library(ggpubr)
library(car)

setwd("/Users/nkrumahgrant/Desktop/Cell Morphology/Manuscript-Cell-size-full-submission/2h.24h-cell-size-csv/")

#Data for Ar+1 was not collected on the 1/24/2020 and 1/25/2020. To compensate, we ran two replicates of Ara+1 on 2/5/2020 and 2/6/2020.
#To make data analysis easier, I changed the date on the folder for one of the replicates collected on 2/5/2020 and 2/6/2020 to the 1/24/2020
#and 1/25/2020. The metadata will maintain the actual date collected. 
CoultClones <-  list.files(pattern = "*.CSV")

CoultCol.namesClones <- c( "x", "seconds", "fL", "um", "population", "generation", "block", "time")

#Here I import the data placing lower and upper bouds on the data at the onset
Coulter.dataClones <- as.data.frame(setNames(do.call(rbind,Map('cbind', lapply(CoultClones, read.csv, skip = 274), V4= CoultClones, V5=CoultClones, V6=CoultClones, V7=CoultClones)), CoultCol.namesClones))  %>%  
  subset(fL > 0.200 & fL < 6.00)

#Extrating generation informaiton from the file name
Coulter.dataClones$generation <- as.factor(gsub("_.*","",Coulter.dataClones$generation))

#Extracting population information from the file name and change data order 
Coulter.dataClones$population <- as.factor(gsub(".*K_|_2020.*","",Coulter.dataClones$population))
Coulter.dataClones$population <- factor(Coulter.dataClones$population, levels = c("606", "M1", "M2", "M3", "M4", "M5", "M6", "607","P1","P2","P3","P4","P5","P6"))

#Change how populations are named
Coulter.dataClones$population <- gsub("M", "Ara-",Coulter.dataClones$population)
Coulter.dataClones$population <- gsub("P", "Ara+",Coulter.dataClones$population)
Coulter.dataClones$population <- gsub("606", "REL606",Coulter.dataClones$population)
Coulter.dataClones$population <- gsub("607", "REL607",Coulter.dataClones$population)

#Assigning a block number based upon date completed
Coulter.dataClones$block <- as.factor(gsub("^.*2020-01-24.*$","1",Coulter.dataClones$block)) #2h block 1 
Coulter.dataClones$block <- as.factor(gsub("^.*2020-01-25.*$","1",Coulter.dataClones$block)) #24h block 1
Coulter.dataClones$block <- as.factor(gsub("^.*2020-02-05.*$","2",Coulter.dataClones$block)) #2h block 2
Coulter.dataClones$block <- as.factor(gsub("^.*2020-02-06.*$","2",Coulter.dataClones$block)) #24h block 2
Coulter.dataClones$block <- as.factor(gsub("^.*2020-02-14.*$","3",Coulter.dataClones$block)) #2h block 3
Coulter.dataClones$block <- as.factor(gsub("^.*2020-02-15.*$","3",Coulter.dataClones$block)) #24h block 3

#Assigning a time based upon date completed
Coulter.dataClones$time <- as.factor(gsub("^.*2020-01-24.*$","Mid-exponential",Coulter.dataClones$time)) #2h block 1 
Coulter.dataClones$time <- as.factor(gsub("^.*2020-01-25.*$","Stationary phase",Coulter.dataClones$time)) #24h block 1
Coulter.dataClones$time <- as.factor(gsub("^.*2020-02-05.*$","Mid-exponential",Coulter.dataClones$time)) #2h block 2
Coulter.dataClones$time <- as.factor(gsub("^.*2020-02-06.*$","Stationary phase",Coulter.dataClones$time)) #24h block 2
Coulter.dataClones$time <- as.factor(gsub("^.*2020-02-14.*$","Mid-exponential",Coulter.dataClones$time)) #2h block 3
Coulter.dataClones$time <- as.factor(gsub("^.*2020-02-15.*$","Stationary phase",Coulter.dataClones$time)) #24h block 3

#What are the fifth and 95th percentiles of the data between the initial bounds of input data?

#5th?
Coulter.dataClones <- Coulter.dataClones %>%  group_by(population, generation, block, time) %>% 
  mutate(fifth = quantile(fL, probs = (.05)))

min(Coulter.dataClones$fifth) #.207

#95th?
Coulter.dataClones <- Coulter.dataClones %>%  group_by(population, generation, block, time) %>% 
  mutate(ninety_fifth = quantile(fL, probs = (.95)))

max(Coulter.dataClones$ninety_fifth) #5.6594

#Now filter dataset keeping only the middle 90% of cell size data  (0.21 <= fL <= 5.66)
Coulter.dataClones <- as.data.frame(Coulter.dataClones %>%  group_by(population, generation, block, time) %>%
  subset(fL > 0.207 & fL < 5.6594))

#Plots the histograms for 2h
Coulter.dataClones %>% filter(time=="Mid-exponential") %>% 
  ggplot(aes(fL, fill = block)) +
  geom_histogram(bins = 100) +
  facet_wrap(~population, nrow = 2)

#Plots the histograms for 24h
Coulter.dataClones %>% filter(time=="Stationary phase") %>% 
  ggplot(aes(fL, fill = block)) +
  geom_histogram(bins = 100) +
  facet_wrap(~population, nrow = 2)


#Code for calculating standard deviation of a population
#https://www.dummies.com/education/math/statistics/standard-deviation-r/

sd.p=function(x){sd(x)*sqrt((length(x)-1)/length(x))}

d1 <- as.data.frame(Coulter.dataClones %>% group_by(population,block,time) %>% 
  summarise(n=n(), median.fL = median(fL), mean.fL = mean(fL), sd.p.fL = sd.p(fL), sd.fL = sd(fL)))

#Order the populations 
d1$population <- factor(d1$population, levels = c("REL606", "Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara-5", "Ara-6", "REL607","Ara+1","Ara+2","Ara+3","Ara+4","Ara+5","Ara+6"))

#True minus between population levels
levels(d1$population) <- c("REL606", "Ara−1", "Ara−2", "Ara−3", "Ara−4", "Ara−5", "Ara−6", "REL607","Ara+1","Ara+2","Ara+3","Ara+4","Ara+5","Ara+6")
d1$time <- factor(d1$time,levels = c("Mid-exponential", "Stationary phase"))
levels(d1$time) <- c("Exponential", "Stationary")

#Plot median estimate of cell size for each replicate and the confidence intervals 

d1 %>% ggplot(aes(x = population, y = median.fL, colour = time, group=time)) + 
  geom_point(size = 2) +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", position="identity", color = "red", width = .5, size=.3) +
  labs(x= "Population", y= "Median cell volume (fL)") +
  scale_y_continuous(limits=c(0,3.2), breaks = seq(.2,3.2,.6)) +
  scale_color_manual(values = c("Exponential" = "Black", "Stationary" = "Grey"), name = "") +
  theme_classic()+
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 3, b=5)), angle = (90), vjust = .5)) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(axis.text.y = element_text(colour = "black",  size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Need to spread data based upon the "time" coloumn. 
d2 <- as.data.frame(spread(d1, time, median.fL))
d2.1 <- na.omit(d2[,c(1,2,7)])
d2.2 <- na.omit(d2[,c(1,2,8)])

#View(d2.1)
#View(d2.2)
#View(d2.3)

d2.3 <- left_join(d2.1,d2.2)

d2.3 <- d2.3 %>% mutate (GSRatio = `Exponential`/`Stationary`) #Absolute fold change

d2.3$descent <- ifelse(d2.3$population == "REL606", "Ancestor", "Evolved")
d2.3$descent <- ifelse(d2.3$population == "REL607", "Ancestor", d2.3$descent)

#d2.3 <- d2.3 %>% mutate (GSRatio = (`Exponential`-`Stationary`)/ (`Exponential`)) 

d2.3 %>% 
  ggplot(aes(x = population, y = GSRatio)) +
  geom_point(size = 2) + 
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", position="identity", color = "red", width = .5, size=.3) +
  labs(x= "Population", y= "Proporiontal difference in cell volume") +
  scale_y_continuous(breaks = c(1.6,1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2)) +
  theme_classic()+
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)), angle = (90))) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(axis.text.y = element_text(colour = "black", size = 12, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Summary data 
#Here I am calculating the confidence interval treating each of the replicates as an entry.
d2.3.1 <- d2.3 %>% 
  group_by(descent) %>%
  summarise(grand.mean.fL = mean(GSRatio), grand.sd.fL = sd(GSRatio), n = n()) %>% 
  mutate(CIlow = grand.mean.fL - abs(qt(0.05/2, n-1)*(grand.sd.fL/sqrt(n)))) %>% 
  mutate(CIHigh = grand.mean.fL + abs(qt(0.05/2, n-1)*(grand.sd.fL/sqrt(n))))

#Plot showing absolute cell size change between ancestor and evolved strain during each growth stage
#Error bars are 95% confidence intervals
# d2.3.1 %>% 
#   group_by(descent) %>% 
#   ggplot(aes(x=descent, y = grand.mean.fL)) + 
#   geom_point() + 
#   geom_errorbar(aes(ymin=`CIlow`, ymax=`CIHigh`), width=.3) +
#   labs(y = "Absolute fold change in median cell volume", x = "Line of descent") +
#   geom_bracket(xmin = "Ancestor", xmax= "Evolved", label = "XXXX", y.position = 2.5)+
#   theme_classic()+
#   theme(axis.text.y = element_text(colour = "black",  size = 14, margin = (margin(l = 5, r=5)))) +
#   theme(axis.title.y = element_text( size = 14)) +
#   theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
#   theme(axis.title.x = element_text(size = 14)) +
#   theme(legend.text = element_text(size = 14)) +
#   theme(legend.title = element_text(size = 14)) +
#   theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#This line of code does the same as above. That is generates absolute cell size change for ancestor and evolved clone
#but shows all points. Error bars 95% confidence intervals. 
d2.3 %>% group_by(descent) %>%
  ggplot(aes(x = descent, y = GSRatio)) +
  scale_y_continuous(breaks = c(1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8)) +
  scale_x_discrete(labels = c("Ancestor", "Evolved")) +
  labs(y = "Proportional difference in cell volume (exponential/stationary)", x = "Strain category") +
  geom_jitter(aes(x = descent, y=GSRatio), width = .25)+
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.3)+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "point", position="identity",size=2.5, color="red") +
  #geom_errorbar(aes(ymin=d2.3.1$CIlow, ymax=d2.3.1$CIHigh), width=.3) +
  geom_bracket(xmin = "Ancestor", xmax= "Evolved", label = "<0.0001", y.position = 2.9)+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black",  size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Filter dataset above as to check norma lity assumptions
d2.3.2 <- d2.3 %>%  filter(descent == "Ancestor") 
d2.3.3 <- d2.3 %>%  filter(descent == "Evolved")

shapiro.test(d2.3.2$GSRatio) #W = 0.91051, p-value = 0.4398 
shapiro.test(d2.3.3$GSRatio) #W = 0.97121, p-value = 0.4596 ## Data are normally distributed

bartlett.test(GSRatio ~ block, data = d2.3.2) #Bartlett's K-squared = 1.377, df = 2, p-value = 0.5023
bartlett.test(GSRatio ~ block, data = d2.3.3) #Bartlett's K-squared = 2.1451, df = 2, p-value = 0.3421 ##Variances homogenous

t.test(GSRatio ~ descent, data = d2.3, alternative = "less", var.equal = T, paired = F) #t = -4.5541, df = 40, p-value = 2.421e-05

#Correlation plot between mid-exponential and stationary phase.
d2.4 <- d2.3 %>% group_by(population) %>%  summarise(mean.exponential = mean(Exponential), 
                                             mean.stationary = mean(Stationary), 
                                             mean.gsratio = mean(GSRatio))

ggplot(d2.4, aes(y=`mean.exponential`, x=`mean.stationary`)) +
  geom_smooth(method = lm, se = F, colour = "red", lty = 2) +
  geom_point(size= 2) +
  scale_y_continuous(breaks = seq(.25,5,.5)) +
  labs(y= "Exponential cell volume (fL)", x= "Stationary cell volume (fL)") +
  theme(panel.background = element_rect(fill = "white"), plot.margin = margin(-.3, 1, 0, 1, "cm"),
  plot.background = element_rect(
  fill = "white",
  size = 1)) +
  theme(axis.text.y = element_text(colour = "black",  size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text( size = 14)) +
  theme(panel.border = element_blank(), axis.line = element_line())

#size = 6.33 x 5.03
#Visualize data to determine wheter normally distributed
ggqqplot(d2.4$`mean.exponential`)
ggqqplot(d2.4$`mean.stationary`)

#Shapiro test to test normality assumptions
shapiro.test(d2.4$`mean.exponential`) #W = 0.90458, p-value = 0.1315
shapiro.test(d2.4$`mean.stationary`) #W = 0.92453, p-value = 0.2553

#Homogenity of variacne
#Calculate mean of medians
#this is the same estimates as in dataframe d2.4.
#Long format.

d4 <- as.data.frame(d1 %>% group_by(population,time) %>% 
                      summarise(mean.fL = mean(median.fL)))

bartlett.test(mean.fL ~ time, data = d4) #p-value = 0.008229 variainces are not equal. Use kendall's 

cor.test(x=d2.4$mean.stationary, y=d2.4$mean.exponential, method = "kendall") 
#tau = 0.7582418, T = 80, p-value = 3.948e-05

#Test whether there is a difference in the mean cell volume of stationary and mid-exponentially growing cells
#Data is normally distributed but the variances are different. Established above. 
#Use wilcoxon tests 

wilcox.test(mean.fL ~ time, data = d4, alternative = "greater") #W = 172, p-value = 0.0001688

#Summary data 
d4 %>% 
      group_by(time) %>%
      summarise(grand.mean.fL = mean(mean.fL), grand.sd.fL = sd(mean.fL), n = n()) %>% 
                      mutate(CIlow = grand.mean.fL - abs(qt(0.05/2, n-1)*(grand.sd.fL/sqrt(n)))) %>% 
                      mutate(CIHigh = grand.mean.fL + abs(qt(0.05/2, n-1)*(grand.sd.fL/sqrt(n))))

d4 %>% group_by(population,time) %>% 
  summarise(mean.fL = mean(mean.fL)) %>% 
  group_by(time) %>% 
  ggplot(aes(x = time, y = mean.fL)) +
  scale_y_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4)) +
  scale_x_discrete(labels = c("Exponential", "Stationary")) +
  labs(y = "Cell volume (fL)", x = "Growth stage") +
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "point", position="identity",size=2) +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width=.3) +
  geom_bracket(xmin = "Exponential", xmax= "Stationary", label = "0.0002", y.position = 2.4)+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black",  size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))




