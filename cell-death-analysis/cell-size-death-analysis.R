#Clear environment
rm(list=ls()) 

#Load libraries
library(tidyverse)
library(cowplot)
library(dplyr)
library(devtools)
library(gridExtra)
library(reshape2)
library(grid)
library(ggplot2)
library(ggpubr)
library(PairedData)
library(car)


#mutate(RelFluor = Fluor2Mean/Fluor1Mean) %>% #dead/alive
#mutate(PropFluor = Fluor2Mean/(Fluor2Mean + Fluor1Mean)) #dead/dead+alive

D1<- read.csv("/Users/nkrumahgrant/Desktop/Cell Morphology/Manuscript-Cell-size-full-submission/cell-death-analysis/cell-death-raw-data.csv")

#str(D1)

D1$Clone <- factor(D1$Clone, levels = c("REL606", "REL11364"))

#View(D1)

D2 <- D1 %>%  
  group_by(Clone, Treatment, Block) %>% 
  count(RelFluorPrime) %>% 
  mutate(PropPrime = n / sum(n))


#Plot 1: Distribution showing relative fluoresences (Dead/Alive)
D1 %>% 
  ggplot(aes(x=as.factor(Block), y=log(RelFluor))) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Clone) +
  labs(y="Log(ratio dead cells)", x = "Block") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), axis.line = element_line())

#Plot 2 : Boxplots showing the distribution of proportions (Dead/Dead+Alive) 
D1 %>% 
  ggplot(aes(x=as.factor(Block), y=log(PropFluor))) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~Clone) +
  labs(y="Log(proportion dead cells)", x = "Block")+
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), axis.line = element_line())

# Here I plot a 2D density plot to look at whether I can can score cells as dead 
# or alive given some value of fluoresence. 

#Plot 3 : density plots showing seperation 
#log-value(x) > 1 dead. 
D1 %>%
  ggplot(aes(x=(log(RelFluor)), y=log(PropFluor), color = RelFluorPrime)) +
  geom_density2d() +
  facet_wrap(~Clone) +
  labs(x="Log(proportion dead cells)", y="Log(ratio dead cells)") +
  guides(color=guide_legend(title="Cell viability"))

#Plot 4: Association between cell size and death?

ggplot(D1, aes(ShortAxis, LongAxis, color = RelFluorPrime)) + 
  geom_point(size=2) +
  scale_color_manual(name = "Viability", values = c("Green", "Red")) +
  facet_grid(~Clone + RelFluorPrime) +
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = margin(l = 5, r = 5))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = margin(b = 5, t =5))) +
  theme(strip.background = element_rect(color = "black", fill = "white")) +
  theme(strip.text =element_text(size = 14, color = "black")) +
  xlab(expression(paste("Cell width" ," (",mu,"m) "))) +
  ylab(expression(paste("Cell length" ," (",mu,"m) "))) +
  theme(legend.position = "none")

ggplot(D1, aes(x=Clone, y=LongAxis, fill = RelFluorPrime)) +
  geom_boxplot() + 
  scale_fill_manual(name = "Viability", values = c("Green", "Red")) +
  theme(legend.position = "none") +  
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = margin(l = 5, r = 5))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = margin(b = 5, t =5))) +
  theme(strip.background = element_rect(color = "black", fill = "white")) +
  theme(strip.text =element_text(size = 14, color = "black")) +
  ylab(expression(paste("Cell length" ," (",mu,"m) "))) +
  theme(legend.position = "none")

#Cell width and cell viability plot
ggplot(D1, aes(x=Clone, y=ShortAxis, fill = RelFluorPrime)) +
  geom_boxplot() + 
  scale_fill_manual(name = "Viability", values = c("Green", "Red")) +
  theme(legend.position = "none") +  
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = margin(l = 5, r = 5))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = margin(b = 5, t =5))) +
  theme(strip.background = element_rect(color = "black", fill = "white")) +
  theme(strip.text =element_text(size = 14, color = "black")) +
  ylab(expression(paste("Cell width" ," (",mu,"m) "))) +
  theme(legend.position = "none")

#Cell length and cell viability plot
ggplot(D1, aes(x=Clone, y=LongAxis, fill = RelFluorPrime)) +
  geom_boxplot() + 
  scale_fill_manual(name = "Viability", values = c("Green", "Red")) +
  theme(legend.position = "none") +  
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = margin(l = 5, r = 5))) +
  theme(axis.title.y = element_text( size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text( size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = margin(b = 5, t =5))) +
  theme(strip.background = element_rect(color = "black", fill = "white")) +
  theme(strip.text =element_text(size = 14, color = "black")) +
  ylab(expression(paste("Cell length" ," (",mu,"m) "))) +
  theme(legend.position = "none")


D1.606 <- D1 %>%  filter(Clone == "REL606")
D1.11364 <- D1 %>%  filter(Clone == "REL11364")

D1.dead <- D1 %>% filter(RelFluorPrime == "Dead")


#Two sample t-tests 

#F-test for homogenity in variances
var.test(LongAxis ~ RelFluorPrime, data = D1.606) #p < 2.2e-16 
var.test(ShortAxis ~ RelFluorPrime, data = D1.606) #p = 1.252e-07
var.test(LongAxis ~ RelFluorPrime, data = D1.11364) #p = 0.6405, variances are equal 
var.test(ShortAxis ~ RelFluorPrime, data = D1.11364) #p < 2.2e-16

var.test(ShortAxis ~ Clone, data = D1.dead) #p <2.2e-16
var.test(LongAxis ~ Clone, data = D1.dead) #p <2.2e-16

#Is the cell length and width of dead ancestral cells smaller than than that of the cit+ clone?
#This is an obvious question. Cells have increased in cell size, so of course dead cells will be 
#bigger between populations.It makes more sense to ask whether there are differences between 
#the cell lengths and widths of living and dead cells WITHIN populaitons. 
t.test(LongAxis ~ Clone, data = D1.dead, var.equal = FALSE, alternative = "l") #p < 2.2e-16
t.test(ShortAxis ~ Clone, data = D1.dead, var.equal = FALSE, alternative = "l") #p < 2.2e-16


#Is the cell size of living cells less than that of dead cells within populations? No, all p > .99. 
t.test(LongAxis ~ RelFluorPrime, data = D1.606, var.equal = FALSE, alternative = "l") 
t.test(ShortAxis ~ RelFluorPrime, data = D1.606, var.equal = FALSE, alternative = "l")

t.test(LongAxis ~ RelFluorPrime, data = D1.11364, var.equal = TRUE, alternative = "l") 
t.test(ShortAxis ~ RelFluorPrime, data = D1.11364, var.equal = FALSE, alternative = "l") 


#Plot 5 : proportion of dead:alive cells

#Figure for cell size paper
#https://www.r-graph-gallery.com/48-grouped-barplot-with-ggplot2.html


D2$RelFluorPrime <- factor(D2$RelFluorPrime, levels = c("Dead", "Alive"))

D2 %>% 
  ggplot(aes(x=Clone, y = n, fill = RelFluorPrime)) +
  geom_bar(stat = "identity", position = "fill") + 
  ylab("Proportion of cells") +
  scale_fill_manual(breaks = c("", ""), values=c("Red", "Green"))+
  labs(fill = "Cell viability") + 
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 14, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(colour = "black", size = 14, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(panel.border = element_blank(), axis.line = element_line())

#Plot 6: plot of the proportions for indepent replicates. 
#First must construct the data frame piece wise.

D2.2 <- D2 %>%  filter(RelFluorPrime == "Dead") %>% rename(P.Dead = n) %>% rename(PropPrime.Dead = PropPrime)
D2.3 <- D2 %>%  filter(RelFluorPrime == "Alive") %>% rename(P.Alive = n) %>%  rename(PropPrime.Alive = PropPrime)

D2.2$P.Alive <- D2.3$P.Alive
D2.2$PropPrime.Alive <- D2.3$PropPrime.Alive

D2.2 <- D2.2 %>%  
  mutate(N = P.Dead + P.Alive)

D2.2 <- D2.2[,-4] #removes RelFluorPrime column 

#Statistical analysis 
#Is the proportion of death higher in REL11364 than the ancestor?
#H0 = There is no difference in the amount of death (null = 0)
#HA = death is higher in clone REL11364, use less referring to the ref pop being less than.
 
#Assumptions:independent samples, data for the two groups follow normal distribution, homogenity in variances

#Check Assumptions

#Shapiro-Wilk normality test for proportion of dead cells in REL606
with(D2.2, shapiro.test(PropPrime.Dead[Clone == "REL606"])) #p = 0.3719

#Shapiro-Wilk normality test for proportion of dead cells in REL11364
with(D2.2, shapiro.test(PropPrime.Dead[Clone == "REL11364"])) #p = 0.2345 

#Shapiro-Wilk test indicates that the distribution of the data is not significantly different
#from the normal distribution. 

#F-test for homogenity in variances
res.ftest <- var.test(PropPrime.Dead ~ Clone, data = D2.2)
res.ftest #p = 0.07019

#Variances are homogenous. All three assumptions met for the classic t-test 

#T.test
t.test(PropPrime.Dead ~ Clone, data = D2.2, var.equal = TRUE, alternative = "less") 
#t = -2.9304, df = 8, p-value = 0.009494

#Percentage dead cells 
D2 %>% group_by(Clone, RelFluorPrime) %>% 
  summarise(mean.prop = mean(PropPrime))
