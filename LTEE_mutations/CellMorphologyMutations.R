rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyverse)
library(grid)
library(reshape2)

setwd("/Users/nkrumahgrant/Desktop/Cell Morphology/Cell-size-full-submission/LTEE_mutations/")

ShinyData <- list.files(pattern = "*.csv") 
Data <- as.data.frame(do.call(rbind,Map('cbind', lapply(ShinyData, read.csv)))) %>% 
  select(-starts_with("html"), -gene_name, -snp_type, -X)

RodMain <- Data %>% filter(gene_list %in% c("mreB", "mreC", "mreD", "mrdA", "mrdB"))

#View(RodMain)

#Arc <- Data %>% filter(grepl(pattern="arc", tolower(gene_list))) #related to anaerobic project 

CloneA <- RodMain %>% filter(clone == "A") 

#NonMutators <- CloneA  %>%  filter(mutator_status == "IS-mutator"| mutator_status == "non-mutator" )
##by_population <- CloneA %>%  group_by(population) %>% arrange(population) 

CloneA$population <- factor(CloneA$population, c("Ara-1","Ara-2","Ara-3","Ara-4","Ara-5","Ara-6","Ara+1","Ara+2","Ara+3","Ara+4","Ara+5","Ara+6"))

#Using true minus symbol
levels(CloneA$population) <- c("Ara−1","Ara−2","Ara−3","Ara−4","Ara−5","Ara−6","Ara+1","Ara+2","Ara+3","Ara+4","Ara+5","Ara+6")

# Manual levels
gene_list_table <- table(CloneA$gene_list)
gene_list_levels <- names(gene_list_table)[order(gene_list_table)]
CloneA$gene_list_2 <- factor(CloneA$gene_list, levels = gene_list_levels)

#gene_list_table <- table(NonMutators$gene_list)
#NonMutators$gene_list_2 <- factor(NonMutators$gene_list, levels = gene_list_levels)

#Change labels for facet
CloneA$time <- factor(CloneA$time)
levels(CloneA$time) <- c("Generation 2k", "Generation 10k", "Generation 50k")

#7.85 x 5.35
ggplot(CloneA, aes(x=population, y=gene_list, height = .5, width = 1)) + 
 geom_tile(aes(fill = mutation_category)) + 
  scale_x_discrete(drop = F)+
  facet_grid (~time) +
  theme_bw() + 
  scale_fill_discrete(name = "Mutation", labels = c("Indel","Nonsynonymous", "Synonymous")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  theme(legend.justification=c(0,1), legend.position=c(0.005,.995)) +
  theme(axis.text.y = element_text(face = "italic", colour = "black", size = 12, margin = (margin(l = 5, r=5)))) +
  theme(axis.title.y = element_blank()) +
  theme(axis.ticks = element_blank())+
  theme(legend.position = "bottom",
        legend.margin = margin (10,0,0,100), 
        legend.box.margin = margin(-20,0,0,10))+
  theme(axis.text.x = element_text(colour = "black", size = 12, margin = (margin(t = 5, b=5)))) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 12)) +
  theme(strip.text = element_text(size=12)) +
  coord_fixed(ratio=1.5)

#binomial tests. Treating as two families, Mre and Mrd. 
# binom.test(5,8,p = 0.40, alternative = "g") #p=0.1737
# binom.test(6,12,p=0.40, alternative = "g") #p=0.3348
# binom.test(8,20,p=0.40,  alternative = "g") #p=0.5841


