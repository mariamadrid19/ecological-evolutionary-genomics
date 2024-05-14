setwd("/Users/mariamadrid/Documents/Eco-Evo-Genomics/lab8")
beetle <- read.delim("~/Documents/Eco-Evo-Genomics/lab8/WG_table.txt")
boxplot(beetle$Fst~beetle$comp)

library(ggplot2)
library(RColorBrewer)

comp_colors <- brewer.pal(6, "Set3")

#plot Fst as histogram to determine the distribution of differences between the populations 
ggplot(data=beetle)+geom_histogram(aes(x=Fst,col=comp))+facet_grid(comp~.)+scale_y_log10()

#plot a scatterplot of the Fst values for each population across the entire genome, locate positional data
ggplot(data=beetle)+geom_point(aes(x=genome_pos,y=Fst,col=comp))+facet_grid(comp~.)

beetleL <- subset(beetle, beetle$Fst < 0.2)
beetleH <- subset(beetle, beetle$Fst > 0.2)

ggplot(data=beetleH)+geom_boxplot(aes(x=comp,y=TajDPop1,col=comp))+theme_minimal()
ggplot(data=beetleH)+geom_boxplot(aes(x=comp,y=TajDPop2,col=comp))+theme_minimal()

beetle_BELGIUM <- subset(beetleH, beetleH$comp=="Belgium")
beetle_FRANCE <- subset(beetleH, beetleH$comp=="France")

merged_beetles <- merge(beetle_BELGIUM,beetle_FRANCE,by="genome_pos")
