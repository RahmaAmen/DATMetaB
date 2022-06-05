# DATMetaB R script calculates :
#         1) Diet overlap indexes 
#         2) The related statistics
#         3) NMDS
#
# Note: Standard libraries vegan, ggplot2, and spaa are required
# Note: External libraries pairwise.adonis2 and chull are required
#
# Rahma Amen 25.05.2022 Potsdam University

## 1) clean environment and load library

rm(list=ls())

setwd("/path/to folder/of DATMetaB matlab program/")

library(vegan)
library(ggplot2)
library(spaa)

source("/path/to/pairwise.adonis2.R")
source("clacshannon.R")
source("chull.R")

set.seed(1000) #reproducible results

## 2) Load data (pooled and individual)

# pooled
data.classp<-read.csv("RRA_CLASS.csv")
data.orderp<-read.csv("RRA_ORDER.csv")
data.familyp<-read.csv("RRA_FAMILY.csv")
data.genusp<-read.csv("RRA_GENUS.csv")
data.speciesp<-read.csv("RRA_SPECIES.csv")

# Individual
data.classi<-read.csv("RRA_CLASS.csv")
data.orderi<-read.csv("RRA_ORDER.csv")
data.familyi<-read.csv("RRA_FAMILY.csv")
data.genusi<-read.csv("RRA_GENUS.csv")
data.speciesi<-read.csv("RRA_SPECIES.csv")

# Pooled data in columns form
rra.classp<-t(as.matrix(data.classp[,2:ncol(data.classp)]))
rra.orderp<-t(as.matrix(data.orderp[,2:ncol(data.orderp)]))
rra.familyp<-t(as.matrix(data.familyp[,2:ncol(data.familyp)]))
rra.genusp<-t(as.matrix(data.genusp[,2:ncol(data.genusp)]))
rra.speciesp<-t(as.matrix(data.speciesp[,2:ncol(data.speciesp)]))

# Individual data in rows form
rra.classi<-as.matrix(data.classi[,4:ncol(data.classi)])
rra.orderi<-as.matrix(data.orderi[,4:ncol(data.orderi)])
rra.familyi<-as.matrix(data.familyi[,4:ncol(data.familyi)])
rra.genusi<-as.matrix(data.genusi[,4:ncol(data.genusi)])
rra.speciesi<-as.matrix(data.speciesi[,4:ncol(data.speciesi)])

## 4) Diet indices of the pooled samples (Pianka and Schöner)

# calculate Pianka indices
class.pianka<-niche.overlap(rra.classp, method = "pianka")
order.pianka<-niche.overlap(rra.orderp, method = "pianka")
family.pianka<-niche.overlap(rra.familyp, method = "pianka")
genus.pianka<-niche.overlap(rra.genusp, method = "pianka")
species.pianka<-niche.overlap(rra.speciesp, method = "pianka")

# calculate Schoener indices
class.schoener<-niche.overlap(rra.classp, method = "schoener")
order.schoener<-niche.overlap(rra.orderp, method = "schoener")
family.schoener<-niche.overlap(rra.familyp, method = "schoener")
genus.schoener<-niche.overlap(rra.genusp, method = "schoener")
species.schoener<-niche.overlap(rra.speciesp, method = "schoener")

## 5) Dispersion Dietary niche width and composition

# Shannon index mean and SD (we have to transpose the data)
class.shanon<-calcshannon(t(rra.classi))
order.shanon<-calcshannon(t(rra.orderi))
family.shanon<-calcshannon(t(rra.familyi))
genus.shanon<-calcshannon(t(rra.genusi))
species.shanon<-calcshannon(t(rra.speciesi))

# Shannon index figure (boxplot) at order level
a<-as.vector(
  unlist(
    niche.width(t(rra.orderi), method = "shannon")
  )
)
a1<-data.frame(group = "C. alces", value = a[1:2])
a2<-data.frame(group = "C. compress.", value = a[3:12])
a3<-data.frame(group = "C. curvirostris", value = a[13:14])
a4<-data.frame(group = "C. numenius", value = a[15])
a5<-data.frame(group = "G. petersii", value = a[16:18])
a6<-data.frame(group = "C. tshokwe", value = a[19:27])

plot.data<-rbind(a1, a2, a3, a4, a5, a6)
ggplot(plot.data, aes(x=group, y=value)) +
  geom_boxplot(colour = "black", fill="gray") +
  ylab("Shannon diversity index") +
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color="black",fill=NA, size =0.75 ),
        plot.title = element_blank(),
        axis.title.y = element_text(color="black", size=18, family="Times"),
        axis.title.x = element_blank(),
        axis.text.x=element_text(colour="black", size = 14,face="italic", family="Times"),
        axis.text.y=element_text(colour="black", size = 14, family="Times")) 

ggsave("/path/to/output/folder/Shannon_order.pdf")
ggsave("/path/to/output/folder/Shannon_order.png")

## 6) Run perMANOVA statistics based on Bray–Curtis dissimilarity

# Bray–Curtis dissimilarity
classp.bc<-vegdist(sqrt(t(rra.classp)), method='bray')
orderp.bc<-vegdist(sqrt(t(rra.orderp)), method='bray')
familyp.bc<-vegdist(sqrt(t(rra.familyp)), method='bray')
genusp.bc<-vegdist(sqrt(t(rra.genusp)), method='bray')
speciesp.bc<-vegdist(sqrt(t(rra.speciesp)), method='bray')

# perMANOVA based on species groups (food taxa at order level)
order.species.div <- 
  adonis2(rra.orderi~SPECIES, data=data.orderi,
          permutations = 999, method="bray", sqrt.dist = TRUE)

# perMANOVA based on EOD groups (food taxa at order level)
order.eod.div <- 
  adonis2(rra.orderi~EOD, data=data.orderi,
          permutations = 999, method="bray", sqrt.dist = TRUE)

# perMANOVA based on SNOUT groups (food taxa at order level)
order.snout.div <- 
  adonis2(rra.orderi~SNOUT, data=data.orderi,
          permutations = 999, method="bray", sqrt.dist = TRUE)

# pairwise perMANOVA based on species groups (food taxa at order level)
order.species.pwdiv <- 
  pairwise.adonis2(rra.orderi~SPECIES, data=data.orderi)

# pairwise perMANOVA based on eod groups (food taxa at order level)
order.eod.pwdiv <- 
  pairwise.adonis2(rra.orderi~EOD, data=data.orderi)

# pairwise perMANOVA based on snout
order.snout.pwdiv <- 
  pairwise.adonis2(rra.orderi~SNOUT, data=data.orderi)

## 7) NMDS

# Run metaMDS and extract data
order.MDS <- metaMDS(sqrt(rra.orderi), distance="bray", k=2, trymax=35, autotransform=TRUE)
NMDS1 <- order.MDS$points[,1]
NMDS2 <- order.MDS$points[,2]
order_NMDS_DATA <- as.data.frame(NMDS1)
order_NMDS_DATA$NMDS2<- NMDS2
order_NMDS_DATA$SPECIES <- data.orderi$SPECIES
order_NMDS_DATA$EOD <- data.orderi$EOD
order_NMDS_DATA$SNOUT <- data.orderi$SNOUT

# Plot and save stress
pdf("/path/to/output/folder/Stress_order.pdf")
stressplot(order.MDS)
dev.off()

# Plot NMDS
ggplot(order_NMDS_DATA, aes(x=NMDS1, y=NMDS2)) +
  geom_point(size = 3, aes( shape=SPECIES, colour=SNOUT)) +
  scale_shape_manual(values=c(15, 16, 17, 18, 10, 12))+
  stat_chull(aes(colour=EOD, fill=EOD), data=order_NMDS_DATA, alpha=0.2) +
  annotate("text", x=0.75, y=min(NMDS2), label=paste('Stress =',round(order.MDS$stress,3))) +
  theme(axis.text.y = element_text(colour="black", size=12, family="Times"), 
        axis.text.x = element_text(colour="black", family="Times", size=12), 
        legend.text = element_text(size = 12, face="italic",family="Times", colour ="black"), 
        legend.position = "right",
        axis.title.y = element_text(family="Times", size=14), 
        axis.title.x = element_text(family="Times", size=14, colour="black"), 
        legend.title = element_text(size=14), 
        panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA, size=1),
        legend.key=element_blank()) 

ggsave("/path/to/output/folder/NMDS_order.pdf")
ggsave("/path/to/output/folder/NMDS_order.png")

