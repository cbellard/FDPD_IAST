rm(list=ls())

#####################################################################################################################################
################"Looming extinctions due to invasive species: Irreversible loss of ecological strategy and evolutionary history"#####
#####################################################################################################################################
#author: Céline Bellard

###############################################
###############################################
############ Load package and data ############
###############################################
###############################################
library(dplyr)
library(stringr)
library(varhandle)
library(gtools)

setwd("./script/Functions")
source("species_to_FE.R")
source("FE_metrics.R")
source("quality_funct_space.R")
source("plot_funct_space.R")
source("multidimFD.R")

traits_f<-readRDS(file = "./data/traits_f.rds")
presence<-readRDS(file = "./data/presenceMammalsEcosyt.rds")

table(presence)
###############################################
###############################################
####### Computing FUNCTIONAL ENTITIES #########
###############################################
###############################################

# Folder where outputs will be saved
setwd("D:/Documents/Projets/ImpactTPF2/")

# Grouping species into FE
species_to_FE_mammals<-species_to_FE(traits_f)

# Number of species per FE
apply(species_to_FE_mammals$FE_sp_01,1,sum)             # from 1 to 82 species per FE

# FE to which each species belongs
FE<-species_to_FE_mammals$FE
FE_mammals_01<-species_to_FE_mammals$FE_sp_01

# Trait values of the FE
FE_traits<-species_to_FE_mammals$FE_traits

# Matrix of FE biomass in assemblages
assemblage_FE_weight<-matrix(0, nrow(presence), nrow(FE_mammals_01),
                             dimnames=list( row.names(presence), row.names(FE_mammals_01) ) ) # empty matrix

for (k in row.names(FE_mammals_01) ) # loop on FE
{
  sp_k<-names(which(FE_mammals_01[k,]==1))
  if( length(sp_k)==1 ) {assemblage_FE_weight[,k]<-presence[,sp_k] } # if only one species in FE k
  if(length(sp_k)>1 ) {assemblage_FE_weight[,k]<-apply(presence[,sp_k],1,sum)  } # if more than 1 species in FE k
}# end of k

sum(presence)==sum(assemblage_FE_weight) # check total biomass kept constant

# Matrix of FE occurence (0/1) in assemblages
assemblage_FE_occ<-assemblage_FE_weight
assemblage_FE_occ[which(assemblage_FE_occ>0)]<-1

# Computing diversity metrics based on Funct ent for the set of fruits baskets studied
assemblages_FE_metrics<-FE_metrics(species_to_FE_mammals, presence, check_species_pool=FALSE, folder_plot="D:/Documents/Projets/ImpactTPF/Results", nm_asb_plot=row.names(presence))
round(assemblages_FE_metrics,3)

# plots illustrating distribution of species into Funct entities are in the subfolder in your working directory.

##############################################################################################
# Computing MULTIDIMENSIONAL FUNCTIONAL DIVERISTY INDICES based on FE position in a functional space to assess
# how weight is distributed in the functional space independently from packing of species into FE (which is assessed by metrics presented above)

# Computing Functional space based on trait value of Functional entities (not based on species trait values
# because we want to represent distance between combinations of trait values independently from their frequency among species, i.e. to give same weight to each FE whatever its number of species)
qual_funct_space_FE<-quality_funct_space(FE_traits, traits_weights=NULL, nbdim=5, metric="Gower", dendro=FALSE, plot="quality_funct_space_FE")
qual_funct_space_FE$meanSD # 
 #mat_dissim <- daisy(FE_traits, metric = "gower", weights = 1)
 #mat_pcoa <- pcoa(mat_dissim)
 #mat_pcoa$values

# FE coordinates in the best space
coord_FE_6D<-qual_funct_space_FE$details_funct_space$mat_coord[,1:5]

# Species coordinates in the space according to those of FE
coord_sp_6D<-coord_FE_6D[FE,]
row.names(coord_sp_6D)<-names(FE)

# Computing FD indices according to species weights
presence <- presence[ , order(colnames(presence))]
coord_sp_6D <- coord_sp_6D[order(rownames(coord_sp_6D)) , ]

FD_assemblage_sp <- multidimFD(coord_sp_6D, presence, check_species_pool=FALSE, verb=TRUE  )
FD_assemblage_FE <- multidimFD(coord_FE_6D, assemblage_FE_weight, check_species_pool=FALSE, verb=TRUE  )

# Comparing values between species- and FE-based FD indices
round(rbind( FRic_sp=FD_assemblage_sp[,"FRic"], FRic_FE= FD_assemblage_FE[,"FRic"]) ,3)

# GENERAL COMMENT:
# computing FD indices based on species or FE weights are both correct, choice should be done according to question addressed :
#  - species-based indices are relevant to detect assembly rules since processes such as dispersion and competition act on species
#  - FE-based indices are relevant to assess the links between FD and ecosystem processes since redundant species are expected to have same ecological roles

###############################################
###############################################
####### Computing NULL MODELS  ################
###############################################
###############################################

rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec <= 0)) {
    negs <- vec <= 0
    pos  <- vec > 1
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}
head(FD_assemblage_FE)
#################################################################################################
nbrep<-99 # number of replicates
nbFE <- FD_assemblage_sp["TRUE","Nb_sp"]
nbFE # 43 FEs
nbSp <- FD_assemblage_sp["TRUE","Tot_weight"]
nbSp # 43 species
#Empty matrix to store simulated species occurences
BRU_H0<-matrix(0, nbrep, ncol(assemblage_FE_weight), dimnames=list(1:nbrep, colnames(assemblage_FE_weight)))
for (k in 1:nbrep)
{
  BRU_H0[k, sample(colnames(assemblage_FE_weight), nbFE) ] <- 1                                     # random sorting of FEs
  BRU_H0[k,] <- replace(BRU_H0 [k,], BRU_H0 [k,] != 0, rand_vect(nbFE, nbSp, pos.only = TRUE))    # random sorting of Species
  print(sum(BRU_H0[k,]))                #Number of species
  print(ncol(BRU_H0[,BRU_H0[k,]!= 0]))  #Number of FEs
}# end of k
#Computing FD indices on these assemblages, check_species_pool=FALSE since by chance some species could be never picked but this is not an issue
FD_assemblage_FE_BRU_H0 <- multidimFD(coord_FE_6D, BRU_H0, check_species_pool=FALSE, verb=TRUE  )

summary(FD_assemblage_FE_BRU_H0)


##Comparing observed and expected values under H0 using SES and p-value metrics
SES_FRic_AA <- (FD_assemblage_FE["TRUE","FRic"]-mean(FD_assemblage_FE_BRU_H0[,"FRic"]) ) / sd(FD_assemblage_FE_BRU_H0[,"FRic"])
##Comparing observed and expected values under H0 using SES and p-value metrics
SES_FRic_AA <- (FD_assemblage_FE["TRUE","FRic"]-mean(FD_assemblage_FE_BRU_H0[,"FRic"]) ) / sd(FD_assemblage_FE_BRU_H0[,"FRic"])
SES_FRic_AA # SES<(-1) means that observed FRic is lower than expected
pvalue_FRic_AA <- length(which(FD_assemblage_FE["TRUE","FRic"]<=FD_assemblage_FE_BRU_H0[,"FRic"]))/ ( length(FD_assemblage_FE_BRU_H0[,"FRic"]) +1 )
pvalue_FRic_AA # p-value >0.975 => FRic is significantly lower than expected under H0


Fric<- cbind(FD_assemblage_FE["TRUE","FRic"],SES_FRic_AA,pvalue_FRic_AA)
colnames(Fric)<-c("Obs","SES","Pvalue")





