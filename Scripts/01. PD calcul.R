rm(list=ls())

#####################################################################################################################################
################"Looming extinctions due to invasive species: Irreversible loss of ecological strategy and evolutionary history"#####
#####################################################################################################################################
#author: Celine Bellard

###############################################
###############################################
############ Load package and data ############
###############################################
###############################################

library(dplyr)
library(ape)
library(permute)
library(lattice)
library(vegan)
library(nlme)
library(picante)

#phylogeny from Faurby and Svenning - Complete phylogeny of mammals
phylo2 <- read.nexus("./data/mammals/Complete_phylogeny_Phylacine.nex")
sp<-readRDS(file = "./data/spPhyloMammals.rds")
spIAS<-readRDS(file = "./data/spIASPhyloMammals.rds")

#############################################################################
############### Remove species missing in the Phylogeny #####################
#############################################################################

list_spPhylo <- as.data.frame(phylo2[[1]]$tip.label)
colnames(list_spPhylo)<-"sp"
ListToremove<-anti_join(list_spPhylo, sp)

phylo2_all <-phylo2
class(phylo2_all)
phylo2_all <- .uncompressTipLabel(phylo2_all)

for (i in 1:length(phylo2_all)){
  phylo2_all[[i]] <- drop.tip(phylo2_all[[i]], as.character(ListToremove$sp))
}


ListToremove<-anti_join(list_spPhylo, spIAS)

phylo2_IAST <-phylo2
class(phylo2_IAST)
phylo2_IAST <- .uncompressTipLabel(phylo2_IAST)

for (i in 1:length(phylo2_IAST)){
  phylo2_IAST[[i]] <- drop.tip(phylo2_IAST[[i]], as.character(ListToremove$sp))
}



#############################################################################
############### Build matrix ################################################
#############################################################################


nb<-length(phylo2_all[[1]]$tip.label)
sp<-phylo2_all[[1]]$tip.label
comm_matrix <- matrix(nrow = 1,ncol = nb) 
labels_columns <- sp
colnames(comm_matrix)<-labels_columns

row.names(comm_matrix)<-"bioinv"

comm_matrix[1,1:25]
for (j in 1:nb){
  for (i in 1){
    comm_matrix[i,j]=labels_columns[j]
  }
}

labels_thr_sp_bioinv<-phylo2_IAST[[1]]$tip.label

for (i in 1:length(labels_thr_sp_bioinv)){
  for (j in 1:nb){
    if (comm_matrix[1,j]==labels_thr_sp_bioinv[i]){
      comm_matrix[1,j]=1
    }
  }
}


for (j in 1:nb){
  if (comm_matrix[1,j] != 1){
    comm_matrix[1,j]=0
  }
}

table(comm_matrix)

#############################################################################
############### PD metrics ##################################################
#############################################################################


#Try for one tree
arbre1 <- phylo2_all[[1]]
results_pd <- pd(comm_matrix, arbre1, include.root = T)
results_pd


comm_matrix_tot <- matrix(nrow = 1,ncol = nb) 
str(sp)
labels_columns <- sp
colnames(comm_matrix_tot)<-labels_columns
row.names(comm_matrix_tot)<-"Total"

for (i in (1:ncol(comm_matrix_tot))){
  comm_matrix_tot[1,i] <- 1
}




sum=0
count=0
L=vector()

for (i in 1:99){
  print(i)
  results_pd <- pd(comm_matrix, phylo2_all[[i]], include.root=T)
  target_pd <- results_pd$PD[1]
  tot_pd <-  pd(comm_matrix_tot, phylo2[[i]], include.root = T)
  tot_pd <- tot_pd$PD[1]
  ratio_pd <- target_pd/tot_pd
  L <- c(L,ratio_pd)
  sum=sum+ratio_pd
  count=count+1
}




PD_bioinv_thr_sp <- sum/count
PD_bioinv_thr_sp
med_PD_bioinv_thr_sp <- median(L)
med_PD_bioinv_thr_sp
min_PD_bioinv_thr_sp <- min(L)
min_PD_bioinv_thr_sp
max_PD_bioinv_thr_sp <- max(L)
max_PD_bioinv_thr_sp

L<-as.data.frame(L)


#############################################################################
############### PD metrics - Null model #####################################
#############################################################################

sample(1:99,1)

comm_matrix[1,]<-as.numeric(comm_matrix[1,])

#############################################################################
#Randomization function from Swenson book p128
rand.pd.fun <- function (x){
  tmp.phylo <-tipShuffle(phylo2_all[[1]])#randomize sp names on the phylogeny 
  pd_sub <- pd(comm_matrix, tmp.phylo, include.root = TRUE)[,1]
  pd_tot <- pd(comm_matrix_tot, tmp.phylo,include.root = TRUE)[,1]
  #tot <- pd_tot[1] + pd_tot[2] + pd_tot[3]+pd_tot[4] + pd_tot[5] + pd_tot[6] + pd_tot[7] + pd_tot[8] + pd_tot[9] + pd_tot[10] + pd_tot[11]
  pd_sub/pd_tot
  
}

null.output <-replicate(99, rand.pd.fun(tmp.phylo))
hist(null.output)
abline(v=pd(comm_matrix, phylo2_all[[1]])[1,1]/pd(comm_matrix_tot, phylo2_all[[1]])[1,1],col="red",lwd=2)


pd_sub_obs <- pd(comm_matrix,phylo2_all[[369]], include.root = TRUE)[,1]
pd_tot_obs <- pd(comm_matrix_tot, phylo2_all[[369]], include.root = TRUE)[,1]
observed.value <- pd_sub_obs/pd_tot_obs

#SES & pvalue
ses.all <- (observed.value-mean(null.output))/sd(null.output)
pvalue_bioinv <- length(which(observed.value<=null.output))/ (length(null.output) +1 )

