
## Elena Martínez Castillejo & H. Christoph Liedtke

##HYPOTHESIS TESTING FOR VOCAL SAC EVOLUTION.
  ##Continuous time Markov models for discrete character evolution --- fitMk()
  ##Seven states vocal sac coding

##Library---------------------------------------------------------------------------

library(ape)
library(tidyverse)
library(readxl)
library(qgraph)
library(phytools)


##1- DATA LOADING AND CLEANING---------------------------------------------------------

#setwd("~/EEDcloud/TFM/Elena/Anc.Rec and Stoch.Mapping/")

 #Load data
tree <- read.tree("../../data/clean_data/portik_phy_asw2026.tre")
data <- read_xlsx ("../../data/clean_data/vocal_sac_coding.xlsx", sheet= "coding") 

 #Clean data

data$seven_states_code %>% unique()

trait.df <- data %>%
  distinct(Species_asw, .keep_all = TRUE) %>%
  filter(Species_asw %in% tree$tip.label) %>%
  arrange(match(Species_asw, tree$tip.label)) %>%
  mutate(
    state_probs =case_when(
      seven_states_code == "?" ~ list(setNames(rep(1/7,7), 0:6)),
      seven_states_code == "[0 1]" ~ list(c("0" = 0.5, "1" = 0.5)),
      seven_states_code == "[3 1]" ~ list(c("3" = 0.5, "1" = 0.5)),
      seven_states_code == "[1 3]" ~ list(c("1" = 0.5, "3" = 0.5)),
      seven_states_code == "[3 4]" ~ list(c("3" = 0.5, "4" = 0.5)),
      seven_states_code == "[0 5]" ~ list(c("0" = 0.5, "5" = 0.5)),
      seven_states_code == "[0 3]" ~ list(c("0" = 0.5, "3" = 0.5)),
      seven_states_code == "[0 4]" ~ list(c("0" = 0.5, "4" = 0.5)),
      TRUE ~ lapply(seven_states_code, function(x) setNames(1.0, x))
    )) %>%
  filter(!is.na(seven_states_code))

states <- as.character(0:6)

   ##Making a matrix for character state probabilities (uncertainties)
trait <- do.call(rbind, lapply(trait.df$state_probs, function(p) {
  res <- setNames(rep(0, length(states)), states)
  res[names(p)] <- p
  return (res)
}))

rownames(trait) <- trait.df$Species_asw

    #colSums(trait)

    #trait.df %>%
   #  count(seven_states_code, name = "num_species")

## prune tree to match data
 clean_tree<-keep.tip(tree, rownames(trait))
 clean_tree
 
 all(rownames(trait)==clean_tree$tip.label)

states <- colnames(trait) 
n <- length (states)

q<- list()

# DEFINING TRANSITION RATE MATRIX FOR ALL MODELS-------------
 #2.1 FULL MODELS-------------------------------------------------

## ER
q$er<-matrix(nrow=n, ncol=n,1) # set all transitions to be the same (index=1)
diag(q$er) <- 0 # make diagonals 0

## SYM
q$sym<-q$er
q$sym[lower.tri(q$sym)] <- 1:(n*3)
q$sym[upper.tri(q$sym)] <- t(q$sym)[upper.tri(q$sym)]

## ARD
q$ard<-q$er
q$ard[lower.tri(q$ard)] <- 1:(n*3)
q$ard[upper.tri(q$ard)] <- 22:(n*3+21)

  #2.2 RADIAL MODELS---------------------------------------------------------

## ER
q$er_radial<-matrix(nrow=n, ncol=n,0) # set all transitions to be 0
q$er_radial[2:n,1]<-1
q$er_radial[1,2:n]<-1
q$er_radial

## SYM
q$sym_radial<-matrix(nrow=n, ncol=n,0) # set all transitions to be 0
q$sym_radial[2:n,1]<-1:6
q$sym_radial[1,2:n]<-1:6
q$sym_radial

## ARD
q$ard_radial<-matrix(nrow=n, ncol=n,0) # set all transitions to be 0
q$ard_radial[2:n,1]<-1:6
q$ard_radial[1,2:n]<-7:12
q$ard_radial

 #2.3 HYPOTHESIS BASED MODELS-------------------------------------------------
 #Making a function to build matrix easily
T.rate <- function (matrix, from, to, value) {
  matrix[from +1, to +1] <- value
  return(matrix)
}

##Absent to Single Subgular
  #In this hypothesis we state that Single subgular (1) is always the first state to appear
  #and vocal sacs are also lost from a Single subgular state only.

#ER
q$ab.to.Ssub_er <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0, to = 1, value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 2, value = 1) %>%
  T.rate (from = 2, to = 1, value = 1) %>%
  T.rate (from = 1, to = 3, value = 1) %>%
  T.rate (from = 3, to = 1, value = 1) %>%
  T.rate (from = 1, to = 4, value = 1) %>%
  T.rate (from = 4, to = 1, value = 1) %>%
  T.rate (from = 1, to = 5, value = 1) %>%
  T.rate (from = 5, to = 1, value = 1) %>%
  T.rate (from = 1, to = 6, value = 1) %>%
  T.rate (from = 6, to = 1, value = 1)

rownames(q$ab.to.Ssub_er) <- colnames(q$ab.to.Ssub_er) <- states
q$ab.to.Ssub_er

#SYM
q$ab.to.Ssub_sym <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0, to = 1, value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 2, value = 2) %>%
  T.rate (from = 2, to = 1, value = 2) %>%
  T.rate (from = 1, to = 3, value = 3) %>%
  T.rate (from = 3, to = 1, value = 3) %>%
  T.rate (from = 1, to = 4, value = 4) %>%
  T.rate (from = 4, to = 1, value = 4) %>%
  T.rate (from = 1, to = 5, value = 5) %>%
  T.rate (from = 5, to = 1, value = 5) %>%
  T.rate (from = 1, to = 6, value = 6) %>%
  T.rate (from = 6, to = 1, value = 6)

rownames(q$ab.to.Ssub_sym) <- colnames(q$ab.to.Ssub_sym) <- states
q$ab.to.Ssub_sym

#ARD
q$ab.to.Ssub_ard <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0, to = 1, value = 1) %>%
  T.rate (from = 1, to = 0, value = 2) %>%
  T.rate (from = 1, to = 2, value = 3) %>%
  T.rate (from = 2, to = 1, value = 4) %>%
  T.rate (from = 1, to = 3, value = 5) %>%
  T.rate (from = 3, to = 1, value = 6) %>%
  T.rate (from = 1, to = 4, value = 7) %>%
  T.rate (from = 4, to = 1, value = 8) %>%
  T.rate (from = 1, to = 5, value = 9) %>%
  T.rate (from = 5, to = 1, value = 10) %>%
  T.rate (from = 1, to = 6, value = 11) %>%
  T.rate (from = 6, to = 1, value = 12)

rownames(q$ab.to.Ssub_ard) <- colnames(q$ab.to.Ssub_ard) <- states
q$ab.to.Ssub_ard

##Paired individual (Paired_ind)
 ##We test if paired vocal sacs evolved individualy, starting from paired subgular (4) as basal state for paire vocal sacs.

#ER
q$paired_ind_er <- matrix (0, nrow= n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 3, value = 1) %>%
  T.rate (from = 3, to = 1, value = 1) %>%
  T.rate (from = 0, to = 4, value = 1) %>%
  T.rate (from = 4, to = 0, value = 1) %>%
  T.rate (from = 1, to = 2, value = 1) %>%
  T.rate (from = 2, to = 1, value = 1) %>%
  T.rate (from = 4, to = 5, value = 1) %>%
  T.rate (from = 5, to = 4, value = 1) %>%
  T.rate (from = 5, to = 6, value = 1) %>%
  T.rate (from = 6, to = 5, value = 1) %>%
  T.rate (from = 4, to = 6, value = 1) %>%
  T.rate (from = 6, to = 4, value = 1) 
rownames(q$paired_ind_er) <- colnames(q$paired_ind_er) <- states
q$paired_ind_er

#SYM
q$paired_ind_sym <- matrix (0, nrow= n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 3, value = 2) %>%
  T.rate (from = 3, to = 1, value = 2) %>%
  T.rate (from = 0, to = 4, value = 3) %>%
  T.rate (from = 4, to = 0, value = 3) %>%
  T.rate (from = 1, to = 2, value = 5) %>%
  T.rate (from = 2, to = 1, value = 5) %>%
  T.rate (from = 4, to = 5, value = 6) %>%
  T.rate (from = 5, to = 4, value = 6) %>%
  T.rate (from = 5, to = 6, value = 7) %>%
  T.rate (from = 6, to = 5, value = 7) %>%
  T.rate (from = 4, to = 6, value = 8) %>%
  T.rate (from = 6, to = 4, value = 8) 
rownames(q$paired_ind_sym) <- colnames(q$paired_ind_sym) <- states
q$paired_ind_sym

##ARD
q$paired_ind_ard <- matrix (0, nrow= n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 2) %>%
  T.rate (from = 1, to = 3, value = 3) %>%
  T.rate (from = 3, to = 1, value = 4) %>%
  T.rate (from = 0, to = 4, value = 5) %>%
  T.rate (from = 4, to = 0, value = 6) %>%
  T.rate (from = 1, to = 2, value = 7) %>%
  T.rate (from = 2, to = 1, value = 8) %>%
  T.rate (from = 4, to = 5, value = 9) %>%
  T.rate (from = 5, to = 4, value = 10) %>%
  T.rate (from = 5, to = 6, value = 11) %>%
  T.rate (from = 6, to = 5, value = 12) %>%
  T.rate (from = 4, to = 6, value = 13) %>%
  T.rate (from = 6, to = 4, value = 14) 
rownames(q$paired_ind_ard) <- colnames(q$paired_ind_ard) <- states
q$paired_ind_ard

##Absent to single subgular with paired vocal sacs evolving individually.
  #A mix hypothesis of absent to single subgular and Paired independen. Vocal sac always appear and dissapear from a single
  # subgular state and paired vocal sacs form an independent evolutionary path.

#ER
q$ab.to.Ssub.Pind_er <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 2, value = 1) %>%
  T.rate (from = 2, to = 1, value = 1) %>%
  T.rate (from = 1, to = 3, value = 1) %>%
  T.rate (from = 3, to = 1, value = 1) %>%
  T.rate (from = 1, to = 4, value = 1) %>%
  T.rate (from = 4, to = 1, value = 1) %>%
  T.rate (from = 4, to = 5, value = 1) %>%
  T.rate (from = 5, to = 4, value = 1) %>%
  T.rate (from = 5, to = 6, value = 1) %>%
  T.rate (from = 6, to = 5, value = 1) %>%
  T.rate (from = 4, to = 6, value = 1) %>% 
  T.rate (from = 6, to = 4, value = 1)
rownames(q$ab.to.Ssub.Pind_er) <- colnames(q$ab.to.Ssub.Pind_er) <- states
q$ab.to.Ssub.Pind_er

#SYM 
q$ab.to.Ssub.Pind_sym <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 2, value = 2) %>%
  T.rate (from = 2, to = 1, value = 2) %>%
  T.rate (from = 1, to = 3, value = 3) %>%
  T.rate (from = 3, to = 1, value = 3) %>%
  T.rate (from = 1, to = 4, value = 4) %>%
  T.rate (from = 4, to = 1, value = 4) %>%
  T.rate (from = 4, to = 5, value = 5) %>%
  T.rate (from = 5, to = 4, value = 5) %>%
  T.rate (from = 5, to = 6, value = 6) %>%
  T.rate (from = 6, to = 5, value = 6) %>%
  T.rate (from = 4, to = 6, value = 7) %>% 
  T.rate (from = 6, to = 4, value = 7)
rownames(q$ab.to.Ssub.Pind_sym) <- colnames(q$ab.to.Ssub.Pind_sym) <- states
q$ab.to.Ssub.Pind_sym

#ARD
q$ab.to.Ssub.Pind_ard <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 2, value = 2) %>%
  T.rate (from = 2, to = 1, value = 3) %>%
  T.rate (from = 1, to = 3, value = 4) %>%
  T.rate (from = 3, to = 1, value = 5) %>%
  T.rate (from = 1, to = 4, value = 6) %>%
  T.rate (from = 4, to = 1, value = 7) %>%
  T.rate (from = 4, to = 5, value = 8) %>%
  T.rate (from = 5, to = 4, value = 9) %>%
  T.rate (from = 5, to = 6, value = 10) %>%
  T.rate (from = 6, to = 5, value = 11) %>%
  T.rate (from = 4, to = 6, value = 12) %>% 
  T.rate (from = 6, to = 4, value = 13)
rownames(q$ab.to.Ssub.Pind_ard) <- colnames(q$ab.to.Ssub.Pind_ard) <- states
q$ab.to.Ssub.Pind_ard

## Origin of triple vocal sac hypothesis
 #triple from bilobate

#ER
q$triple.from.bilobate_er <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 3, value = 1) %>%
  T.rate (from = 3, to = 1, value = 1) %>%
  T.rate (from = 2, to = 3, value = 1) %>%
  T.rate (from = 3, to = 2, value = 1) %>%
  T.rate (from = 1, to = 4, value = 1) %>%
  T.rate (from = 4, to = 1, value = 1) %>%
  T.rate (from = 1, to = 5, value = 1) %>%
  T.rate (from = 5, to = 1, value = 1) %>%
  T.rate (from = 1, to = 6, value = 1) %>%
  T.rate (from = 6, to = 1, value = 1) 
rownames(q$triple.from.bilobate_er) <- colnames(q$triple.from.bilobate_er) <- states
q$triple.from.bilobate_er

#SYM
q$triple.from.bilobate_sym <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 3, value = 2) %>%
  T.rate (from = 3, to = 1, value = 2) %>%
  T.rate (from = 2, to = 3, value = 3) %>%
  T.rate (from = 3, to = 2, value = 3) %>%
  T.rate (from = 1, to = 4, value = 4) %>%
  T.rate (from = 4, to = 1, value = 4) %>%
  T.rate (from = 1, to = 5, value = 5) %>%
  T.rate (from = 5, to = 1, value = 5) %>%
  T.rate (from = 1, to = 6, value = 6) %>%
  T.rate (from = 6, to = 1, value = 6) 
rownames(q$triple.from.bilobate_sym) <- colnames(q$triple.from.bilobate_sym) <- states
q$triple.from.bilobate_sym

#ARD
q$triple.from.bilobate_ard <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 3, value = 2) %>%
  T.rate (from = 3, to = 1, value = 3) %>%
  T.rate (from = 2, to = 3, value = 4) %>%
  T.rate (from = 3, to = 2, value = 5) %>%
  T.rate (from = 1, to = 4, value = 6) %>%
  T.rate (from = 4, to = 1, value = 7) %>%
  T.rate (from = 1, to = 5, value = 8) %>%
  T.rate (from = 5, to = 1, value = 9) %>%
  T.rate (from = 1, to = 6, value = 10) %>%
  T.rate (from = 6, to = 1, value = 11) 
rownames(q$triple.from.bilobate_ard) <- colnames(q$triple.from.bilobate_ard) <- states
q$triple.from.bilobate_ard

#Origin of paired subgular 
#Paired subgular from bilobate

#ER
q$Psub.from.bilobate_er <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 2, value = 1) %>%
  T.rate (from = 2, to = 1, value = 1) %>%
  T.rate (from = 1, to = 3, value = 1) %>%
  T.rate (from = 3, to = 1, value = 1) %>%
  T.rate (from = 3, to = 4, value = 1) %>%
  T.rate (from = 4, to = 3, value = 1) %>%
  T.rate (from = 1, to = 5, value = 1) %>%
  T.rate (from = 5, to = 1, value = 1) %>%
  T.rate (from = 1, to = 6, value = 1) %>%
  T.rate (from = 6, to = 1, value = 1)
rownames(q$Psub.from.bilobate_er) <- colnames(q$Psub.from.bilobate_er) <- states
q$Psub.from.bilobate_er

#SYM
q$Psub.from.bilobate_sym <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 1) %>%
  T.rate (from = 1, to = 2, value = 2) %>%
  T.rate (from = 2, to = 1, value = 2) %>%
  T.rate (from = 1, to = 3, value = 3) %>%
  T.rate (from = 3, to = 1, value = 3) %>%
  T.rate (from = 3, to = 4, value = 4) %>%
  T.rate (from = 4, to = 3, value = 4) %>%
  T.rate (from = 1, to = 5, value = 5) %>%
  T.rate (from = 5, to = 1, value = 5) %>%
  T.rate (from = 1, to = 6, value = 6) %>%
  T.rate (from = 6, to = 1, value = 6)
rownames(q$Psub.from.bilobate_sym) <- colnames(q$Psub.from.bilobate_sym) <- states
q$Psub.from.bilobate_sym

#ARD
q$Psub.from.bilobate_ard <- matrix (0, nrow = n, ncol = n) %>%
  T.rate (from = 0 , to = 1 , value = 1) %>%
  T.rate (from = 1, to = 0, value = 2) %>%
  T.rate (from = 1, to = 2, value = 3) %>%
  T.rate (from = 2, to = 1, value = 4) %>%
  T.rate (from = 1, to = 3, value = 5) %>%
  T.rate (from = 3, to = 1, value = 6) %>%
  T.rate (from = 3, to = 4, value = 7) %>%
  T.rate (from = 4, to = 3, value = 8) %>%
  T.rate (from = 1, to = 5, value = 9) %>%
  T.rate (from = 5, to = 1, value = 10) %>%
  T.rate (from = 1, to = 6, value = 11) %>%
  T.rate (from = 6, to = 1, value = 12)
rownames(q$Psub.from.bilobate_ard) <- colnames(q$Psub.from.bilobate_ard) <- states
q$Psub.from.bilobate_ard

##Changing colnames in all q matrix to be the same

q <- lapply(q, function(m) {
  if (is.matrix(m)) {
    rownames(m) <- colnames(m) <- states
  }
  return(m)
})

save.image ("Seven_states_models_fitMk.RData")

##3 FIT MODELS USING FITMK()------------------------------------------

#define root state probability
fixed_root<- setNames(c(1, rep(0, ncol(trait) - 1)), colnames(trait))

roots<-list()
for(i in 1:length(q)){
  if(names(q)[i] %in% c("er","sym","ard")){
    roots[[i]]="estimated"
  }else roots[[i]]=fixed_root
}
roots


fit<- list()

for(i in 1:length(q)){
  print(paste("Running model", i , "of", length(q)))
  
  fit[[names(q)[i]]] <-fitMk(
    tree   = clean_tree,
    x      = trait,
    type   = "discrete",
    model  = q[[i]],
    pi     = roots[[i]],
    ncores = 20
  )
  
}

#4 COMPARE MODELS WITH AIC TEST------------------------------

## we will use AIC to compare the model fit
aic_table<-
  lapply(X=fit, FUN=AIC) %>% unlist() %>%
  enframe(name = "model", value = "AIC") %>%
  mutate(
    logLik=lapply(X=fit, FUN=logLik) %>% unlist()
  ) %>%
  arrange(AIC)

aic_table

#5. PLOT TRANSITIONS

#Transition plots
par(mfrow=c(4,4))
for(i in 1:length(fit)){
  plot(fit[[i]], color=TRUE,lwd=3, main=names(fit)[i])
}
