
#### Elena Martínez Castillejo & H. Christoph Liedtke

  ##HYPOTHESIS TESTING FOR VOCAL SAC EVOLUTION.
  ##Continuous time Markov models for discrete character evolution --- fitMk()
  ##Three states coding: absent/buccal floor expanded/present

##Library-------------------------------------------

library(ape)
library(tidyverse)
library(readxl)
library(qgraph)
library(phytools)

##1- DATA LOADING AND CLEANING---------------------------------------------------------

#load Data
setwd("C:/Users/elena/Nextcloud/Elena/data/clean_data")

tree <- read.tree("C:/Users/elena/Nextcloud/Elena/data/clean_data/portik_phy_asw2026.tre")
data <- read_xlsx ("../../data/clean_data/vocal_sac_coding.xlsx", sheet= "coding") 

#Cleaning
data <- data %>%
rename(twenty_states_code = `20 states-code`)
 data$twenty_states_code %>% unique
   
 simple.trait.df <- data %>%
  # remove duplicates
    distinct(Species_asw, .keep_all = TRUE) %>%
  # filter and arrange data to match tree
    filter(Species_asw %in% tree$tip.label) %>%
    arrange(match(Species_asw, tree$tip.label)) %>%
   rowwise() %>%
  # recode uncertainties
     mutate(
       state_probs = case_when(
     twenty_states_code == "?" ~ list(setNames(rep(1/20,20), 0:19)),
     twenty_states_code == "[0 12]" ~ list(c("0" = 0.5, "12" = 0.5)),
     twenty_states_code == "[0 15]" ~ list(c("0" = 0.5, "15" = 0.5)),
     twenty_states_code == "[0 3]" ~  list(c("0" = 0.5, "3" = 0.5)),
     twenty_states_code == "[0 9]" ~  list(c("0" = 0.5, "9" = 0.5)),
     twenty_states_code == "[1 2]" ~  list(c("1" = 0.5, "2" = 0.5)),
     twenty_states_code == "[3 11]" ~ list(c("3" = 0.5, "11" = 0.5)),
     twenty_states_code == "[3 4]" ~  list(c("3" = 0.5, "4" = 0.5)),
     twenty_states_code == "[3 6]" ~  list(c("3" = 0.5, "6" = 0.5)),
     twenty_states_code == "[3 9]" ~  list(c("3" = 0.5, "9" = 0.5)),
     twenty_states_code == "[9 12]" ~ list(c("9" = 0.5, "12" = 0.5)),
     TRUE ~ lapply(twenty_states_code, function(x) setNames(1.0, x))
       )) %>%
   filter(!is.na(seven_states_code)) %>%
  ungroup()

 unk <- which(simple.trait.df$twenty_states_code == "?")


 ##Making a matrix for vocal sac state probabilities
 
twenty_states <- as.character (0:19)

Mat_20 <- do.call(rbind, lapply(simple.trait.df$state_probs,function(p) {
  res <- setNames(rep(0, 20), twenty_states)
  res[names(p)] <- p
  return(res)
}))   

trait <- matrix (0,
                 nrow = nrow(Mat_20),
                 ncol = 3,
                 dimnames = list(simple.trait.df$Species_asw, c("0", "1", "2")))

#Sum probabilities
 trait [,"0"] <- Mat_20[,"0"]
 trait [,"1"] <- rowSums(Mat_20[,c("1", "2"), drop = FALSE])
 trait [,"2"] <- rowSums(Mat_20[, as.character(3:19), drop = FALSE])
 trait [unk,] <- 1/3
head(trait) 
all(rowSums(trait) == 1)  


## prune tree to match data
 clean_tree<-keep.tip(tree, rownames(trait))
 clean_tree
 
all(rownames(trait)==clean_tree$tip.label)

states <- (0:2)
n <- length(states)

# DEFINING TRANSITION RATE MATRIX FOR ALL MODELS-------------

q <- list()

T.rate <- function (matrix, from, to, value) {
 matrix[from +1, to +1] <- value
 return(matrix)
 }

#2.1 FULL MODELS-------------------------------------------------

# ER
 q$er<-matrix(nrow=n, ncol=n,1) 
 diag(q$er) <- 0 
 rownames(q$er) <- colnames(q$er) <- states
 q$er

 #SYM
 q$sym<-matrix(nrow=n, ncol=n,0) %>% 
   T.rate(from = 0, to = 1, value = 1) %>% 
   T.rate(from = 1, to = 0, value = 1) %>% 
   T.rate(from = 0, to = 2, value = 2) %>% 
   T.rate(from = 2, to = 0, value = 2) %>%
   T.rate(from = 1, to = 2, value = 3) %>%
   T.rate(from = 2, to = 1, value = 3) 
  rownames(q$sym) <- colnames(q$sym) <- states
  q$sym

  #ARD
  q$ard<-matrix(nrow=n, ncol=n,0) %>% 
    T.rate(from = 0, to = 1, value = 1) %>% 
    T.rate(from = 1, to = 0, value = 2) %>% 
    T.rate(from = 0, to = 2, value = 3) %>% 
    T.rate(from = 2, to = 0, value = 4) %>%
    T.rate(from = 1, to = 2, value = 5) %>%
    T.rate(from = 2, to = 1, value = 6) 
  rownames(q$ard) <- colnames(q$ard) <- states
  q$ard
  
#2.2 RADIAL MODELS---------------------------------------------------------
 
  #ER
  q$er_radial <- matrix(0, nrow = n, ncol = n) %>%
     T.rate(from = 0, to = 1, value = 1) %>% 
     T.rate(from = 1, to = 0, value = 1) %>% 
     T.rate(from = 0, to = 2, value = 1) %>% 
     T.rate(from = 2, to = 0, value = 1) 
   rownames(q$er_radial) <- colnames(q$er_radial) <- states
   q$er_radial 

   #SYM
   q$sym_radial <- matrix(0, nrow = n, ncol = n) %>%
     T.rate(from = 0, to = 1, value = 1) %>% 
     T.rate(from = 1, to = 0, value = 1) %>% 
     T.rate(from = 0, to = 2, value = 2) %>% 
     T.rate(from = 2, to = 0, value = 2) 
   rownames(q$sym_radial) <- colnames(q$sym_radial) <- states
   q$sym_radial 

   # ARD
   q$ard_radial <- matrix(0, nrow = n, ncol = n) %>%
     T.rate(from = 0, to = 1, value = 1) %>% 
     T.rate(from = 1, to = 0, value = 2) %>% 
     T.rate(from = 0, to = 2, value = 3) %>% 
     T.rate(from = 2, to = 0, value = 4) 
   rownames(q$ard_radial) <- colnames(q$ard_radial) <- states
   q$ard_radial   

#2.3 HYPOTHESIS BASED MODELS-------------------------------------------------
   
 #Serial hypothesis
   #Where we test if state 1 (absent vocal sac but with buccal floor expansion) is an intermidiate state between absence and presence.
   
   q$serial_ard <- matrix(0, nrow = n, ncol = n) %>%
     T.rate(from = 0, to = 1, value = 1) %>% 
     T.rate(from = 1, to = 0, value = 2) %>% 
     T.rate(from = 1, to = 2, value = 3) %>% 
     T.rate(from = 2, to = 1, value = 4) 
   rownames(q$serial_ard) <- colnames(q$serial_ard) <- states
   q$serial_ard

 ##False loss hypothesis
   #Where we want to see wether vocal sacs are lost completely or not
   
#Case 1: Once vocal sacs are lost they can't be regained
   
  q$false.loss <- matrix(0, nrow = n, ncol = n) %>%
    T.rate(from = 0, to = 2, value = 1) %>% 
    T.rate(from = 2, to = 1, value = 2) %>%
    T.rate(from = 0, to = 1, value = 3) 
   rownames(q$false.loss) <- colnames(q$false.loss) <- states
   q$false.loss
   
#Case 2: Regaining of presence is allowed
   
   q$false.loss.2 <- matrix(0, nrow = n, ncol = n) %>%
     T.rate(from = 0, to = 2, value = 1) %>% 
     T.rate(from = 2, to = 1, value = 2) %>%
     T.rate(from = 0, to = 1, value = 3) %>%
     T.rate(from = 1, to = 2, value = 4) 
   rownames(q$false.loss.2) <- colnames(q$false.loss.2) <- states
   q$false.loss.2
   
save.image("Three_states_LossVS.RData")   

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

fit <- list()   

save.image("Three_states_LossVS.RData")  

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


