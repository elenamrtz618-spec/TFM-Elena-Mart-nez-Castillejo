## Elena Martínez Castillejo & H. Christoph Liedtke

##HYPOTHESIS TESTING FOR VOCAL SAC EVOLUTION.
  ##Continuous time Markov models for discrete character evolution --- fitMk()
   ##THREE STATES - ABSENT/ SINGLE VS/ PAIRED VS

##Library--------------------------------------------------------------

library(ape)
library(tidyverse)
library(readxl)
library(qgraph)
library(phytools)

##1- DATA LOADING AND CLEANING---------------------------------------------------------

 #Load data
 tree <- read.tree("C:/Users/elena/Nextcloud/Elena/data/clean_data/portik_phy_asw2026.tre")
 data <- read_xlsx ("vocal_sac_coding.xlsx", sheet= "coding") 
 
 data$seven_states_code %>% unique()
 
 simple.trait.df <- data %>%
   distinct(Species_asw, .keep_all = TRUE) %>%
   filter(Species_asw %in% tree$tip.label) %>%
   arrange(match(Species_asw, tree$tip.label)) %>%
   rowwise() %>%
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
   filter(!is.na(seven_states_code)) %>%
   ungroup()
 
 unk <- which(simple.trait.df$seven_states_code == "?")
 

 ##Making a matrix for vocal sac state probabilities
 
 states <- as.character (0:6)
 
 Mat_20 <- do.call(rbind, lapply(simple.trait.df$state_probs,function(p) {
   res <- setNames(rep(0, length(states)), states)
   res[names(p)] <- p
   return(res)
 }))   
 
 trait <- matrix (0,
                   nrow = nrow(Mat_20),
                   ncol = 3,
                   dimnames = list(simple.trait.df$Species_asw, c("0", "1", "2")))
 
 #Sum probabilities
 trait [,"0"] <- Mat_20[,"0"]
 trait [,"1"] <- rowSums(Mat_20[,c("1","2","3"), drop = FALSE])
 trait [,"2"] <- rowSums(Mat_20[,c("4","5","6"), drop = FALSE])
 trait [unk,] <- 1/3
 colSums(trait)
 

 #----------------------------------------------
 states <-as.character(0:2)
 n <- length(states)
 
 ## prune tree to match data
 clean_tree<-keep.tip(tree, rownames(trait))
 clean_tree
 
 all(rownames(trait)==clean_tree$tip.label)
 
 # DEFINING TRANSITION RATE MATRIX FOR ALL MODELS-------------
 q<-list()
 
 #Making a function to build matrix easily
 T.rate <- function (matrix, from, to, value) {
   matrix[from +1, to +1] <- value
   return(matrix) 
 }
 
 #2.1 FULL MODELS-------------------------------------------------
 
 ##ARD
  q$ard<-matrix(nrow=n, ncol=n,0) %>% # set all transitions to be 0
    T.rate(from = 0, to = 1, value = 1) %>% 
    T.rate(from = 1, to = 0, value = 2) %>% 
    T.rate(from = 0, to = 2, value = 3) %>% 
    T.rate(from = 2, to = 0, value = 4) %>%
    T.rate(from = 1, to = 2, value = 5) %>%
    T.rate(from = 2, to = 1, value = 6) 
  rownames(q$ard) <- colnames(q$ard) <- states
  q$ard

  ##SYM model
 
 q$sym<-matrix(nrow=n, ncol=n,0) %>% # set all transitions to be 0
   T.rate(from = 0, to = 1, value = 1) %>% 
   T.rate(from = 1, to = 0, value = 1) %>% 
   T.rate(from = 0, to = 2, value = 2) %>% 
   T.rate(from = 2, to = 0, value = 2) %>%
   T.rate(from = 1, to = 2, value = 3) %>%
   T.rate(from = 2, to = 1, value = 3) 
 rownames(q$sym) <- colnames(q$sym) <- states
 q$sym

 ## ER model
 q$er<-matrix(nrow=n, ncol=n,1) # set all transitions to be the same (index=1)
  diag(q$er) <- 0 # make diagonals 0
  rownames(q$er) <- colnames(q$er) <- states
  q$er 

# 2.2 RADIAL MODELS------------------------------------
  ##ER Radial
  q$er_radial <- matrix(0, nrow = n, ncol = n) %>%
     T.rate(from = 0, to = 1, value = 1) %>% 
     T.rate(from = 1, to = 0, value = 1) %>% 
     T.rate(from = 0, to = 2, value = 1) %>% 
     T.rate(from = 2, to = 0, value = 1) 
   rownames(q$er_radial) <- colnames(q$er_radial) <- states
   q$er_radial  

   # ARD_Radial
 q$ard_radial <- matrix(0, nrow = n, ncol = n) %>%
   T.rate(from = 0, to = 1, value = 1) %>% 
   T.rate(from = 1, to = 0, value = 2) %>% 
   T.rate(from = 0, to = 2, value = 3) %>% 
   T.rate(from = 2, to = 0, value = 4) 
 rownames(q$ard_radial) <- colnames(q$ard_radial) <- states
 q$ard_radial   

 #sym_radial
 q$sym_radial <- matrix(0, nrow = n, ncol = n) %>%
   T.rate(from = 0, to = 1, value = 1) %>% 
   T.rate(from = 1, to = 0, value = 1) %>% 
   T.rate(from = 0, to = 2, value = 2) %>% 
   T.rate(from = 2, to = 0, value = 2) 
 rownames(q$sym_radial) <- colnames(q$sym_radial) <- states
 q$sym_radial 

 #2.3 HYPOTHESIS BASED MODELS-------------------------------------------------
 
 ##Secuential from simple to paired (allowing reverse transitions)
 
#ER
 q$seq_er <- matrix(0, nrow = n, ncol = n) %>%
   T.rate(from = 0, to = 1, value = 1) %>%
   T.rate(from = 1, to = 0, value = 1) %>% 
   T.rate(from = 1, to = 2, value = 1) %>% 
   T.rate(from = 2, to = 1, value = 1) 
 rownames(q$seq_er) <- colnames(q$seq_er) <- states
 q$seq_er
 
#SYM 
  q$seq_sym <- matrix(0, nrow = n, ncol = n) %>%
   T.rate(from = 0, to = 1, value = 1) %>%
   T.rate(from = 1, to = 0, value = 1) %>% 
   T.rate(from = 1, to = 2, value = 2) %>% 
   T.rate(from = 2, to = 1, value = 2) 
  rownames(q$seq_sym) <- colnames(q$seq_sym) <- states
  q$seq_sym

#ARD     
  q$seq_ard <- matrix(0, nrow = n, ncol = n) %>%
    T.rate(from = 0, to = 1, value = 1) %>%
    T.rate(from = 1, to = 0, value = 2) %>% 
    T.rate(from = 1, to = 2, value = 3) %>% 
    T.rate(from = 2, to = 1, value = 4) 
  rownames(q$seq_ard) <- colnames(q$seq_ard) <- states
  q$seq_ard

##Sequential from simple to paired but can only go back to absent

#ER  
  q$seq.only.loss_er <- matrix(0, nrow = n, ncol = n) %>%
    T.rate(from = 0, to = 1, value = 1) %>%
    T.rate(from = 1, to = 0, value = 1) %>% 
    T.rate(from = 1, to = 2, value = 1) %>% 
    T.rate(from = 2, to = 0, value = 1) 
  rownames(q$seq_er) <- colnames(q$seq.only.loss_er) <- states
  q$seq.only.loss_er
  
#SYM
  q$seq.only.loss_sym <- matrix(0, nrow = n, ncol = n) %>%
    T.rate(from = 0, to = 1, value = 1) %>%
    T.rate(from = 1, to = 0, value = 1) %>% 
    T.rate(from = 1, to = 2, value = 2) %>% 
    T.rate(from = 2, to = 0, value = 2) 
  rownames(q$seq.only.loss_sym) <- colnames(q$seq.only.loss_sym) <- states
  q$seq.only.loss_sym
  
#ARD
  q$seq.only.loss_ard <- matrix(0, nrow = n, ncol = n) %>%
    T.rate(from = 0, to = 1, value = 1) %>%
    T.rate(from = 1, to = 0, value = 2) %>% 
    T.rate(from = 1, to = 2, value = 3) %>% 
    T.rate(from = 2, to = 0, value = 4) 
  rownames(q$seq.only.loss_ard) <- colnames(q$seq.only.loss_ard) <- states
  q$seq.only.loss_ard
  
##Changing colnames in all q matrix to be the same  
  q <- lapply(q, function(m) {
    if (is.matrix(m)) {
      rownames(m) <- colnames(m) <- states
    }
    return(m)
  })
  
  
save.image("Three_states_SinglePaired.RData")

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

save.image("Three_states_SinglePaired.RData")

##FitMk

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


