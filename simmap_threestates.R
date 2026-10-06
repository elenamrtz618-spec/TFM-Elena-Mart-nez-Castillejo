
##Plot tree for three states models

##Library-----------------------------------------------

library(phytools)
library(ggtree)
library(tidyverse)
library(ape)

##Load model data------------------------------------------

setwd(this.path::here())
load("Three_states_SinglePaired.RData") 

##Stochastic character mapping----------------------------

fitMk_simmap <- simmap(fit$ard,nsim=500)

##Summarize all simulations

fitMk_simmap_desc <- describe.simmap(fitMk_simmap, plot=F)

save.image("singlepaired_simmaploaded.RData")

anc_states<-fitMk_simmap_desc$ace
head(anc_states)

##Plot node pies

source("C:/Users/elena/Nextcloud/Elena/Anc.Rec and Stoch.Mapping/tree_making/find_transitions.R")

# extract posterior probabilities of nodes
anc_states[1,]
node_states_max<-colnames(anc_states)[apply(X = anc_states, M=1, FUN = which.max)]
head(node_states_max)

#EXTRACT TIP STATES
#Max probability extracted from data matrix
max_prob_obs <- apply(trait,1,max)

#Extract max state for tip state observed from trait data
tip_states_obs <-colnames(trait)[apply(trait,1,which.max)]

#Extract max state predicted by simmap posterior probabilities
tip_states_simmap <- colnames(fitMk_simmap_desc$tips)[apply(fitMk_simmap_desc$tips,1,which.max)]

#Make a hybrid vector for tip states. Simmap prediction for uncertain states and max prob observed for known states
tip_states <- ifelse(max_prob_obs < 1, tip_states_simmap, tip_states_obs)

tip_states <- tip_states[clean_tree$tip.label]



# make a simmap with painted branches according to the maximum posterior probabilities.
phy_painted<-find_transitions(phy=clean_tree,
                              tip_states = tip_states,
                              node_states = node_states_max,
                              simmap = T,
                              stem_prop = 0.001)

save.image("singlepaired_simmaploaded.RData")

##PLOT---------------

#Define colors
cols_three <- viridisLite::viridis(3, option = "A")
cols_three [1] <- "grey60"
cols_three [2] <- "violetred1"
cols_three [3] <- "gold" 
names(cols_three)<-as.character(0:2)
cols_three

#plot tree

par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
plotSimmap(phy_painted, type= "fan", 
           ftype= "off",     
           lwd= 0.5,
           colors = cols_three)
tiplabels(pch=15, col=cols_three[tip_states[phy_painted$tip.label]],cex=0.25)

## PLOT ALL NODES---------------------------------------------------------------
nodelabels(pie=anc_states, piecol= cols_three, cex = 0.2)
#-------------------------------------------------------------------------------

##PLOT SHIFT NODES ONLY---------------------------------------------------------

state_transitions<-find_transitions(phy=clean_tree,
                                    tip_states = tip_states,
                                    node_states = node_states_max,
                                    simmap = F,
                                    stem_prop = 0.001)
state_transitions

# Find node numbers where states are different from the previous node and lets also exclude shifts occuring at tips 
shift_nodes=state_transitions %>%
  filter(shifts) %>%
  filter(end_node>Ntip(clean_tree))
shift_nodes

# plot 
nodelabels(pie=anc_states[rownames(anc_states) %in% shift_nodes$end_node,], 
           node = shift_nodes$end_node, piecol= cols_three, cex = 0.2)

#HEATMAP, FAMILY LABELS AND LEGEND----------------------------------------------

library(plotrix)

max_tree_depth <- max(nodeHeights(clean_tree))
lim_val        <- max_tree_depth * 1.45  # Expand limits so outer labels & legend fit

# Set device margins
par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), xpd = TRUE)

plotSimmap(
  phy_painted, 
  type = "fan", 
  ftype = "off", 
  lwd = 0.5, 
  colors = cols_three,
  xlim = c(-lim_val, lim_val), 
  ylim = c(-lim_val, lim_val)
)

# Plot shift node pies--------------------------------------------
nodelabels(
  pie = anc_states[rownames(anc_states) %in% shift_nodes$end_node, ], 
  node = shift_nodes$end_node, 
  piecol = cols_three, 
  cex = 0.2
)

#Plot all nodes-------------------------------------------
nodelabels(pie=anc_states, piecol= cols_three, cex = 0.2)

# ==============================================================================
# HEATMAP RING, FAMILY LABELS, AND LEGEND
# ==============================================================================
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
n_tips <- Ntip(clean_tree)

# Define radial distances relative to tree height
ring_inner_r <- max_tree_depth * 1.01
ring_outer_r <- max_tree_depth * 1.06
arc_r        <- max_tree_depth * 1.08
label_r      <- max_tree_depth * 1.12

# ------------------------------------------------------------------------------
# A. HEATMAP RING
# ------------------------------------------------------------------------------
  ##tip state vector dedicated to heatmap. Make uncertain states transparent
tip_states_heatmap <- tip_states_obs
tip_states_heatmap[max_prob_obs[clean_tree$tip.label] < 1] <- NA #Set uncertainities as NAs

for (i in 1:n_tips) {
  state_val <- tip_states_heatmap[i]
  if (is.na(state_val)) next  #this skips drawing for NAs (uncertai)
  
  color_val <- cols_three[as.character(state_val)]
  angle_rad <- atan2(lastPP$yy[i], lastPP$xx[i])
  
  lines(
    c(ring_inner_r * cos(angle_rad), ring_outer_r * cos(angle_rad)),
    c(ring_inner_r * sin(angle_rad), ring_outer_r * sin(angle_rad)),
    col = color_val,
    lwd = 1.8
  )
}

# ------------------------------------------------------------------------------
# B. FAMILY CLADE LABELS & ARCS (Continuous Arcs)
# ------------------------------------------------------------------------------
family_nodes <- simple.trait.df %>%
  filter(Species_asw %in% clean_tree$tip.label) %>%
  arrange(match(Species_asw, clean_tree$tip.label)) %>%
  filter(!is.na(Family_asw)) %>%
  group_by(Family_asw) %>%
  summarise(node = ifelse(n() < 2, 
                          which(clean_tree$tip.label == Species_asw), 
                          ape::getMRCA(clean_tree, Species_asw))) %>% 
  filter(!is.na(node))

for (i in seq_len(nrow(family_nodes))) {
  f_node  <- family_nodes$node[i]
  f_label <- family_nodes$Family_asw[i]
  
  if (f_node <= n_tips) {
    clade_tips <- f_node
  } else {
    clade_tips <- phytools::getDescendants(clean_tree, node = f_node)
    clade_tips <- clade_tips[clade_tips <= n_tips]
  }
  
  angles <- atan2(lastPP$yy[clade_tips], lastPP$xx[clade_tips])
  min_a  <- min(angles)
  max_a  <- max(angles)
  mid_a  <- (min_a + max_a) / 2
  
  # Draw continuous arc
  draw.arc(
    x = 0, y = 0, 
    radius = arc_r, 
    deg1 = min_a * (180 / pi), 
    deg2 = max_a * (180 / pi), 
    col = "black", 
    lwd = 1.2
  )
  
  # Text rotation & orientation adjustment
  deg_rot <- mid_a * (180 / pi)
  if (cos(mid_a) < 0) deg_rot <- deg_rot + 180 
  
  text(
    x = label_r * cos(mid_a),
    y = label_r * sin(mid_a),
    labels = f_label,
    srt = deg_rot,
    adj = ifelse(cos(mid_a) < 0, 1, 0),
    cex = 0.35
  )
}

