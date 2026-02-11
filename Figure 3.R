setwd("")
library(tidyverse)
library(ape)
library(caper)
library(phytools)
library(ggtree)
library(ggtreeExtra)

# read the phylogeny 
aostree <- read.tree('T400F_AOS_Clements_sppnames.tre')
# read the species data
ana=read.csv('B-Species level data.csv')
# align species names in the two files
ana$Species <- gsub(' ', '_', ana$Species)

# add node number to the dataset
ana$node <- nodeid(aostree, ana$Species)
treedat <-  full_join(aostree, ana, by = "node")
    
      
# figure 3a
ggtree(treedat, layout = 'fan', size = 0.5, aes(col = Sexual_selection_score)) +
              scale_color_gradient(low = '#1d4b96', high = 'orange') +   
              labs(col = 'Sexual Selection') +
              theme(legend.position="bottom",
                    legend.box ="horizontal") +
              guides(col = guide_colourbar(title.position="top", title.hjust = 0.5)) +
              geom_fruit(geom = geom_bar,
                         stat = 'identity', fill = 'black',
                         mapping = aes(y=node,  x=Number_of_note_types), 
                         width = 1.2, pwidth = 0.5, 
                         offset = 0.05) 

# figure 3b
ggtree(treedat, layout = 'fan', size = 0.5, aes(col = Territoriality)) +
              scale_color_gradient(low = '#1d4b96', high = 'orange',
                                   breaks = c(0,1,2),
                                   labels = c(0,1,2)) +
              labs(col = 'Territoriality') +
              theme(legend.position="bottom",
                    legend.box ="horizontal") +
              guides(col = guide_colourbar(title.position="top", title.hjust = 0.5)) +
              geom_fruit(geom = geom_bar,
                        stat = 'identity', fill = 'black', 
                        mapping = aes(y=node, x=sqrt(Note_count)),
                        width = 1.8, pwidth = 0.6,
                        offset = 0.05)


# test for the strength of phylogenetic signals
anatree <- comparative.data(phy = aostree, data = ana,
                           names.col = 'Species', 
                           vcv = TRUE,
                           warn.dropped = TRUE,
                           na.omit = FALSE)

phylosig(anatree$phy, anatree$data$Note_count, method = "lambda")  # 0.40
phylosig(anatree$phy, anatree$data$Number_of_note_types, method = "lambda")  # 0.56

