library(phytools)
library(ape)
library(caper)
library(tidyverse)
setwd('')

########## 1. comparing method against phylogenetic PCA (this can only be done at the species level) #########

ana = read.csv('B - Species level data.csv') 
ana = ana[,c('Species', 'Peak_frequency_SD', 'Bandwidth_SD', 'Note_slope_SD',
             'Frequency_variation')]

ana$Peak_frequency_SD = scale(ana$Peak_frequency_SD)
ana$Bandwidth_SD = scale(ana$Bandwidth_SD)
ana$Note_slope_SD = scale(ana$Note_slope_SD)
ana = na.omit(ana)

aostree <- read.tree('T400F_AOS_Clements_sppnames.tre')
ana$Species = gsub(' ', '_', ana$Species)
anatree <- comparative.data(phy = aostree, data = ana,
                            names.col = 'Species', 
                            vcv = TRUE,
                            warn.dropped = TRUE,
                            na.omit = FALSE)
ppca <- phyl.pca(anatree$phy,
                 anatree$data[, -4],
                 mode = "corr")

# check pPCA results
summary(ppca) # pPC1 alone accounts for ~70% variation
ppca$Evec[, 1] # each variable has almost identical contribution to PC1, all negatively

# compare results
anatree$data$pPC1 = ppca$S[,1] # add pPC1 to data
cor.test(anatree$data$pPC1, anatree$data$Frequency_variation)   # correlation is nearly -1.


# Figure S1A
figs1a = ggplot(anatree$data, aes(y=Frequency_variation, x=pPC1)) + theme_classic() +
    geom_point(size = 1, alpha = 0.5) +
    theme(axis.title = element_text(size=18), axis.text = element_text(size=18)) +
    labs(y='Frequency variation', x = 'Phylogenetic PC1 (69.1%)') +
    annotate('text', label=paste0("Pearson's correlation = -1.00",
                                  '\np < 0.001 \nn = ', nrow(anatree$data), ' species'), 
             hjust=0, x = -30.5, y = 0, size = 5.5)



########## 2. comparing method against normal PCA - species level #############

pca <- prcomp(ana[,2:4], center = TRUE, scale. = TRUE)

# check PCA results
summary(pca) # PC1 alone accounts for ~75% variation.
pca$rotation[,1]  # each variable has almost identical contribution to PC1, all negatively.

# compare results
ana$pc1 = pca$x[,1] # add PC1 to data
cor.test(ana$pc1, ana$Frequency_variation)   # correlation is nearly -1.

figs1b = ggplot(ana, aes(y=Frequency_variation, x=pc1)) + theme_classic() +
    geom_point(size = 1, alpha = 0.5) +
    theme(axis.title = element_text(size=18), axis.text = element_text(size=18),
          axis.title.y = element_blank()) +
    labs(x = 'PC1 (74.9%)') +
    annotate('text', label=paste0("Pearson's correlation = -1.00",
                                  '\np < 0.001 \nn = ', nrow(anatree$data), ' species'), 
             hjust=0, x = -12, y = 0, size = 5.5)



########## 2. comparing method against normal PCA - song level #############

song <- read.csv('C - Song level data.csv')

# select relevant variales
song = song[, c('Peak_frequency_SD', 'Bandwidth_SD', 'Note_slope_SD',
                'Frequency_variation')]

song$Peak_frequency_SD = scale(song$Peak_frequency_SD)
song$Bandwidth_SD = scale(song$Bandwidth_SD)
song$Note_slope_SD = scale(song$Note_slope_SD)
song = na.omit(song)

pca <- prcomp(song[,-4], center = TRUE, scale. = TRUE)

# check PCA results
summary(pca) # PC1 alone accounts for >70% variation.
pca$rotation[,1]  # each variable has almost identical contribution to PC1, all negatively.

# compare results
song$pc1 = pca$x[,1] # add PC1 to data
cor.test(song$pc1, song$Frequency_variation)   # correlation is nearly -1.

figs1c = ggplot(song, aes(y=Frequency_variation, x=pc1)) + theme_classic() +
    geom_point(size = 1, alpha = 0.5) +
    theme(axis.title = element_text(size=18), axis.text = element_text(size=18),
          axis.title.y = element_blank()) +
    labs(x = 'PC1 (71.2%)') +
    annotate('text', label=paste0("Pearson's correlation = -1.00",
                                  '\np < 0.001 \nn = ', nrow(song), ' songs'), 
             hjust=0, x = -13.8, y = 0.4, size = 5.5)


### Figure S1. Comparing frequency variation against PCA methods
cowplot::plot_grid(figs1a, figs1b, figs1c, nrow = 1, rel_widths = c(1.05, 1, 1))
ggsave('Figures/figure S1.svg', width = 15, height = 5)


