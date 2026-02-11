setwd('')

library(ggdist)
library(tidyverse)
library(brms)

# add y-axis labels
ylabs <- c('Species richness', 'Habitat density', 'Beak size', 'Mass',
           'Social group size', 'Territoriality', 'Sexual selection score')

# function for summarising and visualising model results
distplot <- function(brms_model, sigCI=0.95, LL, RR, cus.title=brms_model$formula$resp, xlabs = '', ylabels = ylabs){
  
    # extract all predictor estimates
    plotdata <- data.frame()
    for (i in 1:length(brms_model[["fit"]]@sim[["samples"]])){
        plotdata1 <- as.data.frame(brms_model[["fit"]]@sim[["samples"]][[i]]) 
        plotdata <- rbind(plotdata, plotdata1) 
    }
    rep = nrow(plotdata)
    plotdata <- plotdata %>% dplyr::select(all_of(starts_with(c('b_', 'bsp_')))) %>% t() %>% as.data.frame()
    
    # calculate the 95 CI for each predictor
    plotdata$lowerCI <- apply(plotdata[,1:rep], 1, quantile, (1-sigCI)/2)
    plotdata$upperCI <- apply(plotdata[,1:rep], 1, quantile, (1+sigCI)/2)
    # use paler colours if the 95CI spans zero
    plotdata$shade <- ifelse(plotdata$lowerCI * plotdata$upperCI > 0, 1, 0.4)
    # re-order predictors
    plotdata$Predictors <- row.names(plotdata) %>% 
      factor(levels = c('b_Species.richness', 'b_Habitat_density', 'b_relative.beak.size', 'b_Mass',
                        'b_Sociality', 'b_territorialityterritorial', 'bsp_moSexual_selection_score', 
                        'b_ave.freq', 'b_ave.dur', 'b_Intercept'))
    
    # remove irrelevant predictors
    plotdata <- filter(plotdata, ! Predictors %in% c('b_Intercept', 'b_ave.freq', 'b_ave.dur'))
    
    plotdata <- pivot_longer(plotdata, cols = 1:rep, names_to = 'run')
    ggplot(plotdata, aes(y = Predictors, x = value, fill = Predictors, alpha=I(shade)) )+
                  stat_slab(scale=0.9)+
                  theme_classic() +
                  theme(legend.position = 'none') +
                  geom_vline(xintercept = 0, col='grey70') +
      ggtitle(cus.title) + labs(x=xlabs) +
      scale_y_discrete(labels = ylabels) +
      scale_x_continuous(limits = c(LL, RR)) +
      scale_fill_manual(values = c(rep('#293f8f',4),
                                 rep('#c25f19',3))) +
      theme(plot.title = element_text(size=16, hjust = 0.5),
            axis.title.y = element_blank(), axis.title.x = element_text(size=13),
            axis.text.y = element_text(size=14), axis.text.x = element_text(size=12))
     
}

# Main results (Figure 4; Table S2 top rows)
nc <- readRDS("NoteCountrd negbinomial __( main models , species level).rds")
nt <- readRDS("NumberofNoteTypesrd negbinomial __( main models , species level).rds")
sf <- readRDS("sdfreq gaussian __( main models , species level).rds")

p1=distplot(nc,  0.95, LL=-0.7, RR=1.6,  'Note count      \n')
p2=distplot(nt,  0.95, LL=-0.75, RR=0.5, '  Note type\n', ylabels = NULL)
p3=distplot(sf,  0.95, LL=-0.7, RR=0.7,  '   Frequency variation\n', ylabels = NULL)

# Figure 4
pdf('Figure 4.pdf', width = 9.1, height = 4.5)
cowplot::plot_grid(p1,p2,p3,ncol = 3, rel_widths = c(1.8, 1, 1), scale = 1) 
dev.off()



# SI: high certainty data (Table S2 middle rows)
nchc <- readRDS("NoteCountrd negbinomial __( high certainty data , species level).rds")
nthc <- readRDS("NumberofNoteTypesrd negbinomial __( high certainty data , species level).rds")
sfhc <- readRDS("sdfreq gaussian __( high certainty data , species level).rds")

# SI: removing weakly/seasonally territorial species (Table S2 bottom rows)
nct <- readRDS('NoteCountrd negbinomial __( strong vs none territoriality , species level).rds')
ntt <- readRDS('NumberofNoteTypesrd negbinomial __( strong vs none territoriality , species level).rds')
sft <- readRDS('sdfreq gaussian __( strong vs none territoriality , species level).rds')


# function to summarise and format model results
tidy.table = function(mod1, mod2, mod3){
  # top rows show the main model results using all species data
  t1 = summary(mod1)$fixed %>% mutate(mod = 'all species')
  t1$Predictors = row.names(t1)
  # middle rows show the supplementary model results using high-certainty data
  t2 = summary(mod2)$fixed %>% mutate(mod = 'high certainty')
  t2$Predictors = row.names(t2)
  # bottom rows show the supplementary model results using non-/strongly territorial species only
  t3 = summary(mod3)$fixed %>% mutate(mod = 'strong vs none territoriality')
  t3$Predictors = row.names(t3)

  t = rbind(t1,t2, t3) %>% 
  # re-order the variables
  mutate(Predictors = factor(Predictors, levels = c('Intercept','moSexual_selection_score','territorialityterritorial', 'Sociality',
                    'Mass','relative.beak.size', 'Habitat_density', 'Species.richness', 
                    'ave.freq', 'ave.dur'), 
                    labels = c('(Intercept)', rev(ylabs), 'Ave. Frequency', 'Ave. note length'))) %>% arrange(Predictors)
  t = t[, c(9,1:8)]
  return(t)
}


# Table S2: main results and supplementary analyses (species-level)
tidy.table(nc, nchc, nct) %>% write.csv("Table S2/Note count.csv", row.names = F)
tidy.table(nt, nthc, ntt) %>% write.csv("Table S2/Note type.csv", row.names = F)
tidy.table(sf, sfhc, sft) %>% write.csv("Table S2/Freq variation.csv", row.names = F)


# Table S3: supplementary analyses (song-level)
snc <- readRDS('NoteCount negbinomial __( main models , song level).rds')
snt <- readRDS('NumberofNoteTypes negbinomial __( main models , song level).rds')
ssf <- readRDS('sdfreq gaussian __( main models , song level).rds')

tidy.table.songs = function(mod){
  # top rows show the main model results using all species data
  t1 = summary(mod)$fixed %>% mutate(mod = 'all species')
  t1$Predictors = row.names(t1)
  t = t1 %>% 
    # re-order the variables
    mutate(Predictors = factor(Predictors, levels = c('Intercept','moSexual_selection_score','territorialityterritorial', 'Sociality',
                                                      'Mass','relative.beak.size', 'Habitat_density', 'Species.richness', 
                                                      'ave.freq', 'ave.dur'), 
                               labels = c('(Intercept)', rev(ylabs), 'Ave. Frequency', 'Ave. note length'))) %>% arrange(Predictors)
  t = t[, c(9,1:8)]
  return(t)
}

tidy.table.songs(snc) %>% write.csv("Table S3/Note count (song level).csv", row.names = F)
tidy.table.songs(snt) %>% write.csv("Table S3/Note type (song level).csv", row.names = F)
tidy.table.songs(ssf) %>% write.csv("Table S3/Freq variation (song level).csv", row.names = F)



###### Figure S3: additional complexity metrics ######

pk <- readRDS("PeakFreqsd gaussian __( SI , species level).rds")
bw <- readRDS("Bandwidthsd gaussian __( SI , species level).rds")
ns <- readRDS("Slopesd gaussian __( SI , species level).rds")
nl <- readRDS("NoteLengthsd gaussian __( SI , species level).rds")
sdu <- readRDS("SongDuration gaussian __( SI , species level).rds")
snr <- readRDS("NoteRate gaussian __( SI , species level).rds")

p1.1 = distplot(pk, 0.95, LL=-4, RR=4, 'Peak frequency SD\n')
p1.2 = distplot(bw, 0.95, LL=-3, RR=3, 'Bandwidth SD\n', ylabels = NULL)
p1.3 = distplot(ns, 0.95, LL=-1, RR=1.2, 'Note slope SD\n', ylabels = NULL)
p1.4 = distplot(nl, 0.95, LL=-3, RR=2.5, 'Note length SD\n')
p1.5 = distplot(sdu, 0.95, LL=-0.5, RR=1, 'Song duration\n', ylabels = NULL)
p1.6 = distplot(snr, 0.95, LL=-1.5, RR=1.5, 'Note rate\n', ylabels = NULL)

pdf('Figure S3.pdf', width = 9.1, height = 8.5)
cowplot::plot_grid(p1.1, p1.2, p1.3, 
                   p1.4, p1.5, p1.6, ncol = 3, rel_widths = c(1.8, 1, 1)) 
dev.off()

