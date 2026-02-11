setwd('')

library(dplyr)
library(ape)
library(caper)
library(brms)

# read the suboscine phylogeny
aostree <- read.tree('T400F_AOS_Clements_sppnames.tre')
# read the species level data
ana <- read.csv('B - Species level data.csv')
# align the species names in the two files
ana$Species <- gsub(' ', '_', ana$Species)

# round the count data (species-level median values) to integers to use the negative-binomial family
ana$Note.Count_rd <- round(ana$Note_count)
ana$Number_of_Note_Types.rd <- round(ana$Number_of_note_types)
# normalise the response variable to use the Gaussian family
ana$sdfreq = sqrt(ana$Frequency_variation)

ana$territoriality <- ifelse(ana$Territoriality == 0, 'non-territorial', 'territorial')
# standardise explanatory variables
ana$Sociality <- scale(ana$Sociality)/2
ana$Mass <- scale(log(ana$Mass))/2
ana$relative.beak.size <- scale(ana$Relative_beak_size)/2
ana$Species.richness = scale(ana$Species_richness)/2
ana$ave.freq = scale(ana$Average_peak_frequency)/2
ana$ave.dur = scale(ana$Average_note_length)/2

# add species data to the phylogeny
anatree <- comparative.data(phy = aostree, data = ana,
                              names.col = 'Species', 
                              vcv = TRUE,
                              warn.dropped = TRUE,
                              na.omit = FALSE)

phy_cov <- ape::vcv.phylo(anatree$phy)
anatree$data$sp <- row.names(anatree$data)

# set weak priors for Brms models
weak_priors <- c(prior(normal(0,2), class="Intercept"),  # for the intercept
                 prior(normal(0,1), class="b"))          # for the slopes


mybrms <- function(brms_formula, file.tag){
  brm(brms_formula,
        data = anatree$data,
        data2 = list(A=phy_cov),
      prior = weak_priors,
      init = 0,
        iter = 2000,
        warmup = 1000,
        chains = 4,
        thin = 1,
      seed = 1,
      cores = 4,
      file = paste(brms_formula$resp, brms_formula$family$family, '__(', file.tag, ', species level).rds'),
      backend = "cmdstanr")
}



###### Main analyses ######

brms_formula_note_count <- brmsformula(Note.Count_rd ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness + 
                            (1|gr(sp, cov=A)), 
                            family = negbinomial())
mybrms(brms_formula_note_count, 'main models')


brms_formula_note_type <- brmsformula(Number_of_Note_Types.rd ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness + 
                            (1|gr(sp, cov=A)), 
                            family = negbinomial())
mybrms(brms_formula_note_type, 'main models')  


brms_formula_freq_var <- brmsformula(sdfreq ~ ave.freq + mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness +  
                            (1|gr(sp, cov=A)), 
                            family = gaussian())
mybrms(brms_formula_freq_var, 'main models') 




###### Sensitivity Analysis 1: using alternative complexity metrics ######

# normalise the response variable to use the Gaussian family
ana$Peak.Freq.sd = sqrt(ana$Peak_frequency_SD)
ana$Bandwidth.sd = sqrt(ana$Bandwidth_SD)
ana$Slope.sd = sqrt(ana$Note_slope_SD)
ana$Note.Length.sd = sqrt(ana$Note_length_SD)
ana$Song.Duration = sqrt(ana$Song_duration)
ana$Note.Rate = sqrt(ana$Note_rate)

# add species data to the phylogeny
anatree <- comparative.data(phy = aostree, data = ana,
                              names.col = 'Species', 
                              vcv = TRUE,
                              warn.dropped = TRUE,
                              na.omit = FALSE)

phy_cov <- ape::vcv.phylo(anatree$phy)
anatree$data$sp <- row.names(anatree$data)


# additional spectral metrics
brms_formula_spec1 <- brmsformula(Peak.Freq.sd ~ ave.freq + mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness +   
                            (1|gr(sp, cov=A)),
                            family = gaussian())
mybrms(brms_formula_spec1, 'SI')


brms_formula_spec2 <- brmsformula(Bandwidth.sd ~ ave.freq + mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness + 
                            (1|gr(sp, cov=A)),
                           family = gaussian())
mybrms(brms_formula_spec2, 'SI')


brms_formula_spec3 <- brmsformula(Slope.sd ~ ave.freq + mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness +  
                            (1|gr(sp, cov=A)),
                           family = gaussian())
mybrms(brms_formula_spec3, 'SI')


# additional temporal metrics
brms_formula_temp1 <- brmsformula(Note.Length.sd ~ ave.dur + mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness +   
                            (1|gr(sp, cov=A)),
                            family = gaussian())
mybrms(brms_formula_temp1, 'SI')


brms_formula_song1 <- brmsformula(Song.Duration ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness +    
                            (1|gr(sp, cov=A)),
                            family = gaussian())
mybrms(brms_formula_song1, 'SI')


brms_formula_song2 <- brmsformula(Note.Rate ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness +    
                            (1|gr(sp, cov=A)),
                           family = gaussian())
mybrms(brms_formula_song2, 'SI')



###### Sensitivity Analysis 2: using species with high certainty data ######

ana2 <- filter(ana, Sexual_selection_certainty %in% c('A', 'B') & ana$Territoriality_certainty %in% c('A', 'B') & ana$Sociality_certaint %in% c('A', 'B'))

anatree <- comparative.data(phy = aostree, data = ana2,
                              names.col = 'Species', 
                              vcv = TRUE,
                              warn.dropped = TRUE,
                              na.omit = FALSE)

phy_cov <- ape::vcv.phylo(anatree$phy)
anatree$data$sp <- row.names(anatree$data)


brms_formula_note_count2 <- brmsformula(Note.Count_rd ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness + 
                            (1|gr(sp, cov=A)), 
                            family = negbinomial())
mybrms(brms_formula_note_count2, 'high certainty data')


brms_formula_note_type2 <- brmsformula(Number_of_Note_Types.rd ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness + 
                            (1|gr(sp, cov=A)), 
                            family = negbinomial())
mybrms(brms_formula_note_type2, 'high certainty data') 



brms_formula_freq_var2 <- brmsformula(sdfreq ~ ave.freq + mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness +  
                            (1|gr(sp, cov=A)), 
                            family = gaussian())
mybrms(brms_formula_freq_var2, 'high certainty data') 



###### Sensitivity Analysis 3: removing weakly and seasonally territorial species ######

ana3 = filter(ana, Territoriality != 1)

anatree <- comparative.data(phy = aostree, data = ana3,
                            names.col = 'Species', 
                            vcv = TRUE,
                            warn.dropped = TRUE,
                            na.omit = FALSE)

phy_cov <- ape::vcv.phylo(anatree$phy)
anatree$data$sp <- row.names(anatree$data)


brms_formula_note_count2 <- brmsformula(Note.Count_rd ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness + 
                                            (1|gr(sp, cov=A)), 
                                        family = negbinomial())
mybrms(brms_formula_note_count2, 'strong vs none territoriality')


brms_formula_note_type2 <- brmsformula(Number_of_Note_Types.rd ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness + 
                                           (1|gr(sp, cov=A)), 
                                       family = negbinomial())
mybrms(brms_formula_note_type2, 'strong vs none territoriality') 



brms_formula_freq_var2 <- brmsformula(sdfreq ~ ave.freq + mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness +  
                                          (1|gr(sp, cov=A)), 
                                      family = gaussian())
mybrms(brms_formula_freq_var2, 'strong vs none territoriality')

