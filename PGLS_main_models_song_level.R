setwd('')

library(dplyr)
library(ape)
library(caper)
library(brms)


# read the suboscine phylogeny
aostree <- read.tree('T400F_AOS_Clements_sppnames.tre')

# read the song level complexity data
song=read.csv('C - Song level data.csv')

# read the species level data to use the explanatory variables
ana <- read.csv('B - Species level data.csv')
# remove all song-related variables
variables =  c("Note_count" ,                "Number_of_note_types"  ,    
               "Peak_frequency_SD" ,         "Bandwidth_SD" ,              "Note_slope_SD",             
               "Frequency_variation" ,       "Note_length_SD" ,            "Song_duration" ,            
               "Note_rate",                  'Average_peak_frequency',     'Average_note_length' )
ana <- dplyr::select(ana, !contains(variables))

# merge the two datasets
ana=left_join(ana, song, by='Species')


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

# get var-cov matrix for Brms
ana$Species <- gsub(' ', '_', ana$Species)
aostree <- keep.tip(aostree, ana$Species) 
phy_cov <- ape::vcv.phylo(aostree)

# add an extra species name column for assigning both phylogenetic relationship and species identity as random effects in Brms.
ana$sp = ana$Species

# set weak priors for Brms models
weak_priors <- c(prior(normal(0,2), class="Intercept"),  # for the intercept
                 prior(normal(0,1), class="b"))          # for the slopes


mybrms <- function(brms_formula, file.tag){
  brm(brms_formula,
        data = ana,
        data2 = list(A=phy_cov),
      prior = weak_priors,
      init = 0,
        iter = 2000, 
        warmup = 1000,
        chains = 4,
        thin = 1, 
      seed = 1,
      cores = 4,
      file = paste(brms_formula$resp, brms_formula$family$family, '__(', file.tag, ', song level).rds'),
      backend = "cmdstanr")
}



###### Sensitivity analysis 4: re-analysing at the song level #########


brms_formula_note_count <- brmsformula(Note_count ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness + 
                            (1|gr(sp, cov=A)) + (1|Species), 
                            family = negbinomial())

mybrms(brms_formula_note_count, 'main models') 



brms_formula_note_type <- brmsformula(Number_of_note_types ~ mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness + 
                            (1|gr(sp, cov=A)) + (1|Species),  
                            family = negbinomial())
mybrms(brms_formula_note_type, 'main models')



brms_formula_freq_var <- brmsformula(sdfreq ~ ave.freq + mo(Sexual_selection_score) + territoriality + Sociality + Mass + relative.beak.size + Habitat_density + Species.richness +
                            (1|gr(sp, cov=A)) + (1|Species),  
                            family = gaussian())
mybrms(brms_formula_freq_var, 'main models') 


