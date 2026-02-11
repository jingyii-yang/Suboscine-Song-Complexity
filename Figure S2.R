library(tidyverse)
ana=read.csv('B - Species level data.csv')

# select data of the main explanatory variables
X_factors <- c('Species',
               'Sexual_selection_score',  'Territoriality', 'Sociality', 
               'Mass', 'Relative_beak_size', 'Habitat_density', 'Species_richness')
               
ana_x <- ana[,X_factors]
# calculate Spearman's correlation coefficients
cor_x = round(cor(ana_x[,2:8], use='complete.obs', method = 'spearman'), 2)
cor_x[upper.tri(cor_x)]<- NA

# add axis labels
Vs = c('Sexual selection score', 'Territoriality', 'Social group size', 'Mass', 'Beak size',
       'Habitat density', 'Species richness')

cor_V = reshape2::melt(cor_x, na.rm = T)  %>%
  # Reverse the order of one of the variables so that the x and y variables have opposing orders
  mutate(Var1 = factor(Var1, levels = rev(levels(.$Var1))))


# Figure S1
pdf('Figures/figure S2.pdf', height = 6, width = 6)
ggplot(cor_V, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#03045e", high = "#780000", mid = "white", 
     midpoint = 0, limit = c(-1,1), space = "Lab", name="Spearman's\nCorrelation") +
  scale_x_discrete(labels = Vs) +
  scale_y_discrete(labels = rev(Vs)) +
  coord_fixed() + # to draw squares
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme_minimal() + 
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        panel.grid.major = element_blank(),
        legend.position = c(0.8, 0.7),
        legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))
dev.off()

