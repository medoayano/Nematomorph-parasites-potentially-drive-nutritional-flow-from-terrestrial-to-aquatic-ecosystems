# Source code
# Estimates of body mass for aquatic invertebrates and crickets
# Calculate EPA (mg) per individual

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(brms)
library(forcats)
library(purrr)

# Data set for body mass
dataWeightAll <- read.csv("SourceData_BodyMass.csv", header = T) 

# Parameters
randomDraws <- 1000
idOrder <- unique(dataEpaAqua$order)

#### Estimate of posterior distribution of body mass for aquatic invertebrates ####

# Objects
posteriorWeightAqua <- list()
posteriorRandomWeightAqua <- list()

# Start estimation
for (i in 1:length(idOrder)) {
  
  dataW <- dataWeightAll %>% filter(order %in% c(idOrder[i]))
    
    print(idOrder[i])
    
    # Data frame
    dataWeight <- data.frame(
      weight_mean = mean(dataW$mass..mg.),
      weight_sd = sd(dataW$mass..mg.),
      weight_n = nrow(dataW)
    )
    
    # Standard error
    dataWeight$weight_se <- dataWeight$weight_sd / sqrt(dataWeight$weight_n)
    
    # Prior
    mu_prior <- mean(dataW$mass..mg.)
    sigma_prior <- sd(dataW$mass..mg.)
    
    # A Bayesian model
    modelWeightAqua <- brm(
      weight_mean | se(weight_se) ~ 1,
      data = dataWeight,
      family = gaussian(),
      prior = prior_string(paste0("normal(", mu_prior, ",", sigma_prior,")"), class = "Intercept"),
      iter = 10000, warmup = 2000, chains = 4, cores = 4, seed = 123,
      control = list(adapt_delta = 0.999, max_treedepth = 15)
    )
    
    # Extract pooled estimation of posterior distribution
    posteriorWeight <- as_draws_df(modelWeightAqua)$b_Intercept
    
    # Save the pooled estimation
    posteriorWeightAqua[[i]] <- data.frame(
      post.weight_per_weight = posteriorWeight,
      order = rep(idOrder[i],length(8000*4))
    )
    
    # Draw 1000 posterior samples
    posteriorRandomWeight <- sample(posteriorWeight, randomDraws)
    
    # Save the random samples
    posteriorRandomWeightAqua[[i]] <- data.frame(
      random.weight_per_weight = posteriorRandomWeight,
      order = rep(idOrder[i],length(randomDraws))
    )
}

# Rename lists
names(posteriorWeightAqua) <- idOrder
names(posteriorRandomWeightAqua) <- idOrder


#### Estimate of posterior distribution of body mass for camel crickets ####
# Data frame
dataWeightCrickets <- dataWeightAll %>% 
  filter(family %in% c("Rhaphidophoridae")) %>% 
  summarise(weight_mean = mean(mass..mg.), 
            weight_sd= sd(mass..mg.), 
            weight_n=n())

# Standard error
dataWeightCrickets$weight_se <- dataWeightCrickets$weight_sd / sqrt(dataWeightCrickets$weight_n)

# Prior
mu_prior <- dataWeightCrickets$weight_mean
sigma_prior <- dataWeightCrickets$weight_sd

# A Bayesian model
modelWeightCricket <- brm(
  weight_mean | se(weight_se)~ 1,
  data = dataWeightCrickets,
  family = gaussian(),
  prior = prior_string(paste0("normal(", mu_prior, ",", sigma_prior,")"), class = "Intercept"),
  iter = 10000, warmup = 2000, chains = 4, cores = 4, seed = 123,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

# Extract pooled estimation of posterior distribution
posteriorWeight <- as_draws_df(modelWeightCricket)$b_Intercept

# Save the pooled estimation
posteriorWeightCricket <- data.frame(
  post.fa_per_weight = posteriorWeight,
  order = rep("Cricket",length(8000*4))
)

# Draw 1000 posterior samples
posteriorRandomWeight <- posterior_predict(modelWeightCricket, ndraws = randomDraws)

# Save the random samples
posteriorRandomWeightCricket <- data.frame(
  random.weight_per_weight = posteriorRandomWeight,
  order = rep("Cricket",length(randomDraws))
)

#### Merged data set of body mass ####
mergedPosteriorRandomWeight <- rbind(do.call(rbind,posteriorRandomWeightAqua),
                          posteriorRandomWeightCricket)

# convert body weight unit (mg -> g)
mergedPosteriorRandomWeight$random.weight_per_weight <- mergedPosteriorRandomWeight$random.weight_per_weight*10^-3


#### Calculate EPA per individual ####
estimatedFaPerIndividual <- data.frame(order = mergedPosteriorRandomWeight$order,
           random.fa_per_individual = mergedPosteriorRandomFA1C$random.fa_per_weight*mergedPosteriorRandomWeight$random.weight_per_weight)


#### Figure 1C ####
# Summary table
tableFaPerIndividual <- map_dfr(unique(estimatedFaPerIndividual$order), ~{
  estimatedFaPerIndividual %>% 
    filter(order == .x) %>% 
    summarise(order = .x,
              mean=mean(random.fa_per_individual),
              lower_CI = quantile(random.fa_per_individual, 0.025, na.rm = TRUE),
              upper_CI = quantile(random.fa_per_individual, 0.975, na.rm = TRUE)) 
})

# Plot
tableFaPerIndividual$order <- fct_relevel(tableFaPerIndividual$order,
                                    "Cricket",
                                    "Ephemeroptera",
                                    "Plecoptera",
                                    "Trichoptera",
                                    "Amphipoda")

col2 <- c("#bd934f","#385686","#5582cc","#9bb2e3","#dee5f5")


p.C <- ggplot(tableFaPerIndividual, aes(x = mean, y = order, colour = order)) +
  geom_point(size = 1) + 
  geom_errorbar(aes(xmin = lower_CI, xmax = upper_CI), width = 0.2) +  
  theme_light() + 
  theme(panel.grid=element_blank(),
        text = element_text(size = 15),
        axis.text.y = element_blank(),
        legend.position = "none") +
  labs(title = "C",
       x = "EPA mg per individual",
       y = "") +
  scale_y_discrete(limits = rev(levels(tableFaPerIndividual$order))) +
  scale_color_manual(values = col2)
