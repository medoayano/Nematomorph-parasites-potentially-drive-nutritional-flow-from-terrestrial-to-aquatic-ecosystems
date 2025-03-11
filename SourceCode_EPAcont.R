# Source code
# Estimations of EPA content for aquatic invertebrates and crickets

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(brms)
library(forcats)
library(purrr)

#setwd("/Users/admin/Open Research")

#### EPA content of aquatic invertebrates ####
### A random effect meta-analysis ###

# Data set for EPA content of aquatic invertebrates
dataEpaAqua <- read.csv("SourceData_EPAcont_Aqua.csv", header = T) 

## Order level ##
# Parameters
randomDraws <- 1000
idOrder <- unique(dataEpaAqua$order)

# Objects
posteriorFaAqua <- list()
posteriorRandomFaAqua <- list()

# Start estimation
for (i in 1:length(idOrder)) {
  
  dataF <- dataEpaAqua %>% filter(order %in% c(idOrder[i]))
  dataF$group <- c(1:length(dataF$Reference))
  
  if(nrow(dataF)>0){
    
    print(idOrder[i])
    
    # Data frame
    dataFA <- data.frame(
      study = dataF$group,
      fa_mean = dataF$EPAmg.g_dw_mean,
      fa_sd = dataF$EPAmg.g_dw_sd,
      fa_n = dataF$Sample.size
    )
    
    # Remove 'NA'
    dataFA <- dataFA %>% filter(!is.na(fa_sd))
    
    # Standard error
    dataFA$fa_se <- dataFA$fa_sd / sqrt(dataFA$fa_n)
    
    # Prior
    mu_prior <- mean(dataFA$fa_mean)
    sigma_prior <- sd(dataFA$fa_mean)
    
    # A Bayesian model
    modelFaAqua <- brm(
      fa_mean | se(fa_se) ~ 1 + (1 | study),
      data = dataFA,
      family = gaussian(),
      prior = c(
        prior_string(paste0("normal(", mu_prior, ",", sigma_prior,")"), class = "Intercept"), # 事前分布（仮）
        prior(cauchy(0,0.3), class = "sd")
      ),
      iter = 10000, warmup = 2000, chains = 4, cores = 4, seed = 123,
      control = list(adapt_delta = 0.999, max_treedepth = 15)
    )
    
    # Check results
    print(summary(modelFaAqua))
    print(plot(modelFaAqua))
    print(pp_check(modelFaAqua))
    
    # Extract pooled estimation of posterior distribution
    posteriorFA <- as_draws_df(modelFaAqua)$b_Intercept
    
    # Save the pooled estimation
     posteriorFaAqua[[i]] <- data.frame(
      post.fa_per_weight = posteriorFA,
      order = rep(idOrder[i],length(8000*4))
    )
    
    # Draw 1000 posterior samples
    posteriorRandomFA <- sample(posteriorFA, randomDraws)
    
    # Save the random samples
    posteriorRandomFaAqua[[i]] <- data.frame(
      random.fa_per_weight = posteriorRandomFA,
      order = rep(idOrder[i],length(randomDraws))
    )
    
  }else{}
}

# Rename lists
names(posteriorFaAqua) <- idOrder
names(posteriorRandomFaAqua) <- idOrder

## Arthropoda level ##
# Parameters
randomDraws <- 1000

# Start estimation
dataF <- dataEpaAqua
dataF$group <- c(1:length(dataF$Reference))
 
  # Data frame
  dataFA <- data.frame(
    study = dataF$group,
    fa_mean = dataF$EPAmg.g_dw_mean,
    fa_sd = dataF$EPAmg.g_dw_sd,
    fa_n = dataF$Sample.size
  )
  
  # Standard error
  dataFA$fa_se <- dataFA$fa_sd / sqrt(dataFA$fa_n)
  
  # Prior
  mu_prior <- mean(dataFA$fa_mean)
  sigma_prior <- sd(dataFA$fa_mean)
  
  # A Bayesian model
  modelFaAquaAll <- brm(
    fa_mean | se(fa_se) ~ 1 + (1 | study),
    data = dataFA,
    family = gaussian(),
    prior = c(
      prior_string(paste0("normal(", mu_prior, ",", sigma_prior,")"), class = "Intercept"), # 事前分布（仮）
      prior(cauchy(0,0.3), class = "sd")
    ),
    iter = 10000, warmup = 2000, chains = 4, cores = 4, seed = 123,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )
  
  # Check results
  print(summary(modelFaAquaAll))
  print(plot(modelFaAquaAll))
  
  # Extract pooled estimation of posterior distribution
  posteriorFA <- as_draws_df(modelFaAquaAll)$b_Intercept
  
  # Save the posterior distribution
  posteriorFaAquaAll <- data.frame(
    post.fa_per_weight = posteriorFA,
    order = rep("Cricket",length(8000*4))
  )
  
  # Draw 1000 posterior samples
  posteriorRandomFA <- sample(posteriorFA, randomDraws)
  
  # Save the random samples
  posteriorRandomFaAquaAll <- data.frame(
    random.fa_per_weight = posteriorRandomFA,
    order = rep("Aqua",length(randomDraws))
  )


#### EPA content of camel crickets ####
# Data set for EPA content of camel crickets
dataEpaCricket <- read.csv("SourceData_EPAcont_Cricket.csv", header = T) 

# Parameters
randomDraws <- 1000

# Data frame
dataFA <- data.frame(
  fa_mean = mean(dataEpaCricket$EPA..mg.g.in.dry.weight.),
  fa_sd = sd(dataEpaCricket$EPA..mg.g.in.dry.weight.),
  fa_n = length(dataEpaCricket$EPA..mg.g.in.dry.weight.)
)

# Standard error
dataFA$fa_se <- dataFA$fa_sd / sqrt(dataFA$fa_n)

# Prior
mu_prior <- mean(dataEpaCricket$EPA..mg.g.in.dry.weight.)
sigma_prior <- sd(dataEpaCricket$EPA..mg.g.in.dry.weight.)

# A Bayesian model
modelFaCricket <- brm(
  fa_mean | se(fa_se) ~ 1,
  data = dataFA,
  family = gaussian(),
  prior = prior_string(paste0("normal(", mu_prior, ",", sigma_prior,")"), class = "Intercept"),
  iter = 10000, warmup = 2000, chains = 4, cores = 4, seed = 123,
  control = list(adapt_delta = 0.999, max_treedepth = 15)
)

# Check results
print(summary(modelFaCricket))
print(plot(modelFaCricket))

# Extract pooled estimation of posterior distribution
posteriorFA <- as_draws_df(modelFaCricket)$b_Intercept

# Save the posterior distribution
posteriorFaCricket <- data.frame(
  post.fa_per_weight = posteriorFA,
  order = rep("Cricket",length(8000*4))
)

# Draw 1000 posterior samples
posteriorRandomFA <- posterior_predict(modelFaCricket, ndraws = randomDraws)

# Save the random samples
posteriorRandomFaCricket <- data.frame(
  random.fa_per_weight = posteriorRandomFA,
  order = rep("Cricket",length(randomDraws))
)


#### Merged data set of EPA content ####
mergedPosteriorFA <- as.data.frame(rbind(do.call(rbind,posteriorFaAqua),
                                posteriorFaCricket))

mergedPosteriorRandomFA1C <- rbind(do.call(rbind,posteriorRandomFaAqua),
                      posteriorRandomFaCricket)

mergedPosteriorRandomFA1D <- data.frame(aqua = posteriorRandomFaAquaAll$random.fa_per_weight,
                                        cricket = posteriorRandomFaCricket$random.fa_per_weight)

#### Figure 1B ####
mergedPosteriorFA$order <- fct_relevel(mergedPosteriorFA$order,
                                     "Cricket",
                                     "Ephemeroptera",
                                     "Plecoptera",
                                     "Trichoptera",
                                     "Amphipoda")

col1 <- c("#e4c69f","#385686","#5582cc","#9bb2e3","#dee5f5")

p.B <- ggplot(mergedPosteriorFA,aes(x=post.fa_per_weight,fill=order)) +
  geom_density(aes(y=after_stat(density)/max(y=after_stat(density))),adjust=1) + 
  theme_light() + 
  theme(panel.grid=element_blank(),
        text = element_text(size = 15),
        legend.position = "none") +
  labs(title = "B",
       x = "EPA mg per g",
       y = "Density") +
  scale_fill_manual(values = col1) + 
  facet_grid(order ~ .)



