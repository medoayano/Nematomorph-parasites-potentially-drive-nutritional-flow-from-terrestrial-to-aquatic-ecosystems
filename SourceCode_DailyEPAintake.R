# Source code
# Calculate daily consumption rate

library(ggplot2)
library(dplyr)
library(gridExtra)
library(purrr)
library(brms)

#### Daily consumption rate ####
# Data set 
dataSCW <- read.csv("SourceData_ConsumptionRate_SCW.csv", header = T) 
dataWT <- read.csv("SourceData_ConsumptionRate_WT.csv", header = T) 

# Parameters
time <- unique(dataSCW$Time)
timing <- unique(dataSCW$Timing)
category <- c("aqua_total","cricket_total")  # <- total or aqua_total or cricket_total
randomDraws <- 1000

# Object
fit <- list()
listnames <- NULL
n <- 1

# Start estimation
for (x in 1:length(category)) {
  
  for (i in 1:length(timing)) {
    
    for (t in 1:(length(time)-1)) {
      
      # Count
      print(paste0(category[x],"_",timing[i],"_",time[t],"_",time[t+1]))
      
      # parameter 1 - initial weight in trout stomach
      dataS0 <- dataSCW %>% 
        filter(Time == time[t], Timing == timing[i]) %>% 
        pull(category[x]) %>% 
        as.numeric()
      
      # Table with mean, SD, and sample size
      tableS0 <- data.frame(
        s0_mean = mean(dataS0),
        s0_sd = sd(dataS0),
        s0_n = length(dataS0)
      )
      
      if(mean(dataS0)>0){
        
        # Standard error
        tableS0$s0_se <- tableS0$s0_sd / sqrt(tableS0$s0_n)
        
        # Prior
        mu_prior <- tableS0$s0_mean
        sigma_prior <- tableS0$s0_sd
        
        # A Bayesian model
        modelS0 <- brm(
          s0_mean | se(s0_se) ~ 1,
          data = tableS0,
          family = gaussian(),
          prior = prior_string(paste0("normal(", mu_prior, ",", sigma_prior,")"), class = "Intercept"),
          iter = 10000, warmup = 2000, chains = 4, cores = 4, seed = 123,
          control = list(adapt_delta = 0.999, max_treedepth = 15)
        )
        
        # Check results
        print(summary(modelS0))
        print(plot(modelS0))
        
        # Extract pooled estimation of posterior distribution
        posteriorS0 <- as_draws_df(modelS0)$b_Intercept
        
        # Draw 1000 posterior samples
        posteriorRandomS0 <- posterior_predict(modelS0, ndraws = randomDraws)
        
      }else{
        posteriorRandomS0 <- rep(0,randomDraws)
      }
      
      ## parameter 2 - weight in trout stomach after t hours
      dataSt <- dataSCW %>% 
        filter(Time == time[t+1], Timing == timing[i]) %>% 
        pull(category[x]) %>% 
        as.numeric()
      
      # Table with mean, SD, and sample size
      tableSt <- data.frame(
        st_mean = mean(dataSt),
        st_sd = sd(dataSt),
        st_n = length(dataSt)
      )
      
      if(mean(dataSt)>0){
        
        # Standard error
        tableSt$st_se <- tableSt$st_sd / sqrt(tableSt$st_n)
        
        # Prior
        mu_prior <- tableSt$st_mean
        sigma_prior <- tableSt$st_sd
        
        # A Bayesian model
        modelSt <- brm(
          st_mean | se(st_se) ~ 1,
          data = tableSt,
          family = gaussian(),
          prior = prior_string(paste0("normal(", mu_prior, ",", sigma_prior,")"), class = "Intercept"),
          iter = 10000, warmup = 2000, chains = 4, cores = 4, seed = 123,
          control = list(adapt_delta = 0.999, max_treedepth = 15)
        )
        
        # Check results
        print(summary(modelSt))
        print(plot(modelSt))
        
        # Extract pooled estimation of posterior distribution
        posteriorSt <- as_draws_df(modelSt)$b_Intercept
        
        # Draw 1000 posterior samples
        posteriorRandomSt <- posterior_predict(modelSt, ndraws = randomDraws)
        
      }else{
        posteriorRandomSt <- rep(0,randomDraws)
      }
      
      # Water temperature at period "t"
      WT <- dataWT %>% 
        filter(Time == time[t+1]) %>% 
        pull(Temperature) %>% 
        as.numeric()
      
      # Calculate Rt
      Rt <- exp(0.224*WT-5.44)
      
      # Fit model
      fit[[n]] <- (posteriorRandomSt - posteriorRandomS0*exp(-Rt))*Rt/(1-exp(-Rt))
      listnames[n] <- paste0(category[x],"_",timing[i],"_",time[t],"_",time[t+1])
      
      # Update counter
      n <- n+1
      
    }
  }
}

# Rename list
names(fit) <- listnames

# Save daily consumption rate
dataCR <- data.frame(
  fit_aqua_temp1 = Reduce("+", fit[1:5]),
  fit_aqua_temp2 = Reduce("+", fit[6:10]),
  fit_cricket_temp1 = Reduce("+", fit[11:15])*0.95,
  fit_cricket_temp2 = Reduce("+", fit[16:20])*0.95
)

#### Merged data set of EPA content ####
mergedPosteriorRandomFA1D # from "SourceCode_EPA_cont"

#### Calculate daily EPA intake (DEI) ####
# Parameters
fishMassPer100m <- 155529.0

# Calculate
dataDEI <- data.frame(dataCR[1:2]*mergedPosteriorRandomFA1D$aqua*fishMassPer100m*10^-3,
                        dataCR[3:4]*mergedPosteriorRandomFA1D$cricket*fishMassPer100m*10^-3)

tableDEI <- sapply(dataDEI, function(x) {
  n <- length(x)
  mean_x <- mean(x)
  se_x <- sd(x) / sqrt(n)
  error_margin <- qt(0.975, df = n - 1) * se_x
  c(mean = mean_x, 
    lower = mean_x - error_margin, 
    upper = mean_x + error_margin)
})

colnames(tableDEI) <- c("Aquatic insects (T1)",
                              "Aquatic insects (T2)",
                              "Cricket (T1)",
                              "Cricket (T2)")

tableDEI <- data.frame(
  order = colnames(tableDEI),
  t(tableDEI)
)

### Figure ###
col3 <- c("#385686","#385686","#bd934f","#bd934f")

p.D <- ggplot(tableDEI, aes(x = mean, y = order, colour = order)) +
  geom_point(size = 1) + 
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +  
  theme_light() + 
  theme(panel.grid=element_blank(),
        text = element_text(size = 15),
        axis.text.y = element_blank(),
        legend.position = "none") +
  labs(title = "D",
       x = "daily EPA intake",
       y = "") +
  scale_y_discrete(limits=c("Aquatic insects (T2)",
                            "Aquatic insects (T1)",
                            "Cricket (T2)",
                            "Cricket (T1)")) +
  scale_color_manual(values = col3)

