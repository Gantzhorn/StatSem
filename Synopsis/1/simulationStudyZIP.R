library(ggplot2)
library(magrittr)
library(tibble)
library(ZIM)
library(dplyr)

theme_set(theme_bw())

sample_mixture <- function(N, omega_t, lambda_t) {
  yt <- numeric(N)
  u <- numeric(N)
  u <- rbinom(n = N, size = 1, prob = omega_t)
  
  yt <- rpois(n = N, lambda = (1-u) * lambda_t)
  return(yt)
}

omega_t <- 0.50
lambda_t <- 1
numberOfSamples <- 5000

simulate_fit_transform <- function(Nsamples, omega, lambda) {
  yt_sample <- sample_mixture(Nsamples, omega, lambda)
  samples <- tibble(x_t = seq_along(yt_sample), y_t = yt_sample)
  
  # Fitting the ZIP model using zim
  zip_model <- ZIM::zim(y_t ~ 1 | 1, data = samples, dist = "zip")
  
  # Extracting and transforming parameters
  lambda_hat <- exp(zip_model$para[1])
  omega_hat <- exp(zip_model$para[2]) / (1 + exp(zip_model$para[2]))
  
  return(c(lambda_hat, omega_hat))
}

# Number of repetitions
n_reps <- 500

params_list <- replicate(n_reps, simulate_fit_transform(numberOfSamples, omega_t, lambda_t),
                         simplify = FALSE)

params_df <- do.call(rbind, params_list)
colnames(params_df) <- c("lambda_hat", "omega_hat")

params_df <- as_tibble(params_df)

# Plotting relative error in parameter estimates
params_df %>% 
  ggplot(aes(x = (lambda_hat - lambda_t)/lambda_t)) + 
  geom_density() + 
  geom_rug() +
  labs(x = "Relative Error",
       y = "Density")

params_df %>% 
  ggplot(aes(x = (omega_hat - omega_t)/omega_t)) + 
  geom_density() + 
  geom_rug() +
  labs(x = "Relative Error",
       y = "Density")

# Values to loop through
omega_values <- c(0.1, 0.25, 0.5, 0.75, 0.9)
Nsamples_values <- c(50, 150, 200, 500, 1000, 1500, 2000)
lambda <- 5

# Data frame to store results
results <- tibble()
omega_values <- 0.8
# Loop through omega and Nsamples values
for(omega in omega_values) {
  print(omega)
  for(Nsamples in Nsamples_values) {
    print(Nsamples)
    params_list <- replicate(n_reps, simulate_fit_transform(Nsamples, omega, lambda),
                             simplify = TRUE)
    
    params_df <- as.data.frame(t(params_list))
    colnames(params_df) <- c("lambda_hat", "omega_hat")
    
    params_df$omega <- omega
    params_df$Nsamples <- Nsamples
    
    results <- rbind(results, params_df)
  }
}

# Calculate median and confidence intervals
summary_results <- results %>%
  group_by(omega, Nsamples) %>%
  summarise(
    lambda_median = mean(lambda_hat),
    omega_median = mean(omega_hat),
    lambda_lower = quantile(lambda_hat, 0.025),
    lambda_upper = quantile(lambda_hat, 0.975),
    omega_lower = quantile(omega_hat, 0.025),
    omega_upper = quantile(omega_hat, 0.975)
  ) %>% 
  mutate(lambda_median = lambda_median - 5,
         lambda_lower = lambda_lower - 5,
         lambda_upper = lambda_upper - 5,
         omega_median = omega_median - omega,
         omega_lower = omega_lower - omega,
         omega_upper = omega_upper - omega)

# Plotting
ggplot(summary_results, aes(x = Nsamples, y = omega_median, color = as.factor(omega))) +
  geom_line()  + 
  geom_ribbon(aes(ymin = omega_lower, ymax = omega_upper, fill = as.factor(omega)), alpha = 0.1) +
  labs(
    x = "Sample Size",
    y = "Bias",
    color = "True Omega"
  ) +
  guides(fill=FALSE) +
  scale_x_log10()

plot1 <- ggplot(summary_results, aes(x = Nsamples, y = lambda_median, color = as.factor(omega))) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(aes(ymin = lambda_lower, ymax = lambda_upper, fill = as.factor(omega)), alpha = 0.25) +
  labs(
    x = "Sample Size",
    y = "Bias of Lambda",
    color = "True Omega"
  ) +
  guides(fill=FALSE,
         color = guide_legend(override.aes = list(fill = NA))) +
  facet_wrap(~as.factor(omega), labeller = label_both) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.title = element_text(face = "bold"),
    legend.key = element_rect(fill = "white"),
    axis.title = element_text(face = "bold")
  )

ggsave("plot1.jpeg", width = 175, height = 75, units = "mm")
