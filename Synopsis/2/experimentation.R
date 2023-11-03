library(tidyverse)
library(Rcpp)
library(microbenchmark)
theme_set(theme_bw())
sourceCpp("simulate_birth_death.cpp")
#-------------------------------- Main simulation method --------------------------------#
simulate_birth_death <- function(lambda, mu, n, max_time) {
  current_time <- 0
  previous_time <- 0
  Z <- n
  times <- c(0)
  population_sizes <- c(n)
  extinction <- 0
  B_t <- 0 
  D_t <- 0
  S_t <- 0
  
  while (current_time < max_time) {
    # Simulate waiting time to next event
    rate <- Z * (lambda + mu)
    wait_time <- rexp(1, rate)
    previous_time <- current_time
    current_time <- current_time + wait_time
    
    # Update S_t
    S_t <- S_t + Z * (current_time - previous_time)
    
    # Decide whether the event is a birth or a death
    birth_prob <- lambda / (lambda + mu)
    is_birth <- runif(1) < birth_prob
    
    if (is_birth) {
      Z <- Z + 1
      B_t <- B_t + 1
    } else {
      Z <- Z - 1
      D_t <- D_t + 1
      if (Z == 0) {
        extinction <- 1
        break
      }
    }
    
    times <- c(times, current_time)
    population_sizes <- c(population_sizes, Z)
  }
  
  return(list(times = times, Z = population_sizes, extinction = extinction, S_t = S_t, B_t = B_t, D_t = D_t))
}

#-------------------------------- Parameters for test --------------------------------#

lambda <- 0.3
mu <- 0.6
n <- 50
max_time <- 10
num_simulations <- 1000
num_repeats <- 1000

# Testing that they give the same output with the same seed.
#testseed <- sample.int(1000, 1)
#set.seed(testseed)
#result1 <- simulate_birth_death(lambda, mu, n, max_time)
#set.seed(testseed)
result2 <- simulate_birth_death_cpp(lambda, mu, n, max_time)

# for(i in seq_along(result1)){
#   print(mean(result1[[i]] == result2[[i]]) == 1)
# }

# Calculate MLE for lambda and mu according to paper
lambda_hat <- result2$B_t / result2$S_t
mu_hat <- result2$D_t / result2$S_t

#-------------------------------- Experiment MLE --------------------------------#
# Experiment Parameters
lambda <- 0.6
mu <- 0.3
num_simulations <- 100
n_values <- c(1, 2, 3, 5, 10, 25)
max_time_values <- c(5, 7.5 ,10, 15, 20, 25, 30, 40)

results <- data.frame()

for(n in n_values){
  print(paste0("Working on n = ", n))
  for(max_time in max_time_values){
    print(paste0("Working on max_time = ", max_time))
    MLE_lambda <- numeric(num_simulations)
    MLE_mu <- numeric(num_simulations)
    time_vec <- numeric(num_simulations)
    for(i in 1:num_simulations){
      if (i %% 5 == 0) {
        print(i)
      }
      time_vec[i] <- microbenchmark(
        {
          result <- simulate_birth_death_cpp(lambda, mu, n, max_time)
          MLE_lambda[i] <- result$B_t / result$S_t
          MLE_mu[i] <- result$D_t / result$S_t
        }, times = 1
      )$time/1e9
    }
    
    ARE_lambda <- abs(MLE_lambda - lambda) / lambda
    ARE_mu <- abs(MLE_mu - mu) / mu
    
    quantiles_lambda <- quantile(ARE_lambda, probs = c(0.025, 0.5, 0.975))
    quantiles_mu <- quantile(ARE_mu, probs = c(0.025, 0.5, 0.975))
    quantiles_time <- quantile(time_vec, probs = c(0.025, 0.5, 0.975))
    
    results <- rbind(results, data.frame(n = n, max_time = max_time, 
                                         median_lambda = quantiles_lambda["50%"], 
                                         lq_lambda = quantiles_lambda["2.5%"], 
                                         uq_lambda = quantiles_lambda["97.5%"],
                                         median_mu = quantiles_mu["50%"], 
                                         lq_mu = quantiles_mu["2.5%"], 
                                         uq_mu = quantiles_mu["97.5%"],
                                         median_time = quantiles_time["50%"],
                                         lq_time = quantiles_time["2.5%"],
                                         uq_time = quantiles_time["97.5%"]
                                         ))
  }
}

# Plotting
ggplot(results, aes(x = max_time, y = median_lambda, group = factor(n), color = factor(n))) +
  geom_line() +
  geom_ribbon(aes(ymin = lq_lambda, ymax = uq_lambda, fill = factor(n)), alpha = 0.3) +
  labs(title = "Bias in MLE of Lambda", y = "Bias", x = "Max Time", color = 'n') + 
  facet_wrap(~factor(n)) + scale_y_log10()

ggplot(results, aes(x = max_time, y = median_mu, group = factor(n), color = factor(n))) +
  geom_line() +
  geom_ribbon(aes(ymin = lq_mu, ymax = uq_mu, fill = factor(n)), alpha = 0.3) +
  labs(title = "Bias in MLE of Lambda", y = "Bias", x = "Max Time", color = 'n') + 
  facet_wrap(~factor(n)) + scale_y_log10()

long_results <- results %>% select(-c(median_time, lq_time, uq_time)) %>% pivot_longer(-c(n, max_time), values_to = "Value", names_to = "Parameter") %>% 
  mutate(metric = ifelse(grepl("lambda", Parameter), "lambda", "mu"),
         quantile = case_when(
           grepl("median", Parameter) ~ "median",
           grepl("lq", Parameter) ~ "lq",
           grepl("uq", Parameter) ~ "uq"
         )) %>% pivot_wider(names_from = "quantile", values_from = "Value") %>% select(-Parameter) %>% 
  group_by(n, max_time, metric) %>% summarise(median = sum(median, na.rm = TRUE), lq = sum(lq, na.rm = TRUE), uq = sum(uq, na.rm = TRUE)) %>% ungroup()


p1 <-  long_results %>% 
  ggplot(aes(x = max_time, y = median, col = metric)) + 
  geom_line(linewidth = 1.75) + 
  geom_ribbon(aes(ymin = lq, ymax = uq, fill = metric), linewidth = 0.7,  alpha = 0.3, show.legend = FALSE) +
  facet_wrap(~n) +
  scale_y_log10(name = "ARE (log scale)") +
  labs(x = "Max time", color = "Parameter") +
  theme(
    text = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 12), 
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 12), 
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12) 
  ) +
  guides(color = guide_legend(override.aes = list(fill = NA))) + 
  scale_color_discrete(name = "Parameter", labels = c(expression(lambda), expression(mu)))

ggsave("ARE1.jpeg", plot = p1, height = 5, width = 8)

long_timing_results <- results %>% select(-c(median_lambda, lq_lambda, uq_lambda, median_mu, lq_mu, uq_mu)) %>% 
  pivot_longer(-c(n, max_time), values_to = "Value", names_to = "Quantile") %>% 
  mutate(Quantile = case_when(
           grepl("median", Quantile) ~ "median",
           grepl("lq", Quantile) ~ "lq",
           grepl("uq", Quantile) ~ "uq"
         )) %>% pivot_wider(names_from = "Quantile", values_from = "Value")


p2 <- long_timing_results %>% 
  ggplot(aes(x = max_time, y = median, col = factor(n))) + 
  geom_line(linewidth = 1.5) + 
  geom_ribbon(aes(ymin = lq, ymax = uq, fill = factor(n)), linewidth = 0.7,  alpha = 0.3, show.legend = FALSE) +
  scale_y_log10(name = "Computation time (s) (log scale)") +
  facet_wrap(~factor(n)) +
  labs(x = "Max time", color = "n") +
  theme(
    text = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 12), 
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 12), 
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(fill = NA)))

ggsave("times1.jpeg", plot = p2, height = 5, width = 8)
#-------------------------------- Experiment with extinction --------------------------------#
# avg_extinctions <- numeric(num_repeats)
# 
# for (i in 1:num_repeats) {
#   # Store results of whether each simulation dies out
#   dies_out <- logical(num_simulations)
#   
#   for (j in 1:num_simulations) {
#     result <- simulate_birth_death(lambda, mu, n, max_time)
#     dies_out[j] <- result$extinction
#   }
#   
#   avg_extinctions[i] <- mean(dies_out)
# }
# 
# # Theoretical extinction probability
# theoretical_prob <- min(1, mu/lambda)
# 
# # Compare the average extinction rate to the theoretical value
# bias <- (avg_extinctions - theoretical_prob) / theoretical_prob
# 
# tibble(bias) %>% ggplot(aes(x = bias)) + geom_histogram(col = "black", fill = "lightblue")


