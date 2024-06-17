# Load necessary libraries
library(lme4)    # For glmer
library(dplyr)   # For data manipulation
library(sandwich) # For clustered standard errors
library(lmtest)   # For significance testing
library(ggplot2)  # For plotting

# Parameters
set.seed(123)        # For reproducibility
N <- 8000            # Number of individuals
T <- 24              # Number of time periods (12 pre-treatment, 12 post-treatment)
P <- 0.5             # Proportion in treatment group
outcome_prob <- 0.2  # Baseline probability for outcome == 1
num_simulations <- 1000  # Number of simulations
effect_size <- 0.5   # Expected effect size (difference-in-difference)
alpha <- 0.05        # Significance level
num_clusters <- 73   # Number of clusters

# Function to simulate data and fit the DiD model with clustered standard errors
simulate_DiD <- function(N, T, P, outcome_prob, effect_size, num_clusters) {
  # Create individuals and time periods
  individuals <- rep(1:N, each = T)
  time <- rep(1:T, N)
  
  # Treatment assignment
  treatment <- rbinom(N, 1, P)
  treatment <- rep(treatment, each = T)
  
  # Post-treatment indicator
  post <- ifelse(time > 12, 1, 0)
  
  # Cluster assignment
  clusters <- rep(1:num_clusters, length.out = N)
  clusters <- rep(clusters, each = T)
  
  # Control variables (5 categorical variables)
  X1 <- rep(sample(1:2, N, replace = TRUE), each = T) # criminal history (binary)
  X2 <- rep(sample(1:5, N, replace = TRUE), each = T) # offense type
  X3 <- rep(sample(1:4, N, replace = TRUE), each = T) # age
  X4 <- rep(sample(1:2, N, replace = TRUE), each = T) # gender
  X5 <- rep(sample(1:4, N, replace = TRUE), each = T) # race/ethnicity
  
  # Linear predictor
  linear_predictor <- log(outcome_prob / (1 - outcome_prob)) + 
    1*treatment + 
    1*post + 
    effect_size*(treatment * post) + 
    0.5*X1 + 
    0.5*X2 + 
    0.5*X3 + 
    0.5*X4 + 
    0.5*X5
  
  # Probability (logistic function is applied to this linear predictor to compute the probability of the outcome)
  prob <- 1 / (1 + exp(-linear_predictor))
  
  # Outcome variable (binary)
  Y <- rbinom(N * T, 1, prob)
  
  # Create data frame
  data <- data.frame(individuals, time, treatment, post, Y, X1, X2, X3, X4, X5, clusters)
  
  return(data)
}

# Function to fit the DiD model and get the p-value for the interaction term
fit_model <- function(data) {
  # Difference-in-Difference model for binary outcome
  model <- glm(Y ~ treatment * post + X1 + X2 + X3 + X4 + X5, data = data, family = binomial)
  
  # Compute clustered standard errors
  cluster_se <- vcovCL(model, cluster = ~clusters)
  
  # Get the p-value for the interaction term (treatment:post)
  p_value <- coeftest(model, vcov. = cluster_se)["treatment:post", 4]
  
  # Return whether the null hypothesis is rejected
  return(p_value < alpha)
}

# Run simulations
results <- replicate(num_simulations, {
  data <- simulate_DiD(N, T, P, outcome_prob, effect_size, num_clusters)
  fit_model(data)
})

# Calculate power
power <- mean(results)
print(paste("Estimated power for Difference-in-Difference model with binary outcome and treatment: ", power))

# Function to calculate power for a given effect size
calculate_power <- function(effect_size) {
  results <- replicate(num_simulations, {
    data <- simulate_DiD(N, T, P, outcome_prob, effect_size, num_clusters)
    fit_model(data)
  })
  
  power <- mean(results)
  return(power)
}

# Calculate power for different effect sizes
effect_sizes <- seq(0.1, 1.0, by = 0.1)
power_values <- sapply(effect_sizes, calculate_power)

# Create a data frame for plotting
power_curve_data <- data.frame(effect_size = effect_sizes, power = power_values)

# Plot power curve
ggplot(power_curve_data, aes(x = effect_size, y = power)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Power Curve for Difference-in-Difference Model",
       x = "Effect Size",
       y = "Power")