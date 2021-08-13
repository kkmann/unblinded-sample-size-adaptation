# Estimating Conditional Power and Friends #####################################



# Setup ------------------------------------------------------------------------
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(patchwork)

set.seed(42)

source("R/util.R")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)



# Define & plot prior ----------------------------------------------------------
theta1 <- 0.4
 prior <- TruncatedNormal(mu = theta1, tau = 0.2, a = -0.5, b = 1.0)
# plot prior PDF
tibble(
      theta = seq(-1, 2, by = 0.005),
        PDF = map_dbl(theta, ~probability_density(prior, .))
  ) %>%
  mutate(
      support = case_when(
              theta < -0.5 ~ 'left',
              theta >  1.0 ~ 'right',
                      TRUE ~ 'support'
          )
  ) %>%
  ggplot() +
      aes(theta, PDF, group = support) +
      annotate('rect', xmin = -10, xmax = 0, ymin = -10, ymax = 10, fill = 'black', alpha = 0.2) +
      geom_line() +
      theme_bw() +
      labs(x = expression(theta)) +
      coord_cartesian(xlim = c(-0.5, 1.5), ylim = c(0, 2)) +
      theme(
          legend.position = 'top',
          panel.grid.minor = element_blank()
      )
ggsave("output/figures/prior_density.pdf", width = 8, height = 3.5)



# Derive required sample size for example --------------------------------------
alpha <- 0.025
 crit <- qnorm(1 - alpha)
 beta <- 0.2
# iterative increase n until expected power is matched
n <- 1
while (expected_power(n, crit, prior) < 1 - beta) {
    n <- n + 1
}
cat(sprintf(
    "required n: %i, expected power: %6.2f\n\r",
    n, 100*expected_power(n, crit, prior)
))



# Plot estimators of conditional power vs interim result for example -----------
# define interim time-point as 1/3 of final sample size
m <- round(n/3)
tibble(
          z = seq(-1, 3.5, by = 0.005),
        ACP = map_dbl(z, ~conditional_power(., m, n, crit, theta1)),
        OCP = map_dbl(z, ~observed_conditional_power(., m, n, crit)),
         PP = map_dbl(z, ~predictive_power(., m, n, crit, prior))
    ) %>%
    pivot_longer(-z, names_to = 'quantity') %>%
    ggplot() +
        aes(z, value, linetype = quantity) +
        geom_line() +
        labs(x = expression(z[m]), y = '') +
        scale_y_continuous("conditional power estimate", breaks = seq(0, 1, by = 0.2)) +
        scale_linetype('') +
    #ylim(c(0, .5)) +
        theme_bw() +
        theme(
            legend.position = 'top'
        )
ggsave("output/figures/acp_ocp_pp_example.pdf", width = 8, height = 3.5)



# Plot bias, MAE, and MSE of estimators vs theta -------------------------------
expand_grid(
        theta = seq(0, 0.75, by = 0.025),
           zm = rnorm(1000)
    ) %>%
    mutate(
           zm = zm + sqrt(m)*theta,
           CP = map2_dbl(zm, theta, ~conditional_power(..1, m, n, crit, ..2)),
          ACP =  map_dbl(zm, ~conditional_power(., m, n, crit, theta1)),
          OCP =  map_dbl(zm, ~observed_conditional_power(., m, n, crit)),
           PP =  map_dbl(zm, ~predictive_power(., m, n, crit, prior))
    ) %>%
    pivot_longer(c(ACP, OCP, PP)) %>%
    group_by(theta, name) %>%
    summarize(
        bias = mean(value - CP),
         MAE = mean(abs(value - CP)),
         MSE = mean((value - CP)^2),
        .groups = "drop"
    ) %>%
    pivot_longer(c(bias, MAE, MSE), names_to = "quantity") %>%
    ggplot() +
        aes(theta, value, linetype = name) +
        geom_line() +
        facet_wrap(~quantity) +
        labs(x = expression(theta), y = '', color = '') +
        scale_linetype('') +
        theme_bw() +
        theme(
             legend.position = 'top',
            panel.grid.minor = element_blank()
        )
ggsave("output/figures/acp_ocp_pp_bias_mae_mse.pdf", width = 8, height = 3.5)



# Plot sampling distribution of estimators for different theta -----------------
expand_grid(
        theta = seq(0, 0.6, by = 0.2),
           zm = rnorm(1000)
    ) %>%
    mutate(
           zm = zm + sqrt(m)*theta,
          ACP =  map_dbl(zm, ~conditional_power(., m, n, crit, theta1)),
          OCP =  map_dbl(zm, ~observed_conditional_power(., m, n, crit)),
           PP =  map_dbl(zm, ~predictive_power(., m, n, crit, prior)),
          PP2 =  map_dbl(zm, ~predictive_power(., m, n, crit, prior))
    ) %>%
    pivot_longer(c(ACP, OCP, PP)) %>%
    ggplot() +
        aes(value) +
        geom_histogram(bins = 25, position = 'identity') +
        facet_grid(name ~ theta, scales = 'free_y') +
        theme_bw() +
        labs(x = "conditional power estimate") +
        theme(
            panel.grid.minor = element_blank(),
            plot.subtitle = element_text(hjust = 0.5)
        ) +
        labs(subtitle = expression(theta))
ggsave("output/figures/acp_ocp_pp_sampling_distributions.pdf", width = 8, height = 3.5)
