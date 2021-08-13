# Naive Pre-planned Recalculation Based on Predictive Power and
# Optimal Pre-planned Recalculation ############################################



# Setup ------------------------------------------------------------------------
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(patchwork)
library(future)
library(adoptr)

plan(multisession)

set.seed(42)
source("R/util.R")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)



# Define planning assumptions --------------------------------------------------
theta1 <- 0.4
 prior <- TruncatedNormal(mu = theta1, tau = 0.2, a = -0.5, b = 1.0)
 alpha <- 0.025
  crit <- qnorm(1 - alpha)
  beta <- 0.2
# compute required sample size based on expected power
n <- 1
while (expected_power(n, crit, prior) < 1 - beta) {
    n <- n + 1
}
cat(sprintf(
    "required n: %i, expected power: %6.2f\n\r",
    n, 100*expected_power(n, crit, prior)
))



# Define recalculation rule ----------------------------------------------------
adapt_naive <- function(m, zm, beta_cond, n_min, n_max, n_old, c_old) {
    # original conditional error
    alpha_bar <- conditional_power(zm, m, n_old, c_old, 0)
    target <- function(x) {
        n_bar <- x[1]; c_bar <- x[2]
        (conditional_power(zm, m, n_bar, c_bar, 0) - alpha_bar)^2 +
            (predictive_power(zm, m, n_bar, c_bar, prior) - (1 - beta_cond))^2
    }
    res <- optim(
        par = c(n, crit),
        fn = target,
        lower = c(n_min, 0),
        upper = c(n_max, 5),
        method = "L-BFGS-B",
        control = list(fnscale = 1e-4, parscale = c((n_max + n_min)/2, 1))
    )
    list(
        n = if (abs(res$par[1] - n_max) < 0.1) m else res$par[1],
        c = if (abs(res$par[1] - n_max) < 0.1) Inf else res$par[2],
        optim = res
    )
}



# Plot 'mandatory' recalculation design ----------------------------------------
# interim time-point m = 26
m <- round(n/3)
# define lower and upper bound for recalculation
n_min <- 30
n_max <- 160
tbl_adapted <- tibble(
            zm = seq(0.0, 3.5, by = 0.01),
        design = "adapted",
           tmp = map(zm, function(zm) {
                res <- adapt_naive(m, zm, beta_cond = beta, n_min = n_min, n_max = n_max, n_old = n, c_old = crit)
                tibble(n = res$n, c = res$c)
            }
        )
    ) %>%
    unnest(tmp) %>%
    bind_rows(
        tibble(
            zm = seq(0.0, 3.5, by = 0.01),
        design = "original",
             n = n,
             c = crit
        )
    ) %>%
    mutate(
             PP = pmap_dbl(list(zm, n, c), ~predictive_power(..1, m, ..2, ..3, prior)),
        section = if_else(c == Inf, "early futility", "continue"),
         design = factor(design, levels = c("original", "adapted"))
    )
# first panel: sample size
plt1 <- ggplot(tbl_adapted) +
    aes(zm, n, group = interaction(section, design)) +
    geom_line(aes(linetype = design)) +
    scale_x_continuous(expression(z[m]), limits = c(0, 3.5)) +
    scale_y_continuous('sample size', limits = c(0, n_max)) +
    scale_linetype('') +
    theme_bw() +
    theme(
        legend.position = "top"
    )
# second panel: critical value
plt2 <- ggplot(tbl_adapted) +
    aes(zm, c, group = interaction(section, design)) +
    geom_line(aes(linetype = design)) +
    scale_x_continuous(expression(z[m]), limits = c(0, 3.5)) +
    scale_y_continuous('critical value', limits = c(0, NA)) +
    scale_linetype('') +
    theme_bw() +
    theme(
        legend.position = "top"
    )
# third panel: predictive power
plt3 <- ggplot(tbl_adapted) +
    aes(zm, PP, group = interaction(section, design)) +
    geom_line(aes(linetype = design)) +
    scale_x_continuous(expression(z[m]), limits = c(0, 3.5)) +
    scale_y_continuous('predictive power') +
    scale_linetype('') +
    theme_bw() +
    theme(
        legend.position = "top"
    )
# compose plot
plt1 + plt2 + plt3 + plot_layout(guides = "collect") & theme(legend.position = 'top')



# Compute properties of 'mandatory' recalculation design -----------------------
# derive early stopping boundaries and corresponding probabilities
a <- tbl_adapted %>%
    filter(design == "adapted", section == "early futility") %>%
    arrange(desc(zm)) %>%
    pull(zm) %>%
    head(1)
pr_a <- predictive_cdf(prior, m, a)
b <- tbl_adapted %>%
    filter(design == "adapted", n == n_min) %>%
    arrange(zm) %>%
    pull(zm) %>%
    head(1)
pr_b <- 1 - predictive_cdf(prior, m, b)
# expected sample size
ess <- pr_a * m + # early futility
    integrate_gl( # continue
        function(zm) {
            predictive_pdf(prior, m, zm) *
            adapt_naive(m, zm, beta_cond = 0.2, n_min = n_min, n_max = n_max, n_old = n, c_old = crit)$n
        },
        a, b
    ) +
    pr_b * n_min # early efficacy
# expected squared sample size (for computation of sample size variance)
esssq <- pr_a * m^2 + # early futility
    integrate_gl( # continue
        function(zm) {
            predictive_pdf(prior, m, zm) *
            adapt_naive(m, zm, beta_cond = 0.2, n_min = n_min, n_max = n_max, n_old = n, c_old = crit)$n^2
        },
        a, b
    ) +
    pr_b * n_min^2 # early efficacy
# maximal type one error rate
mtoer <- integrate_gl(
    function(zm) {
        dnorm(zm, 0, 1) * conditional_power(
            zm, m,
            n = adapt_naive(m, zm, beta_cond = 0.2, n_min = n_min, n_max = n_max, n_old = n, c_old = crit)$n,
            c = adapt_naive(m, zm, beta_cond = 0.2, n_min = n_min, n_max = n_max, n_old = n, c_old = crit)$c,
            theta = 0
        )
    },
    a, 10 # cut off at 10 instead of transforming to [a, inf)
)
# expected power
ep <- integrate_gl(
    function(zm) {
        predictive_pdf(condition(prior, 0), m, zm) * predictive_power(
            zm, m,
            n = adapt_naive(m, zm, beta_cond = 0.2, n_min = n_min, n_max = n_max, n_old = n, c_old = crit)$n,
            c = adapt_naive(m, zm, beta_cond = 0.2, n_min = n_min, n_max = n_max, n_old = n, c_old = crit)$c,
            prior = condition(prior, 0)
        )
    } ,
    a, 10 # cut off at 10 instead of transforming to [a, inf)
)
cat(sprintf(
    "optimal two-stage design\n\r   E[N]: %.1f\n\r  SD[N]: %.1f\n\r  MTOER: %.4f\n\r     EP: %.4f\n\r",
    ess, sqrt(esssq - ess^2), mtoer, ep
))


# Compute & plot optimal two-stage design --------------------------------------
# define assumptions and scores in adoptr
   H0 <- PointMassPrior(.0, 1)
   H1 <- ContinuousPrior(function(theta) dnorm(theta, theta1, 0.2)/(pnorm(1, theta1, 0.2) - pnorm(-0.5, theta1, 0.2)), c(-0.5, 1))
  ESS <- ExpectedSampleSize(Normal(two_armed = FALSE), H1)
    N <- ConditionalSampleSize()
ESSSq <- expected(composite({N^2}), Normal(two_armed = FALSE), H1)
 SDSS <- composite({sqrt(ESSSq - ESS^2)})
   EP <- Power(Normal(two_armed = FALSE), adoptr::condition(H1, c(0, 1)))
 TOER <- Power(Normal(two_armed = FALSE), H0)
   PP <- ConditionalPower(Normal(two_armed = FALSE), adoptr::condition(H1, c(0, 1)))
   CE <- ConditionalPower(Normal(two_armed = FALSE), H0)
# define the optimization problem and solve numerically
optimal_design <- minimize(
    ESS,
    subject_to(
          EP >= 1 - beta,
        TOER <= alpha
    ),
    get_initial_design(
        theta = .4,
        alpha = .025,
        beta  = .2,
        type  = 'two-stage',
        dist  = Normal(two_armed = FALSE),
        order = 11L
    )
)$design
# plot
tbl_optimal_design <- tibble(
             zm = seq(0, 3.5, by = 0.01),
         design = "optimal",
              n = n(optimal_design, zm, round = FALSE),
             c2 = c2(optimal_design, zm),
            # transform stage two critical value to overall critical value
              c = sqrt(n - n1(optimal_design)) / sqrt(n) * c2 + sqrt(n1(optimal_design)/n) * zm,
             PP = evaluate(PP, optimal_design, zm),
        section = case_when(
                c ==  Inf ~ 'early futility',
                c == -Inf ~ 'early efficacy',
                     TRUE ~ 'continue'
            )
    ) %>%
    select(
        -c2
    ) %>%
    mutate(
        PP = pmap_dbl(list(zm, n, c), ~predictive_power(..1, m, ..2, ..3, prior))
    )
tbl_designs <- bind_rows(
        tbl_optimal_design,
        tbl_adapted
    ) %>%
    mutate(
        design = factor(design, levels = c("original", "adapted", "optimal"))
    )
# first panel: sample size
plt1 <- ggplot(tbl_designs) +
    aes(zm, n, group = interaction(section, design)) +
    geom_line(aes(linetype = design)) +
    scale_x_continuous(expression(z[m]), limits = c(0, 3.5)) +
    scale_y_continuous('sample size', limits = c(0, n_max)) +
    scale_linetype('') +
    theme_bw() +
    theme(
        legend.position = "top"
    )
# second panel: critical value
plt2 <- ggplot(tbl_designs) +
    aes(zm, c, group = interaction(section, design)) +
    geom_line(aes(linetype = design)) +
    scale_x_continuous(expression(z[m]), limits = c(0, 3.5)) +
    scale_y_continuous('critical value', limits = c(0, NA)) +
    scale_linetype('') +
    theme_bw() +
    theme(
        legend.position = "top"
    )
# third panel: predictive power
plt3 <- ggplot(tbl_designs) +
    aes(zm, PP, group = interaction(section, design)) +
    geom_line(aes(linetype = design)) +
    scale_x_continuous(expression(z[m]), limits = c(0, 3.5)) +
    scale_y_continuous('predictive power') +
    scale_linetype('') +
    theme_bw() +
    theme(
        legend.position = "top"
    )
# compose plot
plt1 + plt2 + plt3 + plot_layout(guides = "collect") & theme(legend.position = 'top')
ggsave("output/figures/naive_adaptive_and_optimal_designs.pdf", width = 8, height = 3.5)



# Properties of optimal design ------------------------------------------------
show(
    summary(
        optimal_design,
         `E[N]` = ESS,
        `SD[N]` = SDSS,
        `MTOER` = TOER,
        `Expected Power` = EP
    )
)
