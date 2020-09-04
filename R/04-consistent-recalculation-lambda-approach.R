# Consistent Recalculation, lambda-approach ####################################



# Setup ------------------------------------------------------------------------
suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
    library(nloptr)
    library(adoptr)
    library(future)
})

plan(multisession)

set.seed(42)
source("util.R")
dir.create("../output/figures", recursive = TRUE, showWarnings = FALSE)



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



# Compute optimal two-stage design ---------------------------------------------
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



# Determine 'Lagrange-multiplier' lambda ---------------------------------
# function to determine optimal design given lambda
f <- function(lambda) {
    minimize(
        composite({ESS - lambda*EP}),
        subject_to(
            TOER <= alpha
        ),
        get_initial_design(
            theta = theta1,
            alpha = alpha,
            beta  = beta,
            type  = 'two-stage',
            dist  = Normal(two_armed = FALSE),
            order = 5L # reduced precision, we only neeed to get a good proxy for the power to match, faster!
        ),
        opts = list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-05, maxeval = 5e4)
    )$design
}
# evaluate on grid (smoothness unknown!) to pick lambda that gives
# target power of 1 - beta
tbl_lambda <- tibble(
        lambda = seq(200, 300, by = 5),
        design = map(lambda, ~future({f(.)})) # spawn futures
    ) %>%
    mutate(
        design = map(design, value), # collect futures
         delta = map_dbl(design, ~abs(evaluate(EP, .) - (1 - beta)))
    )
# diagnostic plot, not saved
tbl_lambda %>%
    ggplot() +
    aes(lambda, delta) +
    geom_line() +
    geom_point() +
    ylab("power difference")
cat("original optimal design:\n\r")
summary(
        optimal_design,
        `E[N]` = ESS,
        `SD[N]` = SDSS,
        `MTOER` = TOER,
        `Expected Power` = EP
    ) %>%
    show()
# determine optimal lambda
lambda <- tbl_lambda%>% arrange(delta) %>% pull(lambda) %>% head(1)
cat(sprintf("matched-lambda optimal design, lambda=%.2f:\n\r", lambda))
tbl_lambda %>%
    arrange(delta) %>%
    pull(design) %>%
    .[[1]] %>%
    summary(
         `E[N]` = ESS,
        `SD[N]` = SDSS,
        `MTOER` = TOER,
        `Expected Power` = EP
    ) %>%
    show()

# Plot adaptations under lambda rule -------------------------------------------
    m <- 35
n_max <- 160
adapt_lambda <- function(priormean, zm) {
    n_old <- n(optimal_design, zm, round = FALSE)
    c_old <- critical_value(optimal_design, zm)
    if (!is.finite(c_old))
        return(list(
            n = n_old,
            c = c_old,
            optim = NULL
        ))
    # old conditional type-I error
    alpha_bar <- conditional_power(zm, m, n_old, c_old, 0)
    new_prior <- TruncatedNormal(priormean, 0.2, -0.5, 1.0)
    objective <- function(x) {
        n_bar <- x[1]; c_bar <- x[2]
        n_bar - lambda*predictive_power(zm, m = m, n = n_bar, c = c_bar, prior = new_prior)
    }
    conditional_error_constraint <- function(x) {
        n_bar <- x[1]; c_bar <- x[2]
        conditional_power(zm, m = m, n = n_bar, c = c_bar, theta = 0) - alpha_bar # <= 0
    }
    res <- nloptr::nloptr(
                 x0 = c(n_old, c_old),
             eval_f = objective,
        eval_g_ineq = conditional_error_constraint,
                 lb = c(m, 0),
                 ub = c(n_max, 5),
        opts = list(
            algorithm = "NLOPT_LN_COBYLA",
             xtol_rel = 1e-4,
              maxeval = 1e5
        )
    )
    if (res$status != 4) stop(res$message)
    n_bar <- if (round(res$solution[1]) >= n_max)   m else res$solution[1]
    c_bar <- if (round(res$solution[1]) >= n_max) Inf else res$solution[2]
    n_bar <- if (round(res$solution[1]) <= m)   m else res$solution[1]
    c_bar <- if (round(res$solution[1]) <= m) Inf else res$solution[2]
    list(
            n = n_bar,
            c = c_bar,
        optim = res
    )
}
# observed effects to look at
thetas <- c(.3, .35, .4)
tbl_adapted_lambda <- expand_grid(
                zm = sqrt(m)*thetas,
        `new mean` = seq(-0.5, 1, length.out = 100)
    ) %>%
    mutate(
        tmp = map2(`new mean`, zm, ~{tmp <- adapt_lambda(..1, ..2); tibble(n = tmp$n, c = tmp$c)})
    ) %>%
    unnest(tmp) %>%
    mutate(
        section = if_else(c == Inf, "early futility", "continue"),
         effect = factor(zm/sqrt(m)),
             PP = pmap_dbl(
                 list(zm, n, c, `new mean`),
                 ~predictive_power(zm = ..1, m = m, n = ..2, c = ..3, prior = TruncatedNormal(..4, 0.2, -0.5, 1.0))
             )
    )
# original performance
tbl_original <- tibble(
        effect = thetas,
             n = map_dbl(effect, ~n(optimal_design, sqrt(m)*.)),
             c = map_dbl(effect, ~critical_value(optimal_design, sqrt(m)*.)),
            PP = pmap_dbl(
                list(effect, n, c),
                ~predictive_power(sqrt(m)*., m = m, n = ..2, c = ..3, prior = prior)
            )
    ) %>%
    mutate(
        effect = factor(effect)
    )
# first panel: sample size
plt1 <- ggplot(tbl_adapted_lambda) +
    aes(`new mean`, n, group = interaction(section, effect)) +
    geom_hline(aes(yintercept = n, linetype = effect), color = "grey", data = tbl_original) +
    geom_vline(xintercept = 0.4, color = "grey") +
    geom_line(aes(linetype = effect)) +
    scale_x_continuous('new prior mean', expand = c(0, 0)) +
    scale_y_continuous('sample size', limits = c(0, n_max)) +
    scale_linetype('observed effect') +
    theme_bw() +
    theme(
        legend.position = 'top',
        panel.grid = element_blank()
    )
# second panel: critical value
plt2 <- ggplot(tbl_adapted_lambda) +
    aes(`new mean`, c, group = interaction(section, effect)) +
    geom_hline(aes(yintercept = c, linetype = effect), color = "grey", data = tbl_original) +
    geom_vline(xintercept = 0.4, color = "grey") +
    geom_line(aes(linetype = effect)) +
    scale_x_continuous('new prior mean', expand = c(0, 0)) +
    scale_y_continuous('critical value', limits = c(0, 3)) +
    scale_linetype('observed effect') +
    theme_bw() +
    theme(
        legend.position = 'top',
        panel.grid = element_blank()
    )
# third panel: predictive power size
plt3 <- ggplot(tbl_adapted_lambda) +
    aes(`new mean`, PP, group = interaction(section, effect)) +
    geom_hline(aes(yintercept = PP, linetype = effect), color = "grey", data = tbl_original) +
    geom_vline(xintercept = 0.4, color = "grey") +
    geom_line(aes(linetype = effect)) +
    scale_x_continuous('new prior mean', expand = c(0, 0)) +
    scale_y_continuous('predictive power', limits = c(0, 1)) +
    scale_linetype('observed effect') +
    theme_bw() +
    theme(
        legend.position = 'top',
        panel.grid = element_blank()
    )
# compose plot
plt1 + plt2 + plt3 + plot_layout(guides = "collect") & theme(legend.position = 'top')
ggsave("../output/figures/recalculation-lambda-approach-original-m.pdf", width = 8, height = 3.5)
