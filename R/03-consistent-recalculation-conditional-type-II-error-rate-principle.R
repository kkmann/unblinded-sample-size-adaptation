# Consistent Unplanned Recalculation Based on Predictive Power #################



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
        theta = theta1,
        alpha = alpha,
        beta  = beta,
        type  = 'two-stage',
        dist  = Normal(two_armed = FALSE),
        order = 11L
    )
)$design



# Plot adaptation at original interim analysis for various interim results -----
adapt <- function(priormean, zm, nstart = 15) {
    n_old <- adoptr::n(optimal_design, zm, round = FALSE)
    c_old <- critical_value(optimal_design, zm)
    if (!is.finite(c_old))
        return(list(
            n = n_old,
            c = c_old,
            multiroot = NULL
        ))
    # original conditional type-I error
    alpha_bar <- conditional_power(zm, m, n_old, c_old, 0)
    # original conditional type-II error
    beta_bar <- 1 - predictive_power(zm, m, n_old, c_old, prior)
    new_prior <- TruncatedNormal(priormean, 0.2, -0.5, 1.0)
    res <- suppressWarnings(rootSolve::multiroot(
        function(x) {
            # need to reparameterize to make sure n_bar > mbar
            n_bar <- x[1] + m + 1
            c_bar <- x[2]
            alpha <-     conditional_power(zm, m, n_bar, c_bar, 0)
            beta <- 1 -  predictive_power(zm, m, n_bar, c_bar, new_prior)
            c(alpha - alpha_bar, beta - beta_bar)
        },
        start = c(
            nstart + if (priormean < theta1) adoptr::n2(optimal_design, zm) else 0,
            critical_value(optimal_design, zm)
        ),
        maxiter = 50,
        positive = TRUE
    ))
    if (any(abs(res$f.root) > 1e-4) | res$root[1] >= n_max) {
      n_bar <- m
      c_bar <- Inf
    } else {
      n_bar <- res$root[1] + m + 1
      c_bar <- res$root[2]
    }
    list(
      n = n_bar,
      c = c_bar,
      multiroot = res
    )
}
# compute adapted values
m     <- n1(optimal_design)
n_max <-  160
# use these effects to compute zm as sqrt(m)*theta (observed z score)
thetahat <- c(.25, .3, .35)
tbl_adapted <- expand_grid(
                zm = sqrt(m)*thetahat,
        `new mean` = seq(-0.5, 1, length.out = 50)
    ) %>%
    mutate(
        tmp = map2(`new mean`, zm, ~{tmp <- adapt(..1, ..2); tibble(n = tmp$n, c = tmp$c)})
    ) %>%
    unnest(tmp) %>%
    mutate(
        section = if_else(c == Inf, "early futility", "continue"),
         effect = zm/sqrt(m), # recover observed effect
             PP = pmap_dbl(
                     list(effect, n, c, `new mean`),
                     ~predictive_power(sqrt(m)*., m = m, n = ..2, c = ..3, prior = TruncatedNormal(..4, 0.2, -0.5, 1.0))
                 ),
         effect = factor(effect)
    )
# properties of original (optimal) design for plotting
tbl_original <- tibble(
        effect = thetahat,
             n = map_dbl(effect, ~n(optimal_design, sqrt(m)*., round = FALSE)),
             c = map_dbl(effect, ~critical_value(optimal_design, sqrt(m)*.))
    ) %>%
    mutate(
        PP = pmap_dbl(
                list(effect, n, c),
                ~predictive_power(sqrt(m)*., m = m, n = ..2, c = ..3, prior = prior)
            ),
        effect = factor(effect)
    )
# first panel: sample size
plt1 <- ggplot(tbl_adapted) +
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
plt2 <- ggplot(tbl_adapted) +
    aes(`new mean`, c, group = interaction(section, effect)) +
    geom_hline(aes(yintercept = c, linetype = effect), color = "grey", data = tbl_original) +
    geom_vline(xintercept = 0.4, color = "grey") +
    geom_line(aes(linetype = effect)) +
    scale_x_continuous('new prior mean', expand = c(0, 0)) +
    scale_y_continuous('critical value', limits = c(1, 3)) +
    scale_linetype('observed effect') +
    theme_bw() +
    theme(
        legend.position = 'top',
        panel.grid = element_blank()
    )
# third panel: predictive power
plt3 <- ggplot(tbl_adapted) +
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
ggsave("../output/figures/recalculation-type-II-original-m.pdf", width = 8, height = 3.5)



# Plot adaptation at original interim analysis for various interim results -----
# calculate conditional power at earlier time point by averaging over potential
# outcomes
conditional_power_before_interim <- function(design, mbar, zmbar, theta) {
    m <- adoptr::n1(design)
    if ((0 >= mbar) | (mbar >= m)) stop()
    tau <- mbar/m
    # distribution of zm | zmbar is conditional normal with
    mu <- sqrt(m)*theta + sqrt(tau) * (zmbar - sqrt(mbar)*theta)
    sd <- sqrt(1 - tau)
    # everything to the left of the futilty bound is no rejection / does not contribute, start with [c1f, c1e]
      mid <- integrate_gl(
          function(zm) {
              conditional_power(zm, m, adoptr::n(design, zm, round = FALSE), critical_value(design, zm), theta) *
              dnorm(zm, mu, sd)
          },
          design@c1f, design@c1e
    )
    # last bit: [c1e, Inf)
    right <- pnorm(design@c1e, mu, sd, lower.tail = FALSE)
    # combine
    mid + right
}
# same for predictive power, need to integrate over 2-d grid
predictive_power_before_interim <- function(design, mbar, zmbar, prior) {
    m <- adoptr::n1(design)
    if ((0 >= mbar) | (mbar >= m)) stop()
    tau <- mbar/m
    # distribution of zm | zmbar is conditional normal with
    mu <- function(theta) sqrt(m)*theta + sqrt(tau) * (zmbar - sqrt(mbar)*theta)
    sd <- sqrt(1 - tau)
    # everything to the left of the futility bound is no rejection / does not contribute, start with [c1f, c1e]
    mid <- integrate_gl(
        function(theta) {
            integrate_gl(
                function(zm) {
                    # conditional power at z_m and time m (future)
                    conditional_power(zm, m, adoptr::n(design, zm, round = FALSE), critical_value(design, zm), theta) *
                    # distribution of Z_m
                    dnorm(zm, mu(theta), sd) *
                    # distribution of Theta|Z_mbar = z_mbar
                    probability_density(posterior(condition(prior, 0), mbar, zmbar), theta)
                },
                design@c1f, design@c1e
            )
        },
        prior$a, prior$b
    )
    # last bit: [c1e, Inf)
    right <- integrate_gl(
        function(theta) {
            pnorm(design@c1e, mu(theta), sd, lower.tail = FALSE) *
            probability_density(posterior(condition(prior, 0), mbar, zmbar), theta)
        },
        prior$a, prior$b
    )
    # combine
    mid + right
}
# define adaptation rule
adapt_before <- function(priormean, mbar, zmbar, nstart = 10) {
    # original conditional type-I error
    alpha_bar <- conditional_power_before_interim(optimal_design, mbar, zmbar, 0)
    # original conditional type-II error
    beta_bar <- 1 - predictive_power_before_interim(optimal_design, mbar, zmbar, prior)
    new_prior <- TruncatedNormal(priormean, 0.2, -0.5, 1.0)
    res <- suppressWarnings(rootSolve::multiroot(
        function(x) {
            # need to reparameterize to make sure n_bar > mbar
            n_bar <- x[1] + mbar + 1
            c_bar <- x[2]
            alpha <-     conditional_power(zmbar, mbar, n_bar, c_bar, 0)
             beta <- 1 -  predictive_power(zmbar, mbar, n_bar, c_bar, new_prior)
            c(alpha - alpha_bar, beta - beta_bar)
        },
        start = c(
            nstart + if (priormean < theta1) adoptr::n2(optimal_design, zmbar) else 0,
            critical_value(optimal_design, zmbar)
        ),
        maxiter = 50,
        positive = TRUE
    ))
    if (any(abs(res$f.root) > 1e-4) | res$root[1] >= n_max) {
        n_bar <- mbar
        c_bar <- Inf
    } else {
        n_bar <- res$root[1] + mbar + 1
        c_bar <- res$root[2]
    }
    list(
        n = n_bar,
        c = c_bar,
        multiroot = res
    )
}
# computed adapted values
tbl_adapted_before <- expand_grid(
              mbar = c(seq(15, 30, by = 5), 34),
        `new mean` = seq(-0.5, 1, length.out = 50)
    ) %>%
    mutate(
        zmbar = sqrt(mbar)*0.25,
          tmp = pmap(
              list(`new mean`, mbar, zmbar),
              function(mu, mbar, zmbar) future({
                      adapt_before(mu, mbar, zmbar)
                  },
                  # future is struggling with passing methods defined in util.R
                  # specify manually
                  globals = structure(TRUE, add = c(
                      "condition", "condition.TruncatedNormal",
                      "posterior", "posterior.TruncatedNormal",
                      "probability_density", "probability_density.TruncatedNormal"
                  ))
              )
          )
    ) %>%
    mutate(
        tmp = map(
            tmp,
            ~tibble(n = value(.)$n, c = value(.)$c)
        )
    ) %>%
    unnest(tmp) %>%
    mutate(
        section = if_else(c == Inf, "early futility", "continue"),
        effect = zmbar/sqrt(mbar), # recover observed effect
        PP = pmap_dbl(
            list(effect, mbar, n, c, `new mean`),
            ~predictive_power(sqrt(..2)*..1, m = ..2, n = ..3, c = ..4, prior = TruncatedNormal(..5, 0.2, -0.5, 1.0))
        ),
        effect = factor(effect),
          mbar = factor(mbar)
    )
# compute properties of original design at m = 35
tbl_original <- tibble(
    effect = 0.25,
      mbar = 35,
         n = map_dbl(effect, ~n(optimal_design, sqrt(mbar)*.)),
         c = map_dbl(effect, ~critical_value(optimal_design, sqrt(mbar)*.))
    ) %>%
    mutate(
      PP = pmap_dbl(
        list(effect, n, c),
        ~predictive_power(sqrt(35)*..1, m = 35, n = ..2, c = ..3, prior = prior)
      ),
      effect = factor(effect),
        mbar = factor(mbar)
    )
# first panel: sample size
plt1 <- ggplot(tbl_adapted_before) +
    aes(`new mean`, n, group = interaction(section, mbar)) +
    #geom_hline(aes(yintercept = n), linetype = 2, data = tbl_original) +
    geom_vline(xintercept = theta1, linetype = 2) +
    geom_line(aes(color = mbar)) +
    scale_x_continuous('new prior mean', expand = c(0, 0)) +
    scale_y_continuous('sample size', limits = c(0, n_max)) +
    scale_color_grey('interim sample size') +
    theme_bw() +
    theme(
        legend.position = 'top',
        panel.grid = element_blank()
    )
# second panel: critical value
plt2 <- ggplot(tbl_adapted_before) +
    aes(`new mean`, c, group = interaction(section, mbar)) +
    geom_hline(aes(yintercept = c), linetype = 2, data = tbl_original) +
    geom_vline(xintercept = theta1, linetype = 2) +
    geom_line(aes(color = mbar)) +
    scale_x_continuous('new prior mean', expand = c(0, 0)) +
    scale_y_continuous('critical value', limits = c(1, 3)) +
    scale_color_grey('interim sample size') +
    theme_bw() +
    theme(
        legend.position = 'top',
        panel.grid = element_blank()
    )
# third panel: predictive power
plt3 <- ggplot(tbl_adapted_before) +
    aes(`new mean`, PP, group = interaction(section, mbar)) +
    geom_hline(aes(yintercept = PP), linetype = 2, data = tbl_original) +
    geom_vline(xintercept = theta1, linetype = 2) +
    geom_line(aes(color = mbar)) +
    scale_x_continuous('new prior mean', expand = c(0, 0)) +
    scale_y_continuous('predictive power', limits = c(.8, 1)) +
    scale_color_grey('interim sample size') +
    theme_bw() +
    theme(
        legend.position = 'top',
        panel.grid = element_blank()
    )
# compose plot
plt1 + plt2 + plt3 + plot_layout(guides = "collect") & theme(legend.position = 'top')
ggsave("../output/figures/recalculation-type-II-earlier.pdf", width = 8, height = 3.5)
