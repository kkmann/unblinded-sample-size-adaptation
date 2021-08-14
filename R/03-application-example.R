# Application example ##########################################################



# Setup ------------------------------------------------------------------------
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(patchwork)
library(nloptr)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# two-arm version
conditional_power_two_arm <- function(theta_obs, m, n, c, theta) {
    if (any(m > n)) stop("m must be smaller than n")
    if (any(m < 2)) stop("m must be larger than 2")
    zm <- sqrt(m/4)*theta_obs
    tau <- m/n
    mu <- sqrt(n/4)*theta + sqrt(tau)*(zm - sqrt(m/4)*theta)
    sd <- sqrt(1 - tau)
    pnorm(c, mu, sd, lower.tail = FALSE)
}

# Define planning assumptions --------------------------------------------------
theta1 <- 0.4
alpha <- 0.025
beta <- 0.2

# Compute GS design ------------------------------------------------------------
get_optimal_gs_design <- function(
    theta1, theta0 = 0, alpha = .025, beta = .2,
    c1e_ = NULL, c1f_ = NULL, cc_ = NULL, n1_ = NULL, n2_ = NULL, #n/n1 are overall
    nstart = 100, x1_obs = NULL, ce_obs = NULL, nmax = 400
) {
    c1f <- function(x) if (!is.null(c1f_)) c1f_ else x["c1f"]
    c1e <- function(x) if (!is.null(c1e_)) c1e_ else x[["c1e"]]
    cc <- function(x) if (!is.null(cc_)) cc_ else x["cc"]
    n1 <- function(x) if (!is.null(n1_)) n1_ else x["n1"]
    n <- function(x) n1(x) + if (!is.null(n2_)) n2_ else x["n2"]
    named <- function(x) {
        nm <- c()
        if (is.null(c1f_)) nm <- c(nm, "c1f")
        if (is.null(c1e_)) nm <- c(nm, "c1e")
        if (is.null(cc_)) nm <- c(nm, "cc")
        if (is.null(n1_)) nm <- c(nm, "n1")
        if (is.null(n2_)) nm <- c(nm, "n2")
        names(x) <- nm
        x
    }
    cdf <- function(z1, x, theta) {
        x <- named(x)
        pnorm(z1, mean = sqrt(n1(x)/4)*theta, sd = 1)
    }
    pdf <- function(z1, x, theta) {
        x <- named(x)
        dnorm(z1, mean = sqrt(n1(x)/4)*theta, sd = 1)
    }
    ess <- function(x, theta) {
        x <- named(x)
        pcontinue <- cdf(c1e(x), x, theta) - cdf(c1f(x), x, theta)
        pcontinue*n(x) + (1 - pcontinue)*n1(x)
    }
    conditional_power <- function(z1, x, theta) {
        x <- named(x)
        tau <- n1(x)/n(x)
        mu <- sqrt(n(x)/4)*theta + sqrt(tau)*(z1 - sqrt(n1(x)/4)*theta)
        sd <- sqrt(1 - tau)
        pnorm(cc(x), mu, sd, lower.tail = FALSE)
    }
    power <- function(x, theta) {
        x <- named(x)
        integrate(
            function(z1) pdf(z1, x, theta) * conditional_power(z1, x, theta),
            c1f(x),
            c1e(x)
        )$value + (1 - cdf(c1e(x), x, theta))
    }
    x0 <- c()
    if (is.null(c1f_)) x0 <- c(x0, c1f = 0)
    if (is.null(c1e_)) x0 <- c(x0, c1e = 3)
    if (is.null(cc_)) x0 <- c(x0, cc = 2)
    if (is.null(n1_)) x0 <- c(x0, n1 = nstart)
    if (is.null(n2_)) x0 <- c(x0, n2 = nstart)
    lb <- c()
    if (is.null(c1f_)) lb <- c(lb, c1f = -1)
    if (is.null(c1e_)) lb <- c(lb, c1e = 0)
    if (is.null(cc_)) lb <- c(lb, cc = 1)
    if (is.null(n1_)) lb <- c(lb, n1 = 20)
    if (is.null(n2_)) lb <- c(lb, n2 = 20)
    ub <- c()
    if (is.null(c1f_)) ub <- c(ub, c1f = 1)
    if (is.null(c1e_)) ub <- c(ub, c1e = 5)
    if (is.null(cc_)) ub <- c(ub, cc = 4)
    if (is.null(n1_)) ub <- c(ub, n1 = nmax)
    if (is.null(n2_)) ub <- c(ub, n2 = nmax)
    res <- nloptr::nloptr(
                 x0 = x0,
             eval_f = function(x) ess(x, theta1),
        eval_g_ineq = function(x) {
                x <- named(x)
                res <- c(
                    power(x, theta0) - alpha, # toer constraint
                    1 - power(x, theta1) - beta, # power constraint
                    x[1] - x[2], # c1f < c1e
                    n(x) - nmax
                )
                if (!is.null(x1_obs)) { # conditional error constraint
                    res <- c(res, conditional_power(x1_obs, x, theta0) - ce_obs)
                }
                res
            },
        lb = lb,
        ub = ub,
        opts = list(
            algorithm = "NLOPT_LN_COBYLA",
             xtol_rel = 1e-4,
              maxeval = 1e4
        )
    )
    named(res$solution)
}

initial_design <- c(c1f = 0, get_optimal_gs_design(theta1, alpha = alpha, beta = beta, c1f = 0)) # fix conservative futility boundary

gs_design <- list(
    name = "original GS",
    m = as.numeric(initial_design["n1"]),
    lower = as.numeric(initial_design["c1f"] / sqrt(initial_design["n1"]/4)), # convert to treatment effect scale
    upper = as.numeric(initial_design["c1e"] / sqrt(initial_design["n1"]/4)),
    n = function(theta_obs, ...) {
        x1_obs <- theta_obs * sqrt(initial_design["n1"]/4)
        case_when(
            x1_obs < initial_design["c1f"] ~ initial_design["n1"],
            between(x1_obs, initial_design["c1f"], initial_design["c1e"]) ~ initial_design["n1"] + initial_design["n2"],
            x1_obs > initial_design["c1e"] ~ initial_design["n1"]
        ) %>% as.numeric
    },
    c = function(theta_obs, ...) {
        x1_obs <- theta_obs * sqrt(initial_design["n1"]/4)
        case_when(
            x1_obs < initial_design["c1f"] ~ Inf,
            between(x1_obs, initial_design["c1f"], initial_design["c1e"]) ~ initial_design["cc"],
            x1_obs > initial_design["c1e"] ~ -Inf
        ) %>% as.numeric
    },
    decision = function(theta_obs, ...) {
        x1_obs <- theta_obs * sqrt(initial_design["n1"]/4)
        case_when(
            x1_obs < initial_design["c1f"] ~ "early futility",
            between(x1_obs, initial_design["c1f"], initial_design["c1e"]) ~ "continue",
            x1_obs > initial_design["c1e"] ~ "early efficacy"
        )
    }
)


# adapt under fixed predictive power rule --------------------------------------
adapt_original_pp_ <- function(theta_new, theta_obs, design, nstart_factor = 200, ...) {
    n_old <- design$n(theta_obs)
    c_old <- design$c(theta_obs)
    m = design$m
    if (!is.finite(c_old))
        return(list(n = n_old, c = c_old))
    # original conditional type-I error
    alpha_bar <- conditional_power_two_arm(theta_obs, m, n_old, c_old, 0)
    # original conditional type-II error
    beta_bar <- 1 - conditional_power_two_arm(theta_obs, m, n_old, c_old, theta1)
    nstart = max(20, n_old - m + nstart_factor*(theta1 - theta_obs))
    res <- rootSolve::multiroot(
        function(x) {
            # need to reparameterize to make sure n_bar > mbar
            n_bar <- x[1] + m + 2
            c_bar <- x[2] + 1.5
            alpha <-    conditional_power_two_arm(theta_obs, m, n_bar, c_bar, 0)
            beta <- 1 - conditional_power_two_arm(theta_obs, m, n_bar, c_bar, theta_new)
            c(alpha - alpha_bar, beta - beta_bar)
        },
        start = c(nstart, c_old - 1.5),
        maxiter = 50,
        positive = TRUE
    )
    n_bar <- res$root[1] + m + 2
    c_bar <- res$root[2] + 1.5
    if (n_bar > 400 | (n_bar < m + 40 & theta_obs < theta1)) {
        n_bar <- 400
        # resolve with fixed sample size
        c_bar <- uniroot(
          function(c_bar) {
                alpha <- conditional_power_two_arm(theta_obs, m, n_bar, c_bar, 0)
                c(alpha - alpha_bar)
            },
            lower = 0.5*c_old, upper = 3*c_old
        )$root
    }
    if (n_bar < m + 40 | (n_bar >= 400 & theta_obs > theta1)) {
        n_bar <- m + 40
        # resolve with fixed sample size
        c_bar <- uniroot(
            function(c_bar) {
                alpha <- conditional_power_two_arm(theta_obs, m, n_bar, c_bar, 0)
                c(alpha - alpha_bar)
            },
            lower = 0.5*c_old, upper = 3*c_old
        )$root
    }
    list(n = n_bar, c = c_bar)
}
adapt_original_pp <- memoise::memoise(adapt_original_pp_) # expensive, memoise!
fixed_pp_design <- list(
    name = "not-worse",
    m = gs_design$m,
    n = function(theta_obs, theta_new, ...) {
        adapt_original_pp(theta_new, theta_obs, gs_design, ...)$n
    },
    c = function(theta_obs, theta_new, ...) {
        adapt_original_pp(theta_new, theta_obs, gs_design, ...)$c
    },
    decision = function(theta_obs, ...) {
        c_old <- gs_design$c(theta_obs)
        case_when(
            c_old == Inf ~ "early futility",
            is.finite(c_old) ~ "continue",
            c_old == -Inf ~ "early efficacy"
        )
    }
)

# adapt under constant predictive power rule --------------------------------------
adapt_fixed_pp_ <- function(theta_new, theta_obs, design, nstart_factor = 1000) {
    n_old <- design$n(theta_obs)
    c_old <- design$c(theta_obs)
    m = design$m
    if (!is.finite(c_old))
        return(tibble(n = n_old, c = c_old))
    # original conditional type-I error
    alpha_bar <- conditional_power_two_arm(theta_obs, m, n_old, c_old, 0)
    # fixed
    beta_bar <- .2
    nstart = max(20, n_old - m + nstart_factor*(theta1 - theta_obs))
    cstart = c_old + .25*c_old*max(0, theta_obs/.6)
    res <- rootSolve::multiroot(
        function(x) {
            # need to reparameterize to make sure n_bar > mbar
            n_bar <- x[1] + m + 2
            c_bar <- x[2]
            alpha <-    conditional_power_two_arm(theta_obs, m, n_bar, c_bar, 0)
            beta <- 1 - conditional_power_two_arm(theta_obs, m, n_bar, c_bar, theta_new)
            c(alpha - alpha_bar, beta - beta_bar)
        },
        start = c(nstart, cstart),
        maxiter = 100,
        positive = TRUE
    )
    n_bar <- res$root[1] + m + 2
    c_bar <- res$root[2]
    if (n_bar > 400) {
        n_bar <- 400
        # resolve with fixed sample size
        c_bar <- uniroot(
        function(c_bar) {
                alpha <- conditional_power_two_arm(theta_obs, m, n_bar, c_bar, 0)
                c(alpha - alpha_bar)
            },
            lower = c_old/5, upper = 5*c_old,
        )$root
    }
    if (n_bar < m + 40 | (n_bar >= 400 & theta_obs > theta1)) {
        n_bar <- m + 40
        # resolve with fixed sample size
        c_bar <- uniroot(
            function(c_bar) {
                alpha <- conditional_power_two_arm(theta_obs, m, n_bar, c_bar, 0)
                c(alpha - alpha_bar)
            },
            lower = c_old/2, upper = 2*c_old,
        )$root
    }
    list(n = n_bar, c = c_bar)
}
adapt_fixed_pp <- memoise::memoise(adapt_fixed_pp_) # expensive, memoise!
pp_geq_08_design <- list(
    name = "PP >= 0.8",
    m = gs_design$m,
    n = function(theta_obs, theta_new, ...) {
        adapt_fixed_pp(theta_new, theta_obs, gs_design, ...)$n
    },
    c = function(theta_obs, theta_new, ...) {
        adapt_fixed_pp(theta_new, theta_obs, gs_design, ...)$c
    },
    decision = function(theta_obs, ...) {
        c_old <- gs_design$c(theta_obs)
        case_when(
            c_old == Inf ~ "early futility",
            is.finite(c_old) ~ "continue",
            c_old == -Inf ~ "early efficacy"
        )
    }
)


# reoptimize ######
adapt_reoptimize_ <- function(theta_new, theta_obs, design, n_max = 4*design$m, n2min = 20, ...) {
    n_old <- design$n(theta_obs)
    c_old <- design$c(theta_obs)
    if (is.na(c_old))
        browser()
    m <- design$m
    if (!is.finite(c_old))
        return(list(n = n_old, c = c_old))
    # original conditional type-I error
    alpha_bar <- conditional_power_two_arm(theta_obs, m, n_old, c_old, 0)
    # compute optional GS design under CE constraint
    new_design <- get_optimal_gs_design(
        theta_new,
        alpha = alpha, beta = beta,
        c1f_ = initial_design["c1f"],
        c1e_ = initial_design["c1e"],
        n1_ = initial_design["n1"],
        x1_obs = theta_obs * sqrt(design$m/4),
        ce_obs = alpha_bar
    )
    cc <- as.numeric(new_design["cc"])
    nn <- as.numeric(new_design["n2"]) + m
    # corresponding predictive power
    predictive_power_bar <- conditional_power_two_arm(theta_obs, m, nn, cc, theta_new)
    # solve for n/c in case conditional error is not binding
    # conditional power constraint of new design
    g <- function(x) {n <- x[1] + m + n2min; c <- x[2]; conditional_power_two_arm(theta_obs, m, n, c, theta_new) - predictive_power_bar}
    # conditional error constraint of new design
    h <- function(x) {n <- x[1] + m + n2min; c <- x[2]; alpha_bar - conditional_power_two_arm(theta_obs, m, n, c, 0)}
    # solve
    res <- rootSolve::multiroot(
        function(x) c(g(x), h(x)),
        start = c(nn - m, cc),
        maxiter = 100,
        positive = TRUE
    )
    n_bar <- as.numeric(res$root[1] + m + n2min)
    c_bar <- as.numeric(res$root[2])
    if (n_bar > 400) {
        n_bar <- 400
        # resolve with fixed sample size
        c_bar <- uniroot(
        function(c_bar) {
                alpha <- conditional_power_two_arm(theta_obs, m, n_bar, c_bar, 0)
                c(alpha - alpha_bar)
            },
            lower = c_old/5, upper = 5*c_old,
        )$root
    }
    list(n = n_bar, c = c_bar)
}
adapt_reoptimize <- memoise::memoise(adapt_reoptimize_) # expensive, memoise!
reoptimize_design <- list(
    name = "reoptimize",
    m = gs_design$m,
    n = function(theta_obs, theta_new, ...) {
        adapt_reoptimize(theta_new, theta_obs, gs_design)$n
    },
    c = function(theta_obs, theta_new, ...) {
        adapt_reoptimize(theta_new, theta_obs, gs_design)$c
    },
    decision = function(theta_obs, ...) {
        c_old <- gs_design$c(theta_obs)
        case_when(
            c_old == Inf ~ "early futility",
            is.finite(c_old) ~ "continue",
            c_old == -Inf ~ "early efficacy"
        )
    }
)

# change ##########

tbl_design <- expand_grid(
        theta_obs = seq(-.1, .6, by = .01),
        design = list(
            pp_geq_08_design,
            fixed_pp_design,
            reoptimize_design,
            gs_design
        ),
        theta_new = c(.35, .45)
    ) %>%
    mutate(
        nstart_factor = case_when(
            map_chr(design, ~.$name) == "PP >= 0.8" & theta_new == .5 ~ 100,
            map_chr(design, ~.$name) == "not-worse" & theta_new == .5 ~ 100,
            map_chr(design, ~.$name) == "PP >= 0.8" ~ 200,
            map_chr(design, ~.$name) == "not-worse" ~ 100,
            TRUE ~ NA_real_
        ),
        n = pmap_dbl(list(design, theta_obs, theta_new, nstart_factor), ~..1$n(..2, ..3, nstart_factor = ..4)),
        c = pmap_dbl(list(design, theta_obs, theta_new, nstart_factor), ~..1$c(..2, ..3, nstart_factor = ..4)),
        decision = map2_chr(design, theta_obs, ~..1$decision(..2)),
        design = map_chr(design, ~.$name) %>%
            factor(
                levels = c("original GS", "PP >= 0.8", "not-worse", "reoptimize")
            ),
        `predictive power` = pmap_dbl(list(theta_obs, n, c, theta_new), ~conditional_power_two_arm(..1, m = gs_design$m, ..2, ..3, ..4))
    ) %>%
    rename(
        `changed alternative` = theta_new,
        `sample size` = n
    ) %>%
    pivot_longer(c(`sample size`, `predictive power`))

p1 <- tbl_design %>%
    ggplot() +
        aes(theta_obs, value, color = design, group = interaction(design, decision)) +
        geom_line() +
        labs(x = "interim observed effect") +
        scale_color_discrete("") +
        scale_y_continuous("", position = "right", limits = c(0, NA_real_)) +
        facet_grid(
            name ~ `changed alternative`,
            scales = "free_y",
            labeller = labeller(.cols = label_both),
            switch = "y"
        ) +
        theme_bw() +
        theme(
            legend.position = "top"
        )

ess <- function(design, theta, theta_new, ...) {
    gg <- Vectorize(function(x) design$n(x, theta_new, ...))
    pdf <- function(theta_obs) dnorm(theta_obs, mean = theta, sd = sqrt(4/gs_design$m))
    cdf <- function(theta_obs) pnorm(theta_obs, mean = theta, sd = sqrt(4/gs_design$m))
    eps <-  sqrt(.Machine$double.eps)
    left <- gg(gs_design$lower - eps) * cdf(gs_design$lower)
    right <- gg(gs_design$upper + eps) * (1 - cdf(gs_design$upper))
    mid <- integrate(function(x) gg(x) * pdf(x), lower = gs_design$lower + eps, upper = gs_design$upper - eps, rel.tol = .001, subdivisions = 10000L)$value
    left + mid + right
}
ep <- function(design, theta, theta_new, ...) {
    gg <- Vectorize(function(x) {
        conditional_power_two_arm(x, m = gs_design$m, n = design$n(x, theta_new, ...), c = design$c(x, theta_new, ...), theta)
    })
    pdf <- function(theta_obs) dnorm(theta_obs, mean = theta, sd = sqrt(4/gs_design$m))
    cdf <- function(theta_obs) pnorm(theta_obs, mean = theta, sd = sqrt(4/gs_design$m))
    eps <-  sqrt(.Machine$double.eps)
    left <- gg(gs_design$lower - eps) * cdf(gs_design$lower)
    right <- gg(gs_design$upper + eps) * (1 - cdf(gs_design$upper))
    mid <- integrate(function(x) gg(x) * pdf(x), lower = gs_design$lower + eps, upper = gs_design$upper - eps, rel.tol = .001, subdivisions = 10000L)$value
    left + mid + right
}

tbl_unconditional <- expand_grid(
        design = list(
            pp_geq_08_design,
            fixed_pp_design,
            reoptimize_design,
            gs_design
        ),
        theta = c(.35, .45, theta1)
    ) %>%
    mutate(
        nstart_factor = case_when(
            map_chr(design, ~.$name) == "PP >= 0.8" & theta == .45 ~ 100,
            map_chr(design, ~.$name) == "not-worse" & theta == .45 ~ 100,
            map_chr(design, ~.$name) == "PP >= 0.8" ~ 200,
            map_chr(design, ~.$name) == "not-worse" ~ 100,
            TRUE ~ NA_real_
        ),
        ess = pmap_dbl(list(design, theta, theta, nstart_factor), ~ess(..1, ..2, ..3, ..4)),
        ep = pmap_dbl(list(design, theta, theta, nstart_factor), ~ep(..1, ..2, ..3, ..4)),
        design = map_chr(design, ~.$name) %>%
            factor(
                levels = c("original GS", "PP >= 0.8", "not-worse", "reoptimize")
            ),
        theta = factor(theta)
    )

p2 <- tbl_unconditional %>%
    ggplot() +
        aes(ess, ep, shape = theta, color = design) +
        geom_point(size = 2) +
        theme_bw() +
        scale_y_continuous(breaks = seq(0, 1, by = .2), limits = c(0, 1)) +
        #xlim(0, NA_real_) +
        labs(y = "(expected) power", x = "expected sample size") +
        scale_shape_discrete(expression(theta)) +
        scale_color_discrete("") +
        guides(color = FALSE) +
        theme(
            legend.position = "top"
        )

p1 + p2 +
    plot_layout(width = c(1.5, 1)) &
    theme(legend.direction = 'horizontal')


ggsave("output/figures/application-example.pdf", width = 8, height = 4.5)


# response adaptive design

tbl_design <- expand_grid(
        theta_obs = seq(-.1, .6, by = .01),
        design = list(
            pp_geq_08_design,
            fixed_pp_design,
            reoptimize_design,
            gs_design
        )
    ) %>%
    mutate(
        theta_new = (.4 + theta_obs)/2,
        nstart_factor = case_when(
            map_chr(design, ~.$name) == "PP >= 0.8" ~ 500,
            map_chr(design, ~.$name) == "not-worse" ~ 500,
            TRUE ~ NA_real_
        ),
        n = pmap_dbl(list(design, theta_obs, theta_new, nstart_factor), ~..1$n(..2, ..3, nstart_factor = ..4)),
        c = pmap_dbl(list(design, theta_obs, theta_new, nstart_factor), ~..1$c(..2, ..3, nstart_factor = ..4)),
        decision = map2_chr(design, theta_obs, ~..1$decision(..2)),
        design = map_chr(design, ~.$name) %>%
            factor(
                levels = c("original GS", "PP >= 0.8", "not-worse", "reoptimize")
            ),
        `predictive power` = pmap_dbl(list(theta_obs, n, c, theta_new), ~conditional_power_two_arm(..1, m = gs_design$m, ..2, ..3, ..4))
    ) %>%
    rename(
        `changed alternative` = theta_new,
        `sample size` = n
    ) %>%
    pivot_longer(c(`sample size`, `predictive power`))

p1 <- tbl_design %>%
    ggplot() +
        aes(theta_obs, value, color = design, group = interaction(design, decision)) +
        geom_line() +
        labs(x = "interim observed effect") +
        scale_color_discrete("") +
        scale_y_continuous("", limits = c(0, NA_real_), position = "right") +
        facet_wrap(~name, scales = "free_y", strip.position = "left") +
        theme_bw() +
        theme(
            legend.position = "top"
        )

ess <- function(design, theta, ...) {
    gg <- Vectorize(function(x) design$n(x, (.4 + x)/2, ...))
    pdf <- function(theta_obs) dnorm(theta_obs, mean = theta, sd = sqrt(4/gs_design$m))
    cdf <- function(theta_obs) pnorm(theta_obs, mean = theta, sd = sqrt(4/gs_design$m))
    eps <-  sqrt(.Machine$double.eps)
    left <- gg(gs_design$lower - eps) * cdf(gs_design$lower)
    right <- gg(gs_design$upper + eps) * (1 - cdf(gs_design$upper))
    mid <- integrate(function(x) gg(x) * pdf(x), lower = gs_design$lower + eps, upper = gs_design$upper - eps, rel.tol = .001, subdivisions = 10000L)$value
    left + mid + right
}
ep <- function(design, theta, ...) {
    gg <- Vectorize(function(x) {
        conditional_power_two_arm(x, m = gs_design$m, n = design$n(x, (.4 + x)/2, ...), c = design$c(x, (.4 + x)/2, ...), theta)
    })
    pdf <- function(theta_obs) dnorm(theta_obs, mean = theta, sd = sqrt(4/gs_design$m))
    cdf <- function(theta_obs) pnorm(theta_obs, mean = theta, sd = sqrt(4/gs_design$m))
    eps <-  sqrt(.Machine$double.eps)
    left <- gg(gs_design$lower - eps) * cdf(gs_design$lower)
    right <- gg(gs_design$upper + eps) * (1 - cdf(gs_design$upper))
    mid <- integrate(function(x) gg(x) * pdf(x), lower = gs_design$lower + eps, upper = gs_design$upper - eps, rel.tol = .001, subdivisions = 10000L)$value
    left + mid + right
}

tbl_unconditional <- expand_grid(
        design = list(
            pp_geq_08_design,
            fixed_pp_design,
            reoptimize_design,
            gs_design
        ),
        theta = c(.35, .45, theta1)
    ) %>%
    mutate(
        nstart_factor = case_when(
            map_chr(design, ~.$name) == "PP >= 0.8" ~ 500,
            map_chr(design, ~.$name) == "not-worse" ~ 500,
            TRUE ~ NA_real_
        ),
        ess = pmap_dbl(list(design, theta, nstart_factor), ~ess(..1, ..2, ..3)),
        ep = pmap_dbl(list(design, theta, nstart_factor), ~ep(..1, ..2, ..3)),
        design = map_chr(design, ~.$name) %>%
            factor(
                levels = c("original GS", "PP >= 0.8", "not-worse", "reoptimize")
            ),
        theta = factor(theta)
    )

p2 <- tbl_unconditional %>%
    ggplot() +
        aes(ess, ep, shape = theta, color = design) +
        geom_point(size = 2) +
        theme_bw() +
        scale_y_continuous(breaks = seq(0, 1, by = .2), limits = c(0, 1)) +
        #xlim(0, NA_real_) +
        labs(y = "(expected) power", x = "expected sample size") +
        scale_shape_discrete(expression(theta)) +
        scale_color_discrete("") +
        guides(color = FALSE) +
        theme(
            legend.position = "top"
        )

p1 + p2 +
    plot_layout(width = c(2, 1)) &
    theme(legend.direction = 'horizontal')


ggsave("output/figures/application-example-response-adaptive.pdf", width = 8, height = 3.5)
