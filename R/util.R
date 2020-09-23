#
# simple fixed gauss legendre integration
#
get_gauss_legendre_pivots_and_weights <- function(order) {
    order <- as.integer(order)
    if (order < 2) stop("At least two nodes are necessary for integration!")
    j   <- 1:(order - 1)
    mu0 <- 2
    b   <- j / (4 * j^2 - 1)^0.5
    A   <- rep(0, order * order)
    A[(order + 1) * (j - 1) + 2] <- b
    A[(order + 1) * j] <- b
    dim(A) <- c(order, order)
    sd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(sd$vectors[1, ]))
    w <- mu0 * w^2
    x <- rev(sd$values)
    list(pivots = x, weights = w)
}

gl_rule_global <- get_gauss_legendre_pivots_and_weights(100)

integrate_gl <- function(f, low, up, pivots = gl_rule_global$pivots, weights = gl_rule_global$weights) {
    a <- (up - low) / 2
    b <- a + low
    a * sum(weights * purrr::map_dbl(a*pivots + b, f))
}



#
# truncated normal priors and correponding methods for posterior, precictive etc.
#
TruncatedNormal <- function(mu, tau, a, b ) {
    res <- list(mu = mu, tau = tau, a = a, b = b)
    attr(res, "class") <- c("TruncatedNormal", class(res))
    res
}

condition <- function(dist, a = -Inf, b = Inf, ...) UseMethod("condition", dist)

condition.TruncatedNormal <- function(dist, a = -Inf, b = Inf, ...) {
    res <- list(mu = dist$mu, tau = dist$tau, a = max(a, dist$a), b = min(b, dist$b))
    if (res$a >= res$b) stop("conditioning interval has collapsed")
    attr(res, "class") <- c("TruncatedNormal", class(res))
    res
}

probability_density <- function(Theta, theta, ...) {
    UseMethod("probability_density", Theta)
}

probability_density.TruncatedNormal <- function(Theta, theta, ...) {
    truncnorm::dtruncnorm(theta, mean = Theta$mu, sd = Theta$tau, a = Theta$a, b = Theta$b)
}

cdf <- function(Theta, theta, ...) {
  UseMethod("cdf", Theta)
}

cdf.TruncatedNormal <- function(Theta, theta, ...) {
  truncnorm::ptruncnorm(theta, mean = Theta$mu, sd = Theta$tau, a = Theta$a, b = Theta$b)
}

posterior <- function(prior, n, z, ...) UseMethod("posterior", prior)

posterior.TruncatedNormal <- function(prior, n, z) {
     mu_post <- 1 / (1/prior$tau^2 + n) * (prior$mu / prior$tau^2 + sqrt(n) * z)
    tau_post <- 1 / (1/prior$tau^2 + n)
    res <- list(mu = mu_post, tau = tau_post, a = prior$a, b = prior$b)
    attr(res, "class") <- c("TruncatedNormal", class(res))
    res
}

joint_density <- function(Theta, theta, n, z)  {
    k <- length(z)
    if (k != length(n)) stop("n and z must be of same length")
    if (k != length(theta)) stop("z and theta must be of same length")
    dnorm(z, mean = sqrt(n)*theta, sd = 1) * probability_density(Theta, theta)
}

predictive_pdf <- function(Theta, n, z) {
      a <- (Theta$b - Theta$a) / 2
      b <- a + Theta$a
    # scale pivots for theta to (common) support
     pivots <- a*gl_rule_global$pivots + b
    res <- matrix(
        purrr::pmap_dbl(
            tidyr::expand_grid(theta = pivots, tibble::tibble(n = n, z = z)),
            ~joint_density(Theta, ..1, ..2, ..3)
        ),
        ncol = length(n),
        byrow = TRUE
    )
    res <- a*matrix(gl_rule_global$weights, nrow = 1) %*% res
    as.vector(res)
}

.f <- function(Theta, theta, n, z)  {
    k <- length(z)
    if (k != length(n)) stop("n and z must be of same length")
    if (k != length(theta)) stop("z and theta must be of same length")
    pnorm(z, mean = sqrt(n)*theta, sd = 1) * probability_density(Theta, theta)
}

predictive_cdf <- function(Theta, n, z) {
      a <- (Theta$b - Theta$a) / 2
      b <- a + Theta$a
    # scale pivots for theta to (common) support
     pivots <- a*gl_rule_global$pivots + b
    res <- matrix(
        purrr::pmap_dbl(
            tidyr::expand_grid(theta = pivots, tibble::tibble(n = n, z = z)),
            ~.f(Theta, ..1, ..2, ..3)
        ),
        ncol = length(n),
        byrow = TRUE
    )
    res <- a*matrix(gl_rule_global$weights, nrow = 1) %*% res
    as.vector(res)
}


#
# Conditional power and friends
#
conditional_power <- function(zm, m, n, c, theta) {
    if (any(m > n)) stop("m must be smaller than n")
    if (any(m < 1)) stop("m must be larger than 1")
    tau <- m/n
     mu <- sqrt(n)*theta + sqrt(tau)*(zm - sqrt(m)*theta)
     sd <- sqrt(1 - tau)
    pnorm(c, mu, sd, lower.tail = FALSE)
}

observed_conditional_power <- function(zm, m, n, c) conditional_power(zm, m, n, c, zm/sqrt(m))

predictive_power <- function(zm, m, n, c, prior) {
    Theta <- posterior(condition(prior, 0), n = m, z = zm)
    f <- function(theta) conditional_power(zm, m, n, c, theta) * probability_density(Theta, theta)
    integrate_gl(f, Theta$a, Theta$b)
}

expected_power <- function(n, c, prior) {
    Theta <- condition(prior, 0)
    integrate_gl(
        function(theta) pnorm(c, sqrt(n)*theta, 1, lower.tail = FALSE) * probability_density(Theta, theta),
        Theta$a, Theta$b
    )
}

# convert stage-two critical value to overall critical value (for adoptr design object)
critical_value <- function(design, zm) {
    n =  adoptr::n(design, zm, round = FALSE)
    n1 = adoptr::n1(design)
    c2 = adoptr::c2(design, zm)
    if (!is.finite(c2)) return(c2)
    sqrt(n - n1) / sqrt(n) * c2 + sqrt(n1/n) * zm
}
