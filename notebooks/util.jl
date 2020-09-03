# Gauss quadrature pivots and weights
ξ, ω = QuadGK.gauss(50)
# rule for interval [a, b]
integrate(f, a, b; ξ = ξ, ω = ω) = (b-a)/2 * sum(ω.*f.((b-a)/2 .* ξ .+ (a + b)/2))
# rule for interval [a, Inf)
function integrate(f, a; ξ = ξ, ω = ω)
    g(t) = f.(a .+ t ./ (1 .- t)) ./ (1 .- t).^2
    integrate(g, 0, 1)
end
# rule for interval (-Inf, a]
function integrate2(f, a; ξ = ξ, ω = ω)
    g(t) = f.(a .- (1 - t) ./ t) ./ t.^2
    integrate(g, 0, 1)
end
# rule for interval (-Inf, Inf)
function integrate(f; ξ = ξ, ω = ω)
    g(t) = f.(t ./ (1 .- t.^2)) .* (1 .+ t.^2) ./ (1 .- t.^2).^2
    integrate(g, -1, 1)
end




function EP(n, prior) 
         norm = integrate(prior, 0)
    cprior(θ) = prior(θ)/norm
    integrate(θ -> (1 - cdf(Normal(√(n)*θ, 1), crit))*cprior(θ), 0)
end

function CP(zm, n, m, c, θ)
    τ = m/n
    if τ >= 1
        return 0.0
    end
    μ = √(n)*θ + √(τ)*(zm - √(m)*θ)
    σ = √(1 - τ)
    1 - cdf(Normal(μ, σ), c)
end

OCP(zm, n, m, c) = CP(zm, n, m, c, zm/sqrt(m))

function posterior(θ, zm, m, prior)
                 prop(θ) = pdf(Normal(√(m)*θ, 1), zm) * prior(θ)
    normalising_constant = integrate(prop)
    prop(θ) / normalising_constant
end

function cposterior(θ, zm, m, prior)
    # condition prior on effect > 0
    normalizing_constant = integrate(prior, 0)
               cprior(θ) = prior(θ) / normalizing_constant
                 prop(θ) = pdf(Normal(√(m)*θ, 1), zm) * cprior(θ)
    normalising_constant = max(1e-9, integrate(prop, 0))
    prop(θ) / normalising_constant
end

PP(zm, n, m, c, prior) = integrate(θ -> CP.(zm, n, m, c, θ) * cposterior.(θ, zm, m, prior), 0)


marginal_pdf(zm, m, prior) = integrate(θ -> pdf.(Normal(√(m).*θ), zm) .* prior.(θ)) 

function cmarginal_pdf(zm, m, prior) 
    # condition prior on effect > 0
    normalizing_constant = integrate(prior, 0)
               cprior(θ) = prior(θ) / normalizing_constant
    integrate(θ -> pdf(Normal(√(m)*θ, 1), zm) * cprior(θ), 0)
end
