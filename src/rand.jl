"""
    PoissonSampler
For a finitely presented group PoissonSampler returns group elements represented
by words of length at most `R ~ Poisson(λ)` chosen uniformly at random.

For finitely presented groups the Product Replacement Algorithm
(see `PRASampler` from `GroupsCore.jl`) doesn't make much sense due to
overly long words it produces. We therefore resort to a pseudo-random method,
where a word `w` of length `R` is chosen uniformly at random among all
words of length `R` where `R` follows the Poisson distribution.

!!! note
    Due to the choice of the parameters (`λ=8`) and the floating point
    arithmetic the sampler will always return group elements represented by
    words of length at most `42`.
"""
struct PoissonSampler{G,T} <: Random.Sampler{T}
    group::G
    λ::Int
end

function PoissonSampler(G::AbstractFPGroup; λ)
    return PoissonSampler{typeof(G),eltype(G)}(G, λ)
end

function __poisson_invcdf(val; λ)
    # __poisson_pdf(k, λ) = λ^k * ℯ^-λ / factorial(k)
    # pdf = ntuple(k -> __poisson_pdf(k - 1, λ), 21)
    # cdf = accumulate(+, pdf)
    # radius = something(findfirst(>(val), cdf) - 1, 0)
    # this is the iterative version:
    pdf = ℯ^-λ
    cdf = pdf
    k = 0
    while cdf < val
        k += 1
        pdf = pdf * λ / k
        cdf += pdf
    end
    return k
end

function Random.rand(rng::Random.AbstractRNG, sampler::PoissonSampler)
    R = __poisson_invcdf(rand(rng); λ = sampler.λ)

    G = sampler.group
    n = length(alphabet(G))
    W = word_type(G)
    T = eltype(W)

    letters = rand(rng, T(1):T(n), R)
    word = W(letters, false)
    return G(word)
end

function Random.Sampler(
    RNG::Type{<:Random.AbstractRNG},
    G::AbstractFPGroup,
    repetition::Random.Repetition = Val(Inf),
)
    return PoissonSampler(G; λ = 8)
end
