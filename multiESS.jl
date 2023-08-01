
using Distributions
using SpecialFunctions
using Statistics
using LinearAlgebra


function fminESS(p; alpha=0.05, eps=0.05, ess=nothing)
    """
    Minimum effective sample size
    """

    crit = cquantile(Chisq(p), 1 - alpha)
    tmp = 2.0 / p

    val = tmp * log(tmp) + log(π) - tmp * logabsgamma(p / 2)[1] + log(crit)

    if ess === nothing
        logminESS = val - 2.0 * log(eps)
        return round(Int, exp(logminESS))
    else
        logEPS = (val - log(ess)) / 2
        return exp(logEPS)
    end
end


function multiESS(X; b="sqroot", Noffsets=10, Nb=nothing)
    """
    Compute multivariate effective sample size of a single Markov chain X,
    using the multivariate dependence structure of the process.

    X: MCMC samples of shape (n, p)
    n: number of samples
    p: number of parameters

    b: specifies the batch size for estimation of the covariance matrix in
       Markov chain CLT. It can take a numeric value between 1 and n/2, or a
       string value between:

    "sqroot"    b=floor(n^(1/2)) (for chains with slow mixing time; default)
    "cuberoot"  b=floor(n^(1/3)) (for chains with fast mixing time)
    "lESS"      pick the b that produces the lowest effective sample size
                for a number of b ranging from n^(1/4) to n/max(20,p); this
                is a conservative choice

    If n is not divisible by b, Sigma is recomputed for up to Noffsets subsets
    of the data with different offsets, and the output mESS is the average over
    the effective sample sizes obtained for different offsets.

    Nb specifies the number of values of b to test when b=\"lESS\"
    (default Nb=200). This option is unused for other choices of b.

    Original source: https://github.com/lacerbi/multiESS and https://github.com/Gabriel-p/multiESS

    Reference:
    Vats, D., Flegal, J. M., & Jones, G. L. "Multivariate Output Analysis
    for Markov chain Monte Carlo", arXiv preprint arXiv:1512.07713 (2015).

    """

    # MCMC samples and parameters
    n, p = size(X)

    if p > n
        throw(ArgumentError(
            "More dimensions than data points, cannot compute effective sample size.")
        )
    end

    # Input check for batch size B
    if isa(b, String)
        if b ∉ ["sqroot", "cuberoot", "lESS"]
            throw(ArgumentError("Unknown string for batch size. Allowed arguments are \"sqroot\", \"cuberoot\" and \"lESS\".")
            )
        end
        if b != "lESS" && Nb !== nothing
            @warn "Nonempty parameter Nb will be ignored (Nb is used only with \"lESS\" batch size b)."
        end
    else
        if !(1.0 < b < (n / 2))
            throw(ArgumentError(
                "The batch size b needs to be between 1 and N/2.")
            )
        end
    end

    # Compute multiESS for the chain
    mESS = multiESS_chain(X, n, p, b, Noffsets, Nb)

    return mESS
end


function multiESS_chain(Xi, n, p, b, Noffsets, Nb)
    """
    Compute multiESS for a MCMC chain.
    """

    if b == "sqroot"
        b = [floor(Int, n^(1.0 / 2))]
    elseif b == "cuberoot"
        b = [floor(Int, n^(1.0 / 3))]
    elseif b == "lESS"
        b_min = floor(n^(1.0 / 4))
        b_max = max(floor(n / max(p, 20)), floor(sqrt(n)))
        if Nb === nothing
            Nb = 200
        end
        # Try Nb log-spaced values of b from b_min to b_max
        println("b_min is $(b_min), b_max is $(b_max)")
        b = floor.(Int, exp.(range(log(b_min), log(b_max), Nb)))
    end

    # Sample mean
    theta = mean(Xi, dims=1)
    # Determinant of sample covariance matrix
    detLambda = det(cov(Xi))

    # Compute mESS
    mESS_i = []
    for bi in b
        push!(mESS_i, multiESS_batch(Xi, n, p, theta, detLambda, bi, Noffsets))
    end
    # Return lowest mESS
    mESS = minimum(mESS_i)

    return mESS
end


function multiESS_batch(Xi, n, p, theta, detLambda, b, Noffsets)
    """
    Compute multiESS for a given batch size b.
    """

    a = floor(Int, n / b) # batch size
    Sigma = zeros((p, p)) # initiate mBM estimator
    offsets = unique(round.(Int, collect(range(0, n - a * b, Noffsets)))) # Unique values for the offsets

    for j in offsets
        for idx_a in 0:a-1 # Run through the batches

            # get the batch
            Y = Xi[1+idx_a*b+j:(idx_a+1)*b+j, :]

            # compute the batch mean
            Ybar = mean(Y, dims=1)

            # Ingredient for formula (10) in the paper
            Z = Ybar - theta

            # Add to the mBM estimator
            Sigma .+= Z' * Z
        end
    end

    Sigma = (Sigma * b) / (a - 1) / length(offsets)
    mESS = floor(Int, n * (detLambda / det(Sigma))^(1 / p))

    return mESS
end

function main()

    n = rand(200:10000)
    p = rand(1:5)

    MC_chain = randn(n, p)

    println("number of samples: $(n) ; number of parameters : $(p)")
    ESS = multiESS(MC_chain, b="sqroot")
    println("ESS is $(ESS)")
end