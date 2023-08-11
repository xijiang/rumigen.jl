#=
"""
    function _optnid(ebv, sex, A, K)
A `-` in the function name means it is a private function, not to be called
directly by end users.
"""
function _optnid(ebv, sex, A, K)
end

"""
    function optconstraint(ped, A, prop)
The optimum constraint, such that a solution exists. It also meet the required
number of breeding ID. The constraint is searched between ``min(K) =
sum(QAQi)/4`` and ``4min(K)``.
"""
function optconstraint(ped, A, prop)
    ["ebv", "sex"] ⊆ names(ped) || error("Not a proper pedigree")
    issymmetric(A) || error("A is not symmetric")
    size(A, 1) == size(ped, 1) || error("A and ped do not match")

    minK = begin
        Q = [ped.sex .== 1 ped.sex .== 0]
        Ai = inv(A)
        sum(inv(Q'Ai*Q))/4
    end
    maxK = 4minK
    k = (minK + maxK) / 2
end

# a bisect search function
function _bisect(f, a, b, tol)
    fa, fb = f(a), f(b)
    fa * fb > 0 && error("f(a) and f(b) must have opposite signs")
    while abs(b - a) > tol
        c = (a + b) / 2
        fc = f(c)
        if fa * fc < 0
            b, fb = c, fc
        else
            a, fa = c, fc
        end
    end
    (a + b) / 2
end
=#

"""
`K` = is the restriction on average relationships (``Meuwissen, 1997, J Anim
Sci``). Let the constraint on the inbreeding in generation `t`` be:

`k(t) = k(t-1) + DF*(1-k(t-1))`

where `k(t-1)` = constraint on inbreeding in generation `t-1` and `k(0)=0`. Then
the constraint on relationship is `K=2k(t)`.

Based on Theo's function below with some simplifications.
"""
function myopt(ped, A, K; silent = false)
    ["ebv", "sex"] ⊆ names(ped) || error("Not a proper pedigree")
    issymmetric(A) || error("A is not symmetric")
    nid, itr = size(ped, 1), 0
    size(A, 1) == nid || error("A and ped do not match")
    id = 1:nid
    silent || tprintln("{green}Searching solution for K = {/green}$(round(K, digits=3)): ")
    while true
        itr += 1
        u = ped.ebv[id]
        Q = [ped.sex[id] .== 1 ped.sex[id] .== 0]
        Ai = inv(A[id, id])
        QAQ = Q'Ai * Q
        QAQi = inv(QAQ)
        f1 = u' * (Ai - Ai * Q * QAQi * Q'Ai) * u # numerator
        f2 = 4K - sum(QAQi) # denominator
        # if f2 ≤ 0
        #     silent || @info "Cannot meet constraint K = $K"
        #     return 0.5Ai * Q * QAQi * ones(length(id))
        # end
        f2 ≤ 0 && return 0.5Ai * Q * QAQi * ones(length(id))
        f1 = max(f1, 1e-10)
        λ₀ = sqrt(f1/f2)
        λ = QAQi * (Q'Ai * u .- λ₀)
        c = Ai * (u - Q * λ) / (2λ₀)
        ix = findall(c .> 0) # indices of ID of next round
        silent || tprint("  iter: $itr, {cyan}nID{/cyan}: $(length(ix))")
        if length(ix) == length(id)
            silent || tprintln("\n{blue}Solution found{/blue}, 
                                n = $(length(ix))", 
                                "c'Ac = ", round(c'A[id, id] * c, digits = 3))
            rc = zeros(nid)
            rc[id] = c
            return rc
        end
        id = id[ix]
    end
end

"""
    fungencont(dat, A, K)

## Usage
```julia
    C = fungencont(dat, A, K)
```

- ``C`` is the vector of optimal contribution
- ``dat`` is a `n × 3` matrix with columns `[EBV I(male) I(female)]` where
    `I(male)` is an indicator whether the animal is male (0/1)
- ``K`` is the restriction on average relationships (Meuwissen, 1997, J Anim
  Sci)
- ``A`` is the numerator relationship matrix of the selection candidates for the
  next generation.

Let the constraint on the inbreeding in generation `t` be:

``k(t) = k(t - 1) + DF × (1 - k(t - 1))``

where ``k(t - 1)`` is constraint on inbreeding in generation ``t - 1`` and
``k(0) = 0``.

Then the constraint on relationship is ``K = 2k(t)``
"""
function fungencont(dat, A, K)
    #dat[:,1]=EBV
    #dat[:,2:3]=sex indicator (0/1)
    n = size(dat, 1)
    ind = collect(1:n)
    println(" Ncandidates            = ", n)
    println(" contraint relationship = ", K)
    if (all(dat[:, 3] .== 0))
        K = K / 4  #since c adds to .5
    end
    c = zeros(n)
    ierr = 0
    for it = 1:1000
        u = dat[ind, 1]
        if (any(dat[:, 3] .> 0))
            Q = dat[ind, 2:3]
            isex = 2
        else  #one sex
            Q = dat[ind, 2]
            isex = 1
        end #if
        AI = inv(A[ind, ind])
        QAQ = Q' * AI * Q
        QAQI = inv(QAQ)
        denominat = 4 * K - sum(QAQI)
        numerat = u' * (AI - AI * Q * QAQI * Q' * AI) * u
        if (denominat <= 0.0)
            println(" cannot achieve constraint ", K, " MINIMISATION OF RELATIONSHIPS ", size(QAQI))
            ierr = 0
            if (isex == 2)
                c = 0.5 * AI * Q * QAQI * ones(size(QAQI, 1))
            else
                c = 0.5 * AI * Q * QAQI  #*ones(size(QAQI,1))
            end
            #	    return c
        else
            numerat = max(numerat, 1.e-10)
            #            println(" denominat ",denominat)
            #            println(" numerator ",numerat)
            lamb02 = numerat / denominat
            lamb0 = sqrt(lamb02)
            #          println(" lamb0 ",lamb0)    
            lamb = QAQI * (Q' * AI * u .- lamb0)
            c = AI * (u - Q * lamb) / (2 * lamb0)
        end
        ind2 = findall(c .>= 0.0)
        println(" iter ", it, " still in solution ", length(ind2), " (old= ", length(ind))
        if (length(ind2) == length(ind))
            println(" solution found; n=", length(ind))
            break
        end
        ind = ind[ind2]
    end

    if (ierr == 0)
        cc = zeros(n)
        cc[ind] = c
        println("  cAc ", c' * A[ind, ind] * c, " K=  ", K)
        if (any(dat[:, 3] .> 0))
            return cc
        else
            return 2 * cc
        end
    end

end #function
