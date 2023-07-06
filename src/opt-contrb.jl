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
    println(" Ncandidates           = ", n)
    println(" contraint relationship= ", K)
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
            lamb = QAQI * (Q' * AI * u - lamb0)
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
