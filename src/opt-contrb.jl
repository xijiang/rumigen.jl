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
        λ₀ = sqrt(f1 / f2)
        λ = QAQi * (Q'Ai * u .- λ₀)
        c = Ai * (u - Q * λ) / (2λ₀)
        ix = findall(c .> 0) # indices of ID of next round
        silent || tprint("  iter: $itr, {cyan}nID{/cyan}: $(length(ix))")
        if length(ix) == length(id)
            silent || tprintln(
                "\n{blue}Solution found{/blue}, 
                 n = $(length(ix))",
                "c'Ac = ",
                round(c'A[id, id] * c, digits = 3),
            )
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
            println(
                " cannot achieve constraint ",
                K,
                " MINIMISATION OF RELATIONSHIPS ",
                size(QAQI),
            )
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

# DOSc optimises gains with known cost factor lamb0
# find optimal c and optimal numbers of parents to selected (contributions are equal)
# uses Linear Algebra

# DOSc optimises gains with known costfactor lamb0
# find optimal c and optimal numbers of parents to selected (contributions are equal)
# uses LinearAlgebra
function DOSc(lamb0, uhat, A, A12, s, Ktilde, Nx, sex)
    N = size(uhat, 1)
    Nsex = size(Nx, 1)
    (Nsex == 2) ? fact = 0.5 : fact = 1.0

    #step0:
    c = zeros(N, 1)
    nx = zeros(Int, size(Nx, 1))
    lamb0A12 = 2 * lamb0 * (1 - s) * A12
    NS = collect(1:N) #not selected set
    IDS = collect(1:N) #list of IDs

    sumS = sumD = sumSD = sumA12S = sumA12D = Kval = 0.0
    setS = Int64.([])
    setD = Int64.([])

    for ianim = 1:sum(Nx)  #select Nx animals
        #step1:
        crit = uhat
        sum(nx) > 0 ? crit -= 2 * fact * lamb0 * s * A * c : nothing
        (1 - s) > 0 ? crit -= lamb0A12 : nothing

        #step2:
        idsel = 0
        while (idsel == 0)
            (amax, imax) = findmax(crit[NS])
            isex = sex[NS][imax]
            if (nx[isex] < Nx[isex])
                idsel = IDS[NS][imax]
                nx[isex] += 1
                if (isex == 1)
                    push!(setS, idsel)
                    sumS += 2 * sum(A[setS, idsel]) - A[idsel, idsel]
                    sumSD += sum(A[setD, idsel])
                    (1 - s > 0) ? sumA12S += A12[idsel] : nothing
                else
                    push!(setD, idsel)
                    sumD += 2 * sum(A[setD, idsel]) - A[idsel, idsel]
                    sumSD += sum(A[setS, idsel])
                    (1 - s > 0) ? sumA12D += A12[idsel] : nothing
                end
            end
            deleteat!(NS, imax)
        end
        c .= 0
        (nx[1] > 0) ? c[setS] .= 1 / nx[1] : nothing
        (Nsex == 2) && (nx[2] > 0) ? c[setD] .= 1 / nx[2] : nothing

        #Step3
        Kval = 0.0
        (nx[1] > 0) ? Kval += s * fact * fact * sumS / nx[1] / nx[1] : nothing
        (size(nx, 1) == 2) && (nx[2] > 0) ? Kval += s * fact * fact * sumD / nx[2] / nx[2] :
        nothing
        (size(nx, 1) == 2) && (nx[2] * nx[1] > 0) ?
        Kval += s * fact * sumSD / nx[1] / nx[2] : nothing
        if (1 - s > 0)
            (nx[1] > 0) ? Kval += 2 * (1 - s) * fact * sumA12S / nx[1] : nothing
            (nx[2] > 0) ? Kval += 2 * (1 - s) * fact * sumA12D / nx[2] : nothing
        end
        #    println("Kval,cAc ",Kval," ",c'*A*c*fact^2)
        #    println("S:Kval,cAc ",fact*fact*sumS/nx[1]/nx[1]," ",c[sex.==1]'*A[sex.==1,sex.==1]*c[sex.==1]*fact^2)
        #    nx[2]>0 && println("D:Kval,cAc ",fact*fact*sumD/nx[2]/nx[2]," ",c[sex.==2]'*A[sex.==2,sex.==2]*c[sex.==2]*fact^2)
        if (prod(nx) > 0)    #both sexes selected
            Kval <= Ktilde ? break : nothing
            all(nx .== Nx) ? break : nothing
        end
    end
    # println("lamb, Kval,cAc ", lamb0, " ", Kval, " ", c' * A * c * fact^2, " nx ", nx)
    return c, nx
end #function

"""
    function DOSop(uhat, A, A12, s, Ktilde, Nx, sex)
DOSop.jl:  OC selection with equal number of parents

`DOSop` = Discrete optimal contribution selection with optimized number of
parents find optimal cost factor to achieve maximum gain uses Linear Algebra,
Statistics, DOSc.jl


Usage: `c = DOSop(EBV, A, A12, s, Ktilde, Nx, sex)`

Where:

- `EBV = ebvs(N x 1)`, here it is argument `uhat`
- `A = (N x N)` relationship matrix (could be G or IBD)
- `A12`= relationships with already selected animals (=Nx1 vector of zeros in
  your case)
- `s` = contribution yet to be made by current selection = 1 in your case
- `Ktilde` = constraint on relationship (calculated as usual)
- `Nx = (2 x 1)` vector with number of male and female candidates
- `Sex = (N x 1)` vector with 1 or 2 indicating male or female, resp.

``DOSc.jl`` is called by DOSop and does the actual optimization.
"""
function DOSop(uhat, A, A12, s, Ktilde, Nx, sex)
    GOLD = 0.618   #golden bracket search
    R = GOLD
    C = 1 - GOLD

    N = size(uhat, 1)

    #find interval for lamb0
    c = zeros(size(uhat, 1))
    if (s * (Ktilde - s * mean(A) - 2 * (1 - s) * mean(A12)) <= 0) #constraint cannot be achieved; minimisation
        c[sex.==1] .= 1 / Nx[1]
        size(Nx, 1) == 2 ? c[sex.==2] .= 1 / Nx[2] : nothing
        return c
    else
        lamb1 = 0   #lower bound
        lamb2 = sqrt(N * var(uhat) / (s * (Ktilde - s * mean(A) - 2 * (1 - s) * mean(A12))))   #upperbound (guess)
        if (lamb2 <= 0)
            return zeros(size(uhat, 1))
        end
        lamb3 = lamb2 + (1 + GOLD) * (lamb2 - lamb1)
    end

    # golden section search
    x0 = x1 = x2 = x3 = 0.0  #initialise
    x0 = lamb1
    (c, n) = DOSc(x0, uhat, A, A12, s, Ktilde, Nx, sex)
    f0 = sum(c' * uhat)
    x1 = lamb2
    (c, n) = DOSc(x1, uhat, A, A12, s, Ktilde, Nx, sex)
    f1 = sum(c' * uhat)

    if (f0 > f1)  #maximum is within interval x0--x1 => find point that exceeds f0
        xa = fa = 0.0
        for i = 1:50
            xa = 0.5 * (x0 + x1)
            (c, n) = DOSc(xa, uhat, A, A12, s, Ktilde, Nx, sex)
            fa = sum(c' * uhat)
            (fa > f0) ? break : (x1, f1) = (xa, fa)
        end
        (x3, f3) = (x1, f1)
        (x1, f1) = (xa, fa)
        x2 = x1 + C * (x3 - x2)
        (c, n) = DOSc(x2, uhat, A, A12, s, Ktilde, Nx, sex)
        f2 = sum(c' * uhat)
    else
        x3 = lamb3
        x2 = x1 + C * (x3 - x1)
        (c, n) = DOSc(x2, uhat, A, A12, s, Ktilde, Nx, sex)
        f2 = sum(c' * uhat)
        (c, n) = DOSc(x3, uhat, A, A12, s, Ktilde, Nx, sex)
        f3 = sum(c' * uhat)
    end
    if !(f0 <= max(f1, f2) >= f3)
        println(
            "problem is NOT bracketed by $(x0), $(x1), $(x2), $(x3) with gains $(f0), $(f1), $(f2), $(f3) ",
        )
        return zeros(size(uhat, 1))
    end

    while (abs(x3 - x0) > 0.001 * (abs(x1) + abs(x2)))
        if (f2 > f1)   #f2 is best
            x0 = x1
            x1 = x2
            f1 = f2
            x2 = R * x1 + C * x3
            (c, n) = DOSc(x2, uhat, A, A12, s, Ktilde, Nx, sex)
            f2 = sum(c' * uhat)
        else
            x3 = x2
            x2 = x1
            f2 = f1
            x1 = R * x2 + C * x0
            (c, n) = DOSc(x1, uhat, A, A12, s, Ktilde, Nx, sex)
            f1 = sum(c' * uhat)
        end
    end #while

    # println(" costfact = $(x1), $(x2)")
    if (f1 > f2)
        (c, n) = DOSc(x1, uhat, A, A12, s, Ktilde, Nx, sex)
        return vec(c)
    else
        (c, n) = DOSc(x2, uhat, A, A12, s, Ktilde, Nx, sex)
        return vec(c)
    end
end #function
