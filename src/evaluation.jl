"""
    function simpleBLUP(ped, giv, h²)
Calculate the BLUP of a trait with a pedigree, giv and a heritability.
The pedigree should have a column `:grt` to be as fixed effects.
This is a temporary dirty function for slides on May 17, 2023.

This function can server both as GBLUP and PBLUP.
"""
function simpleBLUP(ped, giv, h²)
    λ = (1.0 - h²) / h²
    X = incidence_matrix(ped.grt)
    Y = collect(skipmissing(ped.pht))
    # Z = I here
    lhs = [ X'X X'
            X   I + λ * giv]
    rhs = [X'Y; Y]
    nf = size(X, 2)
    ebv = (lhs \ rhs)[nf + 1 : end]
    ped.ebv = ebv
end

function animalModel(ped, giv, h²; fix = [:grt])
    λ = (1.0 - h²) / h²
    p = .!ismissing.(ped.pht) # only for the phenotyped animals
    Z = Zmat(p)
    X = incidence_matrix(select(ped, fix))[p, :]
    Y = collect(skipmissing(ped.pht))
    lhs = if issparse(giv)
        [X'X X'Z
         Z'X Z'Z + λ * giv]
    else
        Matrix([X'X X'Z
                Z'X Z'Z + λ * giv])
    end
    rhs = vec([X'Y; Z'Y])
    nf = size(X, 2)
    ebv = (lhs \ rhs)[nf + 1 : end]
    # LAPACK.posv!('U', lhs, rhs)
    ped.ebv = ebv
end

"""
    function rrblup_mme(x, z, y, h²; dd = 0.01, norm = false)
SNP effect calculation with rrBLUP and in MME way.

The design matrix `x` and SNP `012` genotype matrix `z` (nSNP by nID) are fed 
separately, with an option to `norm`alize the genotypes, or not. 

The returns are fixed effects, SNP effects and `lhs`.

Returning the `lhs` is for GWAS.  Its `L` is actually the LHS Cholesky factor, 
which can be inversed to calculate test statistics for QTL mapping.

## Note !!!
When `norm` is true, there should be **NO 1** column in `x`. 
As this column is accounted for in the genotypes.

## Not tested.
"""
function rrblup_mme(x, z, y, h²; dd = 0.01, norm = false)
    nlc, nid = size(z)
    λ = (1.0 - h²) / h² * nlc + dd
    x = reshape(x, nid, :) # to garantee `x` a matrix, if it's of only one column
    nf = size(x, 2)
    nb = nf + nlc           # total number of factors (fixed + random)

    mem = memavail() * 99 ÷ 100 # not all available memory
    mlh = nb^2 * 8                  # memory for LHS
    mem < mlh && error("Not enough memory for this calculation: $(mem/1024^3)G < $(mlh/1024^3)G")

    # the left hand side
    lhs = zeros(nb, nb)
    matmul!(view(lhs, 1:nf, 1:nf), x', x)  # block up left
    matmul!(view(lhs, nf+1:nb, 1:nf), z, x)  # block lower left
    matmul!(view(lhs, nf+1:nb, nf+1:nb), z, z') # block lower right
    if norm
        p = mean(z, dims=2)
        q = 1 .- p ./ 2
        v = 1 ./ sqrt.(p .* q)
        s = sum(z, dims=2)
        se = view(lhs, nf+1:nb, nf+1:nb) # south east of lhs, or genotype sub
        Threads.@threads for i in 1:nlc
            for j in 1:i        # ignore the upper triangle
                se[i, j] = se[i, j] - p[i] * s[j] - p[j] * s[i] + nid * p[i] * p[j]
            end
        end
        se .*= v
        se .*= v'
        sw = view(lhs, nf+1:nb, 1:nf)
        nf = size(x)[2]
        s = sum(x, dims=1)
        Threads.@threads for i in 1:nlc
            for j in 1:nf
                sw[i, j] -= p[i]s[j]
            end
        end
        sw .*= v
    end
    for i in nf+1:nb
        lhs[i, i] += λ
    end

    # the right hand side
    rhs = zeros(nb)
    matmul!(view(rhs, 1:nf), x', y)
    matmul!(view(rhs, nf+1:nb), z, y)

    # the solver
    LAPACK.posv!('L', lhs, rhs) # only use data in 'L' triangle of lhs
    (fixed=rhs[1:nf], a=rhs[nf+1:end], lhs=lhs)
end

"""
    function inbreeding(xy, loci)
Calculate inbreeding coefficients on uniquely coded alleles for everybody in `xy`.
"""
function inbreeding(xy, loci)
    hdr = readhdr(xy)
    mt, et, mj, ir, ic = xyhdr(hdr)
    (mt == 'F' && et == UInt16 && mj == 0) || error("Not a properfile.")
    gt = Mmap.mmap(xy, Matrix{UInt16}, (ir, ic), 24)
    nlc, nid = sum(loci), ic ÷ 2
    F = zeros(nid)
    for i in 1:nid
        f = 0.
        for j in 1:nlc
            loci[j] || continue
            gt[j, 2i-1] == gt[j, 2i] && (f += 1.)
        end
        F[i] = f / nlc
    end
    F
end

"""
    function uhp2gp(hp, loc, σₑ)
Calculate TBV and phenotype from haplotypes.
"""
function uhp2gp(hp, loc, σₑ)
    nlc, nhp = size(hp)
    nid = nhp ÷ 2
    nlc == nrow(loc) || error("Length of effects doesn't match number of loci")
    qg = isodd.(hp[loc.qtl, 1:2:end]) + isodd.(hp[loc.qtl, 2:2:end])
    tbv = qg'loc.efct[loc.qtl]
    pht = tbv + randn(nid) * σₑ
    return tbv, pht
end

"""
    function idealPop(xy, grt, lmp)
Given a genotype file `xy`, generation info `grt` and a linkage map `lmp`,
calculate the maximum level of an ideal population.
"""
function idealPop(xy, grt, lmp)
    # requirement check
    hdr = readhdr(xy)
    mt, et, mj, ir, ic = xyhdr(hdr)
    (mt == 'F' && mj == 0) || error("Not a proper file.")
    ir == nrow(lmp) || error("Genotypes and linkage map don't match.")
    ic == 2length(grt) || error("Genotypes and pedigree don't match.")
    (et == Int8 || et == UInt16) || error("Not a proper genotype file.")
    gt = if et == Int8
        Mmap.mmap(xy, Matrix{et}, (ir, ic), 24)
    else
        Int8.(isodd.(Mmap.mmap(xy, Matrix{et}, (ir, ic), 24)))
    end
    qg = gt[lmp.qtl, 1:2:end] + gt[lmp.qtl, 2:2:end]
    efct = lmp.efct[lmp.qtl]
    best = sum(efct[efct .> 0])
    ideal, va, np, nn = Float64[], Float64[], Int64[], Int64[]
    for ig in sort(unique(grt))
        iqg = qg[:, grt .== ig] # QTL genotypes of the ith generation
        plost, nlost = 0, 0  # positive and negative QTL lost
        for (i, v) in enumerate(efct)
            if v > 0
                if all(iqg[i, :] .== 0)
                    best -= v
                    plost += 1
                end
            else
                if all(iqg[i, :] .== 1)
                    best += v
                    nlost += 1
                end
            end
        end
        p = mean(iqg, dims=2) ./ 2
        vv = (2p .* (1 .- p))' * (efct .* efct)
        push!(ideal, best)
        push!(va, vv[1])
        push!(np, plost)    # number of positive loci lost
        push!(nn, nlost)    # number of negative loci lost
    end
    ideal, va, np, nn
end

"""
    function idealID(nqtl, d; shape = 0.75, ϵ = 1e-5)
Simulate effects of `nqtl` QTL of distribution `d`, and their ``Beta(.75, .75)``
distributed frequencies, this function return how many TBV SD an ideal ID is
away from population mean. The calculation stops when the ideal ID value SD
stabilizes.
"""
function idealID(nqtl, d; shape = 0.75, ϵ = 1e-5)
    vi, n = -1., 0
    ideal = Float64[]
    while true
        for _ in 1:10
            a = rand(d, nqtl) .* rand([-1, 1], nqtl) # QTL effects
            p = rand(Beta(shape, shape), nqtl) # QTL frequencies
            m = (2p .- 1)'a   # population mean
            v = sum(2 .* p .* (1 .- p) .* a .* a) # population genetic variance
            b = (sum(a[a .> 0]) - m) / sqrt(v)    # best ID from population mean
            push!(ideal, b)
        end
        n += 10
        t = var(ideal)
        abs(vi - t) < ϵ ? break : vi = t
    end
    mean(ideal), std(ideal), n
end
