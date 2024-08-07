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
    lhs = [
        X'X X'
        X I+λ*giv
    ]
    rhs = [X'Y; Y]
    nf = size(X, 2)
    ebv = (lhs\rhs)[nf+1:end]
    ped.ebv = ebv
end

"""
    function animalModel(ped, giv, h²; fix = [:grt])

Calculate the BLUP of a trait with a pedigree, giv and a heritability. The model
use `:grt` as fixed effects and `:id` as random effects. Note this function is
updating all animals' EBV. You can restore the original EBV by pre-saving them.
"""
function animalModel(ped, giv, h²; fix = [:grt])
    λ = (1.0 - h²) / h²
    p = .!ismissing.(ped.pht) # only for the phenotyped animals
    @show length(p), length(ped.pht)
    Z = Zmat(p)
    X = incidence_matrix(select(ped, fix))[p, :]
    Y = collect(skipmissing(ped.pht))
    lhs = if issparse(giv)
        [
            X'X X'Z
            Z'X Z'Z+λ*giv
        ]
    else
        Matrix([
            X'X X'Z
            Z'X Z'Z+λ*giv
        ])
    end
    rhs = vec([X'Y; Z'Y])
    nf = size(X, 2)
    ebv = (lhs\rhs)[nf+1:end]
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
    mem < mlh &&
        error("Not enough memory for this calculation: $(mem/1024^3)G < $(mlh/1024^3)G")

    # the left hand side
    lhs = zeros(nb, nb)
    matmul!(view(lhs, 1:nf, 1:nf), x', x)  # block up left
    matmul!(view(lhs, nf+1:nb, 1:nf), z, x)  # block lower left
    matmul!(view(lhs, nf+1:nb, nf+1:nb), z, z') # block lower right
    if norm
        p = mean(z, dims = 2)
        q = 1 .- p ./ 2
        v = 1 ./ sqrt.(p .* q)
        s = sum(z, dims = 2)
        se = view(lhs, nf+1:nb, nf+1:nb) # south east of lhs, or genotype sub
        Threads.@threads for i = 1:nlc
            for j = 1:i        # ignore the upper triangle
                se[i, j] = se[i, j] - p[i] * s[j] - p[j] * s[i] + nid * p[i] * p[j]
            end
        end
        se .*= v
        se .*= v'
        sw = view(lhs, nf+1:nb, 1:nf)
        nf = size(x)[2]
        s = sum(x, dims = 1)
        Threads.@threads for i = 1:nlc
            for j = 1:nf
                sw[i, j] -= p[i]s[j]
            end
        end
        sw .*= v
    end
    for i = nf+1:nb
        lhs[i, i] += λ
    end

    # the right hand side
    rhs = zeros(nb)
    matmul!(view(rhs, 1:nf), x', y)
    matmul!(view(rhs, nf+1:nb), z, y)

    # the solver
    LAPACK.posv!('L', lhs, rhs) # only use data in 'L' triangle of lhs
    (fixed = rhs[1:nf], a = rhs[nf+1:end], lhs = lhs)
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
    for i = 1:nid
        F[i] = sum(gt[loci, 2i-1] .== gt[loci, 2i]) / nlc
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
    nlc == size(loc, 1) || error("Length of effects doesn't match number of loci")
    qg = isodd.(hp[loc.qtl, 1:2:end]) + isodd.(hp[loc.qtl, 2:2:end])
    tbv = qg'loc.efct[loc.qtl]
    pht = tbv + randn(nid) * σₑ
    return tbv, pht
end

function qtl_fixed(qhp, hgrt, efct, maf)
    #i₀ = findall(hgrt .== minimum(hgrt))  # indices of the starting generation
    i₀ = findall(hgrt .== 0)  # indices of the starting generation
    frq = mean(qhp[:, i₀], dims = 2)
    lmf = frq .< maf .|| frq .> 1 - maf # low maf loci
    best, nhp = sum(efct[efct.>0]), size(qhp, 2)
    ideal, va, np, nn, nmp, nmn =
        Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]

    for ig in sort(unique(hgrt)) # !!! below can be vectorized !!!
        ihp = view(qhp, :, hgrt .== ig) # QTL haplotypes of the ith generation
        ss = sum(ihp, dims = 2) # sum of haplotypes
        plf = (ss .== 0 .&& efct .> 0)
        prf = (ss .== nhp .&& efct .< 0)
        nlf = (ss .== 0 .&& efct .< 0)
        nrf = (ss .== nhp .&& efct .> 0)
        plost = sum(plf) + sum(prf)
        nlost = sum(nlf) + sum(nrf)
        pmlst = sum(plf .&& lmf) + sum(prf .&& lmf)
        nmlst = sum(nlf .&& lmf) + sum(nrf .&& lmf)
        down = efct'plf - efct'prf
        push!(ideal, best - down[1])
        q = mean(ihp, dims = 2)
        vv = (2q .* (1 .- q))' * (efct .* efct)
        push!(va, vv[1])
        push!(np, plost)    # number of positive loci lost
        push!(nn, nlost)    # number of negative loci lost
        push!(nmp, pmlst)   # number of positive loci lost of maf
        push!(nmn, nmlst)   # number of negative loci lost of maf
    end
    2ideal, va, np, nn, nmp, nmn
end

function snp_fixed(gt, hgrt, maf)
    i₀ = findall(hgrt .== minimum(hgrt))  # indices of the starting generation
    frq = mean(gt[:, i₀], dims = 2)
    tlst, mlst, nhp = Float64[], Float64[], size(gt, 2)

    for ig in sort(unique(hgrt))
        fx = sum(view(gt, :, hgrt .== ig), dims = 2)
        push!(tlst, sum(fx .== 0) + sum(fx .== nhp))
        push!(mlst, sum((fx .== 0 .|| fx .== nhp) .&& ((frq .< maf .|| frq .> 1 - maf))))
    end
    tlst, mlst
end

"""
    function idealPop(xy, grt, lmp; maf = 0.2)
Given a genotype file `xy`, generation info `grt` and a linkage map `lmp`,
calculate the maximum level of an ideal population.

Add QTL change for allele frequencies in (0, maf = 0.2). 2023-07-17.

Ideal needs to be doubled, as the maximum genotype is 2. will modify this later.
2023-10-22, or remember to double the ideal value in summary.

There was an error about the definition of positive and negative alleles.
When effect is positive then allele 1 is positive. When effect is negative
then allele 0 is positive. This version will also count the number of loci that
fixed on 'good' alleles. This is to compare the effects of migration and selection.
2024-01-09.
"""
function idealPop(xy, grt, lmp; maf = 0.2)
    # requirement check
    hdr = readhdr(xy)
    mt, et, mj, ir, ic = xyhdr(hdr)
    (mt == 'F' && mj == 0) || error("Not a proper file.")
    ir == size(lmp, 1) || error("Genotypes and linkage map don't match.")
    ic == 2length(grt) || error("Genotypes and pedigree don't match.")
    hgrt = repeat(grt, inner = 2)
    (et == Int8 || et == UInt16) || error("Not a proper genotype file.")
    gt = if et == Int8
        Mmap.mmap(xy, Matrix{et}, (ir, ic), 24)
    else
        Int8.(isodd.(Mmap.mmap(xy, Matrix{et}, (ir, ic), 24)))
    end
    ideal, va, np, nn, nmp, nmn =
        qtl_fixed(view(gt, lmp.qtl, :), hgrt, lmp.efct[lmp.qtl], maf)
    clst, cmls = snp_fixed(view(gt, lmp.chip, :), hgrt, maf)
    rlst, rmls = snp_fixed(view(gt, lmp.ref, :), hgrt, maf)
    cvr, cvc, cvq = begin
        g0 = findall(hgrt .== 0)
        gn = findall(hgrt .== hgrt[end])
        q0 = mean(gt[:, g0], dims = 2)
        qn = mean(gt[:, gn], dims = 2) - q0
        vld = q0 .≠ 0.0 .&& q0 .≠ 1.0
        σ = sqrt.(q0 .* (1 .- q0))
        q0 .-= 0.5
        q0 ./= σ
        qn ./= σ
        ref = lmp.ref .&& vld
        chp = lmp.chip .&& vld
        qtl = lmp.qtl .&& vld
        cov(q0[ref], qn[ref]), cov(q0[chp], qn[chp]), cov(q0[qtl], qn[qtl])
    end
    (
        ideal = ideal,
        va = va,
        np = np,
        nn = nn,
        nmp = nmp,
        nmn = nmn,
        clst = clst,
        rlst = rlst,
        cmls = cmls,
        rmls = rmls,
        cvr = cvr,
        cvc = cvc,
        cvq = cvq,
    )
end

"""
    function idealID(nqtl, d; shape = 0.75, ϵ = 1e-5)
Simulate effects of `nqtl` QTL of distribution `d`, and their ``Beta(.75, .75)``
distributed frequencies, this function return how many TBV SD an ideal ID is
away from population mean. The calculation stops when the ideal ID value SD
stabilizes.
"""
function idealID(nqtl, d; shape = 0.75, ϵ = 1e-5)
    vi, n = -1.0, 0
    ideal = Float64[]
    while true
        for _ = 1:10
            a = rand(d, nqtl) .* rand([-1, 1], nqtl) # QTL effects
            p = rand(Beta(shape, shape), nqtl) # QTL frequencies
            m = (2p .- 1)'a   # population mean
            v = sum(2 .* p .* (1 .- p) .* a .* a) # population genetic variance
            b = (sum(a[a.>0]) - m) / sqrt(v)    # best ID from population mean
            push!(ideal, b)
        end
        n += 10
        t = var(ideal)
        abs(vi - t) < ϵ ? break : vi = t
    end
    mean(ideal), std(ideal), n
end
