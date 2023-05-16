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
    Y = ped.pht
    # Z = I here
    lhs = [ X'X X'
            X   I + λ * giv]
    rhs = [X'Y; Y]
    nf = size(X, 2)
    ebv = (lhs \ rhs)[nf + 1 : end]
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

function uhp2gp(hp, loc, σₑ)
    nlc, nhp = size(hp)
    nid = nhp ÷ 2
    nlc == nrow(loc) || error("Length of effects doesn't match number of loci")
    qg = isodd.(hp[loc.qtl, 1:2:end]) + isodd.(hp[loc.qtl, 2:2:end])
    tbv = qg'loc.efct[loc.qtl]
    pht = tbv + randn(nid) * σₑ
    return tbv, pht
end