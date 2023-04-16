"""
    function quickgt(nlc::Int, nid::Int; maf = .2, bp = .75)
A quick way to simulate genotypes of `nLoci` by `nID`.
Allele frequencies are sampled from a `Beta(.75, .75)`,
conditioned on `maf`.
It has a U-shaped distribution by default.
"""
function quickgt(nlc::Int, nid::Int; maf = .2, bp = .75)
    (maf ≤ 0 || maf ≥ .5) && error("maf $maf not in (0, 0.5)")
    gt = zeros(Int8, nlc, nid)
    for iic in 1:nlc
        p = 0
        while p <= maf || p >= 1 - maf
            p = rand(Beta(bp, bp))
        end
        rand!(Binomial(2, p), view(gt, iic, :))
    end
    gt
end

"""
    function quickhap(nlc::Int, nid::Int; maf = .2, bp = .75)
A quick way to simulate haplotypes of `nLoci` by `2nID`.  Allele
Frequencies are sampled from a `Beta(.75, .75)`, conditioned on
`maf`.  It has a U-shaped distribution by default.
"""
function quickhap(nlc, nid; maf = 0.2, bp = .75)
    (maf ≤ 0 || maf ≥ .5) && error("maf $maf not in (0, 0.5)")
    nhp = 2nid
    hp = zeros(Int8, nlc, nhp)
    for iic in 1:nlc
        p = 0
        while p <= maf || p >= 1 - maf
            p = rand(Beta(bp, bp))
        end
        rand!(Binomial(1, p), view(hp, iic, :))
    end
    hp
end
