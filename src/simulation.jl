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

"""
    function norm_qtl(Q::Matrix{Int8}, efct, ϵ)
Normalize QTL effect, such that the TBV variance is within `1.0 ± ϵ`.
"""
function norm_qtl(Q::Matrix{Int8}, efct, ϵ)
    nqtl= size(Q, 1)
    bv = Q'efct
    m, s = mean(bv), std(bv)
    while abs(m) > ϵ || abs(s - 1) > ϵ
        efct .-= m/nqtl
        efct ./= s
        bv = Q'efct
        m, s = mean(bv), std(bv)
    end
end

"""
    function simQTL(xy::AbstractString, mmp::AbstractString)
Sample QTL effects from genotypes in `xy` and QTL indice in `mmp`.
The genotypes must be of `Int8` type.
"""
function simQTL(xy::AbstractString, mmp::AbstractString; d = Laplace())
    imp = deserialize(mmp)
    hasproperty(imp, :qtl) || error("No QTL column in $mmp")
    qtn = imp.qtl
    nq = sum(qtn)
    ihdr = readhdr(xy)
    mt, et, mj, ir, ic = xyhdr(ihdr)
    (mt == 'F' && mj == 0 && et == Int8 && ir == length(qtn)) || error("Not a proper input $xy")
    snp = Mmap.mmap(xy, Matrix{Int8}, (ir, ic), 24)
    qtl = view(snp, qtn, 1:2:ic) + view(snp, qtn, 2:2:ic)
    efct = rand(d, nq) .* rand([-1, 1], nq)
    norm_qtl(qtl, efct, 1e-5)
    imp.efct = zeros(ir)
    j = 1
    for i in 1:ir
        qtn[i] && (imp.efct[i] = efct[j]; j += 1)
    end
    serialize(mmp, imp)
end

"""
    function initPedigree(xy, lmp, σₑ; fg = 0)
Initialize a pedigree of generation `fg` from genotypes in `xy`, and QTL
information in `lmp`. The pedigree returned is a DataFrame with column:
- `id`: ID of the individual
- `pa`: sire of the ID, default 0
- `ma`: dam of the ID, default 0
- `sex`: randomly assigned
- `grt`: generation of the ID, default 0
- `tbv`: true breeding value of the ID
- `pht`: phenotype of the ID, which allows missing
- `ebv`: estimated breeding value of the ID, default 0
- `F`: inbreeding coefficient of the ID, default 0
- `Fr`: inbreeding coefficient from reference markers, default 0
- `c`: contribution of the ID to the next generation, default 0
"""
function initPedigree(xy, lmp, σₑ; fg = 0)
    (:qtl ∈ propertynames(lmp) && :efct ∈ propertynames(lmp)) || error("No QTL column in $lmp")
    hdr = readhdr(xy)
    mt, _, mj, ir, ic = xyhdr(hdr)
    (mt == 'F' && mj == 0) || error("Not a proper input $xy")
    nid = ic ÷ 2
    tbv, pht = xy2gp(xy, 1:nid, lmp, σₑ)
    DataFrame(id = 1:nid,
              pa = 0,
              ma = 0,
              sex = rand(0:1, nid),
              grt = fg,
              tbv = tbv,
              pht = allowmissing(pht),
              ebv = 0.,
              F = 0.,
              Fr = 0.,
              c = 0.)
end

"""
    function randPed(nid, ng)
Random generation of a pedigree with `nid` individuals and `ng` generations.
Assuming there are half males and half females in each generation. The next
generation is generated by total random mating.
"""
function randPed(nid, ng)
    nsir = nid ÷ 2
    be = [1, nsir, nsir + 1, nid] # begin and end of sire and dam

    sex = [ones(Int, nsir); zeros(Int, nid - nsir)]
    ped = DataFrame(pa = zeros(Int, nid), ma = zeros(Int, nid), sex = sex, grt = 0)

    for ig in 1:ng
        pm = randomMate(be[1]:be[2], be[3]:be[4], noff=nid)
        append!(ped, DataFrame(pa = pm[:, 1], ma = pm[:, 2], sex = sex, grt = ig))
        be .+= nid
    end
    ped
end
