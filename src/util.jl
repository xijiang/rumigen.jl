"""
    incidence_matrix(levels::AbstractVector)
Create an incidence matrix from a vector of effects.
"""
function incidence_matrix(levels::AbstractVector)
    n, u = length(levels), sort(unique(levels))
    m = length(u)
    x = zeros(n, m)
    for i in 1:m
        x[levels .== u[i], i] .= 1
    end
    x
end

"""
    incidence_matrix(df::DataFrame)
Create an incidence matrix from all columns of a data frame provided.
The first level of each factor is ignored. Instead, a intercept is added.
This is to make the matrix full rank.
"""
function incidence_matrix(df::DataFrame)
    n = size(df, 1)
    u = [sort(unique(df[:, i]))[2:end] for i in eachindex(names(df))] # can also be 1:end-1, which has more codes than 2:end
    m = sum([length(u[i]) for i in eachindex(u)])
    x, a = [ones(n) zeros(n, m)], 2
    for i in eachindex(u)
        for j in eachindex(u[i])
            x[df[:, i] .== u[i][j], a] .= 1
            a += 1
        end
    end
    sparse(x)
end

"""
    function Zmat(nm)
Given a vector of `Bool`sel indicating if a phenotype is not missing, return a
`Z` sparse matrix of `m` phenotypes and `n` ID, for an animal model.
"""
function Zmat(nm)
    n, m = length(nm), sum(nm)
    z = sparse(zeros(m, n))
    a = 1
    for i in 1:n
        if nm[i]
            z[a, i] = 1
            a += 1
        end
    end
    z
end

"""
Return the size of available memory on a Linux alike system.
Or, the free memory size is returned.
"""
function memavail()
    if Sys.isunix()
        for line in eachline("/proc/meminfo")
            if occursin("MemAv", line)
                return parse(Int, split(line)[2]) * 1024
            end
        end
    else
        return Sys.free_memory()
    end
end

"""
    function blksz(n::Int, m::Int)
Group `n` ID as evenly as possible with each group has maximally `m` ID.
The function returns the group size.
"""
function blksz(n::Int, m::Int)
    n % m == 0 ? m : begin
        ng = Int(ceil(n / m))
        n % ng == 0 ? n ÷ ng : n ÷ ng + 1
    end
end

# R, C, V: row, column and value specifying values in a sparse matrix
function pushRCV!(R, C, V, r, c, v)
    push!(R, r)
    push!(C, c)
    push!(V, v)
end

"""
    function histfrq(v::Vector, nb)
Return the frequency of `v` in `nb` bins.
"""
function histfrq(v::Vector, nb)
    l, h = extrema(v)
    w = (h - l) / nb
    tbl = zeros(Int, nb + 1)
    for x in v
        i = Int(floor((x - l) / w)) + 1
        tbl[i] += 1
    end
    tbl[end-1] += tbl[end]
    tbl[1:nb]
end

function binapp!(file, v::Vector)
    open(file, "a") do io
        write(io, v)
    end
end

"""
    function pos_qtl_frq(dir, bar, sls, ppsz)
Count the number of positive QTL alleles of each locus in each generation in
file `dir/bar-sel.xy` for each `sel` in `sls`. When QTL effect is positive,
count ones, or count zeros. QTL loci are defined in `dir/bar-map.ser`. Assume
that the population size is constant of `ppsz`.
"""
function pos_qtl_frq(rst, xps, bar, sel, ppsz)
    nhp = 2ppsz
    lmp = deserialize("$rst/$xps/$bar-map.ser")
    snp = xymap("$rst/$xps/$bar-$sel.xy")
    qgt = isodd.(snp[lmp.qtl, :])
    for i in 1:nhp:size(qgt, 2)
        frq = sum(qgt[:, i:i+nhp-1], dims=2)
        cnt = zeros(Int, nhp + 1)
        for x in frq
            cnt[x+1] += 1
        end
        binapp!("$rst/$xps/fqf.bin", cnt) # frequency of frequency of QTL
    end
    qgt = nothing
    ref = isodd.(snp[lmp.ref, :])
    for i in 1:nhp:size(ref, 2)
        frq = sum(ref[:, i:i+nhp-1], dims=2)
        cnt = zeros(Int, nhp + 1)
        for x in frq
            cnt[x+1] += 1
        end
        binapp!("$rst/$xps/frf.bin", cnt) # frequency of frequency of reference
    end
    ref = nothing
    chp = isodd.(snp[lmp.chip, :])
    for i in 1:nhp:size(chp, 2)
        frq = sum(chp[:, i:i+nhp-1], dims=2)
        cnt = zeros(Int, nhp + 1)
        for x in frq
            cnt[x+1] += 1
        end
        binapp!("$rst/$xps/fcf.bin", cnt) # frequency of frequency of chip
    end
end

"""
    function randomc(sex, nsir, ndam)
Given a list of `sex`, where males are indicated as `1`, females are indicated
as `0`, this function returns a vector of `c` values for each animal. `nsir` and
`ndam` are the number of sires and dams to be randomly selected. The `c` values
are randomly selected from a uniform distribution. The sum of `c` values for
selected sires and dams are normalized to `0.5`, respectively. The rest of the
animals have `c` values of `0`.
"""
function randomc(sex, nsir, ndam)
    c = zeros(length(sex))
    sir = findall(sex .== 1)
    dam = findall(sex .== 0)
    cs = rand(nsir)
    cs /= (2sum(cs))
    cd = rand(ndam)
    cd /= (2sum(cd))
    c[shuffle(sir)[1:nsir]] = cs
    c[shuffle(dam)[1:ndam]] = cd
    c
end

function equalc(sex, nsir, ndam)
    c = zeros(length(sex))
    sir = findall(sex .== 1)
    dam = findall(sex .== 0)
    c[shuffle(sir)[1:nsir]] .= 0.5 / nsir
    c[shuffle(dam)[1:ndam]] .= 0.5 / ndam
    c
end

function sumPed(rst, xps, bar, lmp, sel)
    ped = deserialize("$rst/$xps/$bar-$sel+ped.ser")
    smp = DataFrame(
        mbv = Float64[], # mean breeding value
        vbv = Float64[], # variance of breeding values
        mF  = Float64[], # mean inbreeding coefficient
        mFr = Float64[], # mean inbreeding coefficient from relatives loci
        mFp = Float64[], # mean inbreeding coefficient from pedigree
        bcr = Float64[], # cor(EBV, TBV) both sexes
        scr = Float64[], # cor(EBV, TBV) sires
        dcr = Float64[], # cor(EBV, TBV) dams
        np  = Float64[], # number of sires in each generation
        nm  = Float64[], # number of dams in each generation
        ncp = Float64[], # number of positive sire c values in each generation
        ncm = Float64[]) # number of positive dam c values in each generation
    for grt in groupby(ped, :grt)
        mbv  = mean(grt.tbv)
        vbv  = var(grt.tbv)
        mF   = mean(grt.F)
        mFr  = mean(grt.Fr)
        mFp  = mean(grt.Fp)
        bcr  = cor(grt.ebv, grt.tbv)
        sirs = grt.sex .== 1
        dams = grt.sex .== 0
        scr  = cor(grt.ebv[sirs], grt.tbv[sirs])
        dcr  = cor(grt.ebv[dams], grt.tbv[dams])
        np   = length(unique(grt.pa))
        nm   = length(unique(grt.ma))
        ncp  = sum(grt.c[sirs] .> 0)
        ncm  = sum(grt.c[dams] .> 0)
        push!(smp, (mbv, vbv, mF, mFr, mFp, bcr, scr, dcr, np, nm, ncp, ncm))
    end
    ips = idealPop("$rst/$xps/$bar-$sel.xy", ped.grt, lmp)
    open("$rst/$xps/$xps.bin", "a") do io
        write(io,
            smp.mbv, # 1
            smp.vbv, # as requested by SMS on 2023-05-21, by Theo
            smp.mF,  # 3
            smp.mFr, # 4
            smp.mFp, # 5
            smp.bcr, # 6
            smp.scr, # 7
            smp.dcr, # 8
            smp.np,  # 9
            smp.nm,  # 10
            smp.ncp, # 11
            smp.ncm, # 12
            ips.ideal,   # 13
            ips.va,      # 14
            ips.np,      # 15. no. of positive qtl fixed
            ips.nn,      # 16. no. of negative qtl fixed
            ips.nmp,     # 17. no. of positive qtl fixed of maf 0.2
            ips.nmn,     # 18. no. of negative qtl fixed of maf 0.2
            ips.clst,    # 19. no. of chip snps fixed
            ips.rlst,    # 20. no. of reference snps fixed
            ips.cmls,    # 21. no. of chip snps fixed of maf 0.2
            ips.rmls,    # 22. no. of reference snps fixed of maf 0.2
            ips.cvr,     # 23. covariance begween generation 0 and 20, reference
            ips.cvc,     # 24. covariance begween generation 0 and 20, chip
            ips.cvq,     # 25. covariance begween generation 0 and 20, qtl
            )
    end
end

"""
    function prepFdr(rst::AbstractString, xps::AbstractString, quick=false)
Prepare a founder population. If `quick` is `true`, the founder population is
prepared from the test-suite. Otherwise, the founder population is prepared from
the `founder` folder.
"""
function prepFdr(rst::AbstractString,
    quick,
    ppsz,
    )
    fdr, foo = "$rst/test-suite", "founder"
    if quick
        file = "$rst/test-suite/founder-hap.xy"
        create = !isfile(file)
        if isfile(file)
            a = zeros(Int, 3)
            read!(file, a)
            create = a[2] < 2ppsz
        end
        if create
            baz = cattle_base(ppsz, fdr)
            mv("$rst/test-suite/$baz-map.ser", "$rst/test-suite/founder-map.ser", force=true)
            mv("$rst/test-suite/$baz-hap.xy", "$rst/test-suite/founder-hap.xy", force=true)
        end
    else
        fdr = "$rst/founder"
        isdir(fdr) && rm(fdr, recursive=true, force=true)
        mkpath(fdr)
        foo = cattle_base(ppsz, fdr)
    end
    fdr, foo
end

function updateIBDM(xy, bin, snp, mid, nid)
    tid = mid + nid
    ra, rb = 1:mid, mid+1:tid
    G = zeros(tid, tid)
    read!(bin, view(G, ra, ra))  # read the upper left block
    K = gametemat(xy, snp, ra, rb)
    copyto!(view(G, ra, rb), K)  # copy the upper right block
    copyto!(view(G, rb, ra), K') # copy the lower left block
    K = gametemat(xy, snp, rb)
    copyto!(view(G, rb, rb), K)  # copy the lower right block
    write(bin, G)
end

function commas(num::Integer)
    str = string(num)
    return replace(str, r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
end

function Ceiling(;
    nQTL = 1000,
    nRepeat = 1000,
    ad = Normal(),         # can be of other distributions
    qd = Beta(0.75, 0.75), # allele frequency distribution
)
    C = 0.0        # Ceiling height
    for _ in 1:nRepeat
        a = rand(ad, nQTL) # effects of SNP 1s
        p = rand(qd, nQTL) # U-shaped allele frequencies
        σ = sqrt(sum(2p .* (1 .- p) .* a .^ 2)) # genic std
        C += 2sum(a[a.>0.0]) / σ
    end
    C /= nRepeat
end

function nfix(frq, F, maf)
    n = 0.
    for q in frq
        q < maf && (n += exp(-2q / F))
    end
    Int(floor(n))d
end

"""
    function grtfrq(xy, grt)
Given a `xy` file and a vector `grt` showing haplotype generations, this
function counts the allele frequency of each locus in each generation. The
allele frequency is calculated from the `xy` file. The function returns a matrix
of `nlc x ngrt` matrix.

Note: This is a quick and dirty one. As I know the eltype of `xy` is `Int32`.
This will be changed in the future.
"""
function grtfrq(xy, grt)
    dims = begin
        tmp = zeros(Int, 3)
        read!(xy, tmp)
        (tmp[2], tmp[3])
    end
    gs = unique(grt)
    ng = length(gs)
    frq = Int[]
    gt = Mmap.mmap(xy, Matrix{Int16}, dims, 24)
    for g in gs
        cg = grt .== g
        cgt = isodd.(gt[:, cg])
        append!(frq, sum(cgt, dims=2))
    end
    reshape(frq, dims[1], ng)
end