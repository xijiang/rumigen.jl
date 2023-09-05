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
        binapp!("$dir/fqf.bin", cnt) # frequency of frequency of QTL
    end
    qgt = nothing
    ref = isodd.(snp[lmp.ref, :])
    for i in 1:nhp:size(ref, 2)
        frq = sum(ref[:, i:i+nhp-1], dims=2)
        cnt = zeros(Int, nhp + 1)
        for x in frq
            cnt[x+1] += 1
        end
        binapp!("$dir/frf.bin", cnt) # frequency of frequency of reference
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
        bcr  = cor(grt.ebv, grt.tbv)
        sirs = grt.sex .== 1
        dams = grt.sex .== 0
        scr  = cor(grt.ebv[sirs], grt.tbv[sirs])
        dcr  = cor(grt.ebv[dams], grt.tbv[dams])
        np   = length(unique(grt.pa))
        nm   = length(unique(grt.ma))
        ncp  = sum(grt.c[sirs] .> 0)
        ncm  = sum(grt.c[dams] .> 0)
        push!(smp, (mbv, vbv, mF, mFr, bcr, scr, dcr, np, nm, ncp, ncm))
    end
    ideal, va, plst, nlst, pmls, nmls = idealPop("$rst/$xps/$bar-$sel.xy", ped.grt, lmp)
    open("$rst/$xps/$xps.bin", "a") do io
        write(io,
            smp.mbv, # 1
            smp.vbv, # as requested by SMS on 2023-05-21, by Theo
            smp.mF,  # 3
            smp.mFr, # 4
            smp.bcr, # 5
            smp.scr, # 6
            smp.dcr, # 7
            smp.np,  # 8
            smp.nm,  # 9
            smp.ncp, # 10
            smp.ncm, # 11
            ideal,   # 12
            va,      # 13
            plst,    # 14. n. of positive qtl lost
            nlst,    # 15
            pmls,    # 16. no. of positive qtl lost of maf 0.2
            nmls,    # 17. no. of negative qtl lost of maf 0.2
        )
    end
end
