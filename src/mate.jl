"""
    function randomMate(nsire::Int, ndam::Int; noff = 0)
Random pair sire `1:nsire` and dam `nsire+1:nsire+ndam-1`.
Returns a sorted two-column matrix, `[sires dams]`.
If `noff` is given, it is the number of offspring to be generated.
"""
function randomMate(nsir::Int, ndam::Int; noff = 0)
    noff == 0 && (noff = max(nsir, ndam))
    pa, ma = shuffle(1:nsir), shuffle(1:ndam)
    while length(pa) < noff
        append!(pa, shuffle(1:nsir))
    end
    while length(ma) < noff
        append!(ma, shuffle(1:ndam))
    end
    sortslices([pa[1:noff] ma[1:noff] .+ nsir], dims=1, by=x -> (x[1], x[2]))
end

"""
    function randomMate(sires, dams; noff = 0)
Given names listed in `sires` and `dams`, this function randomly match
them, such that every ID is used as much as possible.  Returns a
sorted two columns name matrix.
"""
function randomMate(sires, dams; noff = 0)
    noff == 0 && (noff = max(length(sires), length(dams)))
    pa, ma = shuffle(sires), shuffle(dams)
    while length(pa) < noff
        append!(pa, shuffle(sires))
    end
    while length(ma) < noff
        append!(ma, shuffle(dams))
    end
    sortslices([pa[1:noff] ma[1:noff]], dims=1, by=x -> (x[1], x[2]))
end

"""
    function randomMate(ped::DataFrame, noff)
Create mate pairs according to ID's contribution stated in `ped.c`.  `ped.c` is
calculated from `myopt` function, where `sum(c[sires]) == 0.5`, and
`sum(c[dams]) == 0.5`.  `noff` is the number of offspring to be generated.
Returns the mating pairs. Here it suppose that a female can mate several males,
for example, with MOET.
"""
function randomMate(ped::DataFrame, noff)
    function nrpt(c, noff)
        m = noff / sum(c)
        nr = Int.(round.(c * m))
        x = sum(nr) - noff
        if x > 0  # to make sure sum(nr) == noff
            nr[argmax(nr)] -= x
        else
            nr[argmin(nr)] -= x
        end
        nr
    end

    ["sex", "c"] ⊆ names(ped) || error("Not a proper pedigree")
    nid = size(ped, 1)
    ped.id = 1:nid
    gd = groupby(ped, :sex) # group 1 for females and 2 for males.
    dam = gd[1].id[findall(gd[1].c .> 0)]
    sir = gd[2].id[findall(gd[2].c .> 0)]
    pa, ma = Int[], Int[]
    nr = nrpt(ped.c[sir], noff)
    for (id, r) ∈ eachrow([sir nr])
        append!(pa, repeat([id], r))
    end
    nr = nrpt(ped.c[dam], noff)
    for (id, r) ∈ eachrow([dam nr])
        append!(ma, repeat([id], r))
    end
    sortslices([pa shuffle(ma)], dims=1, by=x -> (x[1], x[2]))
end

"""
    function crossovers(lms)
Give a linkage map summary `lms`, which can be from `summap`, this
function return a vector of crossover points along the whole genome.

DataFrame `lms` has for columns, chromosome number, its length in
Morgen, number of loci it contains, and the beginning number of its
first locus in the whole genome.

The first number is `1`, or `2`, the starting haplotype.  The vector
is then start from locus `1`, and loci that cross-over happens.  A
"cross-over" may also happen at the first locus of a chromsome.
Otherwise, the segment continues with the previous haplotype.  This
eases the segment copy.  At the beginning of a "cross-over" happens on
`rand(false:true)`.  The number of cross-over occurred on a chromosome
follows a Poisson distribution.

The cross-over location is then sampled from a uniform distribution.
It can can also be U-shaped abouth centromere, which might be
implemented later.

The chromosome, or linkage group, should be in the order of the
genotype file.  Or, unpredicted errors can be resulted.

## Caution!
This is better for uniformly distributed many SNP.  In the future, I
will take into account the different location of SNP, hotspots etc.
"""
function crossovers(lms)
    pts = [rand(1:2)]           # crossover points
    for (_, λ, nlc, bgn) in eachrow(lms)
        bgn > 1 && rand(false:true) && push!(pts, bgn)
        nc = rand(Poisson(λ))
        append!(pts, rand(1:nlc, nc) .+ (bgn - 1))
    end
    push!(pts, last(lms).nlc - 1 + last(lms.bgn))
    pts
end

"""
    function sumMap(lmp; cM = 1e6)
Summary of a linkage map, which has 3 columns, :chr, :pos(in bp), and :frq.
A DataFrame of 4-column for each chromosome is returned:
- numbering (1, 2, ...)
- length in Morgen
- number of loci
- beginning number of the first locus
"""
function sumMap(lmp; cM=1e6)
    df = DataFrame(chr=Int8[],
        len=Float64[], # as λ
        nlc=Int[],
        bgn=Int[])
    bgn = 1
    for grp in groupby(lmp, :chr)
        chr = first(grp).chr
        len = (last(grp).pos - first(grp).pos) / cM / 100
        nlc = size(grp, 1)
        push!(df, (chr, len, nlc, bgn))
        bgn += nlc
    end
    return df
end

"""
    function gamete(prt, hap, lms)
Generate a gamete `hap`, a **vector** view of a child's haplotype,
from `prt`, a view of a parents genotypes,
according crossovers generated from a linkage map summary, `lms`.
"""
function gamete(prt, hap, lms)
    cvs = crossovers(lms)
    h, l = cvs[1], 1            # starting
    for cv in cvs[2:end]
        copyto!(view(hap, l:cv), view(prt, l:cv, h))
        l = cv + 1
        h = 3 - h
    end
end

"""
    function drop(pg::AbstractArray, og::AbstractArray, pm, lms)
Drop haplotypes `pg` of parents into `og`, their offspring genotypes.
Parents of each offspring are defined in `pm`, which are rows of ``pa ma``.
Linkage map summary `lms` is from `summap`.

!!! ``Caution``:
- Merged data matrix from `MaCS` is `n-ID × n-loci`. Treat it with care.
- It is supposed all the genes to drop from are in `pg`.
- This function will be removed in the future.
"""
function drop(pg::AbstractArray, og::AbstractArray, pm, lms)
    nf = size(pm)[1]
    Threads.@threads for id in 1:nf
        ip = pm[id, 1]
        pa = view(pg, :, 2ip-1:2ip)
        zi = vec(view(og, :, 2id - 1))
        gamete(pa, zi, lms)
    end
    Threads.@threads for id in 1:nf
        im = pm[id, 2]
        ma = view(pg, :, 2im-1:2im)
        zi = vec(view(og, :, 2id))
        gamete(ma, zi, lms)
    end
end

"""
    function randrop!(gt, ped)
Drop parents' alleles in `gt` into their offspring according to pedigree `ped`.
Assuming no linkage among loci. The empirical `F` is calculated and returned.
"""
function randrop!(gt, ped)
    function _drop(gt, p, oh)
        pg = view(gt, :, 2p-1:2p)
        for i in eachindex(oh)
            oh[i] = pg[i, rand(1:2)]
        end
    end
    (nlc, nhp), nid = size(gt), size(ped, 1)
    nhp == 2nid || error("Haplotypes do not match pedigree")
    F = zeros(nid)
    for id in 1:nid
        pa, ma = ped.pa[id], ped.ma[id]
        pa > 0 && _drop(gt, pa, view(gt, :, 2id-1))
        ma > 0 && _drop(gt, ma, view(gt, :, 2id))
        pa * ma > 0 && (F[id] = sum(gt[:, 2id-1] .== gt[:, 2id]) / nlc)
    end
    F
end

"""
    function drop(fra::AbstractString, til::AbstractString, pm, lms)
Drop haplotypes from haplotypes in file `fra` of xy format, into 
another xy file `til`, according to sire and dam info in `pm` and 
linkage summary `lms`.
"""
function drop(fra::AbstractString, til::AbstractString, pm, lms)
    ihdr = readhdr(fra)
    mt, et, mj, ir, ic = xyhdr(ihdr)
    mt == 'F' || error("Only support F matrix")
    pg = Mmap.mmap(fra, Matrix{et}, (ir, ic), 24)
    noff = size(pm, 1)
    open(til, "w") do io
        ohdr = xyheader('X','Y',' ',mt,mj,ihdr.e,'\n','\n',ir,2noff)
        write(io, Ref(ohdr))
        og = zeros(et, ir, noff * 2)
        drop(pg, og, pm, lms)
        write(io, og)
    end
    nothing
end

"""
    function drop(pg::AbstractArray, pm, lms)
Drop genotyeps in `pg` of ``n-loci × n-Haplotypes`` into a new matrix
`og`, according to sire and dam info in `pm` (``pa ma`` of n-ID rows),
and linkage summary `lms`.  The new matrix `og` is returned.
"""
function drop(pg::AbstractArray, pm, lms)
    nid, nlc = size(pm, 1), size(pg, 1)
    og = zeros(Int8, nlc, nid * 2)

    Threads.@threads for id in 1:nid
        ip = pm[id, 1]
        pa = view(pg, :, 2ip-1:2ip)
        zi = vec(view(og, :, 2id - 1))
        gamete(pa, zi, lms)
    end
    Threads.@threads for id in 1:nid
        im = pm[id, 2]
        ma = view(pg, :, 2im-1:2im)
        zi = vec(view(og, :, 2id))
        gamete(ma, zi, lms)
    end
    og
end

"""
    function drop_by_chr(fph::String, bar::String, pm, lms; merge = true)
Drop haplotypes of parents in `fph` into `bar-{1..nchr}.bin`,
their offspring genotypes.
Matrix in `fph` should be of dimension `nLoci × nHap`.
Parents of each offspring are defined in `pm`, which are rows of ``pa ma``.
Linkage map summary `lms` is from `summap` of module ``Sim``.

When `merge = true`, the dropped haplotypes will be merged into genotypes.
It is also transposed before written to a file.
"""
function drop_by_chr(fph::String, bar::String, pm, lms; merge=false)
    nlc, nhp = Fio.readdim(fph)
    nof = size(pm)[1]
    mph = Mmap.mmap(fph, Matrix{Int8}, (nlc, nhp), 24)

    tprintln("    - Dropping on chromosome:")
    for (chr, len, cln, fra) in eachrow(lms)
        tprint(" $chr")
        til = cln + fra - 1
        oh = zeros(Int8, cln, 2nof)
        ph = copy(mph[fra:til, :])
        drop(ph, oh, pm, DataFrame(chr=chr, len=len, nlc=cln, bgn=1))
        if merge
            open("$bar-$chr.bin", "w") do io
                # write(io, [cln, nof, Fio.typec(Int8)])
                # for i in 1:nof
                #     write(io, oh[:, 2i-1] + oh[:, 2i])
                # end
                write(io, [nof, cln, Fio.typec(Int8)])
                for i in 1:cln
                    write(io, oh[i, 1:2:2nof] + oh[i, 2:2:2nof])
                end
            end
        else
            Fio.writemat("$bar-$chr.bin", oh)
        end
    end
    println()
end

"""
    function reproduce()
Give a matrix `haps` of haplotypes of nHap (= 2nID) by nLoc, and a
pedigree `ped`, (`[ID pa ma]`) this function drop the parents
haplotype according to the crossovers generated on linkage map summary
`lms`.
"""
function reproduce(haps, ped, lms)
    for (id, ip, im) in eachrow(ped)
        pa = view(haps, 2ip-1:2ip, :)
        zi = vec(view(haps, 2id - 1, :)) # offspring
        gamete(pa, zi, lms)
        ma = view(haps, 2im-1:2im, :)
        zi = vec(view(haps, 2id, :))
        gamete(ma, zi, lms)
    end
end

"""
    function collect_gt(bar, lms, loci; rev = false)
In case when genotypes are stored in different file names started with
`bar` by chromosomes, This function collect the genotypes on `loci` in
total genome order.  Linkage map summary `lms` which file to find
according to number in `loci`.

Note, the genotypes are of `nid × nlc` of each chromosome.  When `rev
= true`, the genotypes are collected in reverse order of chrosome
numbers.
"""
function collect_gt(bar, lms, loci; rev=false)
    nid, _ = Fio.readdim("$bar-$(lms.chr[1]).bin")
    # below commented codes are for matrix of `nlc × nid`.
    #_, nid = Fio.readdim("$bar-$(lms.chr[1]).bin")
    gt = zeros(Int8, length(loci), nid)
    i, _lms = rev ? (length(loci), reverse(lms)) : (0, lms)
    for (chr, _, nlc, fra) in eachrow(_lms)
        chunk = intersect(fra:fra+nlc-1, loci) .+ 1 .- fra
        slc = length(chunk)     # number of selected loci
        rev && (i -= slc)
        open("$bar-$chr.bin", "r") do ig
            cg = Mmap.mmap(ig, Matrix{Int8}, (nid, nlc), 24)
            copyto!(view(gt, i+1:i+slc, :), view(cg, :, chunk)')
            # cg = Mmap.mmap(ig, Matrix{Int8}, (nlc, nid), 24)
            # copyto!(view(gt, i+1:i+slc, :), view(cg, chunk, :))
        end
        rev || (i += slc)
    end
    gt
end

"""
    function drop(php::String, gt::String, pa::AbstractVector{Int}, ma::AbstractVector{Int},
                  nsb::Int, lms, ped, ohp::String)
Drop parent haplotypes in file `php` into file `gt` according to
random mating among (unique) sires `pa` and dams `ma`.  Crossovers
were sampled according to linkage map summary `lms`.  Each full
sibship size is `nsb`.  When `gt` exists, new genotypes will be
appended to it, or a new file `gt` is created.

The function will check if 2length of `pa` and `ma` is less or equal
to number of columns in `php`.  It is assumed haplotypes of `pa` are
packed in the first half columns of the `php`, and in that order.  The
function also check if any intersect of `pa` and `ma` exists.

The newly created haplotypes are written in to `ohp`.
"""
function drop(php::String, gt::String, pa::AbstractVector{Int}, ma::AbstractVector{Int},
    nsb::Int, lms, ped, ohp::String)
    ig = ped.grt[end] + 1
    tprintln("  - Dropping into generation $ig")
    # The drop procedure
    pg = Fio.readmat(php)
    nsr, ndm = length(pa), length(ma)
    2(nsr + ndm) ≤ size(pg)[2] || error("Not enough haplotypes in $hap for list pa and ma")
    length(intersect(pa, ma)) > 0 && error("Some ID are both sire and dam")
    pm = repeat(Sim.random_mate(nsr, ndm), inner=(nsb, 1))
    og = drop(pg, pm, lms)
    Fio.writemat(ohp, og)

    # update pedigree
    nlc, nid = size(og)
    nid ÷= 2
    ndum = size(ped, 2) - 4
    sex = rand(0:1, nid)
    for id in 1:nid
        push!(ped, (ig, pa[pm[id, 1]], ma[pm[id, 2]-nsr], sex[id], zeros(ndum)...))
    end

    # Create genotypes file, or append new genotypes.
    isfile(gt) || write(gt, [nlc, 0, Fio.typec(Int8)])

    # update header
    hd = zeros(Int, 2)
    read!(gt, hd)
    open(gt, "a+") do io
        hd[2] += nid
        seekstart(io)
        write(io, hd)
        seekend(io)
        for i in 1:nid
            write(io, og[:, 2i-1] + og[:, 2i])
        end
    end
end
