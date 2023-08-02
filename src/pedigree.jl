"""
    function kinship(ped, i, j)
---
This function is handy if just to calculate relationship of a few (pairs of) ID.
The first 2 columns of ped must be `pa` and `ma`.
It can also speed up by adding `Thread.@threads` before your pair loop.
"""
function kinship(ped, i, j)
    (i == 0 || j == 0) && return 0
    ipa, ima = ped[i, :]          # used both for below and the last
    i == j && (return 1 + .5kinship(ped, ipa, ima))
    if i < j
        jpa, jma = ped[j, :]
        return .5(kinship(ped, i, jpa) + kinship(ped, i, jma))
    end
    return .5(kinship(ped, j, ipa) + kinship(ped, j, ima))
end

"""
    function kinship(ped, i::Int, j::Int, dic::Dict{Tuple{Int, Int}, Float64})
Recursive kinship calculation with kinship of ID pair `(i, j)`
stored in dictionary `dic`.  The first 2 columns of ped must be `pa` and `ma`.
The memory usage may be bigger than Meuwissen and Luo 1992, or Quaas 1995.
The speed is however standable.
The recursive algorithm is also easy to understand.
"""
function kinship(ped, i::Int, j::Int, dic::Dict{Tuple{Int, Int}, Float64})
    (i == 0 || j == 0) && return 0
    ip, im = ped[i, :]
    if i == j
        haskey(dic, (i, i)) || (dic[(i, i)] = 1 + .5kinship(ped, ip, im, dic))
        return dic[(i, i)]
    end
    if i < j
        jp, jm = ped[j, :]
        haskey(dic, (i, jp)) || (dic[(i, jp)] = kinship(ped, i, jp, dic))
        haskey(dic, (i, jm)) || (dic[(i, jm)] = kinship(ped, i, jm, dic))
        return .5(dic[(i, jp)] + dic[(i, jm)])
    end
    return .5(kinship(ped, j, ip, dic) + kinship(ped, j, im, dic))
end

"""
    function ped_F(ped; force = false)
- When column `:F` is not in DataFrame `ped`, this function calculate
inbreeding coefficients of every ID, and add a `:F` column in `ped`.
- If `:F` exists, and `force = true`, drop `:F` from `ped` and do above.
- Or do nothing.
"""
function ped_F(ped; force = false)
    if "F" ∈ names(ped)
        force ? select!(ped, Not([:F])) : return
    end
    F = [0. for _ in eachrow(ped)]
    dic = Relation()
    for i in eachrow(ped)
        # for the `dic` here, no need to consider full siblings
        F[i] = kinship(ped, i, i, dic) - 1
        i % 200 == 0 && print("\t$i")
    end
    ped.F = F
    nothing
end

"""
    function ped_D(ped; force = false)
Calculate diagonals of `D` for `A` or `A⁻¹` calculation.
It is computed as a vector and a new column of `ped` of name `:D`.
"""
function ped_D(ped; force = false)
    if "D" ∈ names(ped)
        force ? select!(ped, Not([:D])) : return
    end
    #"F" ∈ names(ped) || ped_F(ped)
    N = size(ped, 1)
    "F" ∈ names(ped) || (ped.F = diag(Amat(ped, m = N)))
    D = .5ones(N)
    for i in eachrow(ped)
        pa, ma = ped[i, :]
        vp = (pa == 0) ? -.25 : .25ped.F[pa]
        vm = (ma == 0) ? -.25 : .25ped.F[ma]
        D[i] -= (vp + vm)
    end
    ped.D = D
    nothing
end

"""
    function T4A⁻¹(ped)
Give a pedigree DataFrame, with its first 2 column as `pa`, and `ma`,
this function return the T matrix used for A⁻¹ calculation.
"""
function T4A⁻¹(ped)
    R, C, V = Int[], Int[], Float64[]
    for (id, (pa, ma)) in enumerate(eachrow(ped))
        pa > 0 && pushRCV!(R, C, V, id, pa, -.5)
        ma > 0 && pushRCV!(R, C, V, id, ma, -.5)
        pushRCV!(R, C, V, id, id, 1.)
    end
    return sparse(R, C, V)
end

"""
    function Amat(ped; m = 1000)
Given a pedigree `ped`,
this function returns a full numerical relationship matrix, `A`.
This function is better for small pedigrees, and for demonstration only.
The maximal matrix size is thus limited to 1000.
One can try to set `m` to a bigger value if RAM is enough.
"""
function Amat(ped; m = 30_000)
    N = size(ped, 1)
    N > m && error("Pedigree size ($N > $m) too big")
    A = zeros(N, N) + I(N)
    for (i, (pa, ma)) in enumerate(eachrow(select(ped, :pa, :ma)))
        pa * ma ≠ 0 && (A[i, i] += .5A[pa, ma])
        for j in 1:i-1
            pa ≠ 0 && (A[i, j]  = 0.5A[j, pa])
            ma ≠ 0 && (A[i, j] += 0.5A[j, ma])
            A[j, i] = A[i, j]
        end
    end
    A
end

"""
    function A⁻¹(ped)
Given a pedigree `ped`, this function return a sparse matrix of `A⁻¹`,
where `A` is the numerical relationship matrix.
"""
function A⁻¹(ped)
    tpd = select(ped, :pa, :ma)
    T = T4A⁻¹(tpd)
    "D" ∈ names(tpd) || ped_D(tpd)
    D = Diagonal(1. ./ tpd.D)
    T'D*T
end

"""
    function grm(gt)
Given the genotypes of `Matrix{Int8}`, this function calculate the genomic 
relationship matrix `GRM`. If the matrix is too big, the content will be 
calculated block by block and written to a file, else it will return a 
`Matrix{Float64}`.

By default the matrix needs to be 'nlc by nid' to speed up the calculation.
Such a matrix stores the genotypes continuously, which can speed up the 
matrix multiplication very much.  Another reason is that genotypes are
usually stored continuous for each individual.  They can can be read 
continuously in the `gt` columns.

If a `δ`, e.g., `δ = 0.01`, to the diagonals, you have to do this after this
function.
"""
function grm(gt::AbstractArray)
    p = mean(gt, dims = 2) ./ 2 # allele frequencies
    d = 2(1 .- p)'p             # the denominator
    nlc, nid = size(gt)
    mem = memavail() * 99 ÷ 100 # not all available memory
    gmt = nid^2 * 8             # memory by G
    zmt = nid * nlc * 8         # memory by Z
    if gmt + zmt < mem          # brute force
        # @info "G and Z are stored in memory"
        Z = gt .- 2p
        G = Z'Z ./ d
        return G
    else                        # minimal memory mode
        c1 = 2gt'p
        c2 = 4p'p
        if gmt < mem            # G can still be fit
            @info "only G were stored in memory"
            G = zeros(nid, nid)
            matmul!(G, gt', gt)
            G .-= c1
            G .-= c1'
            G .+= c2
            G ./= d
            return G
        else                            # G is too large
            file = basename(tempname()) # will write the result in pwd.
            @warn "G is too big to fit in memory. It is being writting into $file.
              False will be returned. $file can be read back in to memory, if enough,
              with `QTL.MIO.readmat($file)`"
            # ToDo: check disk space here
            m = mem ÷ 8 ÷ nid
            m = blksz(nid, m) # determine number of ID to be dealed a time
            stops = collect(m:m:nid)
            stops[end] == nid || push!(stops, nid)
            start = 1
            open(file, "w") do io
                write(io, [nid, nid, 13]) # 13 for Float64
                for stop in stops
                    sg = zeros(nid, stop - start + 1)
                    matmul!(sg, gt', gt[:, start:stop])
                    sg .-= c1
                    sg .-= c1[start:stop]'
                    sg .+= c2
                    sg ./= d
                    write(io, sg)
                    start = stop + 1
                end
            end
            return file
        end
    end
end

"""
    function grmiv(xy)
Calculate the genomic relationship matrix from file `xy`.
Returns its inverse.
"""
function grmiv(xy::AbstractString, chip)
    hdr = readhdr(xy)
    mt, et, mj, ir, ic = xyhdr(hdr)
    mt == 'F' || error("File $xy is not a genotype file")
    (et == Int8 || et == UInt16) || error("File $xy is not a valid genotype file")
    hps = if mj == 0
        Mmap.mmap(xy, Matrix{et}, (ir, ic), 24)
    else
        Mmap.mmap(xy, Matrix{et}, (ir, ic), 24)'
    end
    snps = view(hps, chip, :)
    gt = if et == Int8
        snps[:, 1:2:end] + snps[:, 2:2:end]
    else
        Int8.(isodd.(snps[:, 1:2:end])) + Int8.(isodd.(snps[:, 2:2:end]))
    end
    g = grm(gt) + 0.01I
    g isa AbstractString && (@error "Not enough memory to calculate GRM")
    LAPACK.potrf!('L', g)
    LAPACK.potri!('L', g)
    for i in 2:size(g, 1)
        g[i-1, i:end] = g[i:end, i-1]
    end
    g
end

"""
    function grm(xy::AbstractString, chip)
Calculate the genomic relationship matrix from file `xy`, and observable SNP information from `chip`.
"""
function grm(xy::AbstractString, chip)
    hdr = readhdr(xy)
    mt, et, mj, ir, ic = xyhdr(hdr)
    mt == 'F' || error("File $xy is not a genotype file")
    (et == Int8 || et == UInt16) || error("File $xy is not a valid genotype file")
    hps = (mj == 0) ? Mmap.mmap(xy, Matrix{et}, (ir, ic), 24) :
        Mmap.mmap(xy, Matrix{et}, (ir, ic), 24)'
    snps = view(hps, chip, :)
    gt = (et == Int8) ? snps[:, 1:2:end] + snps[:, 2:2:end] :
        Int8.(isodd.(snps[:, 1:2:end])) + Int8.(isodd.(snps[:, 2:2:end]))
    grm(gt) + 0.01I
end

"""
    function A⁻¹(ped, ser)
If `ped` is an expanded pedigree. The `A⁻¹` of `ped`'s previous version was
calculated. This function will reuse the previous results of `D`, `RCV` for `T`,
`F` and relationships dictionary serialized in `ser` to calculate `A⁻¹` of
`ped`. New results will be serialized in `ser` for future use.

## Notes
- Note that in `ped`, the row number is the same as ID number.  Unkown parents
are zeros.
- New ID numbers appended must be larger than the previous ones.
"""
function A⁻¹(ped, ser)
    tpd = select(ped, :pa, :ma)

    R, C, V, D = if isfile(ser)
        deserialize(ser)
    else
        Int64[], Int64[], Float64[], Float64[]
    end
    pid, nid = length(D), size(tpd, 1)

    # New T
    for (id, (pa, ma)) in enumerate(eachrow(tpd))
        id ≤ pid && continue
        pa > 0 && pushRCV!(R, C, V, id, pa, -0.5)
        ma > 0 && pushRCV!(R, C, V, id, ma, -0.5)
        pushRCV!(R, C, V, id, id, 1.0)
    end
    T = sparse(R, C, V)

    F = diag(Amat(ped, m = nid)) .- 1

    # update D
    for i in pid+1:nid
        pa, ma = tpd[i, :]
        vp = (pa == 0) ? -.25 : .25F[pa]
        vm = (ma == 0) ? -.25 : .25F[ma]
        push!(D, .5 - vp - vm)
    end

    di = Diagonal(1. ./ D)

    serialize(ser, (R, C, V, D))
    T'di*T
end

function initUSNP(ped; nlc = 1000)
    tpd, nid = select(ped, :pa, :ma), size(ped, 1)
    nfdr = sum(iszero, tpd.pa) + sum(iszero, tpd.ma)
    et = nfdr < 256 ? UInt8 : UInt16
    gt, usnp = zeros(et, nlc, 2nid), 1
    for id in 1:nid
        pa, ma = tpd[id, :]
        if pa == 0
            gt[:, 2id-1] .= usnp
            usnp += 1
        end
        if ma == 0
            gt[:, 2id]   .= usnp
            usnp += 1
        end
    end
    gt
end

"""
    function fastF(ped; nlc = 1000, ϵ = 1e-5, inc = 10)
This function return the exact `F` for the first 10 generations. For the rest
generations, it simulate 1000 independent loci of unique alleles and drop them
according to `ped`. It then count IBD for each ID. This process is repeated
untill the mean std of the approximate `F` mean of the first 10 generations is
smaller than `ϵ = 1e-5`.

## Notes

- The pedigree must to be sorted such that offspring appear after their parents.
- The pedigree must be coded such that row number is the same as ID number.
- Unkown parents are zeros. They are counted for unique alleles.

! ToDo: not finished. may use threads later.
"""
function fastF(ped; nlc = 1000, ϵ = 1e-3, inc = 10)
    grt = unique(ped.grt)
    ngrt, nid = length(grt), size(ped, 1)
    eid = ngrt ≤ 10 ? nid : findfirst(x -> x == grt[11], ped.grt) - 1

    F, dic = zeros(nid), Relation() # exact F for the first 10 generations
    for i in 1:eid
        F[i] = kinship(ped, i, i, dic) - 1
    end
    ped.F = F

    ngrt ≤ 10 && return ped

    # approximate F for the rest generations
    gt, vf = initUSNP!(ped, nlc = nlc), 1.
    randrop(gt, ped) # check pedigree
    any(gt .== 0) && error("Pedigree has wrong ordered IDs.")
    mf, sf = zeros(eid), zeros(eid)
    while vf > ϵ
        for _ in 1:inc
            tf = randrop!(gt, ped)
            tf = tf[1:eid] - F[1:eid]
            mf += tf
            sf += tf.^2
        end
        vf = (sum(sf) - nid * mean(mf)^2) / (nid - 1)
    end
end
