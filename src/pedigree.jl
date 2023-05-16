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
    N = nrow(ped)
    F = zeros(N)
    dic = Dict{Tuple{Int, Int}, Float64}()
    for i in 1:N
        F[i] = kinship(ped, i, i, dic) - 1
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
    "F" ∈ names(ped) || ped_F(ped)
    N = nrow(ped)
    D = .5ones(N)
    for i in 1:N
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
    # R, C, V: row, column and value specifying values in a sparse matrix
    function pushRCV!(R, C, V, r, c, v)
        push!(R, r)
        push!(C, c)
        push!(V, v)
    end
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
function Amat(ped; m = 1000)
    N = nrow(ped)
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
    inv(g)
end
