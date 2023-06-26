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
    n = nrow(df)
    u = [sort(unique(df[:, i]))[2:end] for i in 1:ncol(df)] # can also be 1:end-1, which has more codes than 2:end
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
Given a vector of `Bool`s indicating if a phenotype is not missing, return a
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
