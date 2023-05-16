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
