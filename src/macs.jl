"""
    function make_macs(; tdir=pwd())
Clone `MaCS` into a temp dir under `tdir` and compile it into `tdir`.
## Note
- designed only for Linux systems.
- `g++`, and `boost-devel` must be installed.
- this function return the absolute path of newly compiled `macs` in the end.
"""
function make_macs(; tdir=pwd())
    macs = joinpath(abspath(tdir), "macs")
    @debug "Making macs"
    isfile(macs) && return macs
    isdir(tdir) || mkpath(tdir)
    wdir = mktempdir(tdir)
    run(`git clone https://github.com/gchen98/macs $wdir`)
    src = joinpath.(wdir, ["simulator.cpp",
        "algorithm.cpp",
        "datastructures.cpp"])
    target = joinpath(tdir, "macs")
    run(`g++ -o $target -O3 -Wall $src`)
    return macs
end

"""
    function read_macs(file)
---
Read genotypes and physical positions from simulation results of `macs`.
Returns
- genotypes of Array{Int8, 2}
- physical positions and
- allele frequencies
"""
function read_macs(file)
    gt, ps = Int8[], Int[]
    open(file, "r") do io
        nbp = parse(Float64, split(readline(io))[4])
        for line in eachline(io)
            (line[1:4] ≠ "SITE") && continue
            _, _, p, _, as = split(line)
            push!(ps, Int(ceil(parse(Float64, p) * nbp)))
            a = parse.(Int8, collect(as))
            append!(gt, a)
        end
    end
    nlc = length(ps)
    gt = reshape(gt, :, nlc)'
    fq = mean(gt, dims=2)
    return gt, vec(ps), vec(fq)
end
