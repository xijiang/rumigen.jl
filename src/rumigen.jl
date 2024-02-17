module rumigen

using BenchmarkTools
using DataFrames
using Distributions
using Glob
using LinearAlgebra
using Mmap
using Octavian
using Printf
using Random
using Serialization
using SparseArrays
using Statistics
using StatsBase
using Term

include("app/xps.jl")
include("cattle.jl")
include("evaluation.jl")
include("macs.jl")
include("mate.jl")
include("opt-contrb.jl")
include("pedigree.jl")
include("selection.jl")
include("simulation.jl")
include("term.jl")
include("util.jl")
include("vcf.jl")
include("xymat.jl")

Relation = Dict{Tuple{Int, Int}, Float64}
Base.show(io::IO, f::Float64) = @printf(io, "%.3f", f)
end # module rumigen
