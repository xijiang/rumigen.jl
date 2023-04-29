module rumigen

using BenchmarkTools
using DataFrames
using Distributions
using LinearAlgebra
using Mmap
using Random
using Serialization
using Statistics
using Term

include("cattle.jl")
include("macs.jl")
include("xymat.jl")
include("simulation.jl")
include("term.jl")
include("mate.jl")
include("app/xps.jl")

end # module rumigen
