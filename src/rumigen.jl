module rumigen

using BenchmarkTools
using DataFrames
using Distributions
using Glob
using LinearAlgebra
using Mmap
using Octavian
using Random
using Serialization
using SparseArrays
using Statistics
using Term

include("app/xps.jl")
include("cattle.jl")
include("evaluation.jl")
include("macs.jl")
include("mate.jl")
include("pedigree.jl")
include("selection.jl")
include("simulation.jl")
include("term.jl")
include("util.jl")
include("xymat.jl")

end # module rumigen
