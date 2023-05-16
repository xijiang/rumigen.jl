# SPDX-License-Identifier: MIT

include("effb9-sub.jl")

# Table of contents
# 1. effb9_xps: Repeat Quinton et al. 1991, with genomic selection

"""
    function xps_effb9(; describe=true, debug=true)
A simulation to expand Quinton et al. 1991, with genomic selection.
"""
function xps_effb9(; describe=true, debug=true)
    # Parameters
    rst, nsir, ndam, nbull, ncow = "rst", 50, 50, 100, 100
    nlc, nqtl, ngrt, h², σₐ = 50_000, 10_000, 10, 0.25, 1.0
    ns, nd, ngrt, nrpt = 10, 50, 10, 100
    dir = "$rst/quinton"

    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)

    # The working parts
    if debug
        @warn "Debugging"
    else
        describe && tprintln(Term.parse_md(read("docs/effb9.md", String)))
        for irpt in 1:nrpt
            @info "Repeat $irpt of $nrpt"
            foo = base_effb9(nsir, ndam, rst)                     # ==> base
            bar = fdr_effb9(rst, dir, foo, nsir, ndam, nlc, nqtl) # ==> founder
            f0_effb9(dir, bar, nsir, ndam, nbull, ncow)           # ==> f0 of UInt16 in 3 copies
            simQTL("$dir/$bar-f0.xy", "$dir/$bar-map.ser", d = Normal())
            ped_effb9(dir, bar, σₑ)                               # ==> f0-ped, shared
            lmp = deserialize("$dir/$bar-map.ser")
            ped = deserialize("$dir/$bar-f0-ped.ser")
            simpleSelection("$dir/$bar-pht.xy", ped, lmp, ns, nd, ngrt, σₑ)
            ped = deserialize("$dir/$bar-f0-ped.ser")
            simpleSelection("$dir/$bar-ped.xy", ped, lmp, ns, nd, ngrt, σₑ, ebv = true)
            ped = deserialize("$dir/$bar-f0-ped.ser")
            simpleSelection("$dir/$bar--gs.xy", ped, lmp, ns, nd, ngrt, σₑ, ebv = true, gs = true)
            sum_effb9(dir, bar)
            clean_effb9(rst, dir, bar)
        end
    end
end
