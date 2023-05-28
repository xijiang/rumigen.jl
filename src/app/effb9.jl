# SPDX-License-Identifier: MIT

include("effb9-sub.jl")

# Table of contents
# 1. effb9_xps: Repeat Quinton et al. 1991, with genomic selection

"""
    function xps_effb9(; describe=true, debug=true)
A simulation to expand Quinton et al. 1991, with genomic selection. This
function is valid for `rumigen.jl v0.1.8-9`.
"""
function xps_effb9(; describe = true, debug = true)
    # Parameters
    rst, ppsz, nlc, nqtl, h², σₐ = "rst", 200, 50_000, 10_000, 0.25, 1.0
    nsir, ndam, pres, ngrt, nrpt, dist = 20, 50, 5, 15, 1, Normal()
    fdr, dir = "$rst/base", "$rst/quinton"

    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)

    # The working parts
    if debug
        @warn "Debugging"
        foo = "016Im"
        # improve speed below
        bar = fdr_effb9(fdr, dir, foo, ppsz, nlc, nqtl, d = dist)
        lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
        pre_effb9(dir, bar, nsir, ndam, pres, σₑ)  # ==> pres gnrtn's of rand mat'g
        ped = deserialize("$dir/$bar-uhp+ped.ser")
        cp("$dir/$bar-uhp.xy", "$dir/$bar-spd.xy", force = true)
        simpleSelection("$dir/$bar-spd.xy", ped, lmp, nsir, ndam, ngrt, σₑ, ebv = true)
    else
        describe && tprintln(Term.parse_md(read("docs/effb9.md", String)))
        #foo = base_effb9(ppsz, rst)                     # ==> base
        foo = "016Im"
        for irpt in 1:nrpt
            @info "Repeat $irpt of $nrpt"
            bar = fdr_effb9(fdr, dir, foo, ppsz, nlc, nqtl, d = dist)
            lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
            pre_effb9(dir, bar, nsir, ndam, pres, σₑ)  # ==> pres gnrtn's of rand mat'g

            # phenotype selection
            ped = deserialize("$dir/$bar-uhp+ped.ser")
            cp("$dir/$bar-uhp.xy", "$dir/$bar-spt.xy", force = true)
            simpleSelection("$dir/$bar-spt.xy", ped, lmp, nsir, ndam, ngrt, σₑ)

            # Pedigree selection
            ped = deserialize("$dir/$bar-uhp+ped.ser")
            cp("$dir/$bar-uhp.xy", "$dir/$bar-spd.xy", force = true)
            simpleSelection("$dir/$bar-spd.xy", ped, lmp, nsir, ndam, ngrt, σₑ, ebv = true)

            # Genomic selection
            ped = deserialize("$dir/$bar-uhp+ped.ser")
            cp("$dir/$bar-uhp.xy", "$dir/$bar-sgs.xy", force = true)
            simpleSelection("$dir/$bar-sgs.xy", ped, lmp, nsir, ndam, ngrt, σₑ, ebv = true, gs = true)
            #sum_effb9(dir, bar)
            rm.(glob("$dir/$bar-*"))
        end
        stats_effb9(dir, ngrt)
    end
end
