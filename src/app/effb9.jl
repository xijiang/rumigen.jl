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
    rst, nsir, ndam, nbull, ncow = "rst", 500, 500, 500, 10_000
    nlc, nqtl, ngrt, h², foo = 50_000, 10_000, 10, 0.5, "9UxnV"
    dir = "$rst/quinton"
    isdir(dir) || mkpath(dir)

    # The working parts
    if debug
        @warn "Debugging"
        bar = "xOjA0d"

        for igrt in 1:ngrt
            # TBV, phenotypes, PEBV, GEBV of igrt - 1
            # Selection and produce igrt
            @show igrt
        end
    else
        describe && tprintln(Term.parse_md(read("docs/effb9.md", String)))
        foo = base_effb9(nsir, ndam, rst)                     # ==> base
        bar = fdr_effb9(rst, dir, foo, nsir, ndam, nlc, nqtl) # ==> founder
        f0_effb9(dir, bar, nsir, ndam, nbull, ncow)           # ==> f0 of UInt16
        simQTL("$dir/$bar-f0.xy", "$dir/$bar-map.ser", d = Normal())
        cmn_effb9(dir, bar)
    end
end
