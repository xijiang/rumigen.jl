# SPDX-License-Identifier: MIT

include("effb9-sub.jl")

# Table of contents
# 1. effb9_xps: Repeat Quinton et al. 1991, with genomic selection

"""
    effb9_xps(; bar="", describe = true)
Simulation 500 males and 500 females.
"""
function effb9_main(; foo="", describe=true, debug=true)
    # Parameters
    rst, nsir, ndam = "rst", 500, 500
    dir = "$rst/quinton"
    isdir(dir) || mkpath(dir)

    # The working parts
    if debug
        @warn "Debugging"
        bar = "JTYrLH"
        # ToDo: drop genotypes by 29/05/2023
        effb9_f0(dir, bar, nsir, ndam)
    else
        describe && tprintln(Term.parse_md(read("docs/effb9.md", String)))
        foo == "" && (foo = effb9_base(nsir, ndam, rst)) # ==> base
        bar = cattle_fdr(rst, dir, foo, nsir, ndam)      # ==> founder
        effb9_f0(dir, bar, nsir, ndam)
    end
end
