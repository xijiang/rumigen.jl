#=
Theo came with this idea that what if we use the real IBD matrix instead of G.
He predicts that the EBV accuracy will be lower.
I added two more scenarios below to the simulation. That is, IBD -> EBV and IBD, IBD -> EBV.
It may need just a few hours to finish.

echo EBV with IBD matrix | md5sum -> 75aee9e63f68406d9b6259891a5799bb
=#

function xps_75aee9( ;
    nlc = 50_000,
    nqtl = 10_000,
    nref = 10_000,
    ngrt = 20,
    ΔF = 0.011,
    nrpt = 1,
    rst = "rst",
    ppsz = 200,
    h² = 0.25,
    σₐ = 1.0,
    nsir = 20,
    ndam = 50,
    pres = 5,
    dist = Normal(),
    sim = "75aee9",
    quick_test=true,
    keep = false
    )
    
    fdr, dir = "$rst/founder", "$rst/$sim"
    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)
    serialize("$dir/par.bin", (nlc=nlc, nqtl=nqtl, nref=nref, ngrt=ngrt, ΔF=ΔF,
    nrpt=nrpt, rst=rst, ppsz=ppsz, h²=h², σₐ=σₐ, nsir=nsir, ndam=ndam, pres=pres,
    dist=dist, sim=sim, quick_test=quick_test))
    # use Dict(pairs(par)) to reconstruct the par dict, remember to include Distributions
    scheme_1 = ("ran", "spd", "sgs", "sgi")
    scheme_2 = ("oap", "oag", "ogg", "oig", "oii")
    schemes = ("spd", "sgs", "sgi", "oap", "oag", "ogg", "oig", "oii")

    # The working parts
    @info "Simulation begins"
    for irpt in 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        isdir(fdr) && rm(fdr, recursive=true, force=true)
        foo = if quick_test
            isfile("$rst/test-suite/founder-hap.xy") || error("No test-suite")
            fdr = "$rst/test-suite"
            "founder"
        else
            cattle_base(ppsz, fdr)
        end
        # random selectio for a few generations
        bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, nref, d=dist)
        lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
        ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg=-pres)
        simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ, 1)
        # selection without optimum contribution schemes
        for op in 2:4
            pop, sel = copy(ped), scheme_1[op]
            cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force=true)
            simpleSelection("$dir/$bar-$sel.xy", pop, lmp, nsir, ndam, ngrt, σₑ, op, mp=false)
            sumPed(rst, sim, bar, lmp, sel)
            pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        for op in 1:5
            pop, sel = copy(ped), scheme_2[op]
            cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force=true)
            optSelection("$dir/$bar-$sel.xy", pop, lmp, ngrt, σₑ, ΔF, op=op, k₀=0.027)
            sumPed(rst, sim, bar, lmp, sel)
            pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        keep || rm.(glob("$dir/$bar-*"), force=true)
    end
end
