# echo AGBLUP vs random mating | md5sum → b0b1a06726eb3f2de8fef58bf4244a2c
# this experiment is discontinued
# The original idea was to compare the performance of AGBLUP and random mating
# But the procedure is not that clear cut.

function xps_b0b1a0(;
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
    sim = "b0b1a0",
    quick_test=true,
    keep = false
    )
    dir = "$rst/$sim"
    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) && rm(dir, recursive=true, force=true)
    mkpath(dir)
    serialize("$dir/par.ser", (nlc=nlc, nqtl=nqtl, nref=nref, ngrt=ngrt, ΔF=ΔF,
        nrpt=nrpt, rst=rst, ppsz=ppsz, h²=h², σₐ=σₐ, nsir=nsir, ndam=ndam,
        pres=pres, dist=dist, sim=sim, quick_test=quick_test))
    schemes = ("oap", "oag", "ogg", "oig", "oii", "rag")

    # The working parts
    @info "Simulation begins"
    for irpt in 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        fdr, foo = prepFdr(rst, quick_test, ppsz)

        # random selectio for a few generations
        bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, nref, d=dist)
        lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
        ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg=-pres)
        simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ, 1)

        @info "Calculating IBD matrix of initial population"
        G = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:size(ped, 1))
        write("$dir/$bar-oii.bin", G)
        G = nothing

        for op in (2, 5, 6)
            pop, sel = copy(ped), schemes[op]
            cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force=true)
            optSelection("$dir/$bar-$sel.xy", pop, lmp, ngrt, σₑ, ΔF, op=op, k₀=0.027)
            sumPed(rst, sim, bar, lmp, sel)
            pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        keep || rm("$dir/$bar-*", force=true)
    end
end
