function xps_7252bd(;
    nlc = 300_000, # or 500k
    nqtl = 10_000,
    nref = 10_000,
    ngrt = 20,
    ΔF = 0.011, # 1/(200)
    nrpt = 1,
    rst = "rst",
    ppsz = 400,
    h² = 0.25,
    σₐ = 1.0,
    nsir = 50,
    ndam = 100,
    pres = 5,
    dist = Normal(),
    sim = "7252bd",
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

    # Truncation selection
    scheme_tcs = ("ran", "spd", "sgs", "sis", "sms", "seg")
    # Optimum contribution selection
    scheme_ocs = ("oap", "oag", "ogg", "oig", "oii", "otg", "dos",
                  "osg", "oss", "ois")
    @info "Simulation for Oda's paper-II begins"
    for irpt in 1:nrpt
        println()
        bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, nref, d=dist)
        lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
        ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg=-pres)
        simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ, 1)

        @info "Calculating IBD matrix of initial population"
        G = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:size(ped, 1))
        for sel in ("sis", "oii")
            write("$dir/$bar-$sel.bin", G)
        end
        G = nothing
        @info "Repeat $irpt of $nrpt"

        for op in (2, 3, 4, 6)
        end
        for op in (5, 8, 9, 10)
        end
    end
end
