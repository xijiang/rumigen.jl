# echo Selection by optimum contribution | md5sum
function xps_2c6dee(;
    nrpt = 100,
    nlc = 50_000,
    nqtl = 10_000,
    nref = 10_000,
    ngrt = 20,
    ΔF = 0.011,
    keep = false,
    qt = false
    ) # qt: quick test
    rst, ppsz, h², σₐ = "rst", 200, 0.25, 1.0
    nsir, ndam, pres, dist = 20, 50, 5, Normal()
    fdr, dir = "$rst/base", "$rst/2c6dee"

    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)

    # The working parts
    @info "Simulation begins"
    for irpt in 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        isdir(fdr) && rm(fdr, recursive=true, force=true)
        foo = if qt
            fdr = "$rst/test-suite"
            "Xm6k4"
        else
            cattle_base(ppsz, fdr) # ==> base
        end

        # random selectio for a few generations
        bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, nref, d=dist)
        lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
        ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg=-pres)
        simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ,
            random=true, mp=false)

        # phenotype selection not possible as no male phenotypes
        # Pedigree selection
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-spd.xy", force=true)
        simpleSelection("$dir/$bar-spd.xy", pop, lmp, nsir, ndam, ngrt, σₑ,
            ebv=true, mp=false)

        # Genomic selection
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-sgs.xy", force=true)
        simpleSelection("$dir/$bar-sgs.xy", pop, lmp, nsir, ndam, ngrt, σₑ,
            ebv=true, gs=true, mp=false)

        # Optimum contribution with A on pedigree
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-oap.xy", force=true)
        optSelection("$dir/$bar-oap.xy", pop, lmp, ngrt, σₑ, ΔF, op=1, k₀=0.027)

        # Optimum contribution with A on genomic selection
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-oag.xy", force=true)
        optSelection("$dir/$bar-oag.xy", pop, lmp, ngrt, σₑ, ΔF, op=2, k₀=0.027)

        # Optimum contribution with G on genomic selection
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-ogg.xy", force=true)
        optSelection("$dir/$bar-ogg.xy", pop, lmp, ngrt, σₑ, ΔF, op=3, k₀=0.027)

        # Optimum contribution with E on genomic selection, E is A with true IBD
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-oeg.xy", force=true)
        optSelection("$dir/$bar-oeg.xy", pop, lmp, ngrt, σₑ, ΔF, op=4, k₀=0.027)

        @info "Summarizing repeat $irpt of $nrpt"
        sum_2c6dee(dir, bar, lmp)
        pos_qtl_frq(dir, bar, ["sgs", "spd", "oeg", "oag", "ogg", "oap"], ppsz)
        keep || rm.(glob("$dir/$bar-*"))
    end
end

function sum_2c6dee(dir, bar, lmp)
    for sel in ["sgs", "spd", "oeg", "oag", "ogg", "oap"]
        ped = deserialize("$dir/$bar-$sel+ped.ser")
        sp = DataFrame(mbv=Float64[], vbv=Float64[], mF=Float64[], mFr=Float64[],
            bcr=Float64[], scr=Float64[], dcr=Float64[], np=Float64[], nm=Float64[],
            ncp=Float64[], ncm=Float64[])
        for grt in groupby(ped, :grt)
            mbv = mean(grt.tbv)
            vbv = var(grt.tbv)
            mF = mean(grt.F)
            mFr = mean(grt.Fr)
            bcr = cor(grt.ebv, grt.tbv)
            sirs = grt.sex .== 1
            dams = grt.sex .== 0
            scr = cor(grt.ebv[sirs], grt.tbv[sirs])
            dcr = cor(grt.ebv[dams], grt.tbv[dams])
            np = length(unique(grt.pa))
            nm = length(unique(grt.ma))
            ncp = sum(grt.c[sirs] .> 0)
            ncm = sum(grt.c[dams] .> 0)
            push!(sp, (mbv, vbv, mF, mFr, bcr, scr, dcr, np, nm, ncp, ncm))
        end
        ideal, va, plst, nlst, pmls, nmls = idealPop("$dir/$bar-$sel.xy", ped.grt, lmp)
        open("$dir/2c6dee.bin", "a") do io
            write(io,
                sp.mbv,  # 1
                sp.vbv,  # as requested by SMS on 2023-05-21, by Theo
                sp.mF,   # 3
                sp.mFr,  # 4
                sp.bcr,  # 5
                sp.scr,  # 6
                sp.dcr,  # 7
                sp.np,   # 8
                sp.nm,   # 9
                sp.ncp,  # 10
                sp.ncm,  # 11
                ideal,   # 12
                va,      # 13
                plst,    # 14. n. of positive qtl lost
                nlst,    # 15
                pmls,    # 16. no. of positive qtl lost of maf 0.2
                nmls)    # 17. no. of negative qtl lost of maf 0.2
        end
    end
end

function spl_2c6dee(dir)
    test = "$dir/rst/test-suite"
    ped = deserialize("$test/bar-uhp+ped.ser")
    giv = A⁻¹(ped, "$test/bar-mid.ser")
    animalModel(ped, giv, 0.25)
    lst = 1001:1200
    A = Amat(ped, m=size(ped, 1))
    return ped.ebv[lst], ped.sex[lst], A[lst, lst]
end

"""
    function restraint(t, df)
This function return the restraint for generation `t` with `ΔF = df` and `k₀ =
0`.
"""
function restraint(t, df; k=0.0)
    #for i in 1:t
    #    k += df * (1 - k)
    #end
    # above is equavalent to
    2(1 - (1 - k) * (1 - df)^(t - 1))
end
