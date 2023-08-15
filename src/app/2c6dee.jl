# echo Selection by optimum contribution | md5sum
function xps_2c6dee(; nrpt=100, ΔF=0.012, keep=false)
    rst, ppsz, nlc, nqtl, nref, h², σₐ = "rst", 200, 50_000, 10_000, 10_000, 0.25, 1.0
    nsir, ndam, pres, ngrt, dist = 20, 50, 5, 20, Normal()
    fdr, dir = "$rst/base", "$rst/2c6dee"

    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)

    # The working parts
    @info "Simulation begins"
    #isdir(fdr) && rm(fdr, recursive=true, force=true)
    #foo = cattle_base(ppsz, fdr)                    # ==> base
    foo = "jArVr"
    for irpt in 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
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

        # Optimum contribution with pedigree
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-oap.xy", force=true)
        optSelection("$dir/$bar-oap.xy", pop, lmp, ngrt, σₑ, ΔF, op=1, k₀=0.027)

        # Optimum contribution with genomic selection, A constrained
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-oag.xy", force=true)
        optSelection("$dir/$bar-oag.xy", pop, lmp, ngrt, σₑ, ΔF, op=2, k₀=0.027)

        # Optimum contribution with genomic selection, G constrained
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-ogg.xy", force=true)
        optSelection("$dir/$bar-ogg.xy", pop, lmp, ngrt, σₑ, ΔF, op=3, k₀=0.027)

        sum_2c6dee(dir, bar, lmp)
        pos_qtl_frq(dir, bar, ["sgs", "spd", "oag", "ogg", "oap"], ppsz)
        keep || rm.(glob("$dir/$bar-*"))
    end
end

function sum_2c6dee(dir, bar, lmp)
    for sel in ["sgs", "spd", "oag", "ogg", "oap"]
        ped = deserialize("$dir/$bar-$sel+ped.ser")
        sp = DataFrame(mbv=Float64[], vbv=Float64[], mF=Float64[], mFr=Float64[],
            bcr=Float64[], scr=Float64[], dcr=Float64[], np=Float64[], nm=Float64[])
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
            push!(sp, (mbv, vbv, mF, mFr, bcr, scr, dcr, np, nm))
        end
        ideal, va, plst, nlst, pmls, nmls = idealPop("$dir/$bar-$sel.xy", ped.grt, lmp)
        open("$dir/2c6dee.bin", "a") do io
            write(io,
                sp.mbv[2:end],  # 1
                sp.vbv[2:end],  # as requested by SMS on 2023-05-21, by Theo
                sp.mF[2:end],   # 3
                sp.mFr[2:end],  # 4
                sp.bcr[2:end],  # 5
                sp.scr[2:end],  # 6
                sp.dcr[2:end],  # 7
                sp.np[2:end],   # 8
                sp.nm[2:end],   # 9
                ideal[2:end],   # 10
                va[2:end],      # 11
                plst[2:end],  # no. of positive qtl lost
                nlst[2:end],    # 13
                pmls[2:end],  # no. of positive qtl lost of maf 0.2
                nmls[2:end])  # no. of negative qtl lost of maf 0.2
        end
    end
end

function spl_2c6dee(dir)
    test = "$dir/rst/test-suite"
    ped = deserialize("$test/bar-uhp+ped.ser")
    giv = A⁻¹(ped, "$test/bar-mid.ser")
    animalModel(ped, giv, .25)
    lst = 1001:1200
    A = Amat(ped, m = size(ped, 1))
    return ped.ebv[lst], ped.sex[lst], A[lst, lst]
end

"""
    function restraint(t, df)
This function return the restraint for generation `t` with `ΔF = df` and `k₀ =
0`.
"""
function restraint(t, df; k = 0.0)
    #for i in 1:t
    #    k += df * (1 - k)
    #end
    # above is equavalent to
    2(1 - (1 - k) * (1 - df)^(t-1))
end
