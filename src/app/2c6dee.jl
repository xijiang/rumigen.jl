# echo Selection by optimum contribution | md5sum
function xps_2c6dee(; nrpt=200, ΔF = 0.012, keep = false)
    rst, ppsz, nlc, nqtl, h², σₐ = "rst", 200, 50_000, 10_000, 0.25, 1.0
    nsir, ndam, pres, ngrt, dist = 20, 50, 5, 20, Normal()
    fdr, dir = "$rst/base", "$rst/2c6dee"

    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)

    # The working parts
    @info "Simulation begins"
    isdir(fdr) && rm(fdr, recursive=true, force=true)
    foo = cattle_base(ppsz, fdr)                    # ==> base
    for irpt in 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        # random selectio for a few generations
        bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, d=dist)
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
        cp("$dir/$bar-uhp.xy", "$dir/$bar-opd.xy", force=true)
        optSelection("$dir/$bar-opd.xy", pop, lmp, ngrt, σₑ, ΔF, gs=false, k₀=0.027)

        # Optimum contribution with genomic selection
        pop = copy(ped)
        cp("$dir/$bar-uhp.xy", "$dir/$bar-ogs.xy", force=true)
        optSelection("$dir/$bar-ogs.xy", pop, lmp, ngrt, σₑ, ΔF, gs=true, k₀=0.027)

        sum_2c6dee(dir, bar, lmp)
        pos_qtl_frq(dir, bar, ["sgs", "spd", "ogs", "opd"], ppsz)
        keep || rm.(glob("$dir/$bar-*"))
    end
end

function sum_2c6dee(dir, bar, lmp)
    for sel in ["sgs", "spd", "ogs", "opd"]
        ped = deserialize("$dir/$bar-$sel+ped.ser")
        sp = DataFrame(mbv=Float64[], vbv=Float64[], mF=Float64[], bcr=Float64[],
            scr=Float64[], dcr=Float64[])
        for grt in groupby(ped, :grt)
            mbv = mean(grt.ebv)
            vbv = var(grt.ebv)
            mF = mean(grt.F)
            bcr = cor(grt.ebv, grt.tbv)
            sirs = grt.sex .== 1
            dams = grt.sex .== 0
            scr = cor(grt.ebv[sirs], grt.tbv[sirs])
            dcr = cor(grt.ebv[dams], grt.tbv[dams])
            push!(sp, (mbv, vbv, mF, bcr, scr, dcr))
        end
        ideal, va, plst, nlst, pmls, nmls = idealPop("$dir/$bar-$sel.xy", ped.grt, lmp)
        open("$dir/2c6dee.bin", "a") do io
            write(io, sp.mbv[2:end],
                sp.vbv[2:end],  # as requested by SMS on 2023-05-21, by Theo
                sp.mF[2:end],
                sp.bcr[2:end],
                sp.scr[2:end],
                sp.dcr[2:end],
                ideal[2:end],
                va[2:end],
                plst[2:end],  # no. of positive qtl lost
                nlst[2:end],
                pmls[2:end],  # no. of positive qtl lost of maf 0.2
                nmls[2:end])  # no. of negative qtl lost of maf 0.2
        end
    end
end

function cbvdist(ebv, sex, A, npa, nma; nsmpl = 10000, rnd = true)
    mbv, cn = zeros(nsmpl), zeros(nsmpl)
    for i in 1:nsmpl
        c = rnd ? randomc(sex, npa, nma) : equalc(sex, npa, nma)
        mbv[i] = c'ebv
        cn[i] = c'A * c / 2
    end
    return mbv, cn
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

function debug_2c6dee(dir)
    f0, df = 0.0266642175, 0.012
    return fungencont(dat, A[1001:1200, 1001:1200], 0.13)
    cor(ped.pht, ped.ebv)
    dat = [lpd.ebv lpd.sex .== 1 lpd.sex .== 0]
    K = Kt(6, df)
    lmp = deserialize("$test/bar-map.ser")
    gt = xymap("$dir/$test/bar-uhp.xy")
    lms = sumMap(lmp)
    nsir, ndam, nfam, nsib = 20, 50, 50, 4
    agt = xymap("$dir/$test/bar-uhp.xy")
    dat = [ped.ebv ped.sex .== 1 ped.sex .== 0]
    K 
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
