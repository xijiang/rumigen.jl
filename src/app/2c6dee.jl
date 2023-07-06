# echo Selection by optimum contribution | md5sum
function xps_2c6bee(; debug=true, nrpt=200, keep=false)
    rst, ppsz, nlc, nqtl, h², σₐ = "rst", 200, 50_000, 10_000, 0.25, 1.0
    nsir, ndam, pres, ngrt, dist = 20, 50, 5, 20, Normal()
    fdr, dir = "$rst/base", "$rst/2c6bee"

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

        sum_65be72(dir, bar)
        count_ones(dir, bar, ["sgs", "spd"], ppsz)
        keep || rm.(glob("$dir/$bar-*"))
    end
end

function xps_2c6bee_debug()
    test = "rst/test-suite"
    ped = deserialize("$test/bar-uhp+ped.ser")
    lmp = deserialize("$test/bar-map.ser")
    gt = xymap("$test/bar-uhp.xy")
    lms = sumMap(lmp)
    nsir, ndam, nfam, nsib = 20, 50, 50, 4
    agt = xymap("$test/bar-uhp.xy")
    giv = A⁻¹(ped, "$test/bar-mid.ser")
    animalModel(ped, giv, .25)
    cor(ped.pht, ped.ebv)
    dat = [ped.ebv ped.sex .== 1 ped.sex .== 0]
    A = Amat(ped, m = nrow(ped))
    K 
end

function kt(t, ne)
    ΔF = 1 / (2 * ne)
    K = 0
    for i in 1:t
        K += ΔF * (1 - K)
    end
    return K
end
