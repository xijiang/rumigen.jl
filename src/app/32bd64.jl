# echo Examine variance surge for big Ne | md5sum
# 32bd6437416747a24e916e3c0014ff45
function xps_32bd64(;
    nlc = 50_000,
    nqtl = 10_000,
    nref = 10_000,
    ngrt = 5,
    ΔF = 0.011,
    nrpt = 1,
    rst = "rst",
    ppsz = 1000,
    h² = 0.25,
    σₐ = 1.0,
    nsir = 100,
    ndam = 200,
    pres = 5,
    dist = Normal(),
    sim = "32bd64",
    quick_test = true,
    keep = false,
)
    dir = "$rst/$sim"
    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) && rm(dir, recursive = true, force = true)
    mkpath(dir)
    serialize(
        "$dir/par.ser",
        (
            nlc = nlc,
            nqtl = nqtl,
            nref = nref,
            ngrt = ngrt,
            ΔF = ΔF,
            nrpt = nrpt,
            rst = rst,
            ppsz = ppsz,
            h² = h²,
            σₐ = σₐ,
            nsir = nsir,
            ndam = ndam,
            pres = pres,
            dist = dist,
            sim = sim,
            quick_test = quick_test,
        ),
    )
    # use Dict(pairs(par)) to reconstruct the par dict, remember to include Distributions
    scheme = ("ran", "spd", "sgs", "sis", "sms")

    # The working parts
    @info "Simulation begins"
    for irpt = 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        fdr, foo = prepFdr(rst, quick_test, ppsz)

        # random selectio for a few generations
        bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, nref, d = dist)
        lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
        ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg = -pres)
        simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ, 1)

        # selection without optimum contribution schemes
        for op = 1:3
            pop, sel = copy(ped), scheme[op]
            cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force = true)
            simpleSelection(
                "$dir/$bar-$sel.xy",
                pop,
                lmp,
                nsir,
                ndam,
                ngrt,
                σₑ,
                op,
                mp = false,
            )
            sumPed(rst, sim, bar, lmp, sel)
            #pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        keep || rm.(glob("$dir/$bar-*"), force = true)
    end
end
