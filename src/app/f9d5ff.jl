# echo Simulation with CAU chicken data. Genotype only. | md5sum => f9d5ff6946069db569810bdb1ff881f9
# Similar to 75aee9, but with real (chicken) data.

function xps_f9d5ff(;
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
    nsir = 25,
    ndam = 50,
    pres = 5,
    dist = Normal(),
    sim = "f9d5ff",
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
    scheme_1 = ("ran", "spd", "sgs", "sis", "sms")
    scheme_ocs = ("oap", "oag", "ogg", "oig", "oii")

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

        @info "Calculating IBD matrix of initial population"
        G = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:size(ped, 1))
        for sel in ("sis", "oii")
            write("$dir/$bar-$sel.bin", G)
        end
        G = nothing
        # selection without optimum contribution schemes
        for op = 2:4
            pop, sel = copy(ped), scheme_1[op]
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
            pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        for op = 1:5
            pop, sel = copy(ped), scheme_ocs[op]
            cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force = true)
            optSelection("$dir/$bar-$sel.xy", pop, lmp, ngrt, σₑ, ΔF, op = op, k₀ = 0.027)
            sumPed(rst, sim, bar, lmp, sel)
            pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        keep || rm.(glob("$dir/$bar-*"), force = true)
    end
end
