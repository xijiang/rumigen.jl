#=
Theo came with this idea that what if we use the real IBD matrix instead of G.
He predicts that the EBV accuracy will be lower.
I added two more scenarios below to the simulation. That is, IBD -> EBV and IBD, IBD -> EBV.
It may need just a few hours to finish.

echo EBV with IBD matrix | md5sum -> 75aee9e63f68406d9b6259891a5799bb
=#

function mydebug(;
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
    sim = "debug-bng",
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
    # scheme_1 = ("ran", "spd", "sgs", "sis", "sms")
    # scheme_ocs = ("oap", "oag", "ogg", "oig", "oii", "otg")
    scheme_ocs = ("oap")

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
        op, pop, sel = 1, copy(ped), "oap"
        cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force = true)
        return optSelection(
            "$dir/$bar-$sel.xy",
            pop,
            lmp,
            ngrt,
            σₑ,
            ΔF,
            op = op,
            k₀ = 0.027,
        )
        sumPed(rst, sim, bar, lmp, sel)
        pos_qtl_frq(rst, sim, bar, sel, ppsz)
        keep || rm.(glob("$dir/$bar-*"), force = true)
    end
end
