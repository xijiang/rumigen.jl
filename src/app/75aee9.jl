#=
Theo came with this idea that what if we use the real IBD matrix instead of G.
He predicts that the EBV accuracy will be lower.
I added two more scenarios below to the simulation. That is, IBD -> EBV and IBD, IBD -> EBV.
It may need just a few hours to finish.

echo EBV with IBD matrix | md5sum -> 75aee9e63f68406d9b6259891a5799bb
=#

function xps_75aee9(;
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
    sim = "75aee9",
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
    scheme_ocs = ("oap", "oag", "ogg", "oig", "oii", "otg")

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
        for op = 1:6
            pop, sel = copy(ped), scheme_ocs[op]
            cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force = true)
            optSelection("$dir/$bar-$sel.xy", pop, lmp, ngrt, σₑ, ΔF, op = op, k₀ = 0.027)
            sumPed(rst, sim, bar, lmp, sel)
            pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        keep || rm.(glob("$dir/$bar-*"), force = true)
    end
end

## echo Check random mating and delta F | md5sum -> 49decb5a5aa5bb7bdf7761b54e3fd2aa
## Modified on 2024-03-02. This function uses only the preselection part.
function xps_49decb(;
    nlc = 50_000,
    nqtl = 10_000,
    nref = 10_000,
    nrpt = 1,
    rst = "rst",
    ppsz = 200,
    h² = 0.25,
    σₐ = 1.0,
    nsir = 25,
    ndam = 50,
    pres = 25,
    dist = Normal(),
    sim = "49decb",
    quick_test = true,
    keep = false,
)
    dir = "$rst/$sim"
    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) && rm(dir, recursive = true, force = true)
    mkpath(dir)
    # make sure ppsz is divisible by bigger of nsir and ndam
    ppsz = ppsz ÷ max(nsir, ndam) * max(nsir, ndam)
    serialize(
        "$dir/par.ser",
        (
            nlc = nlc,
            nqtl = nqtl,
            nref = nref,
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
        ped.grt .+= 20
        serialize("$dir/$bar-uhp+ped.ser", ped)
        sumPed(rst, sim, bar, lmp, "uhp")
        pos_qtl_frq(rst, sim, bar, "uhp", ppsz)
        keep || rm.(glob("$dir/$bar-*"), force = true)
    end
end

"""
echo Check relationship by real IBD info | md5sum -> 64d3f8d6cc4ccf8e58a45cc501f2d2d5
"""
function xps_64d3f8()
    dir, fdr, foo = "rst/64d3f8", "rst/test-suite", "founder"
    isdir(dir) || mkpath("rst/64d3f8")
    bar = cattle_founder(fdr, dir, foo, 200, 50_000, 10_000, 10_000, d = Normal())
    lmp = deserialize("$dir/$bar-map.ser")
    ped = initPedigree("$dir/$bar-uhp.xy", lmp, 1.0)
    simpleSelection("$dir/$bar-uhp.xy", ped, lmp, 20, 50, 9, 1.0, 1)
    A = Amat(ped)
    M = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:size(ped, 1))
    serialize("$dir/$bar-A.imt", A)
    serialize("$dir/$bar-M.imt", M)
    M2 = gametemat("$dir/$bar-uhp.xy", lmp.chr .== 1, 1:size(ped, 1))
    serialize("$dir/$bar-M2.imt", M2)
end

"""
echo test function gametemat | md5sum
ddbbc13870afca0c485058c6fe761a1c
"""
function xys_ddbbc1()
    fdr, dir = "rst/test-suite", "rst/ddbbc1"
    isdir(dir) && rm(dir, recursive = true, force = true)
    mkpath(dir)
    bar = cattle_founder(fdr, dir, "founder", 200, 50_000, 10_000, 10_000, d = Normal())
    lmp = deserialize("$dir/$bar-map.ser")
    ped = initPedigree("$dir/$bar-uhp.xy", lmp, 1.0)
    simpleSelection("$dir/$bar-uhp.xy", ped, lmp, 20, 50, 9, 1.0, 1)
    n = size(ped, 1)
    M = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:n, 1:n)
    A = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:n)
    return M, A
end
