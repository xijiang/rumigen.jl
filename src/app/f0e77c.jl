# echo OCS with equal number of parents | md5sum
#-> f0e77cf92515d3e9318085872b868a93
function xps_f0e77c(;
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
    sim = "f0e77c",
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
        ),
    )
    # use Dict(pairs(par)) to reconstruct the par dict, remember to include Distributions
    #scheme_1 = ("ran", "spd", "sgs", "sis", "sms")
    scheme_ocs = ("oap", "oag", "ogg", "oig", "oii", "otg", "dos")

    # The working parts
    @info "Simulation begins"
    for irpt = 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        bar = sampleFdr(
            "$rst/founder/cattle-hap.xy",
            "$rst/founder/cattle-map.ser",
            2ppsz,
            nlc = nlc,
            nqtl = nqtl,
            nref = nref,
            dir = dir,
        )
        simQTL("$dir/$bar-fdr.xy", "$dir/$bar-map.ser", d = dist)
        uniqSNP("$dir/$bar-fdr.xy")
        lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
        ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg = -pres)
        simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ, 1)

        @info "Calculating IBD matrix of initial population"
        G = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:size(ped, 1))
        for sel in ("oii", "dos")
            write("$dir/$bar-$sel.bin", G)
        end
        G = nothing
        for op in [2, 5, 7]
            pop, sel = copy(ped), scheme_ocs[op]
            cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force = true)
            optSelection("$dir/$bar-$sel.xy", pop, lmp, ngrt, σₑ, ΔF, op = op, k₀ = 0.027)
            sumPed(rst, sim, bar, lmp, sel)
            pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        keep || rm.(glob("$dir/$bar-*"), force = true)
    end
end
