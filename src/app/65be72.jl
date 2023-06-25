# echo "No male phenotypes, more intensive phenotypic selection" | md5sum
function xps_65be72(; debug = true, nrpt = 200, keep = false)
    rst, ppsz, nlc, nqtl, h², σₐ = "rst", 200, 50_000, 10_000, 0.25, 1.0
    nsir, ndam, pres, ngrt, dist = 20, 50, 5, 20, Normal()
    fdr, dir = "$rst/base", "$rst/65be72"

    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)

    # The working parts
    if debug
        @warn "Debugging"
    else
        @info "Simulation begins"
        isdir(fdr) && rm(fdr, recursive = true, force = true)
        foo = cattle_base(ppsz, fdr)                    # ==> base
        for irpt in 1:nrpt
            println()
            @info "Repeat $irpt of $nrpt"
            # random selectio for a few generations
            bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, d = dist)
            lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
            ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg = -pres)
            simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ,
                            random = true, mp = false)
            
            
            # phenotype selection not possible as no male phenotypes
            
            # Pedigree selection
            pop = copy(ped)
            cp("$dir/$bar-uhp.xy", "$dir/$bar-spd.xy", force = true)
            simpleSelection("$dir/$bar-spd.xy", pop, lmp, nsir, ndam, ngrt, σₑ,
                            ebv = true, mp = false)

            # Genomic selection
            pop = copy(ped)
            cp("$dir/$bar-uhp.xy", "$dir/$bar-sgs.xy", force = true)
            simpleSelection("$dir/$bar-sgs.xy", pop, lmp, nsir, ndam, ngrt, σₑ,
                            ebv = true, gs = true, mp = false)

            sum_65be72(dir, bar)
            keep || rm.(glob("$dir/$bar-*"))
        end
    end
end

function sum_65be72(dir, bar)
    for sel in ["sgs", "spd"]
        ped = deserialize("$dir/$bar-$sel+ped.ser")
        lmp = deserialize("$dir/$bar-map.ser")
        sp = combine(groupby(ped, :grt), 
                            :tbv => mean => :mbv,
                            :tbv => var => :vbv,
                            :F => mean => :mF,
                            [:ebv, :tbv] => cor => :cor)
        ideal, va, plst, nlst, bns = idealPop("$dir/$bar-$sel.xy", ped.grt, lmp)
        open("$dir/65be72.bin", "a") do io
            write(io, sp.mbv[2:end],
                      sp.vbv[2:end],  # as requested by SMS on 2023-05-21, by Theo
                      sp.mF[2:end],
                      sp.cor[2:end],
                      ideal[2:end],
                      va[2:end],
                      plst[2:end],  # no. of positive qtl lost
                      nlst[2:end])
        end
        open("$dir/65be72.bns", "a") do io
            write(io, bns)
        end
    end
end

"""
```bash
echo Three phenotype scenarios plus pedgree and gs | md5sum
```
-> c140add3300cbfeddb2c1da06f8eb31a
"""
function xps_c140ad(; debug = true, nrpt = 200, keep = false)
    rst, ppsz, nlc, nqtl, h², σₐ = "rst", 200, 50_000, 10_000, 0.25, 1.0
    nsir, ndam, pres, ngrt, dist = 20, 50, 5, 20, Normal()
    fdr, dir = "$rst/base", "$rst/c140ad"

    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)

    # The working parts
    if debug
        @warn "Debugging"
    else
        # base and founders
        isdir(fdr) && rm(fdr, recursive=true, force=true)
        foo = cattle_base(ppsz, fdr)                    # ==> base
        for irpt in 1:nrpt
            println()
            @info "Repeat $irpt of $nrpt"
            # random selectio for a few generations
            bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, d = dist)
            lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
            ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg = -pres)
            simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ, random = true)
            
            # phenotype selection not possible as no male phenotypes
            pop = copy(ped)
            cp("$dir/$bar-uhp.xy", "$dir/$bar-pta.xy", force = true)
            simpleSelection("$dir/$bar-pta.xy", pop, lmp, nsir, ndam, ngrt, σₑ)

            pop = copy(ped)
            cp("$dir/$bar-uhp.xy", "$dir/$bar-ptb.xy", force = true)
            simpleSelection("$dir/$bar-ptb.xy", pop, lmp, nsir÷2, ndam, ngrt, σₑ)

            pop = copy(ped)
            cp("$dir/$bar-uhp.xy", "$dir/$bar-ptc.xy", force = true)
            simpleSelection("$dir/$bar-ptc.xy", pop, lmp, nsir÷4, ndam, ngrt, σₑ)

            # Pedigree selection
            pop = copy(ped)
            cp("$dir/$bar-uhp.xy", "$dir/$bar-spd.xy", force = true)
            simpleSelection("$dir/$bar-spd.xy", pop, lmp, nsir, ndam, ngrt, σₑ,
                            ebv = true)

            # Genomic selection
            pop = copy(ped)
            cp("$dir/$bar-uhp.xy", "$dir/$bar-sgs.xy", force = true)
            simpleSelection("$dir/$bar-sgs.xy", pop, lmp, nsir, ndam, ngrt, σₑ,
                            ebv = true, gs = true)

            sum_c140ad(dir, bar)
            keep || rm.(glob("$dir/$bar-*"))
        end
    end
end

function sum_c140ad(dir, bar)
    for sel in ["pta", "ptb", "ptc", "sgs", "spd"]
        ped = deserialize("$dir/$bar-$sel+ped.ser")
        lmp = deserialize("$dir/$bar-map.ser")
        sp = combine(groupby(ped, :grt), 
                            :tbv => mean => :mbv,
                            :tbv => var => :vbv,
                            :F => mean => :mF,
                            [:ebv, :tbv] => cor => :cor)
        ideal, va, plst, nlst, bns = idealPop("$dir/$bar-$sel.xy", ped.grt, lmp)
        open("$dir/c140ad.bin", "a") do io
            write(io, sp.mbv[2:end],
                      sp.vbv[2:end],  # as requested by SMS on 2023-05-21, by Theo
                      sp.mF[2:end],
                      sp.cor[2:end],
                      ideal[2:end],
                      va[2:end],
                      plst[2:end],  # no. of positive qtl lost
                      nlst[2:end])
        end
        open("$dir/c140ad.bns", "a") do io
            write(io, bns)
        end
    end
end

function overnight(; nrepeats = 100)
    xps_65be72(debug = false, nrpt = nrepeats)
    xps_c140ad(debug = false, nrpt = nrepeats)
end
