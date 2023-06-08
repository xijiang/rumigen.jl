# echo "No male phenotypes, more intensive phenotypic selection" | md5sum
function xps_65be72(; debug = true)
    rst, ppsz, nlc, nqtl, h², σₐ = "rst", 200, 50_000, 10_000, 0.25, 1.0
    nsir, ndam, pres, ngrt, nrpt, dist = 20, 50, 5, 20, 200, Normal()
    fdr, dir = "$rst/base", "$rst/65be72"

    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)

    # The working parts
    if debug
        @warn "Debugging"
    else
        @info "Simulation begins"
        #foo = cattle_base(ppsz, rst)                    # ==> base
        foo = "fRpPu"
        for irpt in 1:nrpt
            @info "\nRepeat $irpt of $nrpt"
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
            rm.(glob("$dir/$bar-*"))
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
        open("$dir/65be72.bin", "a") do io
            ideal, va, _, _ = idealPop("$dir/$bar-$sel.xy", ped.grt, lmp)
            write(io, sp.mbv[2:end],
                      sp.vbv[2:end],  # as requested by SMS on 2023-05-21, by Theo
                      sp.mF[2:end],
                      sp.cor[2:end],
                      ideal[2:end],
                      va[2:end])
            # write(io, plst[2:end])
            # write(io, nlst[2:end])
        end
    end
end

"""
```bash
echo Three phenotype scenarios plus pedgree and gs | md5sum
```
-> c140add3300cbfeddb2c1da06f8eb31a
"""
function xps_c140ad(; debug = true)
    rst, ppsz, nlc, nqtl, h², σₐ = "rst", 200, 50_000, 10_000, 0.25, 1.0
    nsir, ndam, pres, ngrt, nrpt, dist = 20, 50, 5, 20, 200, Normal()
    fdr, dir = "$rst/base", "$rst/65be72"

    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) || mkpath(dir)

    # The working parts
    if debug
        @warn "Debugging"
    else
        # base and founders
        #foo = cattle_base(ppsz, rst)                    # ==> base
        foo = "fRpPu"
        for irpt in 1:nrpt
            @info "\nRepeat $irpt of $nrpt"
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

            sum_65be72(dir, bar)
            rm.(glob("$dir/$bar-*"))
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
        open("$dir/65be72.bin", "a") do io
            ideal, va, _, _ = idealPop("$dir/$bar-$sel.xy", ped.grt, lmp)
            write(io, sp.mbv[2:end],
                      sp.vbv[2:end],  # as requested by SMS on 2023-05-21, by Theo
                      sp.mF[2:end],
                      sp.cor[2:end],
                      ideal[2:end],
                      va[2:end])
        end
    end
end
