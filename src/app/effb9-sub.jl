"""
    function base_effb9(nsir, ndam, rst)
Create `2nid` haplotypes with MaCS, and convert them to the `xy` format.
"""
function base_effb9(nid, rst)
    tprintln("Generating a base population with MaCS")
    macs = make_macs(tdir = rst)
    tmp = cattle_genome(macs, nid, dir = "$rst/base")
    macs2xy(tmp)
end

"""
    function fdr_effb9(fdr, dir, foo, nsir, ndam, ppsz, nlc, nqtl)
Sample founders from the base population and mate them very randomly into F1.
"""
function fdr_effb9(fdr, dir, foo, ppsz, nlc, nqtl; d = Normal())
    @info "Sample founders from the base population created with MaCS"

    # sample haplotypes for the founder population
    nhp = 2ppsz
    bar = sampleFdr("$fdr/$foo-hap.xy", "$fdr/$foo-map.ser",
                    nhp, nlc = nlc, nqtl = nqtl, dir = dir)
    simQTL("$dir/$bar-fdr.xy", "$dir/$bar-map.ser", d = d)
    uniqSNP("$dir/$bar-fdr.xy")
    return bar
end

function pre_effb9(dir, bar, nsir, ndam, ngrt, σₑ)
    @info "Random mate the first $ngrt generations to form a common base population"
    lmp = deserialize("$dir/$bar-map.ser")

    # expand founders to a constant sized population
    ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ)
    ped.grt .= -ngrt
    simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, ngrt, σₑ, random = true)
end

"""
    function sum_effb9(dir, bar)
Calculate mean `tbv`, `F`, of each generation, and number
of fixed loci on chip and QTL in the end.
"""
function sum_effb9(dir, bar)
    for sel in ["sgs", "spt", "spd"]
        ped = deserialize("$dir/$bar-$sel+ped.ser")
        lmp = deserialize("$dir/$bar-map.ser")
        sp = combine(groupby(ped, :grt), 
                            :tbv => mean => :mbv,
                            :tbv => var => :vbv,
                            :F => mean => :mF,
                            [:ebv, :tbv] => cor => :cor)
        open("$dir/effb9.bin", "a") do io
            ideal, _, _ = idealPop("$dir/$bar-$sel.xy", ped.grt, lmp)
            write(io, sp.mbv[2:end],
                      sp.vbv[2:end],  # as requested by SMS on 2023-05-21, by Theo
                      sp.mF[2:end],
                      sp.cor[2:end],
                      ideal[2:end])
            # write(io, plst[2:end])
            # write(io, nlst[2:end])
        end
    end
end
