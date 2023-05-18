"""
    function base_effb9(nsir, ndam, rst)
Create `2nid` haplotypes with MaCS, and convert them to the `xy` format.
"""
function base_effb9(nsir, ndam, rst)
    tprintln("Generating a base population with MaCS")
    macs = make_macs(tdir = rst)
    nid = nsir + ndam
    tmp = cattle_genome(macs, nid, dir = "$rst/base")
    macs2xy(tmp)
end

"""
    function fdr_effb9(rst, dir, foo, nsir, ndam, nlc, nqtl)
Sample founders from the base population and mate them very randomly into F1.
"""
function fdr_effb9(rst, dir, foo, nsir, ndam, nlc, nqtl)
    @info "Sample founders from the base population created with MaCS"
    nhp = 2(nsir + ndam)
    bar = sampleFdr("$rst/base/$foo-hap.xy", "$rst/base/$foo-map.ser",
                    nhp, nlc = nlc, nqtl = nqtl, dir = dir)
    return bar
end

function f0_effb9(dir, bar, nsir, ndam, nbull, ncow)
    @info "Mate founders to generate the first generation"
    lmp = deserialize("$dir/$bar-map.ser")
    lms = sumMap(lmp)
    prt = randomMate(nsir, ndam, noff = nbull + ncow)
    drop("$dir/$bar-hap.xy", "$dir/$bar-f0.xy", prt, lms)
    uniqSNP("$dir/$bar-f0.xy")
    cp("$dir/$bar-uhp.xy", "$dir/$bar-pht.xy", force = true)
    cp("$dir/$bar-uhp.xy", "$dir/$bar-ped.xy", force = true)
    mv("$dir/$bar-uhp.xy", "$dir/$bar--gs.xy", force = true)
end

"""
    function ped_effb9(dir, bar, σₑ)
Create a pedigree file for the first generation. Calculate the TBV and phenotypes.
"""
function ped_effb9(dir, bar, σₑ)
    lmp = deserialize("$dir/$bar-map.ser")
    ped = initPedigree("$dir/$bar-f0.xy", lmp, σₑ)
    serialize("$dir/$bar-f0-ped.ser", ped)
    nothing
end

"""
    function sum_effb9(dir, bar)
Calculate mean `tbv`, `F`, of each generation, and number
of fixed loci on chip and QTL in the end.
"""
function sum_effb9(dir, bar)
    for sel in ["-gs", "pht", "ped"]
        ped = deserialize("$dir/$bar-$sel+ped.ser")
        lmp = deserialize("$dir/$bar-map.ser")
        sp = combine(groupby(ped, :grt), :tbv => mean => :tbv, :F => mean => :F)
        open("$dir/mbvf.bin", "a") do io
            write(io, sp.tbv[2:end])
            write(io, sp.F[2:end])
            ideal, plst, nlst = idealPop("$rst/$bar-$sel.xy", ped.grt, lmp)
            write(io, ideal[2:end])
            write(io, plst[2:end])
            write(io, nlst[2:end])
        end
    end
end

function stats_effb9(dir, ngrt, nlc, nqtl)
end
