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
    nid = begin
        tmp = zeros(Int, 3)
        read!("$dir/$bar--gs.xy",tmp)
        tmp[3] ÷ 2
    end

    tbv, pht = xy2gp("$dir/$bar--gs.xy", 1:nid, lmp, σₑ)
    sex = repeat([1, 1, 0, 0], outer = nid ÷ 4)
    ped = DataFrame(id = 1:nid, pa = 0, ma = 0, sex = sex, grt = 0,
                    tbv = tbv, pht = pht, ebv = 0., F = 0.)
    serialize("$dir/$bar-f0-ped.ser", ped)
    nothing
end

function sum_effb9()
    @info "Under construction"
end

