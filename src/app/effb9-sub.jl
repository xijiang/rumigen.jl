"""
    function base_effb9(nsir, ndam, rst)
Create `2nid` haplotypes with MaCS, and convert them to the `xy` format.
"""
function base_effb9(nsir, ndam, rst)
    tprintln("Generating a base population with MaCS")
    macs = make_macs(tdir = rst)
    nid = nsir + ndam
    cattle_genome(macs, nid, dir = "$rst/base")
    macs2xy("$rst/base")
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
    prt = randomMate(nsir, ndam, noff = nbull + ncow) # parents of 10,500 calves
    drop("$dir/$bar-hap.xy", "$dir/$bar-f0.xy", prt, lms)
    uniqSNP("$dir/$bar-f0.xy")
end

function cmn_effb9(dir, bar)
    # Genotypes of F0, the starting point for downstream selection
    f0, oo = "$dir/$bar-uhp.xy", "$dir/$bar"
    isfile(f0) || error("$f0 doesn't exist")
    cp(f0, "$oo-pht.xy")
    cp(f0, "$oo-ped.xy")
    mv(f0, "$oo-gs.xy")
    #ToDo: Create a pedigree file with ID, Pa, Ma, TBV, phenotype
    # and make two copies of it.
end
