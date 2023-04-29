"""
    function cattle_base(nid, rst)
Create `2nid` haplotypes with MaCS, and convert them to the `xy` format.
"""
function effb9_base(nsir, ndam, rst)
    tprintln("Generating a base population with MaCS")
    macs = make_macs(tdir = rst)
    nid = nsir + ndam
    cattle_genome(macs, nid, dir = "$rst/base")
    macs2xy("$rst/base")
end

"""
    function cattle_fdr(rst, dir, foo; nlc = 50_000)
Sample founders from the base population and mate them very randomly into F1.
"""
function effb9_fdr(rst, dir, foo, nsir, ndam; nlc = 50_000, nqtl = 5_000)
    @info "Sample founders from the base population created with MaCS"
    nhp = 2(nsir + ndam)
    bar = sampleFdr("$rst/base/$foo-hap.xy", "$rst/base/$foo-map.ser",
                    nhp, nlc = nlc, nqtl = nqtl, dir = dir)
    return bar
end

function effb9_f0(dir, bar, nsir, ndam; noff = 21)
    @info "Mate founders to generate the first generation"
    lmp = deserialize("$dir/$bar-map.ser")
    lms = sumMap(lmp)
    prt = randomMate(repeat(1:nsir, noff), repeat(nsir+1:nsir+ndam, noff)) # parents of 10,500 calves
    pg = xymap("$dir/$bar-hap.xy")
    open("$dir/$bar-f0.xy", "w+") do io
        ihdr = xyhdr("$dir/bar-hap.xy")
    end
    #drop(pg, og, prt, lms)
end