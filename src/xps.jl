# Table of contents
# 1. _gsQuinton1991: Repeat Quinton et al. 1991, with genomic selection

"""
    function _gsQuinton1991(; rst="rst", nsir=500, ndam=500, bar = "", describe = true)
Simulation 500 males and 500 females.
"""
function _gsQuinton1991(;
    rst="rst",
    nsir=500,
    ndam=500,
    bar="", # can be "85R6b"
    describe=true)
    describe && tprintln(Term.parse_md(read("docs/effb9.md", String)))

    # 1. A base population from MaCS. Can be from SLiM later
    nid, dir = nsir + ndam, "$rst/Quinton"
    isdir(dir) || mkpath(dir)
    bar == "" && begin
        tprintln("Generating a base population with MaCS")
        macs = make_macs(tdir=rst)
        cattle_genome(macs, nid, dir="$rst/base")
        macs2xy("$rst/base")
    end

    # Select 50k SNP from the base population
    tprintln("Sampling 50k SNPs from the base population")
    br2 = sampleHap("$rst/base/$bar-hap.xy", "$rst/base/$bar-map.ser", nlc=50_000, dir=dir)

    # First generation
    lmp = deserialize("$dir/$br2-map.ser")
    lms = sumMap(lmp)
    prt = begin
        tmp = randomMate(nsir, ndam)
        repeat(tmp, inner=(21, 1))
    end

end

"""
    function _debug(; b2 = "0QNlmg", dir = "rst/Quinton")
A step by step function to add codes to other `xps` functions.
"""
function _debug(; b2="0QNlmg", dir="rst/Quinton")
    nsir, ndam = 500, 500
    lmp = deserialize("$dir/$b2-map.ser")
    lms = sumMap(lmp)
    
    # First generation
    prt = begin
        tmp = randomMate(nsir, ndam)
        repeat(tmp, inner=(21, 1))
    end
    ihdr = readhdr("$dir/$b2-hap.xy")
    mt, et, mj, rs, cs = xyhdr(ihdr)
    oxy = "$dir/$b2-f0.xy"
    open(oxy, "w") do io
        ohdr = xyheader('X', 'Y', ' ', mt, 0, ihdr.e, '\n', '\n', rs, cs)
        write(io, ohdr)
        for i in 1:nsir
            write(io, prt[i, :])
        end
    end
end
