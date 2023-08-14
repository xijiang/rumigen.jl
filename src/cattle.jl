"""
    function cattle_genome(macs, nid::Int; dir=pwd())
Simulate SNP on cattle autosomes.
- `macs` is the path of `macs` executable.
- `nid` is the number of individuals in a base population.
- `dir` is the directory to store the results.

The chromosome sizes are from https://www.ncbi.nlm.nih.gov/assembly/?term=bos+taurus.
"""
function cattle_genome(macs, nid::Int; dir=pwd())
    nchr = 29
    # https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_642d5d40ceff2e2c64293c60&QueryKey=1&ReleaseType=RefSeq&FileType=GENOME_FASTA&Flat=true
    chrs = [158534110, 136231102, 121005158, 120000601, 120089316, 117806340,
        110682743, 113319770, 105454467, 103308737, 106982474, 87216183,
        83472345, 82403003, 85007780, 81013979, 73167244, 65820629, 63449741,
        71974595, 69862954, 60773035, 52498615, 62317253, 42350435, 51992305,
        45612108, 45940150, 51098607]

    isdir(dir) || mkpath(dir)
    wdir = mktempdir(dir)
    tprintln("  - Simulating founder cattle data in parallele")

    seed = rand(Int32, nchr)
    Threads.@threads for chr in 1:nchr
        cmd = `$macs $(2nid) $(chrs[chr]) -s $(seed[chr]) -t 9E-6 -r 3.6E-6
                      -eN 0.011 1.33 -eN 0.019 2.78 -eN 0.036 3.89 -eN 0.053 11.11
                      -eN 0.069 16.67 -eN 0.431 22.22 -eN 1.264 27.78 -eN 1.819 38.89
                      -eN 4.875 77.78 -eN 6.542 111.11 -eN 9.319 188.89 -eN 92.097 688.89
                      -eN 2592.097 688.89`
        tprint(" $chr")
        run(pipeline(cmd,
            stderr=joinpath(wdir, "info.$chr"),
            stdout=joinpath(wdir, "chr.$chr")))
    end
    println()
    wdir
end

"""
    function cattle_base(nid, rst)
Generate a cattle base population with MaCS.
"""
function cattle_base(nid, rst)
    tprintln("Generating a base population with MaCS")
    macs = make_macs(tdir = rst)
    tmp = cattle_genome(macs, nid, dir = "$rst")
    macs2xy(tmp, swap = true)  # returns foo
end

"""
    function cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, nref; d = Normal())
Sample founders from the base population and mate them very randomly into F1.
"""
function cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, nref; d = Normal())
    @info "Sample founders from the base population created with MaCS"

    # sample haplotypes for the founder population
    nhp = 2ppsz
    bar = sampleFdr("$fdr/$foo-hap.xy", "$fdr/$foo-map.ser",
                    nhp, nlc=nlc, nqtl=nqtl, nref=nref, dir=dir)
    simQTL("$dir/$bar-fdr.xy", "$dir/$bar-map.ser", d = d)
    uniqSNP("$dir/$bar-fdr.xy")
    return bar
end
