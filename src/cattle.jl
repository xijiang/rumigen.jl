"""
    function cattle_genome(macs, nid::Int, dir=pwd())
Simulate SNP on cattle autosomes.
- `macs` is the path of `macs` executable.
- `nid` is the number of individuals in a base population.
- `dir` is the directory to store the results.

The chromosome sizes are from https://www.ncbi.nlm.nih.gov/assembly/?term=bos+taurus.
"""
function cattle_genome(macs, nid::Int, dir=pwd())
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
    function macs_2_hap(raw)
Merge and convert `MaCS` results to `01` allele types of `nHap` by `nLoc`.
Each column in the result file is a locus.
Path `raw` stores the `MaCS` results.
Results are written in the parent directory of `raw`.
File names are of format `chr.1`, `.2`, etc, for genotypes.
And `info.1`, etc for error message from `macs`.
The linkage map is a DataFrame,
and serialized to `map.ser` in the parent dir of `raw` also.

## Binary file `macs.gt`:
- nhap, nlc, 1: the first 3 × 8 bytes. 1 is for `Int8`.
- then `01` allele types

Note, elsewhere I use `nlc×nid`, or `nhap×nid` dimensions.

This function has a very low memory footprint. You can use `Fio.transmat` function to transpose the file
resulted from `macs_2_hap`.
"""
function macs_2_hap(raw)
    # @info "this results in a nHap×nloci matrix."
    bar = randstring(5)         # barcode of this simulation
    tprintln("  - Collecting founder data {cyan}$bar{/cyan} from macs of chromosome: ")
    isdir(raw) || error("$raw not exists")
    raw[end] == '/' && (raw = raw[1:end-1])
    chrs = Int8[]
    for f in readdir(raw)
        occursin.(r"^chr", f) && push!(chrs, parse(Int8, split(f, '.')[2]))
    end
    sort!(chrs)           # chromosome number in integer, and in order
    parent = dirname(raw) # path
    fgt = joinpath(parent, "$bar-hap.bin")
    tmap = DataFrame(chr=Int8[], pos=Int64[], frq=Float64[])
    open(fgt, "w") do io
        write(io, [0, 0, 1]) # places for `nhap, nloc` and type, overwitten later
        nlc, as = 0, nothing # to count nID. `as` is for alleles
        for ic in chrs
            this_chr = joinpath(raw, "chr.$ic")
            tprint(" $ic")
            nbp = parse(Float64, split(readline(this_chr))[4])
            for line in eachline(this_chr)
                line[1:4] ≠ "SITE" && continue
                _, _, pos, _, as = split(line)
                pos = Int(round(parse(Float64, pos) * nbp))
                as = parse.(Int8, collect(as))
                frq = mean(as)
                write(io, as)
                push!(tmap, (ic, pos, frq))
                nlc += 1
            end
        end
        nhp = length(as)
        seekstart(io)
        write(io, [nhp, nlc])
    end
    println()
    serialize(joinpath(parent, "$bar-map.ser"), tmap)
    return bar
end
