"""
    function vcfdim(vcf::AbstractString)
Find the dimensions of a VCF file.
Returns the number of loci and the number of individuals.
"""
function vcfdim(vcf::AbstractString)
    @info "Counting loci and individuals in $vcf"
    nlc, nid = 0, 0
    open(vcf, "r") do io
        for line in eachline(io)
            line[2] == '#' && continue
            if line[1] == '#'
                nid = length(split(line)) - 9
                continue
            end
            nlc += 1
            nlc % 100_000 == 0 && print("\r\tn_ID = $nid; n_Loci = $nlc")
        end
    end
    println("\r\tn_ID = $nid; n_Loci = $nlc\n")
    nlc, nid
end

"""
    function findnth(s::AbstractString, c::Union{AbstractString, AbstractChar}, n::Int)
Return the index of the first character of the nth occurrence of `c` in `s` and
the number of occurrences before the return. If the nth `c` is not found, return
`(nothing, i)`, where `i` is the number of accurances.
"""
function findnth(s::AbstractString, c::Union{AbstractString, AbstractChar}, n::Int)
    i, f = 0, 0
    while i < n
        t = findnext(c, s, f + 1)
        isnothing(t) && break
        f = t[1]
        i += 1
    end
    i < n && error("Only $i occurrences of $c found in $s")
    f, i
end

"""
    function vcfln2av(line::AbstractString, av::Vector{Int8}; sep = '\t')
Convert a line of VCF file into a tuple of chromosome, position, frequency of allele 1.
Write the alternative alleles into a vector of Int8.
"""
function vcfln2av(line::AbstractString, av::AbstractVector{Int8}; sep = '\t')
    f, _ = findnth(line, sep, 1)
    chr = parse(Int8, line[1:f-1])
    f1, _ = findnth(line, sep, 2)
    pos = parse(Int32, line[f+1:f1-1])

    f, _ = findnth(line, sep, 9)
    n, k = length(line), 0
    for i in f+1:4:n
        k += 1
        av[2k - 1] = line[i]
        av[2k] = line[i + 2]
    end
    av .-= 48  # '0' = 48

    return chr, pos, mean(av)
end

"""
    function vcf2xy(vcf::AbstractString, xy::AbstractString)
Convert a VCF file into a `xy-hap.xy` and `xy-map.ser`.
It deals with 10k loci parallelly at a time.
"""
function vcf2xy(vcf::AbstractString, xy::AbstractString; nln = 10000)
    nlc, nid = vcfdim(vcf)
    hdr = xyheader(nlc, 2nlc)

    mmp = DataFrame(chr=zeros(Int8, nlc), pos=zeros(Int32, nlc), frq=zeros(nlc)) # map
    write(xy * "-hap.xy", Ref(hdr))
    
    # Create blocks for parallel processing
    bs = blksz(nlc, nln)
    blks = collect(bs:bs:nlc)
    blks[end] < nlc && push!(blks, nlc)
    
    @info "Processing the genotypes:"
    open(`pigz -dc vcf`, "r") do ii
        for line in eachline(ii)
            line[2] ≠ '#' && break
        end     # skip header
        ilc, jlc, ibk, buf = 0, 0, 1, String[]
        open(xy * "-hap.xy", "r+") do oo
            gt = Mmap.mmap(oo, Matrix{Int8}, (nlc, 2nid), 24)
            for line in eachline(ii)
                jlc += 1
                push!(buf, line)
                if jlc == blks[ibk]
                    Threads.@threads for i in eachindex(buf)
                        mmp[ilc + i, :] = vcfln2av(buf[i], view(gt, i + ilc, :))
                    end
                    print("\r\tProcessed loci: $jlc / $nlc")
                    ibk += 1
                    ilc = jlc
                    empty!(buf)
                end
            end
        end
    end
    println('\n')
    serialize("$xy-map.ser", mmp)
    first(mmp)
end

"""
    function vcf1stln(vcf::AbstractString)
A temporary function to find the first line of a VCF file that is not a header.
"""
function vcf1stln(vcf::AbstractString)
    open(vcf, "r") do io
        for line in eachline(io)
            line[1] ≠ '#' && return line
        end
    end
end
