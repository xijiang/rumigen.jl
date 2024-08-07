abstract type Variation end

struct SNP <: Variation
    a::Int8
end

struct USNP <: Variation
    a::Int8
    id::UInt16
end

# the types that might be used in matrices for breeding
const vldtypes = (
    Bool,
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    UInt128,
    Float16,
    Float32,
    Float64,
)

"""
This struct is designed for long term usage.  The idea magic header
was from plink bed file.  Matrix I/O is necessary as data nowadays can
be huge.  Disk I/O is not avoidable.

This header is usually for genotype storage. Whether the genotypes are
loci-majored or ID-majored is determined by the `t` field.
When loci-majored, the genotypes is ease for dropping operations, and
G-BLUP.  When individual-majored, the genotypes is ease for SNP-BLUP.

To initialize a dummy header:
```
hdr = xyheader(zeros(Int, 10)...)
```
"""
struct xyheader
    x::Int8                     # 'X'
    y::Int8                     # 'Y'
    s::Int8                     # a space 0x20
    f::Int8                     # F, L, U, or S: matrix type
    t::Int8                     # 0 for loci-majored, 1 for ID-majored, or else
    e::Int8                     # eltype
    r::Int8                     # '\n', reserved for future use
    u::Int8                     # '\n', reserved for future use
    m::Int64                    # nrow, seek(_, 8) to reach here
    n::Int64                    # ncol
end

"""
For two ID pairs in a pedigree.
Using `Int32` should be enough for record pedigree.
"""
struct parent
    pa::Int32
    ma::Int32
end

# \number is octal, default eltype is Int8
xyheader(nrow::Int64, ncol::Int64) = xyheader(collect("XY F\0\2\n\n")..., nrow, ncol)
# otherwise use explicit xyheader constructor
# use `code = findfirst(x -> x == type, vldtype)` to get the code integer

# header I/O can use the base.read and write
# write("file.bin", Ref(hdr))
# read!("file.bin", Ref(hdr))

"""
    function readhdr(io::IOStream)
Read header from an IOStream in XY format.  Returns an `xyheader`.
"""
function readhdr(io::IOStream)
    x = zeros(Int8, 24)
    read!(io, x)
    reinterpret(xyheader, x)[1]
end

"""
    function readhdr(io::AbstractString)
Read header from a file in XY format.  Returns an `xyheader`.
"""
function readhdr(io::AbstractString)
    isfile(io) || error("File $io doesn't exist")
    filesize(io) ≤ 24 && error("File $io too small to be an 'XY' matrix")
    open(io, "r") do f
        readhdr(f)
    end
end

"""
    function xyhdr(hdr::xyheader)
Interpret the header and return matrix storage style (``FLUS``), eltype,
loci/ID-Majored and dimensions.
"""
function xyhdr(hdr::xyheader)
    join(Char.([hdr.x, hdr.y, hdr.s])) == "XY " || error("Not an XY file")
    mt = Char(hdr.f)
    mt ∈ "FLUS" || error("Matrix type $mt not defined")
    et = vldtypes[hdr.e]
    mt, et, hdr.t, hdr.m, hdr.n
end

"""
    function readxy(file::AbstractString)
Read a matrix `mat` from `file` according to the specifications in its
header.
"""
function readxy(file::AbstractString)
    mat = nothing               # declare a variable
    open(file, "r") do io
        hdr = readhdr(io)
        mt, et, _, row, col = xyhdr(hdr)
        if mt == Bool
            mat = BitArray(undef, row, col)
            read!(io, mat)
        else
            mat = zeros(et, row, col)
            if mt == 'L' || mt == 'S' # matrix type
                for i = 1:col
                    read!(io, view(mat, i:row, i))
                end
                if mt == 'S'
                    for i = 1:col
                        copyto!(view(mat, i, i+1:col), view(mat, i+1:row, i)')
                    end
                end
            elseif mt == 'U'
                for i = 1:col
                    read!(io, view(mat, 1:i, i))
                end
            else
                read!(io, mat)
            end
        end
    end
    return mat
end

"""
    function writexy(file, mat::AbstractMatrix; mattp='F', major=0)
Write a matrix `mat` into `file`, with specified saving type.  Note a
symmetric matrix is written of its lower triangle.
"""
function writexy(file, mat::AbstractMatrix; mattp = 'F', major = 0)
    mattp ∈ "FLUS" || error("Matrix type $mattp not defined")
    et = findfirst(x -> x == eltype(mat), vldtypes)
    et === nothing && error("Element type $(eltype(mat)) not supported")

    mattp ∉ "LUS" && issymmetric(mat) && (mattp = 'S') # force write half of a symm. mat.
    nrow, ncol = size(mat)
    mattp ∈ "LUS" && (nrow ≠ ncol) && error("Matrix type $mattp requires square matrix")

    hdr = xyheader('X', 'Y', ' ', mattp, major, et, '\n', '\n', nrow, ncol)
    open(file, "w") do io
        write(io, Ref(hdr))
        if mattp == 'L' || mattp == 'S'
            for i = 1:ncol
                write(io, mat[i:nrow, i])
            end
        elseif mattp == 'U'
            for i = 1:ncol
                write(io, mat[1:i, i])
            end
        else
            write(io, mat)
        end
    end
    nothing
end

function samplexy(XY:AbstractString, loci, ids)
    isfile(XY) || error("File $XY doesn't exist")
    hdr = readhdr(XY)
    mt, et, mj, ir, ic = xyhdr(hdr)
    (mt == 'F' && et == Int8) || error("Only support F matrix of Int8")
    (mj == 0 && length(loci) ≤ ir && length(ids) ≤ ic ÷ 2) || error("Too many loci or IDs")
    hap = Mmap.mmap(XY, Matrix{et}, (ir, ic), 24)
    hap[loci, 2ids.-1]
end

"""
    function readbed(bed, nid)
Read genotypes in file `bed` of `nid` samples and return a 2-bit matrix.
For reference, the two-bit values and their meaning:

| Value | Meaning |
| :---: | :------ |
| 00 | Homozygous for first allele in .bim file |
| 01 | Missing genotype |
| 10 | Heterozygous |
| 11 | Homozygous for second allele in .bim file |

"""
function readbed(bed, nid)
    isfile(bed) || error("File $bed doesn't exist")
    fs = filesize(bed)
    bs = Int(ceil(nid / 4))
    @show bs
    nlc = (fs - 3) ÷ bs
    nlc * bs + 3 == fs || error("File $bed has wrong size")
    mat = nothing
    open(bed, "r") do io
        magic = Int8[0x6c, 0x1b, 0x01]
        tst = copy(magic)
        read!(io, magic)
        magic == tst || error("File $bed has no magic header")
        mat = zeros(Int8, bs, nlc)
        read!(io, mat)
    end
    mat
end

"""
    function writebed(bed, mat)
Write Int8 matrix `mat` into plink bed format.
"""
function writebed(bed, mat)
    magic = Int8[0x6c, 0x1b, 0x01]
    eltype(mat) == Int8 || error("Can't write non Int8 matrix to plink bed")
    open(bed, "w") do io
        write(io, magic)
        write(io, mat)
    end
end

"""
    function sampleFdr(ixy::AbstractString, imp::AbstractString, nhp; nlc=50_000, nqtl=5_000, dir="")
Sample founder genotypes from base haplotypes in `ixy` and `imp`.
Write the results into 3 files, the genotype file, the map file and 
a binary file indicate which loci are QTLs, and which loci are chip SNP.
"""
function sampleFdr(
    ixy::AbstractString,
    imp::AbstractString,
    nhp;
    nlc = 50_000,
    nqtl = 5_000,
    nref = 10_000,
    dir = "",
)
    # Prepare files
    (isfile(ixy) && isfile(imp)) || error("File $ixy or $imp doesn't exist")
    iseven(nhp) || error("Number of haplotypes $nhp must be even")
    dir == "" && (dir = dirname(ixy))
    bar = randstring(6)
    oxy, omp = joinpath.(dir, bar .* ["-fdr.XY", "-map.ser"])

    # Parameter check
    ihdr = readhdr(ixy)
    mt, et, mj, rs, cs = xyhdr(ihdr)
    (mt == 'F' && mj ∈ (0, 1)) || error("Not a haplotype file: $mj")
    tlc, thp = (mj == 1) ? (cs, rs) : (rs, cs)
    # to make simpler codes
    ((0 < nlc ≤ tlc) && (0 < nhp ≤ thp) && (0 ≤ nqtl < tlc) && (iseven(nhp))) || error(
        "Parameter nlc($nlc), nhp($nhp), or nqtl($nqtl) not right for tlc($tlc), and thp($thp)",
    )

    # Sample loci
    lqtl = sort(shuffle(1:tlc)[1:nqtl]) # QTL loci
    slc = sort(shuffle(1:tlc)[1:nlc])   # chip loci
    rlc = sort(shuffle(1:tlc)[1:nref])  # reference loci
    shp = begin # modified 2024-04-04, to sample 2 haps of an ID together
        nid = nhp ÷ 2
        sid = sort(shuffle(1:thp÷2)[1:nid])
        sort([2sid .- 1; 2sid])
    end
    flc = sort(unique([lqtl; slc; rlc])) # final sampled loci

    # sampled linkage map
    imp = deserialize(imp)
    mmp = imp[flc, :]
    mmp.chip = map(x -> x ∈ slc, flc)
    mmp.qtl = map(x -> x ∈ lqtl, flc)
    mmp.ref = map(x -> x ∈ rlc, flc)
    serialize(omp, select(mmp, [:chr, :pos, :chip, :qtl, :ref]))

    # sampled genotypes
    open(oxy, "w") do io # always write loci majored
        igt =
            (mj == 0) ? Mmap.mmap(ixy, Matrix{et}, (rs, cs), 24) :
            Mmap.mmap(ixy, Matrix{et}, (rs, cs), 24)'
        ohdr = xyheader('X', 'Y', ' ', mt, 0, ihdr.e, '\n', '\n', length(flc), nhp)
        write(io, Ref(ohdr))
        write(io, igt[flc, shp])
    end
    bar
end

"""
    function codesnp(iv::AbstractArray{Int8}, ov::AbstractArray{UInt16}, inc = 2)
Code SNP alleles of `0` and `1` in `iv` uniquely into `ov`.
"""
function codesnp(iv::AbstractArray{Int8}, ov::AbstractArray{UInt16}, inc = 2)
    (length(iv) == length(ov)) || error("Arrays have different sizes")
    x, y = 0, 1
    for i in eachindex(iv)
        if iv[i] == 0
            ov[i] = x
            x += inc
        else
            ov[i] = y
            y += inc
        end
    end
    nothing
end

"""
    function uniqSNP(ixy::AbstractString; LociMajored = true)
Uniquely number SNP alleles and write the results to a new file with the same
`bar` name in the same directory of `ixy`. The result file is usually for
breeding, hence it is by default loci majored to benefit from sequential I/O.
"""
function uniqSNP(ixy::AbstractString; LociMajored = true, inc = 2)
    (isfile(ixy)) || error("File $ixy doesn't exist")
    bar = begin
        bn = basename(ixy)
        ix = findfirst(x -> x == '-', bn)
        isnothing(ix) ? randstring(6) : bn[1:ix-1]
    end
    dir = dirname(ixy)
    ihdr = readhdr(ixy)
    mt, et, mj, ir, ic = xyhdr(ihdr)
    (mt ≠ 'F' || et ≠ Int8) && error("Only support F matrix of Int8")
    oxy = joinpath(dir, bar * "-uhp.XY")
    ei = findfirst(x -> x == UInt16, vldtypes) # to store unique codes for SNP alleles
    cc = LociMajored ? 0 : 1
    or, oc = (mj == cc) ? (ir, ic) : (ic, ir)
    ohdr = xyheader('X', 'Y', ' ', mt, cc, ei, '\n', '\n', or, oc)
    nhp = LociMajored ? oc : or
    nhp > 2^15 && error("Maybe too many loci to be uniquely coded with UInt16")

    gt = Mmap.mmap(ixy, Matrix{et}, (ir, ic), 24) # input genotype matrix
    open(oxy, "w+") do io
        write(io, Ref(ohdr))
        usnp = Mmap.mmap(io, Matrix{UInt16}, (or, oc), 24)
        if LociMajored && mj == cc
            for i = 1:or
                codesnp(view(gt, i, :), view(usnp, i, :), inc)
            end
        elseif !LociMajored && mj == cc
            for i = 1:oc
                codesnp(gt[:, i], usnp[:, i], inc)
            end
        end
        Mmap.sync!(usnp)
    end
    bar
end

"""
    function transxy(ixy::AbstractString, oxy::AbstractString)
Transpose `ixy` of `XY` format and write the results to `oxy`.
"""
function transxy(ixy::AbstractString, oxy::AbstractString)
    (isfile(ixy)) || error("File $ixy doesn't exist")
    ihdr = readhdr(ixy)
    mt, et, mj, ir, ic = xyhdr(ihdr)
    mj == 1 ? (mj = 0) : (mj == 0 && (mj = 1))
    ohdr = xyheader('X', 'Y', ' ', mt, mj, ihdr.e, '\n', '\n', ic, ir)
    imat = Mmap.mmap(ixy, Matrix{et}, (ir, ic), 24)
    open(oxy, "w+") do io
        write(io, Ref(ohdr))
        write(io, imat')
    end
    nothing
end

"""
    function xy2gp(ixy::AbstractString, ids, loc, σₑ)
Generate TBV and phenotypes of `ids` in `ixy`, with loci info and σₑ.
- `ixy` is the filename of XY format of haplotypes, loci majored,
- `ids` specifies the ID to be phenotyped,
- `loc` is a DataFrame, specifies the QTL, markers and their effects,
- `σₑ` specify the residual standard error.
"""
function xy2gp(ixy::AbstractString, ids, loc, σₑ)
    isfile(ixy) || error("File $ixy doesn't exist")
    ihdr = readhdr(ixy)
    _, et, mj, ir, ic = xyhdr(ihdr)
    mj == 0 || error("Only support loci majored XY format")
    nlc, nid = size(loc, 1), length(ids)
    (nlc ≤ ir && nid ≤ ic ÷ 2) || error("Too many loci or IDs")
    hap = Mmap.mmap(ixy, Matrix{et}, (ir, ic), 24)
    uhp2gp(hap, loc, σₑ)
end

function appendxy!(XY::AbstractString, mat::AbstractArray)
    hdr = readhdr(XY)
    mt, et, _, ir, ic = xyhdr(hdr)
    jr, jc = size(mat)
    (mt == 'F' && et == eltype(mat) && ir == jr) || error("Incompatible matrix")
    open(XY, "a+") do io
        write(io, mat)
        seek(io, 16)
        write(io, [ic + jc])
    end
end

"""
    function xymap(XY::AbstractString)
Map the entire matrix in `XY` and return the Matrix.
"""
function xymap(XY::AbstractString)
    mt, et, mj, ir, ic = xyhdr(readhdr(XY))
    mt == 'F' || error("Only support F matrix")
    Mmap.mmap(XY, Matrix{et}, (ir, ic), 24)
end

"""
    function _mIBD(a::T, b::T, c::T, d::T) where T<:AbstractVector
Given the 4 haplotypes of a pair of individuals, calculate the mean IBD.
This function ignores bounder check, and is set internal.
"""
function _mIBD(a::T, b::T, c::T, d::T) where {T<:AbstractVector}
    r = 0.0
    r += sum(a .== c)
    r += sum(a .== d)
    r += sum(b .== c)
    r += sum(b .== d)
    r /= 2length(a)
end

"""
    function gametemat(XY::AbstractString, loc, ids)
Mean gamete matrix (mgm) uss the IBD info stored in the `XY` file to generate a
gamete matrix. It them average the 4 cells to calculate the relationships
between ID based on IBD information. In reality, these IBD information can be
obtained from dense genotypes very accurately.

- `XY` is the filename of XY format of haplotypes, loci majored,
- `loc` is a Bool vector, specifies which loci to be used,
- `ids` specifies the ID whose gamete matrix is to be generated.
"""
function gametemat(XY::AbstractString, loc::AbstractVector{Bool}, ids::AbstractVector{Int})
    issorted(ids) || sort!(ids)
    mat, nid = xymap(XY), length(ids)
    (length(loc) ≠ size(mat, 1) || 2ids[end] > size(mat, 2)) &&
        error("Not number of loci or IDs")
    ihp = sort([2ids .- 1; 2ids])
    gt = copy(mat[loc, ihp])
    IBD = zeros(nid, nid)
    Threads.@threads for (i, j) in [(i, j) for i = 1:nid for j = 1:i]
        IBD[j, i] =
            IBD[i, j] = _mIBD(
                view(gt, :, 2i - 1),
                view(gt, :, 2i),
                view(gt, :, 2j - 1),
                view(gt, :, 2j),
            )
    end
    IBD
    #=
    Threads.@threads for i in 1:nid
        for j in 1:i-1
            A[i, j] += sum(gt[:, 2i-1] .== gt[:, 2j-1])
            A[i, j] += sum(gt[:, 2i-1] .== gt[:, 2j])
            A[i, j] += sum(gt[:, 2i] .== gt[:, 2j-1])
            A[i, j] += sum(gt[:, 2i] .== gt[:, 2j])
            A[j, i] = A[i, j]
        end
        A[i, i] = 2nlc + 2sum(gt[:, 2i-1] .== gt[:, 2i])
    end
    A ./= 2nlc
    =#
end

"""
    function gametemat(XY::AbstractString, loc, ids, jds)
Calculate mean IBD relationships between `ids` and `jds` using the IBD
information stored in the `XY` file. It them average the 4 cells to calculate
the relationships between ID based on IBD information. In reality, these IBD
information can be obtained from dense genotypes very accurately.
"""
function gametemat(
    XY::AbstractString, # uniquely coded genotypes
    loc::AbstractVector{Bool}, # specify which loci to be used
    ids::AbstractVector{Int}, # specify which IDs to be used
    jds::AbstractVector{Int}, # specify which IDs to be used
)
    issorted(ids) || sort!(ids)
    issorted(jds) || sort!(jds)
    mat, mid, nid = xymap(XY), length(ids), length(jds)
    (length(loc) ≠ size(mat, 1) || 2ids[end] > size(mat, 2) || 2jds[end] > size(mat, 2)) &&
        error("Not number of loci or IDs")
    ihp = sort([2ids .- 1; 2ids])
    jhp = sort([2jds .- 1; 2jds])
    igt = copy(mat[loc, ihp])
    jgt = copy(mat[loc, jhp])
    IBD = zeros(mid, nid)
    Threads.@threads for (i, j) in [(i, j) for i = 1:mid for j = 1:nid]
        IBD[i, j] = _mIBD(
            view(igt, :, 2i - 1),
            view(igt, :, 2i),
            view(jgt, :, 2j - 1),
            view(jgt, :, 2j),
        )
    end
    IBD
    #=
    A = zeros(mid, nid)
    Threads.@threads for i in 1:mid
        for j in 1:nid
            A[i, j] += sum(igt[:, 2i-1] .== jgt[:, 2j-1])
            A[i, j] += sum(igt[:, 2i-1] .== jgt[:, 2j])
            A[i, j] += sum(igt[:, 2i] .== jgt[:, 2j-1])
            A[i, j] += sum(igt[:, 2i] .== jgt[:, 2j])
        end
    end
    A ./= 2nlc
    =#
end

"""
    function _vcfstr2allele(str, loc, sep)
Interpret a line of VCF file and return the chromosome number, position, and
allele frequency. It also fills the `loc` array with the parsed genotypes.
"""
function _vcfstr2allele(str, loc, sep)
    v = split(str, sep)
    id = 0
    for gt in v[10:end]
        id += 1
        loc[2id-1] = parse(Int8, gt[1])
        loc[2id] = parse(Int8, gt[3])
    end
    parse(Int8, v[1]), parse(Int32, v[2]), sum(loc) / 2id
end
