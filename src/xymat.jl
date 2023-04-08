# the types that might be used in matrices for breeding
const vldtypes = (Bool, Int8, Int16, Int32, Int64, Int128, UInt8, UInt16,
                  UInt32, UInt64, UInt128, Float16, Float32, Float64)

"""
This struct is designed for long term usage.  The idea magic header
was from plink bed file.  Matrix I/O is necessary as data nowadays can
be huge.  Disk I/O is not avoidable.

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
    t::Int8                     # N, T: transpose or not
    e::Int8                     # eltype
    r::Int8                     # 0, loci majored
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

# \number is octal, default is 8 bits / Int8
xyheader(nrow::Int64, ncol::Int64) = xyheader(collect("XY FN\2\0\12")..., nrow, ncol)
# otherwise use explicit xyheader constructor
# use `code = findfirst(x -> x == type, vldtype)` to get the code integer

# header I/O can use the base.read and write
# write("file.bin", Ref(hdr))
# read!("file.bin", Ref(hdr))

"""
    function readhdr(io::IOStream)
Read header from an IOStream in XY format.  Returns matrix storage
style (``FLUS``), if transposed (``NT``), eltype, if LociMajored and dimensions.
"""
function readhdr(io::IOStream)
    tmp = zeros(Int8, 8)
    read!(io, tmp)
    join(Char.(tmp[1:3])) == "XY " || error("Not in XY format")
    dim = zeros(Int64, 2)
    read!(io, dim)
    nrow, ncol = dim
    itp = tmp[6]
    LociMajored = tmp[7] == 0
    return Char(tmp[4]), Char(tmp[5]), vldtypes[itp], LociMajored, nrow, ncol
end

function readhdr(file::AbstractString)
    open(file, "r") do io
        readhdr(io)
    end
end

"""
    function readxyhdr(xy::AbstractString)
Read header, as is, from a file in XY format.
"""
function readxyhdr(xy::AbstractString)
    isfile(xy) || error("File $xy doesn't exist")
    x = zeros(UInt8, 24)
    read!(xy, x)
    reinterpret(xyheader, x)
end

"""
    function readxy(file::AbstractString)
Read a matrix `mat` from `file` according to the specifications in its
header.
"""
function readxy(file::AbstractString)
    isfile(file) || error("File $file doesn't exist")
    filesize(file) ≤ 24 &&
        error("File $file too small to be an 'XY' matrix")
    mat = nothing               # to return if success
    open(file, "r") do io
        mt, _, et, _, row, col = readhdr(io)
        if mt == Bool
            mat = BitArray(undef, row, col)
            read!(io, mat)
        else
            mat = zeros(et, row, col)
            if mt == 'L' || mt == 'S' # matrix type
                for i in 1:col
                    read!(io, view(mat, i:row, i))
                end
                if mt == 'S'
                    for i in 1:col
                        copyto!(view(mat, i, i+1:col), view(mat, i+1:row, i)')
                    end
                end
            elseif mt == 'U'
                for i in 1:col
                    read!(io, view(mat, 1:i, i))
                end
            else
                read!(io, mat)
            end
        end
    end
    return mat
end

# below can be exported

"""
    function writexy(file, mat; mattp = 'F', trans = 'N')
Write a matrix `mat` into `file`, with specified saving type.  Note a
symmetric matrix is written of its lower triangle.
"""
function writexy(file, mat::AbstractMatrix; mattp='F', trans='N')
    mattp ∈ "FLUS" || error("Matrix type $mattp not defined")
    trans ∈ "NT" || error("ID locus type $trans not defined")
    et = findfirst(x -> x == eltype(mat), vldtypes)
    et === nothing && error("Element type $(eltype(mat)) not supported")

    mattp ∉ "LUS" && issymmetric(mat) && (mattp = 'S') # force write half of a symm. mat.
    nrow, ncol = size(mat)
    mattp ∈ "LUS" && (nrow ≠ ncol) && error("Matrix type $mattp requires square matrix")

    hdr = xyheader('X', 'Y', ' ', mattp, trans, et, '\n', '\n', nrow, ncol)
    open(file, "w") do io
        write(io, Ref(hdr))
        if mattp == 'L' || mattp == 'S'
            for i in 1:ncol
                write(io, mat[i:nrow, i])
            end
        elseif mattp == 'U'
            for i in 1:ncol
                write(io, mat[1:i, i])
            end
        else
            write(io, mat)
        end
    end
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
    xymerge(file::AbstractString, mat::Matrix; trans='N')

Append genotypes in `mat` to `file`.
"""
function xymerge(file::AbstractString, mat::AbstractMatrix; trans='N')
    isfile(file) || error("File $file doesn't exist")
    open(file, "r+") do io
        mt, _, et, _, nrow, ncol = readhdr(io)
        mt ∈ "LUS" && error("Can't merge triangles")
        et ≠ eltype(mat) && error("Matrices have different eltypes")

        seekend(io)
        if trans == 'N'
            ncl2 = size(mat, 2)
            nrow ≠ size(mat, 1) && error("Matrices don't match")
            write(io, mat)
        else
            ncl2 = size(mat, 1)
            nrow ≠ size(mat, 2) && error("Matrices don't match")
            write(io, mat')
        end

        seek(io, 16)
        write(io, [ncol + ncl2])
    end
end

"""
    xymerge(fa::AbstractString, fb::AbstractString)

Append genotypes in `mb` to file `ma`.  This is only valid if the
matrix type is `F`, and two matrices have same number of ID, or same
number of loci.
"""
function xymerge(fa::AbstractString, fb::AbstractString)
    (isfile(fa) && isfile(fb)) || error("File $fa or $fb doesn't exist")

    mtx, ntx, etx, _, nrwx, nclx = readhdr(fa)
    mty, nty, ety, _, nrwy, ncly = readhdr(fb)
    (mtx == 'F' && mty == 'F') || error("Can't merge triangles")
    (etx == ety) || error("Not of t`he same eltype matrices")
    trans = (ntx == nty) ? 'N' : 'T'

    mat = Mmap.mmap(fb, Matrix{ety}, (nrwy, ncly), 24)
    xymerge(fa, mat, trans=trans)
end

"""
    function sampleHap(ixy::AbstractString, imp::AbstractString; nid = 0, nlc = 0)
Sample `nlc` loci and `nid` ID from genotype file `ixy`, and linkage map 
file `imp`. The results are written to a new file with a random name in 
the same directory of `ixy`.
If `nid` or `nlc` is zero, sample all loci, or ID.
"""
function sampleHap(ixy::AbstractString, imp::AbstractString;
                   nhp = 0, nlc = 0)
    (isfile(ixy) && isfile(imp)) || error("File $ixy or $imp doesn't exist")
    isodd(nhp) && error("Number of haplotypes must be even")
    bar, dir = randstring(6), dirname(ixy)
    oxy, omp = joinpath.(dir, bar .* ["-hap.xy", "-map.ser"])

    mt, nt, et, LociMajored, rs, cs = readhdr(ixy)
    tlc, thp = LociMajored ? (rs, cs) : (cs, rs)
    slc = 0 < nlc < tlc ? sort(shuffle(1:tlc)[1:nlc]) : 1:tlc
    sid = 0 < nhp < thp ? sort(shuffle(1:thp)[1:nhp]) : 1:thp
    mmp = deserialize(imp)
    serialize(omp, mmp[slc, :])
    ei = findfirst(x -> x == et, vldtypes)
    hdr = nothing
    if LociMajored
        hdr = xyheader('X', 'Y', ' ', mt, nt, ei, 0, '\n', length(slc), length(sid))
        gt = Mmap.mmap(ixy, Matrix{et}, (tlc, thp), 24)
        open(oxy, "w") do io
            write(io, Ref(hdr))
            write(io, gt[slc, sid])
        end
    else
        hdr = xyheader('X', 'Y', ' ', mt, nt, ei, 1, '\n', length(sid), length(slc))
        gt = Mmap.mmap(ixy, Matrix{et}, (thp, tlc), 24)
        open(oxy, "w") do io
            write(io, Ref(hdr))
            write(io, gt[sid, slc])
        end
    end 
    bar
end

function uniqSNP(ixy::AbstractString)
    (isfile(ixy)) || error("File $ixy doesn't exist")
    bar = begin
        bn = basename(ixy)
        ix = findfirst(x -> x == '-', bn)
        isnothing(ix) ? randstring(6) : bn[1:ix-1]
    end
    dir = dirname(ixy)
    mt, nt, et, LociMajored, rs, cs = readhdr(ixy)
    (mt ≠ 'F' || et ≠ Int8) && error("Only support F matrix of Int8")
    oxy = joinpath(dir, bar * "-uhp.xy")
    ei = findfirst(x -> x == UInt16, vldtypes) # to store unique codes for SNP alleles
    cc = LociMajored ? 0 : 1
    hdr = xyheader('X', 'Y', ' ', mt, nt, ei, cc, '\n', rs, cs)
    gt = Mmap.mmap(ixy, Matrix{et}, (rs, cs), 24)
    open(oxy, "w+") do io
        write(io, Ref(hdr))
        usnp = Mmap.mmap(io, Matrix{UInt16}, (rs, cs), 24)
        if LociMajored # then count row by row
            for i in 1:rs
                x::UInt16, y::UInt16 = 0, 1
                for j in 1:cs
                    if gt[i, j] == 0
                        usnp[i, j] = x
                        x += 2
                    else
                        usnp[i, j] = y
                        y += 2
                    end
                end
            end
        else
            for i in 1:cs
                x::UInt16, y::UInt16 = 0, 1
                for j in 1:rs
                    if gt[j, i] == 0
                        usnp[j, i] = x
                        x += 2
                    else
                        usnp[j, i] = y
                        y += 2
                    end
                end
            end
        end
        Mmap.sync!(usnp)
    end
    bar
end
