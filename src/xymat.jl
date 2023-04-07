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
    r::Int8                     # '\n', reserved for future use
    u::Int8                     # '\n', reserved for future use
    m::Int64                    # nrow, seek(_, 7) to reach here
    n::Int64                    # ncol, seek(_, 15)
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
xyheader(nrow::Int64, ncol::Int64) = xyheader(collect("XY FN\2\12\12")..., nrow, ncol)
# otherwise use explicit xyheader constructor
# use `code = findfirst(x -> x == type, vldtype)` to get the code integer

# header I/O can use the base.read and write
# write("file.bin", Ref(hdr))
# read!("file.bin", Ref(hdr))

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
        hdr = xyheader(zeros(Int, 10)...)
        read!(io, Ref(hdr))
        mtype = vldtype[hdr.e]
        if mtype == Bool           # eltype
            mat = BitArray(undef, hdr.m, hdr.n)
            read!(io, mat)
        else
            mat = zeros(mtype, hdr.m, hdr.n)
            if mt == 'L' || mt == 'S' # matrix type
                for i in 1:ncol
                    read!(io, view(mat, i:nrow, i))
                end
                if mt == 'S'
                    for i in 1:ncol
                        copyto!(view(mat, i, i+1:ncol), view(mat, i+1:nrow, i)')
                    end
                end
            elseif mt == 'U'
                for i in 1:ncol
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
    et = findfirst(x -> x == eltype(mat), vldtype)
    et === nothing && error("Element type $(eltype(mat)) not supported")

    issymmetric(mat) && (mattp = 'S') # force write half of a symm. mat.
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
        mt, _, et, nrow, ncol = readhdr(io)
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

    mtx, ntx, etx, nrwx, nclx = readhdr(fa)
    mty, nty, ety, nrwy, ncly = readhdr(fb)
    (mtx == 'F' && mty == 'F') || error("Can't merge triangles")
    (etx == ety) || error("Not of t`he same eltype matrices")
    trans = (ntx == nty) ? 'N' : 'T'

    mat = Mmap.mmap(fb, Matrix{ety}, (nrwy, ncly), 24)
    xymerge(fa, mat, trans=trans)
end

function macs2xy(dir, xy)
    isdir(dir) || error("Directory $dir doesn't exist")
end
