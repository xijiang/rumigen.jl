using Test, rumigen

tmp = tempdir()
tdir = "$tmp/rumigen-test"
isdir(tdir) && rm(tdir, recursive=true)
mkpath(tdir)

#=
@testset "macs" begin
    macs = rumigen.make_macs(tdir=tdir)
    @test basename(macs) == "macs"
end
=#

@testset "Matrix I/O" begin
    m, n = 3, 3
    mat = rand(1:10, 3, 3)
    mfile = "$tdir/matrix.xy"
    rumigen.writexy(mfile, mat, mattp = 'L')
    @test isfile(mfile)
    @test filesize(mfile) == 24 + 6*sizeof(eltype(mat))


    #m2 = rumigen.readxy(mfile)
    #@test m == m2
end

rm(tdir, recursive=true)
