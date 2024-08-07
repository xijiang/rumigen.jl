using Test, rumigen

tmp = tempdir()
tdir = "$tmp/rumigen-test"
isdir(tdir) && rm(tdir, recursive = true)
mkpath(tdir)

#=
@testset "macs" begin
    macs = rumigen.make_macs(tdir=tdir)
    @test basename(macs) == "macs"
end
=#

@testset "Matrix I/O" begin
    nid, nlc = 10, 20
    hps = rumigen.quickhap(nlc, nid)'
    hpf = "$tdir/hps.xy"
    rumigen.writexy(hpf, hps, major = 1)
    @test isfile(hpf)
    @test filesize(hpf) == 24 + 2nid * nlc * sizeof(eltype(hps))

    m2 = rumigen.readxy(hpf)
    @test hps == m2
    bar = rumigen.uniqSNP(hpf)
    m3 = rumigen.readxy("$tdir/$bar-uhp.xy")
    m3 .%= 2
    m3 = Int8.(m3)
    @test m3 == hps
end

rm(tdir, recursive = true)
