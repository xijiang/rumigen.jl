#=
Theo came with this idea that what if we use the real IBD matrix instead of G.
He predicts that the EBV accuracy will be lower.
I added two more scenarios below to the simulation. That is, IBD -> EBV and IBD, IBD -> EBV.
It may need just a few hours to finish.

echo EBV with IBD matrix | md5sum -> 75aee9e63f68406d9b6259891a5799bb
=#

function xps_75aee9( ;
    nlc = 50_000,
    nqtl = 10_000,
    nref = 10_000,
    ngrt = 20,
    ΔF = 0.011,
    nrpt = 1,
    rst = "rst",
    ppsz = 200,
    h² = 0.25,
    σₐ = 1.0,
    nsir = 25,
    ndam = 50,
    pres = 5,
    dist = Normal(),
    sim = "75aee9",
    quick_test=true,
    keep = false
    )
    
    dir = "$rst/$sim"
    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) && rm(dir, recursive=true, force=true)
    mkpath(dir)
    serialize("$dir/par.ser", (nlc=nlc, nqtl=nqtl, nref=nref, ngrt=ngrt, ΔF=ΔF,
        nrpt=nrpt, rst=rst, ppsz=ppsz, h²=h², σₐ=σₐ, nsir=nsir, ndam=ndam,
        pres=pres, dist=dist, sim=sim, quick_test=quick_test))
    # use Dict(pairs(par)) to reconstruct the par dict, remember to include Distributions
    scheme_1 = ("ran", "spd", "sgs", "sis", "sms")
    scheme_ocs = ("oap", "oag", "ogg", "oig", "oii")

    # The working parts
    @info "Simulation begins"
    for irpt in 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        fdr, foo = prepFdr(rst, quick_test, ppsz)

        # random selectio for a few generations
        bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, nref, d=dist)
        lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
        ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg=-pres)
        simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ, 1)

        @info "Calculating IBD matrix of initial population"
        G = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:size(ped, 1))
        for sel in ("sis", "oii")
            write("$dir/$bar-$sel.bin", G)
        end
        G = nothing
        # selection without optimum contribution schemes
        for op in 2:4
            pop, sel = copy(ped), scheme_1[op]
            cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force=true)
            simpleSelection("$dir/$bar-$sel.xy", pop, lmp, nsir, ndam, ngrt, σₑ, op, mp=false)
            sumPed(rst, sim, bar, lmp, sel)
            pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        for op in 1:5
            pop, sel = copy(ped), scheme_ocs[op]
            cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force=true)
            optSelection("$dir/$bar-$sel.xy", pop, lmp, ngrt, σₑ, ΔF, op=op, k₀=0.027)
            sumPed(rst, sim, bar, lmp, sel)
            pos_qtl_frq(rst, sim, bar, sel, ppsz)
        end
        keep || rm.(glob("$dir/$bar-*"), force=true)
    end
end

## echo Check random mating and delta F | md5sum -> 49decb5a5aa5bb7bdf7761b54e3fd2aa
function xps_49decb(;
    nlc = 50_000,
    nqtl = 10_000,
    nref = 10_000,
    nrpt = 1,
    rst = "rst",
    ppsz = 200,
    h² = 0.25,
    σₐ = 1.0,
    nsir = 25,
    ndam = 50,
    pres = 5,
    ngrt = 20,
    dist = Normal(),
    sim = "49decb",
    quick_test=true,
    keep = false
)
    dir = "$rst/$sim"
    σₑ = sqrt((1 - h²) / h²) * σₐ
    isdir(dir) && rm(dir, recursive=true, force=true)
    mkpath(dir)
    # make sure ppsz is divisible by bigger of nsir and ndam
    ppsz = ppsz ÷ max(nsir, ndam) * max(nsir, ndam)
    serialize("$dir/par.ser", (nlc=nlc, nqtl=nqtl, nref=nref,
        nrpt=nrpt, rst=rst, ppsz=ppsz, h²=h², σₐ=σₐ, nsir=nsir, ndam=ndam,
        pres=pres, ngrt=ngrt, dist=dist, sim=sim, quick_test=quick_test))
    # The working parts
    @info "Simulation begins"
    for irpt in 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        fdr, foo = prepFdr(rst, quick_test, ppsz)

        # random selectio for a few generations
        bar = cattle_founder(fdr, dir, foo, ppsz, nlc, nqtl, nref, d=dist)
        lmp = deserialize("$dir/$bar-map.ser") # each sample has its own map
        ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ, fg=-pres)
        simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, pres, σₑ, 1)
        pop, sel = copy(ped), "ran"
        cp("$dir/$bar-uhp.xy", "$dir/$bar-$sel.xy", force=true)
        simpleSelection("$dir/$bar-$sel.xy", pop, lmp, nsir, ndam, ngrt, σₑ, 1, mp=false)
        sumPed(rst, sim, bar, lmp, sel)
        pos_qtl_frq(rst, sim, bar, sel, ppsz)
        keep || rm.(glob("$dir/$bar-*"), force=true)
    end
end

# echo Try to find correct effective size | md5sum => a13c84043a9efcc583d989d524748c56
function xps_a13c84(;
    nlc  = 50_000,
    nref = 10_000, # no need of QTL
    nrpt = 1,
    rst  = "rst",
    nsir = 25,
    ndam = 50,
    ngrt = 25,
    sim  = "a13c84",
    quick_test=true,
)
    dir = "$rst/$sim"
    isdir(dir) && rm(dir, recursive=true, force=true)
    mkpath(dir)
    serialize("$dir/par.ser", (nlc=nlc, nref=nref, nrpt=nrpt, rst=rst,
            nsir=nsir, ndam=ndam, ngrt=ngrt, sim=sim, quick_test=quick_test))
    nid = nsir + ndam
    for irpt in 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        ped = randped(nsir, ndam, ngrt)
        bar = sampleFdr("$rst/founder/CuaL5-hap.xy", "$rst/founder/CuaL5-map.ser",
                        2nid, nlc=nlc, nqtl=0, nref=nref, dir=dir)
        uniqSNP("$dir/$bar-fdr.xy")
        lmp = deserialize("$dir/$bar-map.ser")
        lms = sumMap(lmp)
        lmp = deserialize("$dir/$bar-map.ser")
        f0q = begin
            hdr = readhdr("$dir/$bar-fdr.xy")
            _, et, _, nr, nc = xyhdr(hdr)
            gt = Mmap.mmap("$dir/$bar-fdr.xy", Matrix{et}, (nr, nc), 24)
            mean(gt, dims = 2)
        end
        @info "Selection on $bar, for $ngrt generations"
        cax, cmx, rax, rmx = [0.], [0.], [0.], [0.] # chip and reference loci fixed
        for igrt in 1:ngrt
            print(" $igrt")
            hdr = readhdr("$dir/$bar-uhp.xy")
            _, et, _, nr, nc = xyhdr(hdr)
            op = view(ped, ped.grt .== igrt, :)
            agt = Mmap.mmap("$dir/$bar-uhp.xy", Matrix{et}, (nr, nc), 24)
            ogt = zeros(et, nr, 2nid)
            drop(agt, ogt, [op.ma op.pa], lms)
            appendxy!("$dir/$bar-uhp.xy", ogt)
            frq = sum(isodd.(ogt), dims = 2)
            fx = frq .== 0 .|| frq .== 2nid
            mx = fx .&& ((f0q .< 0.2) .|| (f0q .> 0.8))
            push!(cax, sum(fx .&& lmp.chip))
            push!(rax, sum(fx .&& lmp.ref))
            push!(cmx, sum(mx .&& lmp.chip))
            push!(rmx, sum(mx .&& lmp.ref))
        end
        println()
        ped.F = inbreeding("$dir/$bar-uhp.xy", lmp.chip)
        ped.Fr = inbreeding("$dir/$bar-uhp.xy", lmp.ref)
        ped.Fp = diag(Amat(ped)) .- 1
        open("$dir/$sim.bin", "a") do io
            mF = Float64[]
            mFr = Float64[]
            mFp = Float64[]
            for grt in groupby(ped, :grt)
                push!(mF, mean(grt.F))
                push!(mFr, mean(grt.Fr))
                push!(mFp, mean(grt.Fp))
            end
            write(io, mF, mFr, mFp, cax, cmx, rax, rmx)
        end
        rm.(glob("$dir/$bar-*"), force=true)
    end
end

"""
echo Check relationship by real IBD info | md5sum -> 64d3f8d6cc4ccf8e58a45cc501f2d2d5
"""
function xps_64d3f8()
    dir, fdr, foo = "rst/64d3f8", "rst/test-suite", "founder"
    isdir(dir) || mkpath("rst/64d3f8")
    bar = cattle_founder(fdr, dir, foo, 200, 50_000, 10_000, 10_000, d = Normal())
    lmp = deserialize("$dir/$bar-map.ser")
    ped = initPedigree("$dir/$bar-uhp.xy", lmp, 1.0)
    simpleSelection("$dir/$bar-uhp.xy", ped, lmp, 20, 50, 9, 1.0, 1)
    A = Amat(ped)
    M = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:size(ped, 1))
    serialize("$dir/$bar-A.imt", A)
    serialize("$dir/$bar-M.imt", M)
    M2 = gametemat("$dir/$bar-uhp.xy", lmp.chr .== 1, 1:size(ped, 1))
    serialize("$dir/$bar-M2.imt", M2)
end

"""
echo test function gametemat | md5sum
ddbbc13870afca0c485058c6fe761a1c
"""
function xys_ddbbc1()
    fdr, dir = "rst/test-suite", "rst/ddbbc1"
    isdir(dir) && rm(dir, recursive=true, force=true)
    mkpath(dir)
    bar = cattle_founder(fdr, dir, "founder", 200, 50_000, 10_000, 10_000, d = Normal())
    lmp = deserialize("$dir/$bar-map.ser")
    ped = initPedigree("$dir/$bar-uhp.xy", lmp, 1.0)
    simpleSelection("$dir/$bar-uhp.xy", ped, lmp, 20, 50, 9, 1.0, 1)
    n = size(ped, 1)
    M = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:n, 1:n)
    A = gametemat("$dir/$bar-uhp.xy", lmp.chip, 1:n)
    return M, A
end
