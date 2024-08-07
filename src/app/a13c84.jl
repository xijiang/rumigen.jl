# echo Try to find correct effective size | md5sum => a13c84043a9efcc583d989d524748c56
function fullsibs(nsir::T, ndam::T, nsib::T, ngrt::T) where {T<:Int}
    nfam = max(nsir, ndam)
    nid = nfam * nsib
    ped = DataFrame(pa = 0, ma = 0, sex = rand(0:1, nid), grt = 0)
    for igrt = 1:ngrt
        pp = @view(ped[ped.grt.==igrt-1, :])
        pa = begin
            tmp = findall(pp.sex .== 1) .+ (igrt - 1) * nid
            length(tmp) > nsir && (tmp = shuffle(tmp)[1:nsir])
            while length(tmp) < nfam
                append!(tmp, shuffle(tmp))
            end
            tmp[1:nfam]
        end
        ma = begin
            tmp = findall(pp.sex .== 0) .+ (igrt - 1) * nid
            length(tmp) > ndam && (tmp = shuffle(tmp)[1:ndam])
            while length(tmp) < nfam
                append!(tmp, shuffle(tmp))
            end
            tmp[1:nfam]
        end
        append!(
            ped,
            DataFrame(
                pa = repeat(pa, inner = nsib),
                ma = repeat(ma, inner = nsib),
                sex = rand(0:1, nid),
                grt = igrt,
            ),
        )
    end
    sort(ped, [:pa, :ma])
end

function xps_a13c84(;
    nlc = 50_000,
    nref = 10_000, # no need of QTL
    nrpt = 1,
    rst = "rst",
    nsir = 17,
    ndam = 34,
    nsib = 2,
    #nid = 51,
    mfr = 2,  # male:female ratio
    ngrt = 25,
    sim = "a13c84",
    maf = 0.2,
    quick_test = true,
    keep = false,
)
    dir, foo = "$rst/$sim", "8qKws"
    isdir(dir) && rm(dir, recursive = true, force = true)
    mkpath(dir)
    nid = max(nsir, ndam) * nsib
    serialize(
        "$dir/par.ser",
        (
            nlc = nlc,
            nref = nref,
            nrpt = nrpt,
            rst = rst,
            maf = maf,
            nid = nid,
            mfr = mfr,
            ngrt = ngrt,
            sim = sim,
            quick_test = quick_test,
            nsir = nsir,
            ndam = ndam,
            nsib = nsib,
        ),
    )
    for irpt = 1:nrpt
        println()
        @info "Repeat $irpt of $nrpt"
        #ped = quickped(nid, ngrt, ratio=mfr)
        ped = fullsibs(nsir, ndam, nsib, ngrt)
        bar = sampleFdr(
            "$rst/founder/$foo-hap.xy",
            "$rst/founder/$foo-map.ser",
            2nid,
            nlc = nlc,
            nqtl = 0,
            nref = nref,
            dir = dir,
        )
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
        cax, cmx, rax, rmx = [0.0], [0.0], [0.0], [0.0] # chip and reference loci fixed
        for igrt = 1:ngrt
            print(" $igrt")
            hdr = readhdr("$dir/$bar-uhp.xy")
            _, et, _, nr, nc = xyhdr(hdr)
            op = view(ped, ped.grt .== igrt, :)
            agt = Mmap.mmap("$dir/$bar-uhp.xy", Matrix{et}, (nr, nc), 24)
            y = zeros(et, nr, 2nid)
            drop(agt, y, [op.ma op.pa], lms)
            appendxy!("$dir/$bar-uhp.xy", y)
            frq = sum(isodd.(y), dims = 2)
            fx = frq .== 0 .|| frq .== 2nid
            mx = fx .&& ((f0q .< maf) .|| (f0q .> 1 - maf))
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
        keep || rm.(glob("$dir/$bar-*"), force = true)
    end
end

function exp_elfx(; qlm = 0.05) # expected loci fixed
    foo, dir = "CuaL5", "rst/a13c84"
    isdir(dir) && rm(dir, recursive = true, force = true)
    mkpath(dir)
    nlc, nsr, ndm, nsb, ngrt = 10_000, 25, 50, 2, 100
    nid = max(nsr, ndm) * nsb
    ped = fullsibs(nsr, ndm, nsb, ngrt)
    bar = sampleFdr(
        "rst/founder/$foo-hap.xy",
        "rst/founder/$foo-map.ser",
        2nid,
        nlc = nlc,
        nqtl = 0,
        nref = 0,
        dir = dir,
    )
    frq = begin
        hdr = readhdr("$dir/$bar-fdr.xy")
        _, et, _, nr, nc = xyhdr(hdr)
        gt = Mmap.mmap("$dir/$bar-fdr.xy", Matrix{et}, (nr, nc), 24)
        mean(gt, dims = 2)
    end
    ped = fullsibs(nsr, ndm, nsb, ngrt)
    uniqSNP("$dir/$bar-fdr.xy")
    lmp = deserialize("$dir/$bar-map.ser")
    lms = sumMap(lmp)
    n1, n2, fs = Int[], Int[], Float64[]
    for igrt = 1:ngrt
        # the drop part
        hdr = readhdr("$dir/$bar-uhp.xy")
        _, et, _, nr, nc = xyhdr(hdr)
        op = view(ped, ped.grt .== igrt, :)
        agt = Mmap.mmap("$dir/$bar-uhp.xy", Matrix{et}, (nr, nc), 24)
        y = zeros(et, nr, 2nid)
        drop(agt, y, [op.ma op.pa], lms)
        appendxy!("$dir/$bar-uhp.xy", y)
        # the expected loci fixed
        pm = [([[2i - 1, 2i] for i in sort(unique([op.pa; op.ma]))]...)...]
        x = view(agt, :, pm)
        f = 0.0  # inbreeding
        # for i in 1:2:size(x, 2)
        #     f += sum(x[:, i] .== x[:, i+1])
        # end
        # f /= nlc * size(x, 2)
        for i = 1:2:size(y, 2)
            f += sum(y[:, i] .== y[:, i+1])
        end
        f /= nlc * size(y, 2)
        if f > 0.0
            push!(fs, f)
            enf = 0.0
            for i in eachindex(frq)
                frq[i] < qlm && (enf += exp(-frq[i] / f))
            end
            qy = mean(isodd.(y), dims = 2)
            onf = sum(qy .== 0 .&& frq .< qlm) #+ sum(qy .== 2nid)
            @show igrt, f, Int(floor(enf)), onf, mean(exp.(-frq[frq.<qlm] / f))
            push!(n1, Int(floor(enf)))
            push!(n2, onf)
        end
    end
    @show sum(frq .< qlm), length(frq), mean(frq[frq.<qlm])
    frq, n1, n2, fs
end
