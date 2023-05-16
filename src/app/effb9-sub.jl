"""
    function base_effb9(nsir, ndam, rst)
Create `2nid` haplotypes with MaCS, and convert them to the `xy` format.
"""
function base_effb9(nsir, ndam, rst)
    tprintln("Generating a base population with MaCS")
    macs = make_macs(tdir = rst)
    nid = nsir + ndam
    tmp = cattle_genome(macs, nid, dir = "$rst/base")
    macs2xy(tmp)
end

"""
    function fdr_effb9(rst, dir, foo, nsir, ndam, nlc, nqtl)
Sample founders from the base population and mate them very randomly into F1.
"""
function fdr_effb9(rst, dir, foo, nsir, ndam, nlc, nqtl)
    @info "Sample founders from the base population created with MaCS"
    nhp = 2(nsir + ndam)
    bar = sampleFdr("$rst/base/$foo-hap.xy", "$rst/base/$foo-map.ser",
                    nhp, nlc = nlc, nqtl = nqtl, dir = dir)
    return bar
end

function f0_effb9(dir, bar, nsir, ndam, nbull, ncow)
    @info "Mate founders to generate the first generation"
    lmp = deserialize("$dir/$bar-map.ser")
    lms = sumMap(lmp)
    prt = randomMate(nsir, ndam, noff = nbull + ncow)
    drop("$dir/$bar-hap.xy", "$dir/$bar-f0.xy", prt, lms)
    uniqSNP("$dir/$bar-f0.xy")
    cp("$dir/$bar-uhp.xy", "$dir/$bar-pht.xy", force = true)
    cp("$dir/$bar-uhp.xy", "$dir/$bar-ped.xy", force = true)
    mv("$dir/$bar-uhp.xy", "$dir/$bar--gs.xy", force = true)
end

"""
    function ped_effb9(dir, bar, σₑ)
Create a pedigree file for the first generation. Calculate the TBV and phenotypes.
"""
function ped_effb9(dir, bar, σₑ)
    lmp = deserialize("$dir/$bar-map.ser")
    nid = begin
        tmp = zeros(Int, 3)
        read!("$dir/$bar--gs.xy",tmp)
        tmp[3] ÷ 2
    end

    tbv, pht = xy2gp("$dir/$bar--gs.xy", 1:nid, lmp, σₑ)
    sex = repeat([1, 1, 0, 0], outer = nid ÷ 4)
    ped = DataFrame(id = 1:nid, pa = 0, ma = 0, sex = sex, grt = 0,
                    tbv = tbv, pht = pht, ebv = 0., F = 0.)
    serialize("$dir/$bar-f0-ped.ser", ped)
    nothing
end

"""
    function sum_effb9(dir, bar)
Calculate mean `tbv`, `F`, of each generation, and number
of fixed loci on chip and QTL in the end.
"""
function sum_effb9(dir, bar)
    for sel in ["-gs", "pht", "ped"]
        ped = deserialize("$dir/$bar-$sel+ped.ser")
        lmp = deserialize("$dir/$bar-map.ser")
        sp = combine(groupby(ped, :grt), :tbv => mean => :tbv, :F => mean => :F)
        open("$dir/mbvf.bin", "a") do io
            write(io, sp.tbv[2:end])
            write(io, sp.F[2:end])
        end
        bgn = findfirst(x -> x == ped.grt[end], ped.grt)
        hdr = readhdr("$dir/$bar-$sel.xy")
        nlc, nhp = hdr.m, hdr.n
        agt = Mmap.mmap("$dir/$bar-$sel.xy", Matrix{UInt16}, (nlc, nhp), 24)
        lgt = view(agt, :, bgn:nhp) # last generation
        snp = isodd.(lgt)
        fixed = zeros(Int, 3)
        for loc in eachrow(snp[lmp.chip, :])
            length(unique(loc)) == 1 && (fixed[1] += 1)
        end
        for loc in eachrow(snp[lmp.qtl, :])
            length(unique(loc)) == 1 && (fixed[2] += 1)
        end
        for loci in eachrow(snp)
            length(unique(loci)) == 1 && (fixed[3] += 1)
        end
        open("$dir/fixed.bin", "a") do io
            write(io, fixed)
        end
    end
end

function clean_effb9(rst, dir, bar)
    rm.(glob("$rst/base/*"), recursive = true, force = true)
    rm.(glob("$dir/$bar*"), recursive = true, force = true)
end

function stats_effb9(dir, ngrt)
    nf = filesize("$dir/mbvf.bin") ÷ 8
    ni = filesize("$dir/fixed.bin") ÷ 8
    mbvf = reshape(Mmap.mmap("$dir/mbvf.bin", Vector{Float64}, nf), ngrt, :)
    fixd = reshape(Mmap.mmap("$dir/fixed.bin", Vector{Int}, ni), 3, :)
    mbv = zeros(ngrt, 3)
    sbv = zeros(ngrt, 3)
    mf = zeros(ngrt, 3)
    sf = zeros(ngrt, 3)
    mfx = zeros(3, 3)
    sfx = zeros(3, 3)
    for i in 1:3
        mbv[:, i] = mean(mbvf[:, (2i-1):6:end], dims = 2)
        sbv[:, i] =  std(mbvf[:, (2i-1):6:end], dims = 2)
        mf[:, i] = mean(mbvf[:, 2i:6:end], dims = 2)
        sf[:, i] =  std(mbvf[:, 2i:6:end], dims = 2)
        mfx[i, :] = mean(fixd[:, i:3:end], dims = 2)
        sfx[i, :] =  std(fixd[:, i:3:end], dims = 2)
    end

    #=
    p1 = plot(mbv, label=["GS" "Phenotype" "Pedigree"])
    p2 = plot(mf, label=["GS" "Phenotype" "Pedigree"])
    p3 = plot(mbv ./ mf, label=["GS" "Phenotype" "Pedigree"])
    p4 = plot(["Marker", "QTL", "SNP"], mfx, label=["GS" "Phenotype" "Pedigree"], legend=:right)
    plot(p1, p2, p3, p4, dpi=300)
    =#

    mbv, sbv, mf, sf, mfx, sfx
end
