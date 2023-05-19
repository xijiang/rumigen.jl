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
    ped = initPedigree("$dir/$bar-f0.xy", lmp, σₑ)
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
        open("$dir/effb9.bin", "a") do io
            write(io, sp.tbv[2:end])
            write(io, sp.F[2:end])
            ideal, plst, nlst = idealPop("$dir/$bar-$sel.xy", ped.grt, lmp)
            write(io, ideal[2:end])
            write(io, plst[2:end])
            write(io, nlst[2:end])
        end
    end
end

function stats_effb9(file, ngrt)
    rst = Float64[]
    open(file, "r") do io
        a = zeros(3ngrt)
        b = zeros(Int, 2ngrt)
        while !eof(io)
            read!(io, a)
            read!(io, b)
            append!(rst, a)
            append!(rst, b)
        end
    end
    rst = reshape(rst, ngrt, :)
    return rst

    #=
    dg = rst[:, 1:5:end]
    df = rst[:, 2:5:end]
    bp = rst[:, 3:5:end]  # best possible population
    dq = rst[:, 4:5:end]  # positive QTL lost
    for i in 1:3:size(dg, 2)
        dg[:, i:i+2] .-= dg[:, i+1]
        df[:, i:i+2] .-= df[:, i+1]
        bp[:, i:i+2] .-= bp[:, i+1]
        dq[:, i:i+2] .-= dq[:, i+1]
    end
    mdg = [mean(dg[:, 1:3:end], dims = 2) mean(dg[:, 3:3:end], dims = 2)]
    mdf = [mean(df[:, 1:3:end], dims = 2) mean(df[:, 3:3:end], dims = 2)]
    mbp = [mean(bp[:, 1:3:end], dims = 2) mean(bp[:, 3:3:end], dims = 2)]
    mdq = [mean(dq[:, 1:3:end], dims = 2) mean(dq[:, 3:3:end], dims = 2)]

    sdg = [std(dg[:, 1:3:end], dims = 2) std(dg[:, 3:3:end], dims = 2)] ./ sqrt(200)
    sdf = [std(df[:, 1:3:end], dims = 2) std(df[:, 3:3:end], dims = 2)] ./ sqrt(200)
    sbp = [std(bp[:, 1:3:end], dims = 2) std(bp[:, 3:3:end], dims = 2)] ./ sqrt(200)
    sdq = [std(dq[:, 1:3:end], dims = 2) std(dq[:, 3:3:end], dims = 2)] ./ sqrt(200)

    gap = 50
    plot( mean(rst[:,  1:15:end], dims=2), label="GS", xlabel="Generation", ylabel=L"\Delta G", dpi = 300, ribbon = sdg[:, 1], color=1, fillalpha = 0.2)
    plot!(mean(rst[:,  3:15:end], dims=2) .-gap, label="max(GS) - $gap", ls=:dash, color=1)
    plot!(mean(rst[:,  6:15:end], dims=2), label="Phenotype", color=2)
    plot!(mean(rst[:,  8:15:end], dims=2) .-gap, label="max(phenotype) - $gap", ls=:dash, color=2, ribbon=std(rst[:, 8:15:end], dims=2) / sqrt(200), fillalpha = 0.2)
    plot!(mean(rst[:, 11:15:end], dims=2), label="Pedigree", color=4, ribbon = sdg[:, 2], fillalpha = 0.2)
    plot!(mean(rst[:, 13:15:end], dims=2) .-gap, label="max(pedigree) - $gap", ls=:dash, color=4, leg=:left)
    savefig("effb9-dg.png") # yerror can be used instead of ribbon


    plot( mean(rst[:,  2:15:end], dims=2), label="GS", xlabel="Generation", ylabel="Mean inbreeding")
    plot!(mean(rst[:,  7:15:end], dims=2), label="Phenotype")
    plot!(mean(rst[:, 12:15:end], dims=2), label="Pedigree", dpi = 300)
    savefig("effb9-inb.png")

    plot( mean(rst[:,  2:15:end], dims=2), mean( rst[:,  1:15:end], dims=2), label = "GS", xlabel="Mean inbreeding", ylabel="Mean " * L"\Delta G")
    plot!(mean(rst[:,  7:15:end], dims=2), mean( rst[:,  6:15:end], dims=2), label = "Phenotype")
    plot!(mean(rst[:, 12:15:end], dims=2), mean( rst[:, 11:15:end], dims=2), label = "Pedigree", dpi = 300)
    savefig("effb9-inb-dg.png")
    =#
end
