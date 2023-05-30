"""
    function base_effb9(nsir, ndam, rst)
Create `2nid` haplotypes with MaCS, and convert them to the `xy` format.
"""
function base_effb9(nid, rst)
    tprintln("Generating a base population with MaCS")
    macs = make_macs(tdir = rst)
    tmp = cattle_genome(macs, nid, dir = "$rst/base")
    macs2xy(tmp)
end

"""
    function fdr_effb9(fdr, dir, foo, nsir, ndam, ppsz, nlc, nqtl)
Sample founders from the base population and mate them very randomly into F1.
"""
function fdr_effb9(fdr, dir, foo, ppsz, nlc, nqtl; d = Normal())
    @info "Sample founders from the base population created with MaCS"

    # sample haplotypes for the founder population
    nhp = 2ppsz
    bar = sampleFdr("$fdr/$foo-hap.xy", "$fdr/$foo-map.ser",
                    nhp, nlc = nlc, nqtl = nqtl, dir = dir)
    simQTL("$dir/$bar-fdr.xy", "$dir/$bar-map.ser", d = d)
    uniqSNP("$dir/$bar-fdr.xy")
    return bar
end

function pre_effb9(dir, bar, nsir, ndam, ngrt, σₑ)
    @info "Random mate the first $ngrt generations to form a common base population"
    lmp = deserialize("$dir/$bar-map.ser")

    # expand founders to a constant sized population
    ped = initPedigree("$dir/$bar-uhp.xy", lmp, σₑ)
    ped.grt .= -ngrt
    simpleSelection("$dir/$bar-uhp.xy", ped, lmp, nsir, ndam, ngrt, σₑ, random = true)
end

"""
    function sum_effb9(dir, bar)
Calculate mean `tbv`, `F`, of each generation, and number
of fixed loci on chip and QTL in the end.
"""
function sum_effb9(dir, bar)
    for sel in ["sgs", "spt", "spd"]
        ped = deserialize("$dir/$bar-$sel+ped.ser")
        lmp = deserialize("$dir/$bar-map.ser")
        sp = combine(groupby(ped, :grt), :tbv => mean => :mbv, :tbv => var => :vbv, :F => mean => :mF)
        open("$dir/effb9.bin", "a") do io
            write(io, sp.mbv[2:end])
            write(io, sp.vbv[2:end])  # as requested by SMS on 2023-05-21, by Theo
            write(io, sp.mF[2:end])
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
        a = zeros(4ngrt) # tbv, vbv, F, ideal
        b = zeros(Int, 2ngrt)
        while !eof(io)
            read!(io, a)
            read!(io, b)
            append!(rst, a)
            # append!(rst, b) # I ignore fixed QTL numbers here
        end
    end
    rst = reshape(rst, ngrt, :) # I ignore fixed QTL numbers here

    gs_mg  = mean(rst[:, 1:12:end], dims=2)
    gs_sg  = std( rst[:, 1:12:end], dims=2)
    gs_σₐ² = mean(rst[:, 2:12:end], dims=2)
    gs_mF  = mean(rst[:, 3:12:end], dims=2)
    gs_clg = mean(rst[:, 4:12:end], dims=2) # clg: ceiling

    pt_mg  = mean(rst[:, 5:12:end], dims=2)
    pt_sg  = std( rst[:, 5:12:end], dims=2)
    pt_σₐ² = mean(rst[:, 6:12:end], dims=2)
    pt_mF  = mean(rst[:, 7:12:end], dims=2)
    pt_clg = mean(rst[:, 8:12:end], dims=2)

    pd_mg  = mean(rst[:, 9:12:end], dims=2)
    pd_sg  = std( rst[:, 9:12:end], dims=2)
    pd_σₐ² = mean(rst[:,10:12:end], dims=2)
    pd_mF  = mean(rst[:,11:12:end], dims=2)
    pd_clg = mean(rst[:,12:12:end], dims=2)

    plot( -4:ngrt-5, gs_mg, label = "Genomic selection", xlabel = "Generation", ylabel = "Mean TBV", dpi = 300, legend = :left)
    plot!(-4:ngrt-5, pt_mg, label = "Phenotypic selection")
    plot!(-4:ngrt-5, pd_mg, label = "Pedigree selection")
    plot!(-4:ngrt-5, gs_clg .- 40, color = 1, label = "GS ceiling - 40")
    plot!(-4:ngrt-5, pt_clg .- 40, color = 2, label = "Phen. sel. ceiling - 40")
    plot!(-4:ngrt-5, pd_clg .- 40, color = 3, label = "Pedi. sel. ceiling - 40")
    savefig("effb9-mg.png")

    plot( gs_mg[6:end], gs_sg[6:end], label = "GS risks", ylabel = "STD of TBV", xlabel = "Mean TBV", dpi = 300)
    plot!(pt_mg[6:end], pt_sg[6:end], label = "Phe.S. risks")
    plot!(pd_mg[6:end], pd_sg[6:end], label = "Ped.S. risks")
    savefig("effb9-mg-sg.png")

    plot( -4:ngrt-5, gs_σₐ², label = "Genomic selection", xlabel = "Generation", ylabel = "Mean σₐ²", dpi = 300)
    plot!(-4:ngrt-5, pt_σₐ², label = "Phenotypic selection")
    plot!(-4:ngrt-5, pd_σₐ², label = "Pedigree selection")
    savefig("effb9-va.png")

    plot( -4:ngrt-5, gs_mF, label = "Genomic selection", xlabel = "Generation", ylabel = "Mean F", dpi = 300)
    plot!(-4:ngrt-5, pt_mF, label = "Phenotypic selection")
    plot!(-4:ngrt-5, pd_mF, label = "Pedigree selection")
    savefig("effb9-mf.png")

    plot( gs_mF, gs_mg, label = "Genomic selection", xlabel = "Mean F", ylabel = "Mean TBV", dpi = 300)
    plot!(pt_mF, pt_mg, label = "Phenotypic selection")
    plot!(pd_mF, pd_mg, label = "Pedigree selection")
    savefig("effb9-mf-mg.png")
end
