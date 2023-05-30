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
            append!(rst, b)
        end
    end
    rst = reshape(rst, ngrt, :) # I ignore fixed QTL numbers here

    mvggs = mean(rst[:,  2:18:end], dims=2)
    mvgpt = mean(rst[:,  8:18:end], dims=2)
    mvgpd = mean(rst[:, 14:18:end], dims=2)
    mfgs = mean(rst[:,  3:18:end], dims=2)
    mfpt = mean(rst[:,  9:18:end], dims=2)
    mfpd = mean(rst[:, 15:18:end], dims=2)

    plot(mfgs, mvggs, label="GS", xlabel="Inbreeding", ylabel = L"\Delta\sigma_a^2", dpi = 300)
    plot!(mfpt, mvgpt, label="Phenotype")
    plot!(mfpd, mvgpd, label="Pedigree")

    mdggs = mean(rst[:, 1:18:end], dims=2)
    mdgpt = mean(rst[:, 7:18:end], dims=2)
    mdgpd = mean(rst[:,13:18:end], dims=2)

    plot( 1 .- mvggs, mdggs, label="GS", xlabel=L"1 - \sigma_a^2", ylabel = L"\Delta G", dpi = 300)
    plot!(1 .- mvgpt, mdgpt, label="Phenotype")
    plot!(1 .- mvgpd, mdgpd, label="Pedigree")
    savefig("effb9-gvsvg.png")
    return
    #=
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
