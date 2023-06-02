# echo Ideal ID number of SD from population mean | md5sum => c791a754bef3dbb2f081284ab66aac0e

function theoretical_c791a()
    qtls = 100:100:15000
    
    @info "Normal(0, 1)"
    d = Normal(0, 1)
    ideal_normal, sdi_normal = Float64[], Float64[]
    for nqtl in qtls
        a, b, _ = idealID(nqtl, d)
        push!(ideal_normal, a)
        push!(sdi_normal, b)
    end

    @info "Uniform(-1, 1)"
    d = Uniform(-1, 1)
    ideal_uniform, sdi_uniform = Float64[], Float64[]
    for nqtl in qtls
        a, b, _ = idealID(nqtl, d)
        push!(ideal_uniform, a)
        push!(sdi_uniform, b)
    end

    @info "Gamma(1, 1)"
    d = Gamma(1, 1)
    ideal_gamma, sdi_gamma = Float64[], Float64[]
    for nqtl in qtls
        a, b, _ = idealID(nqtl, d)
        push!(ideal_gamma, a)
        push!(sdi_gamma, b)
    end

    @info "Laplace(0, 1)"
    d = Laplace(0, 1)
    ideal_laplace, sdi_laplace = Float64[], Float64[]
    for nqtl in qtls
        a, b, _ = idealID(nqtl, d)
        push!(ideal_laplace, a)
        push!(sdi_laplace, b)
    end

    writexy("c791a.xy", [ideal_normal sdi_normal ideal_uniform sdi_uniform ideal_gamma sdi_gamma ideal_laplace sdi_laplace])

    #=
    @info "Plotting"
    qtls /= 100
    plot(qtls, ideal_normal, ribbon = sdi_normal, fillalpha = 0.3,
         label = "Normal(0, 1)",
         xlabel = L"N_{\mathrm{QTL}}/100",
         ylabel = L"\sigma_a", dpi=300)
    plot!(qtls, ideal_uniform, ribbon = sdi_uniform,
          fillalpha = 0.3, label = "Uniform(-1, 1)")
    plot!(qtls, ideal_gamma, ribbon = sdi_gamma,
          fillalpha = 0.3, label = "Gamma(1, 1)")
    plot!(qtls, ideal_laplace, ribbon = sdi_laplace,
          fillalpha = 0.3, label = "Laplace(0, 1)")
    savefig("c791a.svg")
    =#
end


function empirical_c791a()
    nsir, ndam, nsnp = 500, 500, 0
    foo = base_effb9(nsir, ndam, "rst")
    for nqtl in 100:100:15000
        @info "nqtl = $nqtl"
        bar = fdr_effb9("rst", "rst/ideal", foo, nsir, ndam, 0, nqtl)
    end
end


function xps_c791a()
end
