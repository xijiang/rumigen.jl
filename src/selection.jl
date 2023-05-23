"""
    function pedng(ped, sel::Symbol, nsir, ndam)
Given a pedigree `ped`, select `nsir` sires and `ndam` dams from 
the last generation on `ped` column `sel`. This function returen
a random pairs of `nsir` sires and `ndam` dams for form full sibship
families. One can use `repeat(pm, outer=nsib)` for a full next
generation pedigree after this function.
"""
function pedng(ped, sel::Symbol, nsir, ndam)
    pg = maximum(ped.grt)
    sel ∈ propertynames(ped) || error("Pedigree has no property $sel")
    pp = groupby(sort(ped[ped.grt .== pg, :], sel, rev = true), :sex)
    ma = pp[1].id[1:ndam] # sex 0 for females, 1 for males
    pa = pp[2].id[1:nsir]
    randomMate(pa, ma)
end

"""
    function simpleSelection(xy, ped, lmp, nsir, ndam, ngrt, σₑ; ebv = false, gs = false, random = false)
A simple selection strategy.
- The population is of the same size of the last generation in `ped`.
- Everybody has genotypes (in `xy`) and phenotypes (in `ped`).
- Select `nsir` sires and `ndam` dams from the last generation of `ped`.
- Mate them randomly into full sibship families of equal sizes.
- The first half of the sibs are males. The rest are females.
- Selection is on `:pht` by default.
    - If `random = true`, random select the parents.
    - If `ebv = true`, selection is on `:ebv`.
    - If `gs = true`, `:ebv` is of GEBV
    - By default, `:ebv` is of PEBV if to be calculated.
- `h²` is calculated as `σₐ = 1`.
"""
function simpleSelection(xy, ped, lmp, nsir, ndam, ngrt, σₑ; ebv = false, gs = false, random = false)
    hdr, h² = readhdr(xy), 1 / (1 + σₑ^2)
    mt, et, mj, nr, nc = xyhdr(hdr)
    (mt == 'F' && mj == 0 && nc == 2nrow(ped)) || error("$xy, or pedigree not right")
    lms, nfam = sumMap(lmp), max(nsir, ndam)
    nsib = sum(ped.grt .== ped.grt[end]) ÷ nfam
    # to generate the same number of animals in the last generation
    sibs = [ones(Int, nsib ÷ 2); zeros(Int, nsib - nsib ÷ 2)]
    @info "Selection on $xy, for $ngrt generations"

    for igrt in 1:ngrt
        print(" $igrt")
        agt = Mmap.mmap(xy, Matrix{et}, (nr, nc), 24) # ancestors
        ogt = zeros(et, nr, nfam * nsib * 2)
        pm = if random
            ped.F = rand(nrow(ped)) # borrow :F, which is overwritten in the end
            repeat(pedng(ped, :F, nsir, ndam), outer = nsib)
        elseif ebv
            giv = if gs  # GEBV
                grmiv(xy, lmp.chip)
            else   # PEBV
                A⁻¹(ped)
            end
            animalModel(ped, giv, h²) # default using :grt as fixed effect
            repeat(pedng(ped, :ebv, nsir, ndam), outer = nsib)
        else
            repeat(pedng(ped, :pht, nsir, ndam), outer = nsib)
        end
        drop(agt, ogt, pm, lms)
        appendxy(xy, ogt)
        tbv, pht = uhp2gp(ogt, lmp, σₑ)
        df = DataFrame(id = nrow(ped) + 1 : nrow(ped) + nsib * nfam,
                       pa = pm[:, 1], ma = pm[:, 2], 
                       sex = repeat(sibs, outer = nfam), 
                       grt = ped.grt[end] + 1,
                       tbv = tbv, 
                       pht = pht, ebv = 0., F = 0.)
        append!(ped, df)
        agt = nothing
        nc += nsib * nfam * 2
    end
    println()
    ped.F = inbreeding(xy, lmp.chip)
    serialize("$(xy[1:end-3])+ped.ser", ped)
end
