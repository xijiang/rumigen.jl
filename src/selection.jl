"""
    function pedng(ped, sel::Symbol, nsir, ndam)
Given a pedigree `ped`, select `nsir` sires and `ndam` dams from 
the last generation on `ped` column `sel`. This function returen
a random pairs of `nsir` sires and `ndam` dams for form full sibship
families. One can use `repeat(pm, outer=nsib)` for a full next
generation pedigree alrer this function.

Note: this function is deprecated. Use `prt4ng` instead.
"""
function pedng(ped, sel::Symbol, nsir, ndam)
    pg = maximum(ped.grt)
    sel ∈ propertynames(ped) || error("Pedigree has no property $sel")
    pp = groupby(sort(ped[ped.grt .== pg, :], sel, rev = true), :sex)
    ma = pp[1].id[1:ndam] # sex 0 for females, 1 for males
    pa = pp[2].id[1:nsir]
    randomMate(pa, ma)
end

function prt4ng(ped, nsir, ndam)
    [:ebv, :sex, :grt, :id] ⊆ propertynames(ped) || error("Not a proper pedigree")
    pool = filter(row -> row.grt == ped.grt[end], ped)
    sort!(pool, :ebv, rev = true)
    pp = groupby(pool, :sex)
    ma = pp[1].id[1:ndam]
    pa = pp[2].id[1:nsir]
    randomMate(pa, ma)
end

"""
    function simpleSelection(xy, ped, lmp, nsir, ndam, ngrt, σₑ; ebv = false, gs = false, random = false, mp = true)
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
- `h²` is calculated with `σₐ = 1`, and given `σₑ`.
- When `mp = false`, males' phenotypes are set to `missing`.

## Update 2023-06-13
- only update EBV of the last generation
  - Or, the accuracy of EBV will be overestimated for the earlier generations.
"""
function simpleSelection(xy, ped, lmp, nsir, ndam, ngrt, σₑ;
                         ebv = false,
                         gs = false,
                         random = false,
                         mp = true,  # male phenotypes
                         )
    hdr, h² = readhdr(xy), 1 / (1 + σₑ^2)
    mt, et, mj, nr, nc = xyhdr(hdr)
    (mt == 'F' && mj == 0 && nc == 2size(ped, 1)) || error("$xy, or pedigree not right")
    lms, nfam = sumMap(lmp), max(nsir, ndam)
    nsib = sum(ped.grt .== ped.grt[end]) ÷ nfam
    # to generate the same number of animals in the last generation
    sibs = [ones(Int, nsib ÷ 2); zeros(Int, nsib - nsib ÷ 2)]
    sex = repeat(sibs, outer = nfam)
    ser = splitext(xy)[1] * "-mid.ser" # to store mid-results for A⁻¹
    
    @info "Selection on $xy, for $ngrt generations"
    for igrt in 1:ngrt
        print(" $igrt")
        agt = Mmap.mmap(xy, Matrix{et}, (nr, nc), 24) # ancestors
        ogt = zeros(et, nr, nfam * nsib * 2)
        oebv = ped.ebv[ped.grt .< ped.grt[end]]
        if random # update EBV
            ped.ebv = rand(size(ped, 1))
        elseif ebv
            giv = if gs  # GEBV
                grmiv(xy, lmp.chip)
            else   # PEBV
                A⁻¹(ped, ser)
            end
            animalModel(ped, giv, h²) # default using :grt as fixed effect
        else
            ped.ebv = disallowmissing(ped.pht)
        end
        ped.ebv[ped.grt .< ped.grt[end]] = oebv # restore previously calculated EBV
        pm = repeat(prt4ng(ped, nsir, ndam), outer = nsib)
        drop(agt, ogt, pm, lms)
        appendxy!(xy, ogt)
        tbv, pht = uhp2gp(ogt, lmp, σₑ)
        df = DataFrame(id = size(ped, 1) + 1 : size(ped, 1) + nsib * nfam,
                       pa = pm[:, 1], ma = pm[:, 2], 
                       sex = sex, 
                       grt = ped.grt[end] + 1,
                       tbv = tbv, 
                       pht = pht, ebv = 0., F = 0., Fr = 0.)
        append!(ped, df)

        mp || (ped.pht[ped.sex .== 1] .= missing)
        agt = nothing
        nc += nsib * nfam * 2
    end
    println()
    ped.F = inbreeding(xy, lmp.chip)
    ped.Fr = inbreeding(xy, lmp.ref) # inbreeding by reference loci
    serialize("$(xy[1:end-3])+ped.ser", ped)
end

"""
    function optSelection(xy, ped, lmp, ngrt, σₑ, dF; op=1, k₀=0.)
Note, dF is often not big enough, such that to have a solution for `c`.
- `op = 1` for optimum contribution selection with pedigree
- `op = 2` for optimum contribution selection with genomic selection,
  constrained by `A`.
- `op = 3` for optimum contribution selection with genomic selection,
  constrained by `G`.
"""
function optSelection(xy, ped, lmp, ngrt, σₑ, dF; op=1, k₀=0.)
    hdr, h² = readhdr(xy), 1 / (1 + σₑ^2)
    mt, et, mj, nr, nc = xyhdr(hdr)
    (mt == 'F' && mj == 0 && nc == 2size(ped, 1)) || error("$xy, or pedigree not right")
    lms = sumMap(lmp)
    nid = sum(ped.grt .== ped.grt[end])
    #t = length(unique(ped.grt))
    
    @info "Selection on $(basename(xy)), for $ngrt generations"
    for igrt ∈ 1:ngrt
        print(" $igrt")
        agt = Mmap.mmap(xy, Matrix{et}, (nr, nc), 24) # ancestors
        ogt = zeros(et, nr, nid * 2)
        oebv = ped.ebv[ped.grt .< ped.grt[end]]
        giv = A₂₂ = nothing
        pool = findall(ped.grt .== ped.grt[end]) # ID of current generation
        if op == 1
            A = Amat(ped)
            A₂₂ = copy(A[pool, pool])
            giv = inv(A)
            A = nothing
        elseif op == 2
            G = grm(xy, lmp.chip)
            giv = inv(G)
            G = nothing
            A = Amat(ped)
            A₂₂ = copy(A[pool, pool])
            A = nothing
        elseif op == 3
            G = grm(xy, lmp.chip)
            giv = inv(G)
            A₂₂ = copy(G[pool, pool])
            G = nothing
        end

        animalModel(ped, giv, h²) # default using :grt as fixed effect
        ped.ebv[ped.grt .< ped.grt[end]] = oebv # restore previously calculated EBV
        K = 2(1 - (1 - k₀) * (1 - dF) ^ (igrt+1))
        c = myopt(DataFrame(ebv=ped.ebv[pool], sex=ped.sex[pool]), A₂₂, K, silent=true)
        pm = randomMate(DataFrame(sex=ped.sex[pool], c=c), nid) .+ (size(ped, 1) - nid)
        drop(agt, ogt, pm, lms)
        appendxy!(xy, ogt)
        tbv, pht = uhp2gp(ogt, lmp, σₑ)
        df = DataFrame( id = size(ped, 1) + 1: size(ped, 1) + size(pm, 1),
                        pa = pm[:, 1], ma = pm[:, 2],
                        sex = rand(0:1, nid),
                        grt = ped.grt[end] + 1,
                        tbv = tbv,
                        pht = pht, ebv = 0., F = 0., Fr=0.
                      )
        append!(ped, df)
        ped.pht[ped.sex .==1] .== missing
        agt = nothing
        nc += nid * 2
    end
    println()
    ped.F = inbreeding(xy, lmp.chip)
    ped.Fr = inbreeding(xy, lmp.ref) # inbreeding by reference loci
    serialize("$(xy[1:end-3])+ped.ser", ped)
end
