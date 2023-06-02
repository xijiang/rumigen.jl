function nodic_8ac835(ped)
    for i in 1:nrow(ped)
        kinship(ped, i, i)
    end
end

function mddic_8ac835(ped)
    dic = Relation()
    for i in 1:nrow(ped)
        kinship(ped, i, i, dic)
    end
end

function table_8ac835(ped)
    Amat(ped, m = nrow(ped))
end

function xps_8ac835()
    for ig in 1:15 # pedigree depth
        ped = randPed(200, ig)
        @benchmark nodic_8ac835($ped)
        @benchmark mddic_8ac835($ped)
    end
end

# turned out that the table method `Amat` is the fastest.