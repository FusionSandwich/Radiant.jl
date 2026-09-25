using Radiant
using DelimitedFiles
using Printf

length(ARGS) == 1 || error("usage: julia compare_bethe_pstar.jl <qualification-directory>")
root = ARGS[1]
me = 0.51099895069

function interp(table, E)
    grid = vec(table[:,1])
    i = searchsortedlast(grid,E)
    i == length(grid) && return table[end,2]
    w = (E-grid[i])/(grid[i+1]-grid[i])
    return (1-w)*table[i,2]+w*table[i+1,2]
end

open(joinpath(root,"bethe_vs_pstar_stopping.csv"),"w") do io
    println(io,"energy_MeV,material,pstar_electronic_MeV_cm2_g,radiant_bethe_MeV_cm2_g,ratio_radiant_to_pstar")
    for (id,Z,density) in (("Al",13,2.699),("Cu",29,8.96))
        table = readdlm(joinpath(root,"pstar_$(id).csv"),',',Float64;skipstart=1)
        for E in (1.0,1.2,2.0,5.0,10.0,100.0,499.99)
            nist = interp(table,E)
            bethe = Radiant.bethe([Int64(Z)],[1.0],density,E/me,Proton())*me/density
            @printf(io,"%.5f,%s,%.12g,%.12g,%.12g\n",E,id,nist,bethe,bethe/nist)
        end
    end
end
