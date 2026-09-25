using Radiant
using DelimitedFiles
using Printf

length(ARGS) == 1 || error("usage: julia compare_bethe.jl <geant4_stopping.csv>")
table = readdlm(ARGS[1],',',Float64;skipstart=1)
grid = vec(table[:,1])
function interpolated(values,E)
    index = searchsortedlast(grid,E)
    index == length(grid) && return values[end]
    fraction = (E-grid[index])/(grid[index+1]-grid[index])
    return (1-fraction)*values[index]+fraction*values[index+1]
end
out = joinpath(dirname(ARGS[1]),"bethe_vs_geant4_stopping.csv")
open(out,"w") do io
    println(io,"energy_MeV,material,geant4_hIoni_MeV_cm,radiant_bethe_MeV_cm,ratio_radiant_to_geant4")
    for E in (1.0,2.0,10.0,100.0,499.99), (id,Z,density,column) in (("Al",13,2.699,2),("Cu",29,8.96,3))
        g4 = interpolated(vec(table[:,column]),E)
        # bethe() takes kinetic energy in electron-rest-mass units and returns cm^-1.
        # Multiplication by m_e c^2 produces MeV/cm for this direct model comparison.
        me = 0.51099895069
        radiant = Radiant.bethe([Int64(Z)],[1.0],density,E/me,Proton())*me
        @printf(io,"%.5f,%s,%.12g,%.12g,%.12g\n",E,id,g4,radiant,radiant/g4)
    end
end
println(read(out,String))
