using CSV, DataFrames,Plots,LaTeXStrings

# Define the path to your file
file_path = "all_vectors.txt"

# Define your custom column names
column_names = [
    :τ_grid, :util1, :util2, :util3, :util4, 
    :wavec, :prob1, :prob2, :yvec, :dvec1, :dvec2, 
    :cHvec1, :cHvec2, :cLvec1, :cLvec2
]

# Read the data into a DataFrame with specified column names
data = CSV.read(file_path, DataFrame, delim=',', header=column_names)

τn = size(data)[1]
xvec = ones(τn,1).*2 .- data.yvec

# display(Plots.plot(τ_grid,prob, linecolor=[:blue :red], linewidth=[1 1], label=[L"$G_1$" L"$G_2$"],legend_position=:topleft))
# println("Press Enter to Exit ...")
# readline()
# p8=Plots.plot(τ_grid,[rh_rp1],linecolor=:blue ,linewidth=2,ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,label=L"$i=1$",
#         legend=false,ylabel="Participation Constraint")

# p20=Plots.plot(τ_grid,[rh_rp2],linecolor=:red,linewidth=2, ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,label=L"$i=2$",
#         legend=false,legendfontsize=9,ylabel="Participation Constraint")

p7=Plots.plot(data.τ_grid,[data.util1 data.util2 data.util3 data.util4],linestyle=[:solid :dash :dot :dot],linecolor=[:blue :red :green :orange], linewidth=[2 2 2 2],ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,
        label=["Case 1" "Case 2" "Case 3" "Case 4"],legend_columns=2,legend_position=:bottomleft,legendfontsize=9)

p1=Plots.plot(data.τ_grid,[data.prob1 data.prob2],ylims=(0,1),yticks=collect(0.0:0.2:1.0),linewidth=[2 2],linecolor=[:blue :red], linestyle=[:solid :dash] ,ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,
        label=[L"$i=1$" L"$i=2$"],legend_position=:topleft,legendfontsize=9,ylabel=L"$Pr(BR|s=L)$")

p2=Plots.plot(data.τ_grid,[xvec data.yvec],linecolor=[:blue :red],linewidth=[2 2],ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,
        label=[L"$x$" L"$y$"],legend_columns=2,legend_position=:outertop,legendfontsize=9,ylabel=L"$x$,$y$",ylim=(0,2))

p3=Plots.plot(data.τ_grid,[data.dvec1  data.cLvec1  data.cHvec1],linestyle=[:solid :dash :solid],linecolor=[:red :blue :green],ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,
        label=[L"$d_1$" L"$c_{21L}$" L"$c_{21H}$"],legend_columns=3,legend_position=:outertop,legendfontsize=9,ylabel="Consumption Allocations",ylim=(-0.2,2))

p4=Plots.plot(data.τ_grid,[data.dvec2  data.cLvec2  data.cHvec2 ],linestyle=[:solid :dash :solid],linecolor=[:red :blue :green],ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,
        label=[L"$d_2$" L"$c_{22L}$" L"$c_{22H}$"],legend_columns=3,legend_position=:outertop,legendfontsize=9,ylabel="Consumption Allocations",ylim=(-0.2,2))

p5=Plots.plot(data.τ_grid,data.wavec,ytickfontsize=9,linewidth=2,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,label=L"$W^\ast$",
        legend_position=:bottomright,legendfontsize=9)

# p9=Plots.plot(τ_grid,[dvec[:,1] ./ dvec[1,1], cLvec[:,1] ./ cLvec[1,1], cHvec[:,1] ./ cHvec[1,1]],linestyle=[:solid :dash :solid],linecolor=[:red :blue :green],ytickfontsize=9,xlabel=L"$\tau$",
#         xguidefontsize=12,xtickfontsize=9,label=[L"$d_1$" L"$c_{21L}$" L"$c_{21H}$"],legend_columns=1,legend_position=:topleft,legendfontsize=9)

# p10=Plots.plot(τ_grid,[dvec[:,2] ./ dvec[1,2], cLvec[:,2] ./ cLvec[1,2], cHvec[:,2] ./ cHvec[1,2]],linestyle=[:solid :dash :solid],linecolor=[:red :blue :green],ytickfontsize=9,xlabel=L"$\tau$",
#         xguidefontsize=12,xtickfontsize=9,label=[L"$d_2$" L"$c_{22L}$" L"$c_{22H}$"],legend_columns=1,legend_position=:bottomleft,legendfontsize=9)

# p11=Plots.plot(τ_grid,[(params[3] .* (dvec[:,1] .+ params[2]) .^ (-params[1])) ./ ((1-params[3]) .* ((params[4] .* (cHvec[:,1] .+ params[2]) .^ (-params[1]) + params[5] .* (cLvec[:,1] .+ params[2]) .^ (-params[1]))))
#                 , (params[3] .* (dvec[:,2] .+ params[2]) .^ (-params[1])) ./ ((1-params[3]) .* ((params[4] .* (cHvec[:,2] .+ params[2]) .^ (-params[1]) + params[5] .* (cLvec[:,2] .+ params[2]) .^ (-params[1]))))]
#                 ,linestyle=[:solid :solid],linecolor=[:blue :red],ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,label=[L"MRS$^1_{12}$" L"MRS$^2_{12}$"],legend_columns=2,legend_position=:bottomright,legendfontsize=9)

# p12=Plots.plot(1 .+ τ_grid[2:end], [(dvec[2:end,1] .- dvec[1:end-1,1]) ./ (τ_grid[2:end] .- τ_grid[1:end-1]), params[4] .* ((cHvec[2:end,1] .- cHvec[1:end-1,1]) ./ (τ_grid[2:end] .- τ_grid[1:end-1])) .+ params[5] .* ((cLvec[2:end,1] .- cLvec[1:end-1,1]) ./ (τ_grid[2:end] .- τ_grid[1:end-1]))]
#                 , linestyle=[:solid :solid], linecolor=[:red :blue],linewidth=[2 2],ytickfontsize=9,xlabel=L"$\omega_1$",xguidefontsize=12,xtickfontsize=9,ylabel=L"$MPC_1$" , label=[L"$d_1$" L"$E[c_{21}]$"])
# p13=Plots.plot(1 .+ τ_grid[2:end], [(dvec[2:end,2] .- dvec[1:end-1,2]) ./ (τ_grid[2:end] .- τ_grid[1:end-1]), params[4] .* ((cHvec[2:end,2] .- cHvec[1:end-1,2]) ./ (τ_grid[2:end] .- τ_grid[1:end-1])) .+ params[5] .* ((cLvec[2:end,2] .- cLvec[1:end-1,2]) ./ (τ_grid[2:end] .- τ_grid[1:end-1]))]
#                 , linestyle=[:solid :solid], linecolor=[:red :blue],linewidth=[2 2],ytickfontsize=9,xlabel=L"$\omega_2$",xguidefontsize=12,xtickfontsize=9,ylabel=L"$MPC_2$" , label=[L"$d_1$" L"$E[c_{21}]$"],xflip=true)

# p14=Plots.plot(τ_grid,[params[1] .* dvec[:,1] ./ (dvec[:,1] .- params[2]), params[1] .* cHvec[:,1] ./ (cHvec[:,1] .- params[2]), params[1] .* cLvec[:,1] ./ (cLvec[:,1] .- params[2]) ]
#                 , linestyle=[:solid :dash :solid], linecolor=[:red :blue :green],linewidth=[2 2 2],ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,label=[L"$d_1$" L"$c_{21H}$" L"$c_{21L}$"]
#                 , legend_columns=3,legendposition=:outertop,legendfontsize=9,ylabel=L"$RRA$",ylim=(-0.2,6))

# p15=Plots.plot(τ_grid,[params[1] .* dvec[:,2] ./ (dvec[:,2] .- params[2]), params[1] .* cHvec[:,2] ./ (cHvec[:,2] .- params[2]), params[1] .* cLvec[:,2] ./ (cLvec[:,2] .- params[2]) ]
#                 , linestyle=[:solid :dash :solid], linecolor=[:red :blue :green],linewidth=[2 2 2],ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9,label=[L"$d_2$" L"$c_{22H}$" L"$c_{22L}$"]
#                 , legend_columns=3,legendposition=:outertop,legendfontsize=9,ylabel=L"$RRA$",ylim=(-0.2,6))

# p16=Plots.plot(τ_grid,[params[1] .* dvec[:,1] ./ (dvec[:,1] .+ params[2]), params[4] .* (params[1] .* cHvec[:,1] ./ (cHvec[:,1] .+ params[2])) .+ params[5] .* ( params[1] .* cLvec[:,1] ./ (cLvec[:,1] .+ params[2]))]
#                 , linestyle=[:solid :solid], linecolor=[:red :blue],linewidth=[2 2],ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9, label=[L"$d_1$" L"$E[c_{21}]$" ]
#                 , legend_columns=1,legendposition=:topright,legendfontsize=9,ylabel=L"$RRA$")

# p17=Plots.plot(τ_grid,[params[1] .* dvec[:,2] ./ (dvec[:,2] .+ params[2]), params[4] .* (params[1] .* cHvec[:,2] ./ (cHvec[:,2] .+ params[2])) .+ params[5] .* ( params[1] .* cLvec[:,2] ./ (cLvec[:,2] .+ params[2]))]
#                 , linestyle=[:solid :solid], linecolor=[:red :blue],linewidth=[2 2],ytickfontsize=9,xlabel=L"$\tau$",xguidefontsize=12,xtickfontsize=9, label=[L"$d_2$" L"$E[c_{22}]$" ]
#                 , legend_columns=1,legendposition=:right,legendfontsize=9,ylabel=L"$RRA$")

# p18=Plots.plot(1 .+ τ_grid[2:end], [dvec[2:end,1] .- dvec[1:end-1,1] , cHvec[2:end,1] .- cHvec[1:end-1,1] , cLvec[2:end,1] .- cLvec[1:end-1,1] ]
#                 , linestyle=[:solid :solid :solid], linecolor=[:red :blue :green],linewidth=[2 2 2],ytickfontsize=9,xlabel=L"$\omega_1$",xguidefontsize=12,xtickfontsize=9,ylabel=L"$\Delta c_i$" , label=[L"$d_1$" L"$c_{21H}$" L"$c_{21L}$"])

# p19=Plots.plot(1 .+ τ_grid[2:end] , [dvec[2:end,2] .- dvec[1:end-1,2] , cHvec[2:end,2] .- cHvec[1:end-1,2] , cLvec[2:end,2] .- cLvec[1:end-1,2] ]
#                 , linestyle=[:solid :solid :solid], linecolor=[:red :blue :green],linewidth=[2 2 2],ytickfontsize=9,xlabel=L"$\omega_2 $",xguidefontsize=12,xtickfontsize=9,ylabel=L"$\Delta c_i$" , label=[L"$d_2$" L"$c_{22H}$" L"$c_{22L}$"],xflip=true)

display(Plots.plot(p1))
savefig("p1_DRRA_main.png")
# println("Press Enter to Continue ...")
# readline()
display(Plots.plot(p2))
savefig("p2_DRRA_main.png")
# println("Press Enter to Continue ...")
# readline()
display(Plots.plot(p3))
savefig("p3_DRRA_main.png")
# println("Press Enter to Continue ...")
# readline()
display(Plots.plot(p4))
savefig("p4_DRRA_main.png")
# println("Press Enter to Continue ...")
# readline()
display(Plots.plot(p5))
savefig("p5_DRRA_main.png")
# println("Press Enter to Continue ...")
# readline()
display(Plots.plot(p7))
savefig("p7_DRRA_main.png")
# println("Press Enter to Exit ...")
# readline()
# display(Plots.plot(p8))
# savefig("p8_DRRA_main.png")
# println("Press Enter to Exit ...")
# readline()
# display(Plots.plot(p9))
# savefig("p9_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p10))
# savefig("p10_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p11))
# savefig("p11_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p12))
# savefig("p12_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p13))
# savefig("p13_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p14))
# savefig("p14_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p15))
# savefig("p15_DRRA_main.png")
# println("Press Enter to Exit ...")
# readline()
# display(Plots.plot(p16))
# savefig("p16_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p17))
# savefig("p17_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p18))
# savefig("p18_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p19))
# savefig("p19_DRRA_main.png")
# # println("Press Enter to Exit ...")
# # readline()
# display(Plots.plot(p20))
# savefig("p20_DRRA_main.png")
# println("Press Enter to Exit ...")
# readline()