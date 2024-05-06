using Distributed, SharedArrays, DelimitedFiles, Parameters

# Add 30 worker processes
addprocs(6)

@everywhere begin
    using Parameters, SharedArrays
    include("./definitions.jl")
    # Function to handle the distributed computation
    function compute_for_each_yτ(y, τ,obj_function, params, cn)
        γ, ψ, λ, π_H, π_L, r, ω, RH, RL, β = params
        local max_value = -Inf
        local max_d_1, max_d_2, max_c21H, max_c22H, max_c21L, max_c22L, max_y = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, y

        d_1_range = LinRange(0.01 * y / λ , 0.99 * y / λ, cn)
        c21H_range = LinRange(0.01 * RH * (ω - y) / (1 - λ), 0.99 * RH * (ω - y) / (1 - λ), cn)

        for d_1 in d_1_range, c21H in c21H_range
            if obj_function == objective_function  # Adjust for case 1
                c21L_range = LinRange(0.01 * RL * (ω - y) / (1 - λ) , 0.99 * RL * (ω - y) / (1 - λ), cn)
                for d_1 in d_1_range, c21H in c21H_range, c21L in c21L_range
                    value = obj_function(d_1, (y/λ) - d_1, c21H, (RH * (ω-y)/(1-λ)) - c21H, c21L, (RL * (ω-y)/(1-λ)) - c21L, τ, y, params)
                    if value > max_value
                        max_value = value
                        max_d_1 = d_1
                        max_d_2 = (y/λ) - d_1
                        max_c21H = c21H
                        max_c22H = (RH * (ω-y)/(1-λ)) - c21H
                        max_c21L = c21L
                        max_c22L =  (RL * (ω-y)/(1-λ)) - c21L
                        max_y = y
                    end
                end
            else
                if obj_function == objective_function2                             
                    for d_1 in d_1_range, c21H in c21H_range
                        value = obj_function(d_1, (y/λ) - d_1, c21H, (RH * (ω-y)/(1-λ) - c21H), d_1*(r*(ω-y)+y)/(d_1+(y/λ) - d_1), ((y/λ) - d_1)*(r*(ω-y)+y)/(d_1+(y/λ) - d_1), τ, y, params)
                        if value > max_value
                            max_value = value
                            max_d_1 = d_1
                            max_d_2 = (y/λ) - d_1
                            max_c21H = c21H
                            max_c22H = (RH * (ω-y)/(1-λ)) - c21H
                            max_c21L = d_1*(r*(ω-y)+y)/(d_1+(y/λ) - d_1)
                            max_c22L = ((y/λ) - d_1)*(r*(ω-y)+y)/(d_1+(y/λ) - d_1)
                            max_y = y
                        end
                    end
                elseif obj_function == objective_function3
                    for d_1 in d_1_range, c21H in c21H_range
                        value = obj_function(d_1, (y/λ) - d_1, c21H, (RH * (ω-y)/(1-λ)) - c21H, λ*d_1*(r*(ω-y)+y)/(λ*d_1+(y/λ) - d_1), ((y/λ) - d_1)*(r*(ω-y)+y)/(λ*d_1+(y/λ) - d_1), τ, y, params)
                        if value > max_value
                            max_value = value
                            max_d_1 = d_1
                            max_d_2 = (y/λ) - d_1
                            max_c21H = c21H
                            max_c22H =(RH * (ω-y)/(1-λ)) - c21H
                            max_c21L = λ*d_1*(r*(ω-y)+y)/(λ*d_1+(y/λ) - d_1)
                            max_c22L = ((y/λ) - d_1)*(r*(ω-y)+y)/(λ*d_1+(y/λ) - d_1)
                            max_y = y
                        end
                    end
                elseif obj_function == objective_function4
                    for d_1 in d_1_range, c21H in c21H_range
                        value = obj_function(d_1, (y/λ) - d_1, c21H, (RH * (ω-y)/(1-λ)) - c21H, d_1*(r*(ω-y)+y)/(d_1+λ*((y/λ) - d_1)),λ*((y/λ) - d_1)*(r*(ω-y)+y)/(d_1+λ*((y/λ) - d_1)), τ, y,params)
                        if value > max_value
                            max_value = value
                            max_d_1 = d_1
                            max_d_2 = (y/λ) - d_1
                            max_c21H = c21H
                            max_c22H = (RH * (ω-y)/(1-λ)) - c21H
                            max_c21L = d_1*(r*(ω-y)+y)/(d_1+λ*((y/λ) - d_1))
                            max_c22L = λ*((y/λ) - d_1)*(r*(ω-y)+y)/(d_1+λ*((y/λ) - d_1))
                            max_y = y
                        end
                    end                                
                end
            end
        end
        return (τ, max_value, max_d_1, max_d_2, max_c21H, max_c22H, max_c21L, max_c22L, max_y)
    end
end

function optimize_function(obj_function, params, y_grid, τ_grid, cn)
    # Ensure params are available as a shared array
    shared_params = convert(SharedArray, params)
    
    results = pmap((yt) -> compute_for_each_yτ(yt[1], yt[2],obj_function,shared_params, cn), [(y, τ) for y in y_grid, τ in τ_grid])
    max_values = Dict([(result[1], result[2:end]) for result in results])

    τ_values = sort(collect(keys(max_values)))
    max_utility_values = [value == -Inf ? NaN : value[1] for value in values(max_values)]
    max_y_values = [max_values[τ][8] for τ in τ_values]

    max_d1_values = [max_values[τ][2] for τ in τ_values]
    max_c21H_values = [max_values[τ][4] for τ in τ_values]
    max_c21L_values = [max_values[τ][6] for τ in τ_values]
    

    max_d2_values = [max_values[τ][3] for τ in τ_values]
    max_c22H_values = [max_values[τ][5] for τ in τ_values]
    max_c22L_values = [max_values[τ][7] for τ in τ_values]

    return τ_values, max_utility_values, max_d1_values, max_c21H_values, max_c21L_values, max_d2_values, max_c22H_values, max_c22L_values, max_y_values
end

# Parameters and grids
params = [3.0, 0.40, 0.15, 0.80, 0.20, 0.80, 2.0, 1.5, 1.065, 0.94]
yn = 10
τn = 10
cn = 40

y_grid = collect(range(0.001, stop=0.50, length=yn))
τ_grid = collect(range(0.30, stop=0.90, length=τn))

util=zeros(τn,4)

# Call the optimization functions for different cases
@time begin
    results_case1 = optimize_function(objective_function, params, y_grid, τ_grid, cn);
    results_case2 = optimize_function(objective_function2, params, y_grid, τ_grid, cn);
    results_case3 = optimize_function(objective_function3, params, y_grid, τ_grid, cn);
    results_case4 = optimize_function(objective_function4, params, y_grid, τ_grid, cn);
    util[:,1]=results_case1[2][:]
    util[:,2]=results_case2[2][:]
    util[:,3]=results_case3[2][:]
    util[:,4]=results_case4[2][:]
    util1=replace!(util,NaN=>-Inf)
    prob=zeros(τn,2)
    dvec=zeros(τn,2)
    cLvec=zeros(τn,2)
    cHvec=zeros(τn,2)
    yvec=zeros(τn,1)
    wavec = zeros(τn,1)
    for i in eachindex(τ_grid)
        c1 = findmax(util1[i,:])
        if c1[2]==1
            prob[i,1:2].=0
            dvec[i,1]=results_case1[3][i]
            dvec[i,2]=results_case1[6][i]
            cHvec[i,1]=results_case1[4][i]
            cHvec[i,2]=results_case1[7][i]
            cLvec[i,1]=results_case1[5][i]
            cLvec[i,2]=results_case1[8][i]
            yvec[i,1]=results_case1[9][i]
            wavec[i,1]=results_case1[2][i]
        elseif c1[2]==2
            prob[i,1:2].=1
            dvec[i,1]=results_case2[3][i]
            dvec[i,2]=results_case2[6][i]
            cHvec[i,1]=results_case2[4][i]
            cHvec[i,2]=results_case2[7][i]
            cLvec[i,1]=results_case2[5][i]
            cLvec[i,2]=results_case2[8][i]
            yvec[i,1]=results_case2[9][i]
            wavec[i,1]=results_case2[2][i]
        elseif c1[2]==3
            prob[i,1:2].=[0,1]
            dvec[i,1]=results_case3[3][i]
            dvec[i,2]=results_case3[6][i]
            cHvec[i,1]=results_case3[4][i]
            cHvec[i,2]=results_case3[7][i]
            cLvec[i,1]=results_case3[5][i]
            cLvec[i,2]=results_case3[8][i]
            yvec[i,1]=results_case3[9][i]
            wavec[i,1]=results_case3[2][i]
        elseif c1[2]==4
            prob[i,1:2].=[1,0]
            dvec[i,1]=results_case4[3][i]
            dvec[i,2]=results_case4[6][i]
            cHvec[i,1]=results_case4[4][i]
            cHvec[i,2]=results_case4[7][i]
            cLvec[i,1]=results_case4[5][i]
            cLvec[i,2]=results_case4[8][i]
            yvec[i,1]=results_case4[9][i]
            wavec[i,1]=results_case4[2][i]
        end
    end
end

# display(Plots.plot(τ_grid,util))
xvec = ones(τn,1).*2 .- yvec

# Stack vectors horizontally (ensure all vectors are the same length or transposed correctly)
all_data = hcat(τ_grid, util, prob, xvec, yvec, dvec, cLvec, cHvec, wavec)

# Save the combined array to a file
writedlm("all_vectors.txt", all_data, ',')
rmprocs(workers())