
using Distributed, DelimitedFiles, BenchmarkTools
using LaTeXStrings,  Parameters
# using GR
# gr(label=false)

# Set the number of worker processes
num_workers = 1  # You can adjust this based on your available cores

# Add processes
addprocs(num_workers)


#region: Utility functions
@everywhere begin

    # Define the utility function
    function utility_function(x, γ, ψ)
        return ((x-ψ)^(1-γ))/(1-γ)
        # return ψ-γ*log(x)
    end
    
    # Define the objective function for case 1
    function objective_function(d_1, d_2, c21H, c22H, c21L, c22L, τ, y,params)
        γ, ψ, λ, π_H, π_L,r,ω,RH,RL,β = params 
        if c21L >= d_1 && c22L >= d_2 && c21H > c21L && c22H > c22L && min(d_1 , d_2 ) != ψ && min(d_1 , d_2 ) != 0 &&
            λ * utility_function(d_1, γ, ψ) + (1 - λ) * β * (π_H * utility_function(c21H, γ, ψ) + π_L * utility_function(c21L, γ, ψ)) >= utility_function(1 + τ, γ, ψ) &&
            λ * utility_function(d_2, γ, ψ) + (1 - λ) * β * (π_H * utility_function(c22H, γ, ψ) + π_L * utility_function(c22L, γ, ψ)) >= utility_function(1 - τ , γ, ψ) &&
            λ*(d_1+d_2)==y && (1-λ)*(c21H+c22H)==RH*(ω-y) && (1-λ)*(c21L+c22L)==RL*(ω-y) 
            term_H = λ * (utility_function(d_1, γ, ψ) + utility_function(d_2, γ, ψ)) +
                    (1 - λ) * β * (utility_function(c21H, γ, ψ) + utility_function(c22H, γ, ψ))
            term_L = λ * (utility_function(d_1, γ, ψ) + utility_function(d_2, γ, ψ)) +
                    (1 - λ) * β * (utility_function(c21L, γ, ψ) + utility_function(c22L, γ, ψ))
            return π_H * term_H + π_L * term_L
        else
            return -Inf  # Return a large negative value if conditions are not satisfied
        end
    end

    # Define the objective function
    function objective_function2(d_1, d_2, c21H, c22H, c21L, c22L, τ, y, params)
        γ, ψ, λ, π_H, π_L,r,ω,RH,RL,β = params 
        if  d_1 > c21L && d_2 > c22L && c21H >= d_1 && c22H >=d_2 &&  min(c21L,c22L) != ψ &&  min(c21L,c22L) != 0 &&
            r*(ω-y)+y <= d_1 + d_2 &&
            π_H*(λ * utility_function(d_1, γ, ψ) + (1 - λ) * β * utility_function(c21H, γ, ψ)) + π_L * β * utility_function(c21L, γ, ψ) >= utility_function(1 + τ, γ, ψ) &&
            π_H*(λ * utility_function(d_2, γ, ψ) + (1 - λ) * β * utility_function(c22H, γ, ψ)) + π_L * β * utility_function(c22L, γ, ψ) >= utility_function(1 - τ , γ, ψ) &&
            λ*(d_1+d_2)==y && (1-λ)*(c21H+c22H)==RH*(ω-y) 
            term_H = λ * (utility_function(d_1, γ, ψ) + utility_function(d_2, γ, ψ)) +
                    (1 - λ) * β * (utility_function(c21H, γ, ψ) + utility_function(c22H, γ, ψ))
            term_L = β *(utility_function(c21L, γ, ψ) + utility_function(c22L, γ, ψ))
            return π_H * term_H + π_L * term_L
        else
            return -Inf  # Return a large negative value if conditions are not satisfied
        end
    end

    # Define the objective function
    function objective_function4(d_1, d_2, c21H, c22H, c21L, c22L, τ, y, params)
        γ, ψ, λ, π_H, π_L,r,ω,RH,RL,β = params 
        if c21L >= d_1 && c21H > c21L && c22H >= d_2 && d_2 > c22L && min(d_1,c22L) != ψ && min(d_1,c22L) != 0 &&
            r*(ω-y)+y <= λ*d_1 +d_2 &&
            λ * utility_function(d_1, γ, ψ) + (1 - λ) * β * π_H * utility_function(c21H, γ, ψ) + π_L * β * (1 - λ) * utility_function(c21L, γ, ψ) >= utility_function(1 + τ,γ,ψ) &&
            π_H*(λ * utility_function(d_2, γ, ψ) + (1 - λ) * β * utility_function(c22H, γ, ψ)) + π_L * β * utility_function(c22L, γ, ψ) >= utility_function(1 - τ ,γ,ψ) &&
            λ*(d_1+d_2)==y && (1-λ)*(c21H+c22H)==RH*(ω-y)  
            term_H = λ * (utility_function(d_1, γ, ψ) + utility_function(d_2, γ, ψ)) +
                    (1 - λ) * β * (utility_function(c21H, γ, ψ) + utility_function(c22H, γ, ψ))
            term_L = λ*utility_function(d_1,γ,ψ)+(1-λ) * β * (utility_function(c21L, γ, ψ)) + utility_function(c22L, γ, ψ)
            return π_H * term_H + π_L * term_L 
        else
            return -Inf  # Return a large negative value if conditions are not satisfied
        end
    end
    # Define the objective function
    function objective_function3(d_1, d_2, c21H, c22H, c21L, c22L, τ, y, params)
        γ, ψ, λ, π_H, π_L,r,ω,RH,RL,β = params 
        if  c22L >= d_2 && c22H > c22L && c21H >= d_1 && min(d_2,c21L) != ψ && min(d_2,c21L) != 0 && 
            r*(ω-y)+y <= d_1 + λ*d_2 &&  λ * utility_function(d_2, γ, ψ) + (1 - λ) * π_H * β * utility_function(c22H, γ, ψ) + π_L * β * (1 - λ) * utility_function(c22L, γ, ψ) >= utility_function(1 - τ ,γ,ψ) &&
            π_H*(λ * utility_function(d_1, γ, ψ) + (1 - λ) * β * utility_function(c21H, γ, ψ)) + π_L * β * utility_function(c21L, γ, ψ) >= utility_function(1 + τ,γ,ψ) &&
            λ*(d_1+d_2)==y && (1-λ)*(c21H+c22H)==RH*(ω - y) 
            term_H = λ * (utility_function(d_1, γ, ψ) + utility_function(d_2, γ, ψ)) +
                    (1 - λ) * β * (utility_function(c21H, γ, ψ) + utility_function(c22H, γ, ψ))
            term_L = λ*utility_function(d_2,γ,ψ)+(1-λ) * β *(utility_function(c22L, γ, ψ)) +  β * utility_function(c21L, γ, ψ)
            return π_H * term_H + π_L * term_L
        else
            return -Inf  # Return a large negative value if conditions are not satisfied
        end
    end

    #endregion

    #region: Optimized function
    function optimize_function(obj_function, params, y_grid, τ_grid, cn)
        γ, ψ, λ, π_H, π_L,r,ω,RH,RL,β = params 
        max_values = Dict()

        @sync begin
            worker_results = Vector{Future}()
            for i in 1:num_workers
                push!(worker_results, @spawn begin
                    local results = Dict()
                    for τ in τ_grid
                        max_value = -Inf
                        max_d_1 = 0.0
                        max_d_2 = 0.0
                        max_c21H = 0.0
                        max_c22H = 0.0
                        max_c21L = 0.0
                        max_c22L = 0.0
                        max_y = 0.0
                        for y in y_grid
                            d_1_range = LinRange(0.01 * y / λ , 0.99 * y / λ, cn)
                            # d_2_range = LinRange(max(ψ, 0.01*(y/λ)),max(ψ,0.99*y/λ), cn)
                            #d_2_range = LinRange(0.1, y/λ, cn)
                            c21H_range = LinRange(0.01 * RH * (ω - y) / (1 - λ) , 0.99 * RH * (ω - y) / (1 - λ), cn)
                            # c22H_range = LinRange(max(ψ,0.01*(RH * (ω - y) / (1 - λ))),max(ψ,0.99*RH * (ω - y) / ((1 - λ))), cn)
                            #c22H_range = LinRange(0.1, RH * (ω - y) / ((1 - λ)), cn)
                            if obj_function == objective_function  # Adjust for case 1
                                c21L_range = LinRange(0.01 * RL * (ω - y) / (1 - λ) , 0.99 * RL * (ω - y) / (1 - λ), cn)
                                # c22L_range = LinRange(max(ψ,0.01*(RH * (ω - y) / (1 - λ))),max(ψ,0.99*RH * (ω - y) / ((1 - λ))), cn)
                                # c22L_range = LinRange(0.1,RL * (ω - y) / ((1 - λ)), cn)
                                for d_1 in d_1_range, c21H in c21H_range, c21L in c21L_range
                                    value = obj_function(d_1, (y/λ) - d_1, c21H, (RH * (ω-y)/(1-λ)) - c21H, c21L, (RL * (ω-y)/(1-λ)) - c21L, τ, y, params)
                                    # value = obj_function(d_1, d_2, c21H, c22H, c21L, c22L, τ, y, params)
                                    if value > max_value
                                        max_value = value
                                        max_d_1 = d_1
                                        max_d_2 = (y/λ) - d_1
                                        max_c21H = c21H
                                        max_c22H = (RH * (ω-y)/(1-λ)) - c21H
                                        max_c21L = c21L
                                        max_c22L =  (RL * (ω-y)/(1-λ)) - c21L
                                        max_y = y
                                        # max_value = value
                                        # max_d_1 = d_1
                                        # max_d_2 = d_2
                                        # max_c21H = c21H
                                        # max_c22H = c22H
                                        # max_c21L = c21L
                                        # max_c22L = c22L
                                        # max_y = y
                                    end
                                end
                            else
                                if obj_function == objective_function2                             
                                    for d_1 in d_1_range, c21H in c21H_range
                                        value = obj_function(d_1, (y/λ) - d_1, c21H, (RH * (ω-y)/(1-λ) - c21H), d_1*(r*(ω-y)+y)/(d_1+(y/λ) - d_1), ((y/λ) - d_1)*(r*(ω-y)+y)/(d_1+(y/λ) - d_1), τ, y, params)
                                        # value = obj_function(d_1,d_2, c21H,c22H, d_1*(r*(ω-y)+y)/(d_1+d_2), (d_2)*(r*(ω-y)+y)/(d_1+d_2), τ, y, params)
                                        if value > max_value
                                            max_value = value
                                            max_d_1 = d_1
                                            max_d_2 = (y/λ) - d_1
                                            max_c21H = c21H
                                            max_c22H = (RH * (ω-y)/(1-λ)) - c21H
                                            max_c21L = d_1*(r*(ω-y)+y)/(d_1+(y/λ) - d_1)
                                            max_c22L = ((y/λ) - d_1)*(r*(ω-y)+y)/(d_1+(y/λ) - d_1)
                                            max_y = y
                                            # max_value = value
                                            # max_d_1 = d_1
                                            # max_d_2 = d_2
                                            # max_c21H = c21H
                                            # max_c22H = c22H
                                            # max_c21L = c21L
                                            # max_c22L = c22L
                                            # max_y = y
                                        end
                                    end
                                elseif obj_function == objective_function3
                                    for d_1 in d_1_range, c21H in c21H_range
                                        value = obj_function(d_1, (y/λ) - d_1, c21H, (RH * (ω-y)/(1-λ)) - c21H, λ*d_1*(r*(ω-y)+y)/(λ*d_1+(y/λ) - d_1), ((y/λ) - d_1)*(r*(ω-y)+y)/(λ*d_1+(y/λ) - d_1), τ, y, params)
                                        # value = obj_function(d_1, d_2, c21H, c22H, λ*d_1*(r*(ω-y)+y)/(λ*d_1+d_2), (d_2)*(r*(ω-y)+y)/(λ*d_1+ d_2), τ, y, params)
                                        if value > max_value
                                            max_value = value
                                            max_d_1 = d_1
                                            max_d_2 = (y/λ) - d_1
                                            max_c21H = c21H
                                            max_c22H =(RH * (ω-y)/(1-λ)) - c21H
                                            max_c21L = λ*d_1*(r*(ω-y)+y)/(λ*d_1+(y/λ) - d_1)
                                            max_c22L = ((y/λ) - d_1)*(r*(ω-y)+y)/(λ*d_1+(y/λ) - d_1)
                                            max_y = y
                                            # max_value = value
                                            # max_d_1 = d_1
                                            # max_d_2 = d_2
                                            # max_c21H = c21H
                                            # max_c22H = c22H
                                            # max_c21L = c21L
                                            # max_c22L = c22L
                                            # max_y = y
                                        end
                                    end
                                elseif obj_function == objective_function4
                                    for d_1 in d_1_range, c21H in c21H_range
                                        value = obj_function(d_1, (y/λ) - d_1, c21H, (RH * (ω-y)/(1-λ)) - c21H, d_1*(r*(ω-y)+y)/(d_1+λ*((y/λ) - d_1)),λ*((y/λ) - d_1)*(r*(ω-y)+y)/(d_1+λ*((y/λ) - d_1)), τ, y,params)
                                        # value = obj_function(d_1,d_2, c21H, c22H, d_1*(r*(ω-y)+y)/(d_1+λ*d_2),λ*(d_2)*(r*(ω-y)+y)/(d_1+λ*d_2), τ, y,params)
                                        if value > max_value
                                            max_value = value
                                            max_d_1 = d_1
                                            max_d_2 = (y/λ) - d_1
                                            max_c21H = c21H
                                            max_c22H = (RH * (ω-y)/(1-λ)) - c21H
                                            max_c21L = d_1*(r*(ω-y)+y)/(d_1+λ*((y/λ) - d_1))
                                            max_c22L = λ*((y/λ) - d_1)*(r*(ω-y)+y)/(d_1+λ*((y/λ) - d_1))
                                            max_y = y
                                            # max_value = value
                                            # max_d_1 = d_1
                                            # max_d_2 = d_2
                                            # max_c21H = c21H
                                            # max_c22H = c22H
                                            # max_c21L = c21L
                                            # max_c22L = c22L
                                            # max_y = y
                                        end
                                    end                                
                                end
                            end
                        end
                        results[τ] = (max_value, max_d_1, max_d_2, max_c21H, max_c22H, max_c21L, max_c22L, max_y)
                    end
                    results
                end)
            end
            for future in worker_results
                worker_result = fetch(future)
                merge!(max_values, worker_result)
            end
        end

        τ_values = sort(collect(keys(max_values)))
        max_utility_values = [value == -Inf ? NaN : value for (τ, value) in zip(τ_values, [max_values[τ][1] for τ in τ_values])]
        max_y_values = [max_values[τ][8] for τ in τ_values]

        max_d1_values = [max_values[τ][2] for τ in τ_values]
        max_c21H_values = [max_values[τ][4] for τ in τ_values]
        max_c21L_values = [max_values[τ][6] for τ in τ_values]
        

        max_d2_values = [max_values[τ][3] for τ in τ_values]
        max_c22H_values = [max_values[τ][5] for τ in τ_values]
        max_c22L_values = [max_values[τ][7] for τ in τ_values]

        return τ_values, max_utility_values, max_d1_values, max_c21H_values, max_c21L_values, max_d2_values, max_c22H_values, max_c22L_values, max_y_values
    end
end
#endregion

# Parameters
#   params = [γ, ψ, λ, π_H, π_L, r,ω,RH,RL]  # Set your parameter values here
params = [3.0,0.40,0.15,0.80, 0.20, 0.80, 2.0, 1.5,1.065,0.94] 
yn = 10000
τn = 10000
cn = 200

# Grids
y_grid = collect(range(0.10,stop=0.50,length=yn))

τ_grid = collect(range(0.30, stop=0.90,length=τn))  


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
all_data = hcat(τ_grid, util, wavec, prob, yvec,dvec,cHvec,cLvec)

# Save the combined array to a file
writedlm("all_vectors.txt", all_data, ',')
rmprocs(workers())