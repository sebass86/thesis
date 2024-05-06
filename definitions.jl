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
