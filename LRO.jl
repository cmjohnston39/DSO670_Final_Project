using JuMP
using Ipopt
using CSV, DataFrames, DataStructures, MathOptFormat
using Statistics,Distributions, StatsFuns, StatsBase
using PyPlot

# @nbinclude("utils.ipynb")
# function get_gamma_val(N_i, N)
#     gamma = sum(N_i[:,"Occurence"] .* log.(N_i[:,"Occurence"] / N)) - 77.929
# end

# calculate optimal stocking quantity x_var using LRO framework (convexification of problem (14) in the paper)
# df - historical demand data, gamma_step - amount to vary gamma by (for experments), moment_info - true if using moment information
function LRO_opt(df,gamma_step,moment_info)
    model = Model(with_optimizer(Ipopt.Optimizer,print_level=0))
    
    N_i = get_occurences(df) 
    N = sum(N_i[:,"Occurence"])
    
    gamma_temp = get_gamma_val(N_i, N)
    
    gamma = gamma_temp += gamma_step  
    
#     print("gamma is", gamma, "\n")

    @variable(model, lambda_var >=0)
    @variable(model, x_var)
    @variable(model, y_vars[1:length(N_i[:,"Occurence"])] >=0)
    
    if moment_info
        @variable(model, mu_vars[1:3])

    else
        @variable(model, mu_var)
    end

    b = 1
    h = 1

    z = sum(N_i[:,"Occurence"] .* log.(N_i[:,"Occurence"]))
    
    var_len = length(N_i[:,"Occurence"])
    
    #use sample mean/var information
    if moment_info
        mu_hat = mean(df[:,"Column1"])
        sigma_sq_hat = var(df[:,"Column1"])
        b_vec = [1, mu_hat, mu_hat^2 + sigma_sq_hat]
        
        @constraint(model,con1[i = 1:var_len], b * (N_i[:,"Demand"][i] - x_var) + y_vars[i] <= mu_vars[1]
        + mu_vars[2] * N_i[:,"Demand"][i] + mu_vars[3] * N_i[:,"Demand"][i]^2)

        @constraint(model, con2[i = 1:var_len], h * (x_var - N_i[:,"Demand"][i]) + y_vars[i] <= mu_vars[1]
        + mu_vars[2] * N_i[:,"Demand"][i] + mu_vars[3] * N_i[:,"Demand"][i]^2)
    else
        @constraint(model,con1[i = 1:var_len], b * (N_i[:,"Demand"][i] - x_var) + y_vars[i] <= mu_var)
        @constraint(model, con2[i = 1:var_len], h * (x_var - N_i[:,"Demand"][i]) + y_vars[i] <= mu_var)
    end
    
    if moment_info
        @NLobjective(model, Min, sum(b_vec[i] * mu_vars[i] for i in 1:3) + lambda_var * (z - N - gamma) + N * lambda_var * log(lambda_var)
         - lambda_var * sum(N_i[:,"Occurence"][j] * log(y_vars[j]) for j in 1:var_len))
    else
        @NLobjective(model, Min, mu_var + lambda_var * (z - N - gamma) + N * lambda_var * log(lambda_var)
         - lambda_var * sum(N_i[:,"Occurence"][j] * log(y_vars[j]) for j in 1:var_len))
    end

    write_LP(model,"model")
    
    optimize!(model)
    
    println(" x*: ", value(x_var))
    return value(x_var)
end