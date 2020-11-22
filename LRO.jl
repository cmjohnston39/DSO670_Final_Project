using JuMP
using Ipopt
using CSV, DataFrames, DataStructures, MathOptFormat
using Statistics,Distributions, StatsFuns, StatsBase
using PyPlot


function LRO_opt(df,gamma_step,moment_info, x_bounds=nothing, recover_SAA=nothing)
    # calculate optimal stocking quantity x_var using LRO framework (convexification of problem (14) in the paper)
    # df - dataframe, historical observed demand data
    #gamma_step - int, amount to vary gamma by (for experments)
    #moment_info - boolean, true if using moment information, false else
    #x_bounds - int or nothing, value to fix x to (for experiments), free otherwise
    
    b = 1 #underage cost
    h = 1 #overage cost
    
    N_i = get_occurences(df)  #construct N_i observation data
    N = sum(N_i[:,"Occurence"])
    
    gamma_temp = get_gamma_val(N_i, N)  #get gamma value
    
    gamma = gamma_temp += gamma_step  
    
    model = Model(with_optimizer(Ipopt.Optimizer,print_level=0))
    
    @variable(model, lambda_var >=0)

    
    
    @variable(model, y_vars[1:length(N_i[:,"Occurence"])] >=0)
    
    if isnothing(x_bounds)
        @variable(model, x_var)
    else
        @variable(model, x_bounds <= x_var <= x_bounds)
    end
    
    if moment_info
        @variable(model, mu_vars[1:3])

    else
        @variable(model, mu_var)
    end

    z = sum(N_i[:,"Occurence"] .* log.(N_i[:,"Occurence"])) #to store log likelihood
    
    var_len = length(N_i[:,"Occurence"])
    
    if !isnothing(recover_SAA) #only if testing if we can recover SAA for sanity
        data_size = length(df[:,"Column1"])
        
        @variable(model, z_vars[1:data_size]) #define variable to linearize max behavior of objective function

        @constraint(model,con3[i = 1:data_size], z_vars[i] >= x_var - df[:,"Column1"][i])
        @constraint(model, con4[i = 1:data_size], z_vars[i] >= df[:,"Column1"][i] - x_var )
        
        @constraint(model, lambda_var == 1/N *(sum(z_vars[i]  for i in 1:var_len) - mu_var)) # analytically solve for lambda according to Proposition 1, and p_i is empirical distribution
    end
    
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
   
    if isnothing(x_bounds)
        return value(x_var) # return x* (optimal stocking quantity) for experiments
    else
        return objective_value(model) # return objective value for experiments
    end
end