using JuMP
using Ipopt
using CSV, DataFrames, DataStructures, MathOptFormat
using Statistics,Distributions, StatsFuns, StatsBase
using PyPlot


function critical_fractile(dist_type)
    #closed form of solution for newsvendor for normal distribution
    # dist_type - str, "norm" or "exp" for normal or exponential distribution, respectively
    
    b = 1 #underage cost
    h = 1 #overage cost
    
    if dist_type == "norm"
        mean = 50
        var = 50

        beta = b/(b+h)
        
        b_inv = norminvcdf(beta) #inverse of the cdf for normal distribution 
        
        return mean + b_inv * var
        
    elseif dist_type == "exp"
        # plugged and chugged with cdf of exponential function to get warm start value (solver kept converging to point of infeasibility)
        cdf_func(x) = cdf(Exponential(50),x)
        
        # print(cdf_func(34.65))
        
        model3 = Model(with_optimizer(Ipopt.Optimizer,print_level=0))

        @variable(model3, x_var, start=34) #use warm start to avoid converge to infeasible point (see above)

        JuMP.register(model3, :cdf_func, 1, cdf_func; autodiff = true)  #to use cdf function in constraint

        @NLconstraint(model3,con1, cdf_func(x_var) >= b/(b+h))

        @objective(model3,Min,x_var)
        optimize!(model3)

        return value(x_var)
        
    else
        return nothing 
    end
end

function scarf(df)
    #closed form scarf solution
    # df - dataframe, data of observed demand 
    
    b = 1 #underage cost
    h = 1 #overage cost
    
    mu_hat = mean(df[:,"Column1"]) #sample mean
    sigma_sq_hat = var(df[:,"Column1"]) #sample variance
    
    return mu_hat + sigma_sq_hat/2 * (sqrt(b / h) - sqrt(h / b))
end


function SAA(df)
    # compute SAA solution (empirical risk minimization)
    # df - dataframe, data of observed demand 
    
    b = 1 #underage cost
    h = 1 #overage cost
    
    model = Model(with_optimizer(Ipopt.Optimizer,print_level=0))

    data_size = length(df[:,"Column1"])

    @variable(model, x_var) #stocking quantity 
    @variable(model, y_vars[1:data_size]) #variable to linearize max(d_i - x, x - d_i) in objective function

    @constraint(model,con1[i = 1:data_size], y_vars[i] >= x_var - df[:,"Column1"][i])
    @constraint(model, con2[i = 1:data_size], y_vars[i] >= df[:,"Column1"][i] - x_var )

    @NLobjective(model,Min, 1 / data_size * sum(y_vars[i] for i = 1:data_size))
    
    optimize!(model)
    
    return value(x_var)
end