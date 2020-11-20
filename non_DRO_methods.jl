using JuMP
using Ipopt
using CSV, DataFrames, DataStructures, MathOptFormat
using Statistics,Distributions, StatsFuns, StatsBase
using PyPlot


#closed form of solution for newsvendor for normal distribution
function critical_fractile()
    b = 1
    h = 1
    
    mean = 50
    var = 50
    
    beta = b/ (b+h)
    b_inv = norminvcdf(beta) #inverse of the cdf for normal distribution 
    return mean + b_inv * var,'\n'
end

#closed form scarf solution
function scarf(df)
    b=1
    h=1
    
    mu_hat = mean(df[:,"Column1"])
    sigma_sq_hat = var(df[:,"Column1"])
    
    return mu_hat + sigma_sq_hat/2 * (sqrt(b / h) - sqrt(h / b))
end

#SAA 
function SAA(df)
    model = Model(with_optimizer(Ipopt.Optimizer,print_level=0))

    data_size = length(df[:,"Column1"])

    @variable(model, x_var) #stocking quantity 
    @variable(model, y_vars[1:data_size]) #variable to linearize max(d_i - x, x - d_i) in objective function

    b = 1
    h = 1

    @constraint(model,con1[i = 1:data_size], y_vars[i] >= x_var - df[:,"Column1"][i])
    @constraint(model, con2[i = 1:data_size], y_vars[i] >= df[:,"Column1"][i] - x_var )

    @NLobjective(model,Min, 1 / data_size * sum(y_vars[i] for i = 1:data_size))
    
    optimize!(model)
    
    println("x = ", value(x_var))
    return value(x_var)
end