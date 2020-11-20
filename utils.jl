using JuMP
using Ipopt
using CSV, DataFrames, DataStructures, MathOptFormat
using Statistics,Distributions, StatsFuns, StatsBase
using PyPlot

# write out LP file for model according to 'model_name'
# (will not write out NLP constraints/objective functions)
function write_LP(model,model_name)
    lp_file = MathOptFormat.LP.Model()
    MOI.copy_to(lp_file, backend(model))
    MOI.write_to_file(lp_file, string(model_name, ".lp")) 
end

#calculate gamma value for ambiguity set
function get_gamma_val(N_i,N)
    #TODO: need to change the chi-squared value as a function of the support (?) [0,200]
    gamma = sum(N_i[:,"Occurence"] .* log.(N_i[:,"Occurence"] / N)) - 77.929 #chi-squared, with d.o.f. n-1, 1- alpha = 0.95
end

#returns dataframe with number of occurences of each element in given dataframe
function get_occurences(df) 
    N_i = combine(nrow,groupby(df,:Column1))
    rename!(N_i, ["Demand", "Occurence"])
end

