using JuMP
using Ipopt
using CSV, DataFrames, DataStructures, MathOptFormat
using Statistics,Distributions, StatsFuns, StatsBase
using PyPlot

function write_LP(model,model_name)
    # write out LP file for model according to 'model_name'
    # (will not write out NLP constraints/objective functions)
    # model - JUMP model type, model to save
    # model_name - string, name of model to save file as 
    
    lp_file = MathOptFormat.LP.Model()
    MOI.copy_to(lp_file, backend(model))
    MOI.write_to_file(lp_file, string(model_name, ".lp")) 
end

function get_gamma_val(N_i,N)
    #calculate gamma value for ambiguity set
    #N_i - df, dataframe of number of occurences per demand
    #N - int, total number of observations
    
    gamma = sum(N_i[:,"Occurence"] .* log.(N_i[:,"Occurence"] / N)) - 1/2*cquantile(Chisq(201-1),.95) # Theorem 1, chi-squared with d.o.f. n-1, 1- alpha = 0.95
end

function get_occurences(df) 
    #returns dataframe with number of occurences of each element in given dataframe
    #df - dataframe of observed demand data
    
    N_i = combine(nrow,groupby(df,:Column1))
    rename!(N_i, ["Demand", "Occurence"])
end

