# DSO670_Final_Project

This repo expands upon the experiments in [1] and was completed as a final project for DSO 670: Data-Driven Optimization: Theory, Methods, and Current Themes in November 2020.

Files are organized as follows:

main.ipynb -- Jupyter notebook that runs main experiments
non_DRO_methods.jl -- Julia file for calculating the non-DRO methods for the paper (baselines) -- Scarf, SAA, critical fractile
LRO.jl -- Julia file for implementation of LRO method defined in [1]. Utilizes JuMP and Ipopt
utils.jl -- Julia file for utility functions

Questions? Comments? Concerns? Feel free to contact Caroline Johnston, cmjohnst@usc.edu

[1] Wang  Z,  Glynn  PW,  Ye  Y  (2016)  Likelihood  robust  optimization  for  data-driven  problems.Compu-tational Management Science13(2):241–261,  ISSN  16196988,  URLhttp://dx.doi.org/10.1007/s10287-015-0240-3
