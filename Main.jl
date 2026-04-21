# PHD22 Macro labour  
# Replication: Shimer (2005) & Hagedorn-Manovskii (2008) calibrations of DMP 
# Problem set for Sebastian Graves, Lent Term 2026
# April 2026 

## 1. Packages & load functions 
using Parameters, QuantEcon, Roots, Random, Statistics, Plots, Printf
Threads.nthreads()
include("scripts/ModelInfrastructure.jl")
include("scripts/Functions.jl")

## 2. Calibrate c to hit target steady-state unemployment (Question 3) 
@time cˢˢ = fnCalibrateSS!(UsedParameters, SS)
#fnPrintSteadyState(UsedParameters, SS) 

## 3. Solve the model with aggregate uncertainty (Question 4) 
@time fnSolveAggregate!(cˢˢ, UsedParameters, Agg)
fᵐᵃˣ = maximum(Agg.f⃗)
fᵐⁱⁿ = minimum(Agg.f⃗)
@printf "Highest job-finding rate: %.4f\n" fᵐᵃˣ
@printf "Lowest  job-finding rate: %.4f\n" fᵐⁱⁿ

## 4. Simulate and replicate Hagedorn-Manovskii Table 4 (Question 5) 
@time fnSimulate!(UsedParameters, Agg)
Tbl = fnTable4(UsedParameters, Agg; burn = 100)
#fnPrintTable4(Tbl; path = joinpath("tables", "hm_table4.tex"))

## 5. Hagedorn-Manovskii calibration (Question 6): b = 0.94, γ = 0.052 
ParamsHM    = fnSetUpParameters(; b = 0.94, γ = 0.052)
SSHM        = fnSetUpSteadyState(ParamsHM)
AggHM       = fnSetUpAggregate(ParamsHM)

@time cₕₘ = fnCalibrateSS!(ParamsHM, SSHM)
@time fnSolveAggregate!(cₕₘ, ParamsHM, AggHM)
@time fnSimulate!(ParamsHM, AggHM)
#TblHM = fnTable4(ParamsHM, AggHM; burn = 100)

## 6. Grid over (b, γ) for Question 7 
# [... to be written]