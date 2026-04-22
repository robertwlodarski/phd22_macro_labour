# PHD22 Macro labour  
# Replication: Shimer (2005) & Hagedorn-Manovskii (2008) calibrations of DMP 
# April 2026 
# Author: Rob Włodarski

# 1. Packages & load functions 
using Parameters, QuantEcon, Roots, Random, Statistics, Plots, Printf, LinearAlgebra, LaTeXStrings
Threads.nthreads()
include("scripts/ModelInfrastructure.jl")
include("scripts/Functions.jl")
include("scripts/Question3.jl")
include("scripts/Question4.jl")
include("scripts/Question5.jl")

# 2. Calibrate c to hit target steady-state unemployment (Question 3) 
@time cˢˢ = fnCalibrateSS!(UsedParameters, SS)
fnPrintSteadyState(UsedParameters, SS, cˢˢ)
plt = fnPlotCalibration(UsedParameters, SS, cˢˢ)
savefig(plt, "plots/calibration.pdf")

# 3. Solve the model with aggregate uncertainty (Question 4) 
@time fnSolveAggregate!(cˢˢ, UsedParameters, Agg)
fᵐᵃˣ = maximum(Agg.f⃗)
fᵐⁱⁿ = minimum(Agg.f⃗)
@printf "Highest job-finding rate: %.4f\n" fᵐᵃˣ
@printf "Lowest  job-finding rate: %.4f\n" fᵐⁱⁿ
plt2 = fnPlotProductivity(UsedParameters)
savefig(plt2, "plots/productivity.pdf")
fnPrintEquilibriumRange(UsedParameters,SS,Agg) 

# 4. Simulate and replicate Hagedorn-Manovskii Table 4 (Question 5) 
@time fnSimulate!(UsedParameters, Agg)
TabMoments = fnTable4Moments(Agg; burn = 52 * 20, λ = 1600)
fnPrintTable4(TabMoments)

# Discretisation issue
AlternativeParameters   = fnSetUpParameters(Nₚ = 21)
AggAP                   = fnSetUpAggregate(AlternativeParameters)
@time fnSolveAggregate!(cˢˢ, AlternativeParameters, AggAP)
@time fnSimulate!(AlternativeParameters, AggAP)
TabMomentsAlt           = fnTable4Moments(AggAP)   
fnPrintTable4(TabMomentsAlt; path = joinpath("tables", "hm_table4Alt.tex"), caption = "Results from the calibrated model with more gridpoints")

# More serious treatment of the discretisation issue 
N⃗ₚ      = [7, 9, 11, 13, 15, 17, 19, 21]
sweep   = fnSweepNp(cˢˢ, N⃗ₚ)
plt_σ   = fnPlotDiscretisationStd(sweep)
savefig(plt_σ, "plots/discretisation_std.pdf")

plt_M   = fnPlotDiscretisationCorr(sweep)
savefig(plt_M, "plots/discretisation_corr.pdf")

# 5. Hagedorn-Manovskii calibration (Question 6): b = 0.94, γ = 0.052 
ParamsHM    = fnSetUpParameters(; b = 0.94, γ = 0.052)
SSHM        = fnSetUpSteadyState(ParamsHM)
AggHM       = fnSetUpAggregate(ParamsHM)

@time cₕₘ = fnCalibrateSS!(ParamsHM, SSHM)
@time fnSolveAggregate!(cₕₘ, ParamsHM, AggHM)
@time fnSimulate!(ParamsHM, AggHM)

## 6. Grid over (b, γ) for Question 7 
# [... to be written]