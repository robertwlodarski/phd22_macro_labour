# Model infrastructure: simple DMP model (Shimer 2005 / Hagedorn-Manovskii 2008 style)
# Replication for PhD22 Macro Labour problem set

# Content
# 1. Params           (struct & constructor)
# 2. Steady state     (struct & constructor)
# 3. Aggregate unc.   (struct & constructor)

# ------------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------------
@with_kw struct ModelParameters

    # A. Deep parameters (Shimer 2005 calibration from problem set)
    b::Float64      = 0.4                       # Home production / unemployment flow value
    β::Float64      = 0.99^(1/12)               # Discount factor (monthly)
    r::Float64      = 1/β - 1                   # Implied discount rate
    s::Float64      = 0.0081                    # Separation rate
    γ::Float64      = 0.7                       # Worker's bargaining power
    l::Float64      = 0.407                     # Matching function elasticity (Hagedorn-Manovskii)
    c::Float64      = 0.1                       # Vacancy posting cost (to be calibrated)
    pₛₛ::Float64    = 1.0                       # Steady-state productivity normalisation

    # B. Calibration target
    uₛₛ::Float64    = 0.055                     # Target steady-state unemployment rate

    # C. Aggregate productivity process
    ρ::Float64      = 0.9895                    # Productivity persistence
    σ̃::Float64      = 0.0034                    # Unconditional std of log productivity
    σ::Float64      = σ̃ * sqrt(1 - ρ^2)         # Innovation std
    Nₚ::Int         = 15                        # Number of Rouwenhorst nodes
    P::Matrix{Float64}                          # Transition matrix
    p⃗::Vector{Float64}                          # Productivity grid

    # D. Steady-state solver settings 
    δᵍ::Float64     = 1e-12                     # Tolerance: inner (job-filling rate) loop
    δᶜ::Float64     = 1e-10                     # Tolerance: outer (calibration) loop
    q⁰::Float64     = 0.5                       # Initial guess for job-filling rate
    c̲::Float64      = 1e-4                      # Lower bracket for bisection on c
    c̅::Float64      = 10.0                      # Upper bracket for bisection on c

    # E. Aggregate-uncertainty solver settings
    δᵛ::Float64     = 1e-8                      # Tolerance: value function iteration
    n̅ᵛᶠⁱ::Int       = 2000                      # Max VFI iterations

    # F. Simulation settings
    T::Int          = 52 * 100                  # Number of simulated periods
    seed::Int       = 1997                      # RNG seed
end

# Constructor
function fnSetUpParameters(; ρ = 0.9895, σ̃ = 0.0034, Nₚ = 15, kwargs...)

    # A. Innovation std
    σ = σ̃ * sqrt(1 - ρ^2)

    # B. Rouwenhorst discretisation of log productivity
    ℳ𝒞  = rouwenhorst(Nₚ, ρ, σ)
    P   = ℳ𝒞.p
    p⃗   = exp.(ℳ𝒞.state_values)

    # C. Return struct
    return ModelParameters(
        ρ = ρ, σ̃ = σ̃, σ = σ, Nₚ = Nₚ, P = P, p⃗ = p⃗; kwargs...
    )
end

# ------------------------------------------------------------------
# 2. Steady-state variables
# ------------------------------------------------------------------
@with_kw mutable struct SteadyState

    # A. Equilibrium prices and quantities
    θ::Float64      = 0.0               # Labour market tightness
    w::Float64      = 0.0               # Wage
    q::Float64      = 0.0               # Job-filling rate
    f::Float64      = 0.0               # Job-finding rate
    u::Float64      = 0.0               # Unemployment rate
    J::Float64      = 0.0               # Value of a filled job

    # B. Calibrated parameter
    c::Float64      = 0.0               # Calibrated vacancy posting cost

    # C. Convergence diagnostics
    εᵍ::Float64   = Inf                 # Final inner-loop residual
    εᶜ::Float64   = Inf                 # Final outer-loop residual (u - u⋆)
    nᵍ::Int         = 0                 # Inner-loop iteration count
    nᶜ::Int         = 0                 # Outer-loop iteration count
end

# Constructor
function fnSetUpSteadyState(params::ModelParameters)
    return SteadyState()
end

# ------------------------------------------------------------------
# 3. Aggregate-uncertainty variables
# ------------------------------------------------------------------
# Naming convention:
# - Capital letters with arrows (Θ⃗, N⃗, ...): simulated aggregate series
# - Matrices indexed by productivity state: cross-sectional objects on p-grid
@with_kw mutable struct AggregateUncertainty

    # A. Policy and value functions on the productivity grid
    θ⃗::Vector{Float64}              # Tightness policy θ(p)
    w⃗::Vector{Float64}              # Wage policy w(p)
    q⃗::Vector{Float64}              # Job-filling rate q(p)
    f⃗::Vector{Float64}              # Job-finding rate f(p)
    J⃗::Vector{Float64}              # Firm value J(p)
    𝔼J⃗::Vector{Float64}             # Continuation value E[J(p')|p]

    # B. Simulated exogenous process
    p⃗̂::Vector{Float64}              # Simulated productivity path
    p⃗̂ᵢ::Vector{Int}                 # Indices into p-grid

    # C. Simulated aggregate series
    Θ⃗::Vector{Float64}              # Tightness
    W⃗::Vector{Float64}              # Wage
    Q⃗::Vector{Float64}              # Job-filling rate
    F⃗::Vector{Float64}              # Job-finding rate
    U⃗::Vector{Float64}              # Unemployment rate
    V⃗::Vector{Float64}              # Vacancies
    Y⃗::Vector{Float64}              # Output

    # D. Convergence diagnostics
    εᵛ::Float64     = Inf           # Final VFI residual
    nᵛ::Int         = 0             # VFI iteration count
end

# Constructor
function fnSetUpAggregate(params::ModelParameters; T = nothing, seed = nothing)

    # A. Unpack and override if supplied
    @unpack Nₚ, P, p⃗ = params
    T    = isnothing(T)    ? params.T    : T
    seed = isnothing(seed) ? params.seed : seed

    # B. Simulate the productivity process
    Random.seed!(seed)
    p⃗̂ᵢ       = zeros(Int, T)
    p⃗̂        = zeros(Float64, T)
    cs       = (Nₚ + 1) ÷ 2                 # Start at the median state
    p⃗̂ᵢ[1]    = cs
    p⃗̂[1]     = p⃗[cs]
    𝐹P       = cumsum(P, dims = 2)
    for t in 1:(T-1)
        ns      = searchsortedfirst(𝐹P[cs, :], rand())
        p⃗̂ᵢ[t+1] = ns
        p⃗̂[t+1]  = p⃗[ns]
        cs      = ns
    end

    # C. Pre-allocate series
    return AggregateUncertainty(
        θ⃗   = zeros(Nₚ),
        w⃗   = zeros(Nₚ),
        q⃗   = zeros(Nₚ),
        f⃗   = zeros(Nₚ),
        J⃗   = zeros(Nₚ),
        𝔼J⃗  = zeros(Nₚ),
        p⃗̂   = p⃗̂,
        p⃗̂ᵢ  = p⃗̂ᵢ,
        Θ⃗   = zeros(T),
        W⃗   = zeros(T),
        Q⃗   = zeros(T),
        F⃗   = zeros(T),
        U⃗   = zeros(T),
        V⃗   = zeros(T),
        Y⃗   = zeros(T)
    )
end

# ------------------------------------------------------------------
# Driver: instantiate everything
# ------------------------------------------------------------------
UsedParameters  = fnSetUpParameters()
SS              = fnSetUpSteadyState(UsedParameters)
Agg             = fnSetUpAggregate(UsedParameters)