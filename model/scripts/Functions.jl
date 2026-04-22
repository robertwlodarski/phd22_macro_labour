# Content 
# A. Helper functions 
# 1. Matching function helpers 
# 2. Bargained wage 
# 3. Steady state unemployment rate 
# 4. Unemployment law of motion 
# B. Steady state functions 
# 1. Equilibrium
# 2. Steady state calibration residual 
# 3. Steady state calibration (c)  
# C. Aggregate uncertanty functions 
# 1. Bellman residual
# 2. Aggregate uncertainty solver
# 3. Simulator

# -----------------------------------------------------------------------------
#  A. HELPER FUNCTIONS 
# -----------------------------------------------------------------------------

# 1. Matching function helpers 
# A. Job-filling rate from tightness
function fnq(θ,params)
    @unpack l   = params
    return (1 + θ^l)^(-1/l)
end 

# B. Tightness from job-filling rate (inverse)
function fnqInverse(q,params)
    @unpack l   = params
    return (q^(-l) - 1)^(1/l)
end 

# C. Job-finding rate from tightness
function fnf(θ,params)
    @unpack l   = params
    return (1 + θ^(-l))^(-1/l)
end 

# D. Tightness from job-finding rate (inverse)
function fnfInverse(f,params)
    @unpack l   = params
    return (f^(-l) - 1)^(-1/l)
end

# 2. Bargained wage 
function fnW(p,θ,c,params)
    @unpack γ, b   = params
    return γ * (p + θ * c) + (1 - γ) * b
end 

# 3. Steady state unemployment rate (Beveridge curve)
function fnSteadyStateU(f,params)
    @unpack s   = params
    return s / (s + f)
end 

# 4. Unemployment law of motion 
function fnNextU(u,f,params)
    @unpack s   = params
    return u + s * (1 - u) - f * u
end

# -----------------------------------------------------------------------------
# % B. STEADY STATE FUNCTIONS 
# -----------------------------------------------------------------------------

# 1. Equilibrium
function fnEquilibrium!(p,c,params,SS=nothing)
    # A. Unpacking business 
    @unpack β, s, δᵍ,γ,b    = params 

    # B. Residual little helper 
    function fnResidualq(q)
        θ       = fnqInverse(q,params)
        w       = fnW(p,θ,c,params)
        q̃       = c * (1 - β * (1 - s)) / (β * (p - w))
        return q - q̃  
    end 
    # C. Avoiding explosion (learnt the hard way)
    θ̅       = (1 - γ) * (p - b) / (γ * c)
    q̲       = fnq(θ̅, params) * 1.001             

    # D. Solve 
    q       = find_zero(fnResidualq, (q̲, 1 - 1e-6), Roots.Brent(); xatol = δᵍ)

    # E. Equilibrium objects 
    θ       = fnqInverse(q,params)
    w       = fnW(p,θ,c,params)
    f       = fnf(θ,params)

    # F. Update SS struct 
    if !isnothing(SS)
        SS.q    = q 
        SS.θ    = θ 
        SS.w    = w 
        SS.f    = f 
        SS.u    = fnSteadyStateU(f,params)
        SS.J    = (p - w) / (1 - β * (1 - s))
    end 
    return θ, w, q, f 
end 

# 2. Steady state calibration residual 
function fnResidualSS!(c,params,SS=nothing)
    # A. Unpacking business 
    @unpack pₛₛ, uₛₛ    = params

    # B. Solve inner fixed point at steady-state productivity 
    _, _, _, f          = fnEquilibrium!(pₛₛ,c,params,SS)

    # C. Implied unemployment from Beveridge curve 
    u                   = fnSteadyStateU(f,params)

    # D. Residual 
    εᶜ                  = u - uₛₛ 

    # E. Update 
    if !isnothing(SS)
        SS.c            = c 
        SS.εᶜ           = εᶜ 
    end 
    return εᶜ 
end

# 3. Steady state calibration
function fnCalibrateSS!(params,SS)
    # A. Unpacking business 
    @unpack β, s, γ, b, pₛₛ, uₛₛ   = params

    # B. Implied eq objects
    f       = s * (1 - uₛₛ) / uₛₛ 
    θ       = fnfInverse(f, params)
    q       = fnq(θ, params)

    # C. Closed-form c (from free entry + NB)
    c       = q * β * (1 - γ) * (pₛₛ - b) / ((1 - β * (1 - s)) + q * β * γ * θ)

    # D. Recover other objects
    w       = fnW(pₛₛ, θ, c, params)
    J       = (pₛₛ - w) / (1 - β * (1 - s))

    # E. Update 
    SS.c    = c 
    SS.q    = q 
    SS.θ    = θ 
    SS.f    = f 
    SS.w    = w 
    SS.J    = J 
    SS.u    = uₛₛ 
    SS.εᶜ   = 0.0 
    return c 
end

# -----------------------------------------------------------------------------
# C. AGGREGATE UNCERTAINTY FUNCTIONS 
# -----------------------------------------------------------------------------

# 1. Bellman residual
function fnBellman(J⃗,c,params)
    # A. Unpacking business 
    @unpack β, s, P, p⃗, b, γ    = params

    # B. Continuation value 
    𝔼J⃗      = P * J⃗ 

    # C. Back out tightness 
    q⃗       = c ./ (β .* 𝔼J⃗)                                
    q⃗       = clamp.(q⃗, 1e-8, 1 - 1e-8)                    
    θ⃗       = fnqInverse.(q⃗, Ref(params))

    # D. Wage and updated value 
    w⃗       = fnW.(p⃗, θ⃗, c, Ref(params))
    J⃗ⁿᵉʷ    = p⃗ .- w⃗ .+ β * (1 - s) .* 𝔼J⃗ 

    return J⃗ⁿᵉʷ, θ⃗, w⃗, q⃗, 𝔼J⃗ 
end 

# 2. Aggregate uncertainty solver
function fnSolveAggregate!(c,params,agg)
    # A. Unpacking business 
    @unpack Nₚ, p⃗, δᵛ, n̅ᵛᶠⁱ     = params

    # B. Initialise J guess using the steady-state formula 
    @unpack β, s, pₛₛ   = params
    _, wₛₛ, _, _        = fnEquilibrium!(pₛₛ, c, params)
    J⃗                   = fill((pₛₛ - wₛₛ) / (1 - β * (1 - s)), Nₚ)

    # C. Live, love, iterate
    εᵛ      = Inf 
    nᵛ      = 0
    θ⃗       = zeros(Nₚ)
    w⃗       = zeros(Nₚ)
    q⃗       = zeros(Nₚ)
    𝔼J⃗      = zeros(Nₚ)
    while εᵛ > δᵛ && nᵛ < n̅ᵛᶠⁱ
        J⃗ⁿᵉʷ, θ⃗, w⃗, q⃗, 𝔼J⃗   = fnBellman(J⃗, c, params)
        εᵛ                  = maximum(abs.(J⃗ⁿᵉʷ .- J⃗))
        J⃗                   .= J⃗ⁿᵉʷ 
        nᵛ                  += 1 
    end 

    # D. Recover the job-finding rate 
    f⃗       = fnf.(θ⃗, Ref(params))

    # E. Update
    agg.θ⃗   .= θ⃗ 
    agg.w⃗   .= w⃗ 
    agg.q⃗   .= q⃗ 
    agg.f⃗   .= f⃗ 
    agg.J⃗   .= J⃗ 
    agg.𝔼J⃗  .= 𝔼J⃗ 
    agg.εᵛ  = εᵛ 
    agg.nᵛ  = nᵛ 
end

# 3. Simulator 
function fnSimulate!(params,agg)
    # A. Unpacking business 
    @unpack s, T                = params
    @unpack p⃗̂, p⃗̂ᵢ, θ⃗, w⃗, q⃗, f⃗   = agg 

    # B. Initial unemployment guess 
    agg.U⃗[1]    = fnSteadyStateU(f⃗[p⃗̂ᵢ[1]], params)

    # C. Iterate, iterate iterate  
    for t in 1:T
        i           = p⃗̂ᵢ[t] 
        agg.Θ⃗[t]    = θ⃗[i] 
        agg.W⃗[t]    = w⃗[i] 
        agg.Q⃗[t]    = q⃗[i] 
        agg.F⃗[t]    = f⃗[i]

        # (i) Flow equations 
        u           = agg.U⃗[t] 
        agg.V⃗[t]    = agg.Θ⃗[t] * u                          
        agg.Y⃗[t]    = p⃗̂[t] * (1 - u)                        

        # (ii) Update unemployment for t+1 using the law of motion 
        if t < T 
            agg.U⃗[t+1]  = fnNextU(u, agg.F⃗[t], params)
        end 
    end 
end
# This all will save me 0.0001 seconds compared to doing that in Matlab. 
