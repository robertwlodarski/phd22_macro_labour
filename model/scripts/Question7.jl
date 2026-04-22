# 1. Extended container 
@with_kw mutable struct SweepResults
    b⃗::Vector{Float64}              
    γ⃗::Vector{Float64}              
    
    # A. SS objects 
    c::Matrix{Float64}              
    θ::Matrix{Float64}              
    w::Matrix{Float64}              
    q::Matrix{Float64}              
    f::Matrix{Float64}              
    J::Matrix{Float64}              
    
    # B. Aggregate-uncertainty 
    fᵐᵃˣ::Matrix{Float64}           
    fᵐⁱⁿ::Matrix{Float64}
    qᵐᵃˣ::Matrix{Float64}
    qᵐⁱⁿ::Matrix{Float64}
    wᵐᵃˣ::Matrix{Float64}
    wᵐⁱⁿ::Matrix{Float64}
    
    # C.  Simulated moments (HP-filtered logs) 
    σᵤ::Matrix{Float64}             
    σᵥ::Matrix{Float64}            
    σᵥᵤ::Matrix{Float64}            
    σₚ::Matrix{Float64}             
    
    # D. Wage statistics (levels) 
    w̄::Matrix{Float64}              # Mean wage level 
    εₚ::Matrix{Float64}             # Elasticity of log w w.r.t. log p 
end 

# 2. SS + aggregate solve + simulation sweep 
# It may take 7 years to solve 
function fnSweepBGammaFull(b⃗, γ⃗; burn = 52*20, λ = 1600)
    # A. Pre-allocate 
    Nᵦ, Nᵧ  = length(b⃗), length(γ⃗)
    M       = () -> zeros(Nᵦ, Nᵧ)
    res     = SweepResults(
                b⃗ = b⃗, γ⃗ = γ⃗,
                c = M(), θ = M(), w = M(), q = M(), f = M(), J = M(),
                fᵐᵃˣ = M(), fᵐⁱⁿ = M(), qᵐᵃˣ = M(), qᵐⁱⁿ = M(), wᵐᵃˣ = M(), wᵐⁱⁿ = M(),
                σᵤ = M(), σᵥ = M(), σᵥᵤ = M(), σₚ = M(), w̄ = M(), εₚ = M())

    # B. Loop 
    for (i, b) in enumerate(b⃗), (j, γ) in enumerate(γ⃗)
        try 
            # (i) SS calibration 
            params      = fnSetUpParameters(; b = b, γ = γ)
            SS          = fnSetUpSteadyState(params)
            c           = fnCalibrateSS!(params, SS)
            res.c[i, j]     = c 
            res.θ[i, j]     = SS.θ 
            res.w[i, j]     = SS.w 
            res.q[i, j]     = SS.q 
            res.f[i, j]     = SS.f 
            res.J[i, j]     = SS.J 

            # (ii) Aggregate uncertainty 
            Agg         = fnSetUpAggregate(params)
            fnSolveAggregate!(c, params, Agg)
            res.fᵐᵃˣ[i, j]  = maximum(Agg.f⃗)
            res.fᵐⁱⁿ[i, j]  = minimum(Agg.f⃗)
            res.qᵐᵃˣ[i, j]  = maximum(Agg.q⃗)
            res.qᵐⁱⁿ[i, j]  = minimum(Agg.q⃗)
            res.wᵐᵃˣ[i, j]  = maximum(Agg.w⃗)
            res.wᵐⁱⁿ[i, j]  = minimum(Agg.w⃗)

            # (iii) Simulation 
            fnSimulate!(params, Agg)
            u   = Agg.U⃗[burn+1:end] 
            v   = Agg.V⃗[burn+1:end] 
            p   = Agg.p⃗̂[burn+1:end] 
            θ   = Agg.Θ⃗[burn+1:end] 
            w   = Agg.W⃗[burn+1:end] 

            # (iv) HP-filtered cyclical components  
            _, ũ            = hp_filter(log.(u), λ)
            _, ṽ            = hp_filter(log.(v), λ)
            _, θ̃            = hp_filter(log.(θ), λ)
            _, p̃            = hp_filter(log.(p), λ)
            res.σᵤ[i, j]    = std(ũ)
            res.σᵥ[i, j]    = std(ṽ)
            res.σᵥᵤ[i, j]   = std(θ̃)
            res.σₚ[i, j]    = std(p̃)

            # (v) Wage stats: mean (in levels) and elasticity (log w on log p) 
            res.w̄[i, j]     = mean(w)
            logw            = log.(w)
            logp            = log.(p)
            res.εₚ[i, j]    = cov(logw, logp) / var(logp)
        catch e 
            # Leave zeros; print for debugging 
            @printf "Failed at (b=%.3f, γ=%.3f): %s\n" b γ sprint(showerror, e) 
        end 
    end 

    return res 
end 

# 3. Plotting equilibrium values 
function fnSweepPlotSS(res)
    # A. Extract data 
    @unpack b⃗, γ⃗, c, w, J     = res 
    Nᵧ      = length(γ⃗)

    # B. Colour scheme
    colors  = cgrad(:RdBu, Nᵧ; categorical = true)

    # C. Helper to build one panel 
    function panel(Z, ylabel)
        plt = plot(; xlabel = L"b", ylabel = ylabel, 
                     framestyle = :box, grid = true, gridalpha = 0.3,
                     legend = false)
        for j in 1:Nᵧ 
            plot!(plt, b⃗, Z[:, j]; lw = 1.5, color = colors[j])
        end 
        return plt 
    end 

    # D. Build the substantive panels 
    p1  = panel(c, L"c^\star")
    p2  = panel(w, L"w^\star")
    p3  = panel(J, L"J^\star")

    # E. Legend
    cbar    = heatmap(reshape(γ⃗, 1, :); 
                        c = :RdBu, 
                        colorbar = false, 
                        xticks = (1:3:Nᵧ, [@sprintf("%.2f", γ⃗[k]) for k in 1:3:Nᵧ]),
                        yticks = false, 
                        framestyle = :box, 
                        ylabel = L"\gamma",
                        yguidefontrotation = -90, 
                        titlefont = 10,
                        xrotation = 0)

    # F. Assemble: 3 panels on top, colourbar strip below 
    plt     = plot(p1, p2, p3, cbar;
                    layout = @layout([grid(1,3); a{0.05h}]),
                    size   = (1100, 450))
    return plt 
end

# 4. Sweep: Plot answer to question 4 
function fnSweepPlotRange(res)
    # A. Extract 
    @unpack b⃗, γ⃗, fᵐᵃˣ, fᵐⁱⁿ, qᵐᵃˣ, qᵐⁱⁿ, wᵐᵃˣ, wᵐⁱⁿ     = res 
    Nᵧ      = length(γ⃗)

    # B. Colour scheme 
    colors  = cgrad(:RdBu, Nᵧ; categorical = true)

    # C. Panel helper 
    function panel(Z, ylabel)
        plt = plot(; xlabel = L"b", ylabel = ylabel, 
                     framestyle = :box, grid = true, gridalpha = 0.3,
                     legend = false)
        for j in 1:Nᵧ 
            plot!(plt, b⃗, Z[:, j]; lw = 1.5, color = colors[j])
        end 
        return plt 
    end 

    # D. Build panels 
    p1  = panel(fᵐᵃˣ, L"f^{\max}")
    p2  = panel(fᵐⁱⁿ, L"f^{\min}")
    p3  = panel(qᵐᵃˣ, L"q^{\max}")
    p4  = panel(qᵐⁱⁿ, L"q^{\min}")
    p5  = panel(wᵐᵃˣ, L"w^{\max}")
    p6  = panel(wᵐⁱⁿ, L"w^{\min}")

    # E. Heatmap 
    cbar    = heatmap(reshape(γ⃗, 1, :); 
                        c = :RdBu, 
                        colorbar = false, 
                        xticks = (1:3:Nᵧ, [@sprintf("%.2f", γ⃗[k]) for k in 1:3:Nᵧ]),
                        yticks = false, 
                        framestyle = :box, 
                        ylabel = L"\gamma",
                        yguidefontrotation = -90, 
                        titlefont = 10,
                        xrotation = 0)

    # F. Combine it into so much fun for macro labour people 
    plt     = plot(p1, p2, p3, p4, p5, p6, cbar;
                    layout = @layout([grid(3,2); a{0.05h}]),
                    size   = (900, 900))
    return plt 
end

# 5. Plot simulated moments 
function fnSweepPlotMomentsSimulated(res)
    # A. Extract 
    @unpack b⃗, γ⃗, σᵤ, σᵥᵤ, w̄, εₚ     = res 
    Nᵧ      = length(γ⃗)

    # B. Colour scheme 
    colors  = cgrad(:RdBu, Nᵧ; categorical = true)

    # C. Panel helper 
    function panel(Z, ylabel)
        plt = plot(; xlabel = L"b", ylabel = ylabel, 
                     framestyle = :box, grid = true, gridalpha = 0.3,
                     legend = false)
        for j in 1:Nᵧ 
            plot!(plt, b⃗, Z[:, j]; lw = 1.5, color = colors[j])
        end 
        return plt 
    end 

    # D. Four panels of interest for Q7 
    p1  = panel(σᵤ, L"\sigma_u")
    p2  = panel(σᵥᵤ, L"\sigma_{v/u}")
    p3  = panel(w̄, L"\bar{w}")
    p4  = panel(εₚ, L"\varepsilon_{w,p}")

    # E. Heeatmap 
    cbar    = heatmap(reshape(γ⃗, 1, :); 
                        c = :RdBu, 
                        colorbar = false, 
                        xticks = (1:3:Nᵧ, [@sprintf("%.2f", γ⃗[k]) for k in 1:3:Nᵧ]),
                        yticks = false, 
                        framestyle = :box, 
                        ylabel = L"\gamma",
                        yguidefontrotation = -90, 
                        titlefont = 10,
                        xrotation = 0)

    # F. Assemble 
    plt     = plot(p1, p2, p3, p4, cbar;
                    layout = @layout([grid(2,2); a{0.05h}]),
                    size   = (900, 700))
    return plt 
end