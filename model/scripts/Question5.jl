# 1. Compute HM Table 4 moments 
function fnTable4Moments(Agg; burn = 52 * 20, λ = 1600)
    # A. Strip burn-in from simulated series 
    u      = Agg.U⃗[burn+1:end]
    v      = Agg.V⃗[burn+1:end]
    p      = Agg.p⃗̂[burn+1:end]
    θ      = Agg.Θ⃗[burn+1:end]               

    # B. Take logs 
    lu     = log.(u)
    lv     = log.(v)
    lp     = log.(p)
    lθ     = log.(θ)

    # C. HP filter: keep cyclical component 
    _, ũ    = hp_filter(lu, λ)
    _, ṽ    = hp_filter(lv, λ)
    _, p̃    = hp_filter(lp, λ)
    _, θ̃    = hp_filter(lθ, λ)

    # D. Standard deviations 
    σ       = [std(ũ), std(ṽ), std(θ̃), std(p̃)]

    # E. Autocorrelations (lag 1) 
    ρ       = [cor(ũ[1:end-1], ũ[2:end]),
               cor(ṽ[1:end-1], ṽ[2:end]),
               cor(θ̃[1:end-1], θ̃[2:end]),
               cor(p̃[1:end-1], p̃[2:end])]

    # F. Correlation matrix 
    M       = cor(hcat(ũ, ṽ, θ̃, p̃))

    return (σ = σ, ρ = ρ, M = M) 
end 

# 2. LaTeX table: HM Table 4 replication 
function fnPrintTable4(moments; path = joinpath("tables", "hm_table4.tex"), caption = "Results from the calibrated model")
    @unpack σ, ρ, M     = moments 

    open(path, "w") do io 
        write(io, "\\caption{$caption}\n")
        write(io, "\\begin{tabular}{lcccc}\n")
        write(io, "\\toprule\n")
        write(io, " & \$u\$ & \$v\$ & \$v/u\$ & \$p\$ \\\\\n")
        write(io, "\\midrule\n")

        # Standard deviations 
        @printf io "Standard deviation        & %6.3f & %6.3f & %6.3f & %6.3f \\\\\n" σ[1] σ[2] σ[3] σ[4] 
        @printf io "Autocorrelation           & %6.3f & %6.3f & %6.3f & %6.3f \\\\\n" ρ[1] ρ[2] ρ[3] ρ[4] 
        write(io, "\\midrule\n")

        # Correlation matrix 
        write(io, "\\multicolumn{5}{l}{\\textit{Correlation matrix}} \\\\\n")
        @printf io "\$u\$   & %6.3f & %6.3f & %6.3f & %6.3f \\\\\n"  M[1,1] M[1,2] M[1,3] M[1,4]
        @printf io "\$v\$   &   ---  & %6.3f & %6.3f & %6.3f \\\\\n"        M[2,2] M[2,3] M[2,4]
        @printf io "\$v/u\$ &   ---  &   ---  & %6.3f & %6.3f \\\\\n"              M[3,3] M[3,4]
        @printf io "\$p\$   &   ---  &   ---  &   ---  & %6.3f \\\\\n"                   M[4,4]

        write(io, "\\bottomrule\n")
        write(io, "\\end{tabular}\n")
    end 
    println("Table written to: $path")
end

# 3. Resolve and resimulate
function fnSweepNp(cˢˢ, N⃗ₚ; T = 52 * 100)
    # A. Containers for moments 
    K       = length(N⃗ₚ)
    σs      = zeros(K, 4)         # rows: Nₚ values; cols: u, v, v/u, p
    Ms      = [zeros(4, 4) for _ in 1:K]

    # B. Loop over Nₚ 
    for (k, Nₚ) in enumerate(N⃗ₚ)
        @printf "Solving for Nₚ = %d ...\n" Nₚ

        # Rebuild params and aggregate struct with this Nₚ 
        params  = fnSetUpParameters(Nₚ = Nₚ)
        Agg     = fnSetUpAggregate(params; T = T)

        # Solve policy and simulate 
        fnSolveAggregate!(cˢˢ, params, Agg)
        fnSimulate!(params, Agg)

        # Compute moments 
        moments = fnTable4Moments(Agg)
        σs[k, :]    .= moments.σ 
        Ms[k]       .= moments.M 
    end 

    return (Nₚ = N⃗ₚ, σ = σs, M = Ms)
end 

# 3. Standard deviations
function fnPlotDiscretisationStd(sweep)
    @unpack Nₚ, σ   = sweep 
    labels          = [L"u"  L"v"  L"\frac{v}{u}"  L"p"]

    plt     = plot(; 
                    xlabel      = "Nₚ", 
                    ylabel      = "Standard deviation", 
                    framestyle  = :box, 
                    grid        = true, 
                    gridalpha   = 0.3,
                    legend      = :right)
    for j in 1:4
        plot!(plt, Nₚ, σ[:, j]; 
                    label       = labels[j], 
                    lw          = 2, 
                    marker      = :square, 
                    markersize  = 4)
    end 
    return plt 
end 

# 3. Plot correlations vs Nₚ (6 unique off-diagonal pairs)
function fnPlotDiscretisationCorr(sweep)
    @unpack Nₚ, M   = sweep 
    K               = length(Nₚ)

    # Unique off-diagonal pairs of {u, v, v/u, p}
    pairs   = [(1,2,L"u,\, v"), (1,3,L"u,\, v/u"), (1,4,L"u,\, p"),
           (2,3,L"v,\, v/u"), (2,4,L"v,\, p"), (3,4,L"v/u,\, p")]

    plts    = []
    for (i, j, label) in pairs 
        corrs   = [M[k][i, j] for k in 1:K]
        p       = plot(Nₚ, corrs; 
                        title       = label,
                        lw          = 2, 
                        marker      = :square, 
                        markersize  = 4, 
                        label       = "", 
                        framestyle  = :box, 
                        grid        = true, 
                        gridalpha   = 0.3,
                        color       = :maroon,
                        titlefont   = 12,
                        yformatter  = y -> @sprintf("%.3f", y))
        push!(plts, p)
    end 

    plt     = plot(plts...; 
                    layout  = (2, 3), 
                    size    = (900, 500),
                    xlabel  = "Nₚ",
                    ylabel  = "")
    return plt 
end