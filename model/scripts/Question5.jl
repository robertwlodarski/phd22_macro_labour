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
        write(io, "\\begin{table}[H]\n")
        write(io, "\\centering\n")
        write(io, "\\caption{$caption}\n")
        write(io, "\\label{tab:hm_table4}\n")
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
        write(io, "\\end{table}\n")
    end 
    println("Table written to: $path")
end