# 1. Plot the stationary distribution of the productivity process 
function fnPlotProductivity(params)
    # A. Unpacking business 
    @unpack P, p⃗, Nₚ   = params 

    # B. Compute the stationary distribution 
    eigvals_P, eigvecs_P    = eigen(P')
    i                       = argmin(abs.(eigvals_P .- 1))
    π                       = real.(eigvecs_P[:, i])
    π                       = π ./ sum(π)

    # C. Plot as a discrete PDF 
    plt     = bar(p⃗, π; 
                    xlabel      = "p", 
                    ylabel      = "f(p)",
                    label       = "",
                    framestyle  = :box, 
                    grid        = true, 
                    gridalpha   = 0.3, 
                    color       = :maroon,
                    linecolor   = :match,
                    fillalpha   = 0.6)
    vline!(plt, [params.pₛₛ]; linestyle = :dash, color = :navy, label = "p⋆")
    return plt 
end 

# 2. LaTeX table of equilibrium objects across SS / min / max productivity 
function fnPrintEquilibriumRange(params,SS,Agg; path = joinpath("tables", "equilibrium_range.tex"))
    # A. Find indices of the min and max productivity states 
    imin    = argmin(Agg.p⃗̂)                                # fallback: use simulated path
    imax    = argmax(Agg.p⃗̂)
    # Better: use the p-grid directly, policies are stored on it 
    @unpack p⃗   = params 
    imin    = argmin(p⃗)
    imax    = argmax(p⃗)

    # B. Open the file 
    open(path, "w") do io 
        write(io, "\\centering\n")
        write(io, "\\caption{Equilibrium objects across productivity states}\n")
        write(io, "\\begin{tabular}{llccc}\n")
        write(io, "\\toprule\n")
        write(io, "\\textbf{Description} & \\textbf{Symbol} & \\textbf{SS} & \\textbf{Min \$p\$} & \\textbf{Max \$p\$} \\\\\n")
        write(io, "\\midrule\n")

        # C. Productivity row 
        @printf io "Productivity          & \$p\$       & %8.4f & %8.4f & %8.4f \\\\\n" params.pₛₛ p⃗[imin] p⃗[imax]

        # D. Equilibrium rows 
        @printf io "Labour market tightness & \$\\theta\$ & %8.4f & %8.4f & %8.4f \\\\\n" SS.θ Agg.θ⃗[imin] Agg.θ⃗[imax]
        @printf io "Wage                    & \$w\$       & %8.4f & %8.4f & %8.4f \\\\\n" SS.w Agg.w⃗[imin] Agg.w⃗[imax]
        @printf io "Job-filling rate        & \$q\$       & %8.4f & %8.4f & %8.4f \\\\\n" SS.q Agg.q⃗[imin] Agg.q⃗[imax]
        @printf io "Job-finding rate        & \$f\$       & %8.4f & %8.4f & %8.4f \\\\\n" SS.f Agg.f⃗[imin] Agg.f⃗[imax]
        @printf io "Firm value              & \$J\$       & %8.4f & %8.4f & %8.4f \\\\\n" SS.J Agg.J⃗[imin] Agg.J⃗[imax]
        write(io, "\\bottomrule\n")
        write(io, "\\end{tabular}\n")
    end 
    println("Table written to: $path")
end