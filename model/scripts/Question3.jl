# 1. Table of steady-state equilibrium objects 
# 1. LaTeX table of steady-state equilibrium objects 
function fnPrintSteadyState(params,SS,cˢˢ; path = joinpath("tables", "calibration.tex"))
    # A. Unpacking business 
    @unpack b, β, s, γ, l, pₛₛ, uₛₛ     = params 

    # B. Open the file 
    open(path, "w") do io 

        # C. Preamble 
        write(io, "\\begin{table}[H]\n")
        write(io, "\\centering\n")
        write(io, "\\caption{Steady-state calibration}\n")
        write(io, "\\label{tab:calibration}\n")
        write(io, "\\begin{tabular}{llr}\n")
        write(io, "\\toprule\n")
        write(io, "\\textbf{Description} & \\textbf{Description} & \\textbf{Description} \\\\\n")
        write(io, "\\midrule\n")

        # D. Deep parameters 
        write(io, "\\multicolumn{3}{l}{Deep parameters} \\\\\n")
        @printf io "Home production            & \$b\$       & %8.4f \\\\\n"  b 
        @printf io "Discount factor            & \$\\beta\$  & %8.4f \\\\\n"  β 
        @printf io "Separation rate            & \$s\$       & %8.4f \\\\\n"  s 
        @printf io "Bargaining power           & \$\\gamma\$ & %8.4f \\\\\n"  γ 
        @printf io "Matching elasticity        & \$l\$       & %8.4f \\\\\n"  l 
        @printf io "Steady-state productivity  & \$p^\\star\$ & %8.4f \\\\\n" pₛₛ 
        write(io, "\\midrule\n")

        # E. Calibration 
        write(io, "\\multicolumn{3}{l}{Calibration} \\\\\n")
        @printf io "Target unemployment rate   & \$u^\\star\$ & %8.4f \\\\\n" uₛₛ 
        @printf io "Calibrated vacancy cost    & \$c^\\star\$ & %8.4f \\\\\n" cˢˢ 
        write(io, "\\midrule\n")

        # F. Equilibrium objects 
        write(io, "\\multicolumn{3}{l}{Equilibrium objects} \\\\\n")
        @printf io "Labour market tightness    & \$\\theta^\\star\$ & %8.4f \\\\\n" SS.θ 
        @printf io "Wage                       & \$w^\\star\$       & %8.4f \\\\\n" SS.w 
        @printf io "Job-filling rate           & \$q^\\star\$       & %8.4f \\\\\n" SS.q 
        @printf io "Job-finding rate           & \$f^\\star\$       & %8.4f \\\\\n" SS.f 
        @printf io "Unemployment rate          & \$u^\\star\$       & %8.4f \\\\\n" SS.u 
        @printf io "Firm value                 & \$J^\\star\$       & %8.4f \\\\\n" SS.J 

        # G. Closing 
        write(io, "\\bottomrule\n")
        write(io, "\\end{tabular}\n")
        write(io, "\\end{table}\n")
    end 
    println("Table written to: $path")
end

# 2. Plot u(c) with the calibrated solution marked 
function fnPlotCalibration(params,SS,cˢˢ; N = 200, factor = 3.0) 
    # A. Unpacking business 
    @unpack s, pₛₛ, uₛₛ     = params 

    # B. Build a grid of c values around the solution 
    c⃗       = range(cˢˢ / factor, cˢˢ * factor; length = N)

    # C. Compute u(c) along the grid 
    u⃗       = zeros(N) 
    for (i, c) in enumerate(c⃗)
        _, _, _, f  = fnEquilibrium!(pₛₛ, c, params)
        u⃗[i]        = fnSteadyStateU(f, params)
    end 

    # D. Plot 
    plt     = plot(c⃗, u⃗; 
                    xlabel      = "c", 
                    ylabel      = "u",
                    label       = "u(c)",
                    lw          = 2,
                    legend      = :bottomright,
                    framestyle  = :box,
                    grid        = true)
    hline!(plt, [uₛₛ]; linestyle = :dash, color = :red, label = "u⋆ = $(round(uₛₛ, digits=3))")
    vline!(plt, [cˢˢ]; linestyle = :dash, color = :red, label = "c⋆ = $(round(cˢˢ, digits=3))")
    scatter!(plt, [cˢˢ], [uₛₛ]; color = :red, markersize = 6, label = "")

    return plt 
end