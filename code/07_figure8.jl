##---------------------------------------------------------------
## Figure 8
## Stochastic precipitation
##---------------------------------------------------------------

##---------------------------------------------------------------
## Head
##---------------------------------------------------------------

using Pkg

Pkg.activate("../")

using Plots, DataFrames, Distributions,
    SpecialFunctions, NLsolve, CSV, Random,
    ForwardDiff, ProgressBars,
    QuadGK, JuMP, Ipopt, StatsBase

Ninit = 1.0

## ecological parameters
μ::Float64 = 0.11   ## mortality rate
E::Float64 = 0.5   ## evapotranspiration rate
l::Float64 = 1.5   ## leaf area allometric constant
b::Float64 = 3.0   ## biomass allometric constant
F::Float64 = 100.0  ## fecundity per unit biomass
W₀::Float64 = 0.4  ## initial water content (default)
θ_fc::Float64 = 0.4 ## field capacity
P::Float64 = 10.0

## include function headers
include("simulator/utility_functions.jl")
include("simulator/simulation_functions.jl")
include("simulator/eq_functions.jl")
include("simulator/meta_functions.jl")

## set plotting variables
theme(:default)
my_cgrad = cgrad(:roma)

##---------------------------------------------------------------
## Main
##---------------------------------------------------------------

function plot_simulation_dynamics2(results, save::Bool = false, filename = "", cgrad = my_crad, lw = 3.5)
    pdata = DataFrames.stack(results[2])
    pdata.variable = parse.(Int64, string.(pdata.variable))

    t = zeros(length(results[13]))
    t[1] = results[13][1]
    for i in 2:length(results[13])
        t[i] = t[i-1] + results[13][i]
    end
    t = repeat(t, outer = length(unique(pdata.variable)))

    p = plot(t, log.(pdata.value), group = pdata.variable, line_z = pdata.variable,
             ylim = [round(minimum(log.(pdata.value)))-1.0, round(maximum(log.(pdata.value)))+1.0], xlim = [0, maximum(t)],
             seriescolor = cgrad, seriestype = :line,
             legend = :none, colorbar = :right, frame = :box, grid = false, linewidth = lw,
             xlab = "Time (years)", ylab = "Population density")

    pdata.value

    if save
        savefig(p, filename)
    end

    return p
end

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

##---------------------------------------------------------------
## Simulation parameters
##---------------------------------------------------------------

## simulation parameters
Nspp::Int64 = 4; Nyr::Int64 = 400; Ninit::Float64 = 1.0;

## ecological parameters
E::Float64 = 0.5 ## evapotranspiration rate
l::Float64 = 1.5  ## leaf area allometric constant
b::Float64 = 3.0  ## biomass allometric constant
F::Float64 = 10.0  ## fecundity per unit biomass
W₀::Float64 = 0.4 ## initial water content (default)
θ_fc::Float64 = 0.15 ## field capacity
μ::Float64 = 0.11 ## mortality rate
uf::Float64 = 0.1 ## understory growth reduction factor


## species with minimum feasible investment in stem:

μ = 0.21
F = 50.0
Random.seed!(3)
spp_data = generate_spp_data(25, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 3.0, 5.0)


psi = -10.0 .+ (4.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)

cg = cgrad(["#253494", "#41b6c4", "#a1dab4"], 15, rev = false)


P = 5.0
mp = 0.2
sd1 = 0.0
sd2 = 0.1
sd3 = 0.4

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(1500, P, 1000.0, mp, sd1, true, false)

out1 = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics2(out1, false, "", cg)

## start dialing up variation in mean annual precip

r_2 = generate_rainfall_regime(1500, P, 1000.0, mp, sd2, true, false)
plot_rainfall_regime(r_2, missing)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics2(out, false, "", cg)


r_3 = generate_rainfall_regime(1500, P, 1000.0, mp, sd3, true, false)
plot_rainfall_regime(r_3, missing)

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics2(out, false, "", cg)

r_1[1] .= r_1[1] .* 100
r_2[1] .= r_2[1] .* 100
ymax = maximum(r_2[1]) + (0.1 * maximum(r_2[1]))
rplot_1 = plot_rainfall_regime(r_1, missing, ymax)
rplot_1 = plot!(xlab = "", ylab = "")
dynplot_1 = plot(dynplot_1, ylab = "")
p1 = plot(plot(rplot_1, xlab = "", title = "σ map = " * string(sd1)), plot(dynplot_1, colorbar = :none,  legend = :none, yab = "", xlab = ""), layout = [1,1])

rplot_2 = plot_rainfall_regime(r_2, missing, ymax)
rplot_2 = plot!(xlab = "", ylab = "")
dynplot_2 = plot(dynplot_2, ylab = "")
p2 = plot(plot(rplot_2, xlab = "", title = "σ map = " * string(sd2)), plot(dynplot_2, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

rplot_1
savefig(rplot_1, "../figures/figure8/figure8a.png")
rplot_2
savefig(rplot_2, "../figures/figure8/figure8b.png")
dynplot_1
savefig(dynplot_1,"../figures/figure8/figure8c.png")
dynplot_2
savefig(dynplot_2,"../figures/figure8/figure8d.png")
