##---------------------------------------------------------------
## Figure S3
## Constant storm size, variable storm frequency
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


##---------------------------------------------------------------
## Random interstorm intervals, constant storm size
##---------------------------------------------------------------

##---------------------------------------------------------------
## First for low mean map (all coexisting)
##---------------------------------------------------------------

μ = 0.21
F = 50.0
Random.seed!(3)
spp_data = generate_spp_data(4, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 3.0, 5.0)

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (4.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)
cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 4, rev = true)

P = 5.0
mp = 0.05

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(400, P, 0.1, mp, sd1, true, false)
plot_rainfall_regime(r_1, 20)

out = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics2(out, false, "", cg)

## start dialing up variation in mean annual precip
r_2 = generate_rainfall_regime(400, P, 10.0, mp, sd1, false, true)
plot_rainfall_regime(r_2, 20)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics2(out, false, "", cg)


r_3 = generate_rainfall_regime(400, P, 15.0, mp, sd1, false, true, false)
length(r_3[1])
plot_rainfall_regime(r_3, 20)

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics2(out, false, "", cg)


ymax = maximum(r_3[1]) + (0.1 * maximum(r_3[1]))

rplot_1 = plot_rainfall_regime(r_1, 10, ymax)
rplot_1 = plot!(xlab = "", ylab = "")
dynplot_1 = plot(dynplot_1, ylab = "")
p1 = plot(plot(rplot_1, xlab = "", title = "σ map = " * string(sd1)), plot(dynplot_1, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

rplot_2 = plot_rainfall_regime(r_2, 10, ymax)
rplot_2 = plot!(xlab = "", ylab = "")
dynplot_2 = plot(dynplot_2, ylab = "")
p2 = plot(plot(rplot_2, xlab = "", title = "σ map = " * string(sd2)), plot(dynplot_2, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

rplot_3 = plot_rainfall_regime(r_3, 10, ymax)
rplot_3 = plot!(xlab = "", ylab = "")
dynplot_3 = plot(dynplot_3, ylab = "")
p3 = plot(plot(rplot_3, xlab = "", title = "σ map = " * string(sd3)), plot(dynplot_3, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

panel1 = plot(p1, p2, p3, layout = (1,4), size = (1200, 400))

##---------------------------------------------------------------
## Moderate average rainfall
##---------------------------------------------------------------


P = 5.0
mp = 0.2

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(400, P, 0.1, mp, sd1, true, false)
plot_rainfall_regime(r_1, 10)

out = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics2(out, false, "", cg)

## start dialing up variation in mean annual precip

r_2 = generate_rainfall_regime(400, P, 10.0, mp, sd1, false, true)
plot_rainfall_regime(r_2, 20)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics2(out, false, "", cg)


r_3 = generate_rainfall_regime(400, P, 15.0, mp, sd1, false, true, false, 7000)
plot_rainfall_regime(r_3, 20)
length(r_3[1])

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics2(out, false, "", cg)


ymax = maximum(r_3[1]) + (0.1 * maximum(r_3[1]))

rplot_1 = plot_rainfall_regime(r_1, 10, ymax)
rplot_1 = plot!(xlab = "", ylab = "")
dynplot_1 = plot(dynplot_1, ylab = "")
p1 = plot(plot(rplot_1, xlab = "", title = "σ map = " * string(sd1)), plot(dynplot_1, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

rplot_2 = plot_rainfall_regime(r_2, 10, ymax)
rplot_2 = plot!(xlab = "", ylab = "")
dynplot_2 = plot(dynplot_2, ylab = "")
p2 = plot(plot(rplot_2, xlab = "", title = "σ map = " * string(sd2)), plot(dynplot_2, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

rplot_3 = plot_rainfall_regime(r_3, 10, ymax)
rplot_3 = plot!(xlab = "", ylab = "")
dynplot_3 = plot(dynplot_3, ylab = "")
p3 = plot(plot(rplot_3, xlab = "", title = "σ map = " * string(sd3)), plot(dynplot_3, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

panel2 = plot(p1, p2, p3, layout = (1,4), size = (1200, 400))

##---------------------------------------------------------------
## High average rainfall
##---------------------------------------------------------------

P = 5.0
mp = 0.8

## check that stochastic simulator and standard simulator return same results
r_1 = generate_rainfall_regime(400, P, 0.1, mp, sd1, true, false)
plot_rainfall_regime(r_1, 10)

out = sim_water_ppa_stochastic(spp_data, length(r_1[1]), nrow(spp_data), 1.0, r_1, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)

dynplot_1 = plot_simulation_dynamics2(out, false, "", cg)

## start dialing up variation in mean annual precip

r_2 = generate_rainfall_regime(400, P, 10.0, mp, sd1, false, true)
plot_rainfall_regime(r_2, 20)

out = sim_water_ppa_stochastic(spp_data, length(r_2[1]), nrow(spp_data), 1.0, r_2, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_2 = plot_simulation_dynamics2(out, false, "", cg)


r_3 = generate_rainfall_regime(400, P, 15.0, mp, sd1, false, true, false, 7000)
plot_rainfall_regime(r_3, 20)
length(r_3[1])

out = sim_water_ppa_stochastic(spp_data, length(r_3[1]), nrow(spp_data), 1.0, r_3, F, θ_fc, μ, zeros(1,1),
                               0.4, b, 0.1, true)
dynplot_3 = plot_simulation_dynamics2(out, false, "", cg)


ymax = maximum(r_3[1]) + (0.1 * maximum(r_3[1]))

rplot_1 = plot_rainfall_regime(r_1, 10, ymax)
rplot_1 = plot!(xlab = "", ylab = "")
dynplot_1 = plot(dynplot_1, ylab = "")
p1 = plot(plot(rplot_1, xlab = "", title = "σ map = " * string(sd1)), plot(dynplot_1, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

rplot_2 = plot_rainfall_regime(r_2, 10, ymax)
rplot_2 = plot!(xlab = "", ylab = "")
dynplot_2 = plot(dynplot_2, ylab = "")
p2 = plot(plot(rplot_2, xlab = "", title = "σ map = " * string(sd2)), plot(dynplot_2, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

rplot_3 = plot_rainfall_regime(r_3, 10, ymax)
rplot_3 = plot!(xlab = "", ylab = "")
dynplot_3 = plot(dynplot_3, ylab = "")
p3 = plot(plot(rplot_3, xlab = "", title = "σ map = " * string(sd3)), plot(dynplot_3, colorbar = :none, legend = :none, yab = "", xlab = ""), layout = [1,1])

panel3 = plot(p1, p2, p3, layout = (1,4), size = (1200, 400))

plot(panel1, panel2, panel3, layout = (3,1), size = (1000, 1600))
savefig("../figures/figureS3/figureS3.pdf")
