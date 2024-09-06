using Base: attr_list

##---------------------------------------------------------------
## Figure S4
## Simulations without the step-function approximation
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
Aₘ::Float64 = 0.2
r::Float64 = 0.03
k::Float64 = 0.01
big_gamma::Float64 = 50
γ::Float64 = 0.5
Dₛ::Float64 = 1.5
D₀::Float64 = 1.5
Cₐ::Float64 = 450.0
Vmax::Float64 = 0.55

## include function headers
include("simulator/utility_functions.jl")
include("simulator/simulation_functions.jl")
include("simulator/eq_functions.jl")
include("simulator/meta_functions.jl")
include("simulator/nostep_sim_functions.jl")

## set plotting variables
theme(:default)
my_cgrad = cgrad(["#253494", "#41b6c4", "#a1dab4"], 10, rev = true)
F = 10.0
μ = 0.09

Random.seed!(2)
spp_data = generate_spp_data(10, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 1.0, 30.0)


psi = -10.0 .+ (3.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)
spp_data.Wᵢ_first .= calc_Wᵢ_first.(psi, calc_phi())


ρ = 0.017
k = 0.05
include("nostep_sim_functions.jl")
spp_data.Wᵢ_first .= calc_Wᵢ_first.(psi, calc_phi())
## high k value
outstep = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
pstep = plot_simulation_dynamics2(outstep, false, "", my_cgrad)

outnostep = sim_water_ppa_nostep(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
pnostep1 = plot_simulation_dynamics2(outnostep, false, "", my_cgrad)


## low k value
k = 0.01
include("nostep_sim_functions.jl")
spp_data.Wᵢ_first .= calc_Wᵢ_first.(psi, calc_phi())

outnostep = sim_water_ppa_nostep(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
pnostep2 = plot_simulation_dynamics2(outnostep, false, "", my_cgrad)

## extreme k value
k = 0.003
include("nostep_sim_functions.jl")
spp_data.Wᵢ_first .= calc_Wᵢ_first.(psi, calc_phi())

outnostep = sim_water_ppa_nostep(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
pnostep3 = plot_simulation_dynamics2(outnostep, false, "", my_cgrad)

plot(plot(pstep, colorbar = :false), plot(pnostep1, ylab = "", colorbar = :false),
     plot(pnostep2, colorbar = :false), plot(pnostep3, ylab = "", colorbar = :false), layout = (2,2))
savefig("../figures/figureS4/figureS4full.pdf")
