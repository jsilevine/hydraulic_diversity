##---------------------------------------------------------------
## Figure 5
## Coexistence algorithm plots
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

P = 10.0
μ = 0.09
F = 10.0
Random.seed!(3)
spp_data = generate_spp_data(10, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 5.0, 30.0)

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (4.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi) .+ rand(Normal(0.0, 0.0005), 10)
sort!(spp_data, :Wᵢ, rev = true)
spp_data.spp .= 1:10

## canopy open only last species
ρ = 0.016
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
hs = calc_hstar(spp_data, ρ, F, μ, 1.0 / P, uf, E, 3.0, θ_fc)

indx = findall(spp_data.eqN .> 0.0)
cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 20, rev = true)

p1 = plot(vcat([0.0], spp_data.τ[indx]) ./ 0.1, vcat([ρ + spp_data.Wᵢ[nrow(spp_data)]], spp_data.Wᵢ[indx]),
     linewidth = 2, linestyle = :dash, color = :black)
p1 = plot!(spp_data.τ ./ 0.1, spp_data.Wᵢ, marker_z = spp_data.spp,
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black)
p1 = plot!([0.0], [ρ + spp_data.Wᵢ[nrow(spp_data)]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.08], xlim = [0.0, 2.0], widen = :false, grid = :false, colorbar = :none)

#savefig("../figures/figure3/figure3a.pdf")

out_16 = sim_water_ppa(spp_data, 600, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

p2 = plot_simulation_dynamics(out_16, false, "", cg, 5)
p2 = plot!(legend = :false, grid = :false, size = (950,800), colorbar = :false,
           ylab = "", xlab = "",  xlim = [0.0, 400.0])
#savefig("../figures/figure3/figure3b.pdf")
panel1 = plot(p1, p2, layout = (3,1), size = (400, 600))



## canopy open only last species
ρ = 0.028
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
hs = calc_hstar(spp_data[1:3,:], ρ, F, μ, 1.0 / P, uf, E, 3.0, θ_fc)

indx = findall(spp_data.eqN .> 0.0)
cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 20, rev = true)


τ_2 = calc_τ_cc.(hs, 1.0 / P, spp_data.Aₘ[2:nrow(spp_data)], spp_data.r[2:nrow(spp_data)], spp_data.α[2:nrow(spp_data)], spp_data.δ[2:nrow(spp_data)],
              μ, uf)
τ_2 = vcat([(spp_data.Wᵢ[3] + ρ - spp_data.Wᵢ[1]) / E], τ_2)

calc_g.(τ_2, spp_data.Aₘ, spp_data.r, spp_data.α, spp_data.δ, 0.1)

cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 20, rev = true)

p1 = plot(vcat([0.0], τ_2[indx]) ./ 0.1, vcat([ρ + spp_data.Wᵢ[nrow(spp_data)]], spp_data.Wᵢ[indx]),
     linewidth = 2, linestyle = :dash, color = :black)
#p1 = plot!(spp_data.τ ./ 0.1, spp_data.Wᵢ,  marker_z = spp_data.spp,
#     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
#     frame = :box, legend = :none, color = :black, alpha = 0.4)
p1 = plot!(τ_2 ./ 0.1, spp_data.Wᵢ, marker_z = spp_data.spp,
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black)
p1 = plot!([1.0], seriestype = :vline, color = :gray, linestyle = :dot, linewidth = 2)

p1 = plot!([0.0], [ρ + spp_data.Wᵢ[nrow(spp_data)]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.08], xlim = [0.0, 2.0], widen = :false, grid = :false, colorbar = :none)

#savefig("../figures/figure3/figure3a.pdf")

out_28 = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
p2 = plot_simulation_dynamics(out_28, false, "", cg, 5)
p2 = plot!(legend = :false, grid = :false, size = (950,800), colorbar = :false,
           ylab = "", xlab = "",  xlim = [0.0, 400.0])
#savefig("../figures/figure3/figure3b.pdf")
panel2 = plot(p1, p2, layout = (3,1), size = (400, 600))

plot(panel1, panel2, layout = (1,2), size = (800, 800))
savefig("../figures/figure5/figure5_full.pdf")
