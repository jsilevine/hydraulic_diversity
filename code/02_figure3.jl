##---------------------------------------------------------------
## Figure 3
## Two species competition and coexistence
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

F = 10.0
μ = 0.09

Random.seed!(1)
spp_data = generate_spp_data(2, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 1.0, 30.0)
#spp_data = generate_spp_data(2, 0.7, 1, 1.0 / P, F, μ,
#                             3.0, 0.4, 0.0, 0.06, 0.03, 1.5, 1.5, 0.5, 4.0)

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (3.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)


## canopy open species 2 wins
ρ = 0.008
cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 2, rev = false)
p1 = plot([0.0, spp_data.τ[1], spp_data.τ[2]] .* 10, [ρ + spp_data.Wᵢ[2], spp_data.Wᵢ[1], spp_data.Wᵢ[2]],
     linewidth = 2, linestyle = :dash, color = :black)
p1 = plot!(spp_data.τ .* 10, spp_data.Wᵢ, marker_z = [2.0, 0.0, 0.0],
     markercolor = cg, seriestype = :scatter, markersize = 10,
     frame = :box, legend = :none, color = :black, xlab = "Break-even time (τᵢ*; proportion of interval)",
     ylab = "Critical soil water content (wᵢ*)", markerstrokewidth = 0)
p1 = plot!([0.0], [ρ + spp_data.Wᵢ[2]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.09], xlim = [0.0, 0.7], widen = :false, grid = :false, colorbar = :none)
p1 = plot!(p1, ylab = "", xlab = "")

#savefig("../figures/figure2/figure2a.pdf")

out_08 = sim_water_ppa(spp_data, 300, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
p2 = plot_simulation_dynamics(out_08, false, "", cg, 5)
p2 = plot!(colorbar = :false, legend = :false, grid = :false, size = (950,800),
           ylab = "", xlab = "", ylim = [0.0, 120.0])
#savefig("../figures/figure2/figure2b.pdf")

p3 = plot_canopy_cover(out_08, cg, 5)
p3 = plot!(colorbar = :false, legend = :false, widen = :false,
           grid = :false, ylim = [0.0, 1.1], size = (950,800))
#savefig("../figures/figure2/figure2c.pdf")

panel1 = plot(p1, p2, p3, layout = (3,1), size = (400, 900))
## canopy open, coexistence
ρ = 0.015
p1 = plot([0.0, spp_data.τ[1], spp_data.τ[2]] .* 10, [ρ + spp_data.Wᵢ[2], spp_data.Wᵢ[1], spp_data.Wᵢ[2]],
     linewidth = 2, linestyle = :dash, color = :black)
p1 = plot!(spp_data.τ .* 10, spp_data.Wᵢ, marker_z = [2.0, 0.0, 0.0],
     markercolor = cg, seriestype = :scatter, markersize = 10,
     frame = :box, legend = :none, color = :black, xlab = "Break-even time (τᵢ*; proportion of interval)",
     ylab = "Critical soil water content (wᵢ*)", markerstrokewidth = 0)
p1 = plot!([0.0], [ρ + spp_data.Wᵢ[2]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.09], xlim = [0.0, 0.7], widen = :false, grid = :false, colorbar = :false)
p1 = plot!(p1, ylab = "", xlab = "")
#savefig("../figures/figure2/figure2d.pdf")

out_23 = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
p2 = plot_simulation_dynamics(out_23, false, "", cg, 5)
p2 = plot!(colorbar = :false, legend = :false, grid = :false, size = (950,800),
           ylab = "", xlab = "", ylim = [0.0, 120.0])
#savefig("../figures/figure2/figure2e.pdf")

p3 = plot_canopy_cover(out_23, cg, 5)
p3 = plot!(colorbar = :false, legend = :false, widen = :false, grid = :false, ylim = [0.0, 1.1], size = (950,800))
#savefig("../figures/figure2/figure2f.pdf")

panel2 = plot(p1, p2, p3, layout = (3,1), size = (400, 900))


## canopy closed, coexistence
ρ = 0.017
p1 = plot([0.0, spp_data.τ[1], spp_data.τ[2]] .* 10, [ρ + spp_data.Wᵢ[2], spp_data.Wᵢ[1], spp_data.Wᵢ[2]],
     linewidth = 2, linestyle = :dash, color = :black)
p1 = plot!(spp_data.τ .* 10, spp_data.Wᵢ, marker_z = [2.0, 0.0, 0.0],
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black, xlab = "Break-even time (τᵢ*; proportion of interval)",
     ylab = "Critical soil water content (wᵢ*)")
p1 = plot!([0.0], [ρ + spp_data.Wᵢ[2]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.09], xlim = [0.0, 0.7], widen = :false, grid = :false, colorbar = :false)
p1 = plot!(p1, ylab = "", xlab = "")
#savefig("../figures/figure2/figure2g.pdf")

out_28 = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
p2 = plot_simulation_dynamics(out_28, false, "", cg, 5)
p2 = plot!(colorbar = :false, legend = :false, grid = :false, size = (950,800),
           ylab = "", xlab = "", ylim = [0.0, 120.0])
#savefig("../figures/figure2/figure2h.pdf")

p3 = plot_canopy_cover(out_28, cg, 5)
p3 = plot!(colorbar = :false, legend = :false, widen = :false, grid = :false, ylim = [0.0, 1.1], size = (950,800))
#savefig("../figures/figure2/figure2i.pdf")

panel3 = plot(p1, p2, p3, layout = (3,1), size = (400, 900))

## canopy closed, species 1 wins
ρ = 0.025
hs = calc_hstar(spp_data, ρ, F, μ, 1.0 / P, uf, E, 3.0, θ_fc)
τ_2 = calc_τ_cc(hs, 1.0 / P, spp_data.Aₘ[2], spp_data.τ[2], spp_data.α[2], spp_data.δ[2],
              μ, uf)

p1 = plot([0.0, ρ / E, τ_2] .* 10, [ρ + spp_data.Wᵢ[2], spp_data.Wᵢ[1], spp_data.Wᵢ[2]],
     linewidth = 2, linestyle = :dash, color = :black)
p1 = plot!(spp_data.τ .* 10, spp_data.Wᵢ, marker_z = [2.0, 0.0, 0.0],
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black, xlab = "Break-even time (τᵢ*; proportion of interval)",
     ylab = "Critical soil water content (wᵢ*)")
p1 = plot!(spp_data.τ .* 10, spp_data.Wᵢ, seriestype = :scatter, markersize = 8,
     frame = :box, legend = :none, color = :white, xlab = "Break-even time (τᵢ*; proportion of interval)",
     ylab = "Critical soil water content (wᵢ*)")
p1 = plot!([ρ / E, τ_2] .* 10, spp_data.Wᵢ, marker_z = [2.0, 0.0, 0.0],
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black, xlab = "Break-even time (τᵢ*; proportion of interval)",
     ylab = "Critical soil water content (wᵢ*)")
p1 = plot!([1.0], seriestype = :vline, color = :gray, linestyle = :dot, linewidth = 2)
p1 = plot!([0.0], [ρ + spp_data.Wᵢ[2]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.09], xlim = [0.0, 1.5], widen = :false, grid = :false, colorbar = :false)
p1 = plot!(p1, ylab = "", xlab = "")
#savefig("../figures/figure2/figure2j.pdf")

out_33 = sim_water_ppa(spp_data, 150, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
p2 = plot_simulation_dynamics(out_33, false, "", cg, 5)
p2 = plot!(colorbar = :false, legend = :false, grid = :false, size = (950,800),
           ylab = "", xlab = "", ylim = [0.0, 120.0], xlim = [0.0, 150.0])
plot!(colorbar = :false)
#savefig("../figures/figure2/figure2k.pdf")

p3 = plot_canopy_cover(out_33, cg, 5)
p3 = plot!(colorbar = :false, legend = :false, widen = :false, grid = :false, ylim = [0.0, 1.1], size = (950,800),
           xlim = [0.0, 150.0])
#savefig("../figures/figure2/figure2l.pdf")

panel4 = plot(p1, p2, p3, layout = (3,1), size = (400, 900))

plot(panel1, panel2, panel3, panel4, layout = (1,4), size = (1600, 800))
savefig("../figures/figure3/figure3_full.pdf")
