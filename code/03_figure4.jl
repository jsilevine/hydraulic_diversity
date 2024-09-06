
##---------------------------------------------------------------
## Figure 4
## Coexistence plots for 20 species communities
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

μ = 0.09
F = 10.0
Random.seed!(1)
spp_data = generate_spp_data(20, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 5.0, 30.0)

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (4.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)


## canopy open only last species
ρ = 0.0015
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
hs = calc_hstar(spp_data, ρ, F, μ, 1.0 / P, uf, E, 3.0, θ_fc)

cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 20, rev = true)

p1 = plot(vcat([0.0], spp_data.τ[15:nrow(spp_data)]) ./ 0.1, vcat([ρ + spp_data.Wᵢ[nrow(spp_data)]], spp_data.Wᵢ[15:nrow(spp_data)]),
     linewidth = 2, linestyle = :dash, color = :black)
p1 = plot!(spp_data.τ ./ 0.1, spp_data.Wᵢ, marker_z = spp_data.spp,
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black)
p1 = plot!([0.0], [ρ + spp_data.Wᵢ[nrow(spp_data)]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.08], xlim = [0.0, 1.0], widen = :false, grid = :false, colorbar = :none)

savefig("../figures/figure4/figure4a.png")

out_01 = sim_water_ppa(spp_data, 600, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

p2 = plot_simulation_dynamics(out_01, false, "", cg, 5)
p2 = plot!(legend = :false, grid = :false, size = (950,800), colorbar = :false,
           ylab = "", xlab = "",  xlim = [0.0, 600.0], xtickfontcolor = :white)

savefig("../figures/figure4/figure4b.png")

p3 = plot_canopy_cover(out_01, cg, 5)
p3 = plot!(legend = :false, widen = :false, grid = :false, ylim = [0.0, 1.1], size = (950,800),
           xlim = [0.0, 600.0])

p3 = plot!(spp_data.Wᵢ, Vector(out_01[8][6000,2:21]), seriestype = :scatter,
           inset = (1, bbox(0.05, -0.05, 0.45, 0.45, :center, :right)),
           subplot = 2, markerstrokewidth = 0, xticks = [0.052, 0.054, 0.056],
           marker_z = spp_data.spp, yticks = [0.0, 0.01, 0.02],
           markercolor = cg, markersize = 8, grid = :false,
           frame = :box, legend = :none, color = :black)
savefig("../figures/figure4/figure4c.png")

#savefig("../figures/figure4/figure4d.pdf")

panel1 = plot(p1, p2, p3, layout = (3,1), size = (400, 900))




## canopy open only last species
ρ = 0.017
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
hs = calc_hstar(spp_data, ρ, F, μ, 1.0 / P, uf, E, 3.0, θ_fc)

cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 20, rev = true)

p1 = plot(vcat([0.0], spp_data.τ[1:nrow(spp_data)]) ./ 0.1, vcat([ρ + spp_data.Wᵢ[nrow(spp_data)]], spp_data.Wᵢ[1:nrow(spp_data)]),
     linewidth = 2, linestyle = :dash, color = :black)
p1 = plot!(spp_data.τ ./ 0.1, spp_data.Wᵢ, marker_z = spp_data.spp,
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black)
p1 = plot!([0.0], [ρ + spp_data.Wᵢ[nrow(spp_data)]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.08], xlim = [0.0, 1.0], widen = :false, grid = :false, colorbar = :none,
           ytickfontcolor = :white)

savefig("../figures/figure4/figure4d.png")

out_17 = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
p2 = plot_simulation_dynamics(out_17, false, "", cg, 5)
p2 = plot!(legend = :false, grid = :false, size = (950,800), colorbar = :false,
           ylab = "", xlab = "",  xlim = [0.0, 400.0], xtickfontcolor = :white)
savefig("../figures/figure4/figure4e.png")

p3 = plot_canopy_cover(out_17, cg, 5)
p3 = plot!(legend = :false, widen = :false, grid = :false, ylim = [0.0, 1.1], size = (950,800),
           xlim = [0.0, 400.0], ytickfontcolor = :white)
p3 = plot!(spp_data.Wᵢ, Vector(out_17[8][4000,2:21]), seriestype = :scatter,
           inset = (1, bbox(0.05, -0.05, 0.45, 0.45, :center, :right)),
           subplot = 2, markerstrokewidth = 0, xticks = [0.052, 0.054, 0.056],
           marker_z = spp_data.spp, yticks = [0.05, 0.1],
           markercolor = cg, markersize = 8, grid = :false,
           frame = :box, legend = :none, color = :black)
savefig("../figures/figure4/figure4f.png")

panel2 = plot(p1, p2, p3, layout = (3,1), size = (400, 900))




## canopy closed
ρ = 0.0185
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
hs = calc_hstar(spp_data[1:17,:], ρ, F, μ, 1.0 / P, uf, E, 3.0, θ_fc)

τ_2 = calc_τ_cc.(hs, 1.0 / P, spp_data.Aₘ[2:nrow(spp_data)], spp_data.r[2:nrow(spp_data)], spp_data.α[2:nrow(spp_data)], spp_data.δ[2:nrow(spp_data)],
              μ, uf)
τ_2 = vcat([(spp_data.Wᵢ[17] + ρ - spp_data.Wᵢ[1]) / E], τ_2)

calc_g.(τ_2, spp_data.Aₘ, spp_data.r, spp_data.α, spp_data.δ, 0.1)


cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 20, rev = true)

p1 = plot(vcat([0.0], spp_data.τ[1:nrow(spp_data)]) ./ 0.1, vcat([ρ + spp_data.Wᵢ[nrow(spp_data)]], spp_data.Wᵢ[1:nrow(spp_data)]),
     linewidth = 2, linestyle = :dash, color = :black)
#p1 = plot!(spp_data.τ ./ 0.1, spp_data.Wᵢ,  marker_z = spp_data.spp,
#     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
#     frame = :box, legend = :none, color = :black, alpha = 0.4)
p1 = plot!(τ_2 ./ 0.1, spp_data.Wᵢ, marker_z = spp_data.spp,
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black)
p1 = plot!([1.0], seriestype = :vline, color = :gray, linestyle = :dot, linewidth = 2)

p1 = plot!([0.0], [ρ + spp_data.Wᵢ[nrow(spp_data)]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.08], xlim = [0.0, 1.0], widen = :false, grid = :false, colorbar = :none, ytickfontcolor = :white)

savefig("../figures/figure4/figure4g.png")

out_08 = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
p2 = plot_simulation_dynamics(out_08, false, "", cg, 5)
p2 = plot!(legend = :false, grid = :false, size = (950,800), colorbar = :false,
           ylab = "", xlab = "",  xlim = [0.0, 400.0], xtickfontcolor = :white)
savefig("../figures/figure4/figure4h.png")

p3 = plot_canopy_cover(out_08, cg, 5)
p3 = plot!(legend = :false, widen = :false, grid = :false, ylim = [0.0, 1.1], size = (950,800),
           xlim = [0.0, 400.0], ytickfontcolor = :white)
p3 = plot!(spp_data.Wᵢ, Vector(out_08[8][4000,2:21]), seriestype = :scatter,
           inset = (1, bbox(0.05, -0.05, 0.45, 0.45, :center, :right)),
           subplot = 2, markerstrokewidth = 0, xticks = [0.052, 0.054, 0.056], yticks = [0.2, 0.6, 1.0],
           marker_z = spp_data.spp,
           markercolor = cg, markersize = 8, grid = :false,
           frame = :box, legend = :none, color = :black,
        ylim = [0.0, 1.0])

savefig("../figures/figure4/figure4i.png")

panel3 = plot(p1, p2, p3, layout = (3,1), size = (400, 900))


ρ = 0.025
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
hs = calc_hstar(spp_data, ρ, F, μ, 1.0 / P, uf, E, 3.0, θ_fc)

cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 20, rev = true)


τ_2 = calc_τ_cc.(hs, 1.0 / P, spp_data.Aₘ[2:nrow(spp_data)], spp_data.r[2:nrow(spp_data)], spp_data.α[2:nrow(spp_data)], spp_data.δ[2:nrow(spp_data)],
              μ, uf)
τ_2 = vcat([(spp_data.Wᵢ[17] + ρ - spp_data.Wᵢ[1]) / E], τ_2)

calc_g.(τ_2, spp_data.Aₘ, spp_data.r, spp_data.α, spp_data.δ, 0.1)


cg = cgrad(["#0868ac", "#7bccc4", "#bae4bc"], 20, rev = true)

p1 = plot(vcat([0.0], τ_2) ./ 0.1, vcat([ρ + spp_data.Wᵢ[nrow(spp_data)]], spp_data.Wᵢ[1:nrow(spp_data)]),
     linewidth = 2, linestyle = :dash, color = :black)
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
p1 = plot!(rectangle(0.8, 0.08, 1.0, 0.0), alpha = 0.3, color = :gray)
p1 = plot!(spp_data.τ ./ 0.1, spp_data.Wᵢ,  marker_z = spp_data.spp,
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black, alpha = 0.4)
p1 = plot!(τ_2 ./ 0.1, spp_data.Wᵢ, marker_z = spp_data.spp,
     markercolor = cg, seriestype = :scatter, markersize = 10, markerstrokewidth = 0,
     frame = :box, legend = :none, color = :black)
p1 = plot!([1.0], seriestype = :vline, color = :gray, linestyle = :dot, linewidth = 2)

p1 = plot!([0.0], [ρ + spp_data.Wᵢ[nrow(spp_data)]], seriestype = :scatter, markersize = 10, color = :black,
      ylim = [0.05, 0.08], xlim = [0.0, 1.8], widen = :false, grid = :false, colorbar = :none,
           ytickfontcolor = :white)


savefig("../figures/figure4/figure4j.png")

out_3 = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
p2 = plot_simulation_dynamics(out_3, false, "", cg, 5)
p2 = plot!(legend = :false, grid = :false, size = (950,800), colorbar = :false,
           ylab = "", xlab = "",  xlim = [0.0, 400.0], xtickfontcolor = :white)
savefig("../figures/figure4/figure4k.png")

p3 = plot_canopy_cover(out_3, cg, 5)
p3 = plot!(legend = :false, widen = :false, grid = :false, ylim = [0.0, 1.1], size = (950,800),
           xlim = [0.0, 400.0], ytickfontcolor = :white)
p3 = plot!(spp_data.Wᵢ, Vector(out_3[8][4000,2:21]), seriestype = :scatter,
           inset = (1, bbox(0.05, -0.05, 0.45, 0.45, :center, :right)),
           subplot = 2, markerstrokewidth = 0, xticks = [0.052, 0.054, 0.056], yticks = [0.1, 0.3],
           marker_z = spp_data.spp,
           markercolor = cg, markersize = 8, grid = :false,
           frame = :box, legend = :none, color = :black)
savefig("../figures/figure4/figure4l.png")


panel4 = plot(p1, p2, p3, layout = (3,1), size = (400, 900))

plot(panel1, panel2, panel3, panel4, layout = (1,4), size = (1600,800))
savefig("../figures/figure4/figure4_full.png")
savefig("../figures/figure4/figure4_full.pdf")
