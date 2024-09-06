##---------------------------------------------------------------
## Figure 2
## Single-species dynamics
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

Random.seed!(10)
spp_data = generate_spp_data(1, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 1.0, 30.0)

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (3.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)

p = [-10:0.01:-1;]
plot(p, log.(w_psi.(p)))


## biomass growth plot

t = [1.0:1.0:100.0;]

g1 = calc_g(0.03, spp_data.Aₘ[1], spp_data.r[1], spp_data.α[1], spp_data.δ[1], 0.1)
g2 = calc_g(0.05, spp_data.Aₘ[1], spp_data.r[1], spp_data.α[1], spp_data.δ[1], 0.1)
g3 = calc_g(0.08, spp_data.Aₘ[1], spp_data.r[1], spp_data.α[1], spp_data.δ[1], 0.1)

b1 = (g1 .* t) .^ b
b2 = (g2.* t) .^ b
b3 = (g3 .* t) .^ b

p0 = plot(t ./ 10.0, b1, linewidth = 5, color = "#7a0177")
p0 = plot!(t ./ 10.0, b2, linewidth = 5, color = "#c51b8a")
p0 = plot!(t ./ 10.0, b3, linewidth = 5, color = "#f768a1",
      frame = :box, grid = :false, widen = :false, legend = :false)


ρ = 0.01
out = sim_water_ppa(spp_data, 200, nrow(spp_data),
                    Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

cg = cgrad(["#fed976", "#bd0026"], 2, rev = false)
p1 = plot_simulation_dynamics(out, false, "", cg, 5)
p1 = plot!(colorbar = :false, legend = :false, grid = :false, ylim = [0.0, 75.0], size = (950,800))
#savefig("../figures/figure1/figure1a.pdf")

p2 = plot_canopy_cover(out, cg, 5)
p2 = plot!(colorbar = :false, legend = :false, widen = :false, grid = :false, ylim = [0.0, 1.0], size = (950,800))
#savefig("../figures/figure1/figure1b.pdf")

## keep biomass plot code in case I need it later
#bdata = copy(out[1])
#plot(bdata.rowkey ./ 10, reverse(bdata[:,2]), linewidth = 4, color = "#fed976")
#plot!(bdata.rowkey ./ 10, reverse(bdata[:,3]), linewidth = 4, color = "#bd0026",
#      frame = :box, grid = :false, xlab = "Time (years)", ylab = "Biomass (kg)")

t = collect(1:length(out[5][1:1000]))
p3 = plot(t ./ 10.0, out[5][1:1000], linewidth = 5,
     frame = :box, grid = :false, legend = :false, color = "#3182bd", widen = :false,
     xlab = "Time (years)", ylab = "Volumetric soil water content", ylim = [0.03, 0.22])
p3 = plot!(t ./ 10.0, out[6][1:1000], linewidth = 5,
     frame = :box, grid = :false, legend = :false, color = "#9ecae1", widen = :false,
     xlab = "Time (years)", ylab = "Volumetric soil water content", size = (950,800))
#savefig("../figures/figure1/figure1c.pdf")

p4 = plot(out[14][:,1], out[14][:,2], linewidth = 5,
     frame = :box, grid = :false, legend = :false, color = "#2ca25f", widen = :false,
     xlab = "Time (years)", ylab = "Shutddown time (t₁)", ylim = [0.04, 0.101], size = (950,800))
#savefig("../figures/figure1/figure1d.pdf")


ρ = 0.03
out = sim_water_ppa(spp_data, 200, nrow(spp_data),
                    Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

cg = cgrad(["#fed976", "#bd0026"], 2, rev = false)
p5 = plot_simulation_dynamics(out, false, "", cg, 5)
p5 = plot!(colorbar = :false, legend = :false, grid = :false, size = (950,800))
#savefig("../figures/figure1/figure1e.pdf")

p6 = plot_canopy_cover(out, cg, 5)
p6 = plot!(colorbar = :false, legend = :false, widen = :false, grid = :false, ylim = [0.0, 1.1], size = (950,800))
#savefig("../figures/figure1/figure1f.pdf")


t = collect(1:length(out[5][1:1000]))
p7 = plot(t ./ 10.0, out[5][1:1000], linewidth = 5,
     frame = :box, grid = :false, legend = :false, color = "#3182bd", widen = :false,
     xlab = "Time (years)", ylab = "Volumetric soil water content", ylim = [0.03, 0.22])
p7 = plot!(t ./ 10.0, out[6][1:1000], linewidth = 5,
     frame = :box, grid = :false, legend = :false, color = "#9ecae1", widen = :false,
     xlab = "Time (years)", ylab = "Volumetric soil water content", size = (950,800))
#savefig("../figures/figure1/figure1g.pdf")

p8 = plot(out[14][:,1], out[14][:,2], linewidth = 5,
     frame = :box, grid = :false, legend = :false, color = "#2ca25f", widen = :false,
     xlab = "Time (years)", ylab = "Shutddown time (t₁)", ylim = [0.04, 0.101], size = (950,800))
#savefig("../figures/figure1/figure1h.pdf")

ps = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout = (2,4), size = (1200, 600))
pss = plot(p0, p0)
plot(pss, ps, layout = (2,1), size = (1200, 1000))
savefig("../figures/figure1/figure1_full.pdf")

## keep biomass plot code in case I need it later
#bdata = copy(out[1])
#plot(bdata.rowkey ./ 10, reverse(bdata[:,2]), linewidth = 4, color = "#fed976")
#plot!(bdata.rowkey ./ 10, reverse(bdata[:,3]), linewidth = 4, color = "#bd0026",
#      frame = :box, grid = :false, xlab = "Time (years)", ylab = "Biomass (kg)")

## create water dynamics data
w_data = Vector{Float64}(undef, 0)
t_data = Vector{Float64}(undef, 0)

t = 0.0
for i in 1:nrow(out[14])

    push!(w_data, out[6][i])
    push!(t_data, t)

    if out[14][i,1] < 0.1
        if out[14][i,2] < 0.1
            t = t + out[14][i,1]
            push!(t_data, t)
            push!(w_data, spp_data[1,:Wᵢ])

            t = t + (out[14][i,2] - out[14][i,1])
            push!(t_data, t)
            push!(w_data, spp_data[1,:Wᵢ])
        else
            t = t + out[14][i,1]
            push!(t_data, t)
            push!(w_data, spp_data[1,:Wᵢ])

            t = t + (0.1 - out[14][i,1])
            push!(t_data, t)
            push!(w_data, out[5][i])
        end
    else
        t = t + 0.1
        push!(t_data, t)
        push!(w_data, out[5][i])
    end

end

plot(t_data[1:2000], w_data[1:2000], linewidth = 0.1,
      frame = :box, grid = :false, legend = :false,
     xlab = "Time (years)", ylab = "Volumetric soil water content")
savefig("../figures/figure1/figure1b.pdf")

plot(t_data[400:415], w_data[400:415], linewidth = 2,
      frame = :box, grid = :false, legend = :false,
     xlab = "Time (years)", ylab = "Volumetric soil water content")
savefig("../figures/figure1/figure1c.pdf")
