##---------------------------------------------------------------
## Figure 6
## Misc. eq dynamics plots
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
spp_data = generate_spp_data(100, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 1.0, 30.0, 0.9, 1.1, true)

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (4.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)


ρ = 0.001
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 10.0])

ρ = 0.005
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 10.0])

ρ = 0.01
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 15.0])

ρ = 0.02
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 15.0])

ρ = 0.03
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 25.0])

ρ = 0.04
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 500.0])




μ = 0.09
F = 10.0
Random.seed!(1)
spp_data = generate_spp_data(100, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 1.0, 30.0, 0.9, 1.1)

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (4.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)


ρ = 0.001
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 10.0])

ρ = 0.005
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 10.0])

ρ = 0.01
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)

indx = findall(spp_data.eqN .> 0.0)

p1 = plot(spp_data[indx, :δ], spp_data[indx, :eqN], seriestype = :scatter, ylim = [0.0, 6.0],
     color = :black, markersize = 5, frame = :box, grid = :false, legend = :none)

spp_data[:,:eqN_norm] .= 0.0
for i in 2:nrow(spp_data)-1
    spp_data[i,:eqN_norm] = spp_data[i, :eqN] / (spp_data[i-1, :Wᵢ] - spp_data[i+1, :Wᵢ])
end

#spp_data[nrow(spp_data), :eqN_norm] = spp_data[nrow(spp_data), :eqN] / (spp_data[nrow(spp_data)-1, :Wᵢ] - spp_data[nrow(spp_data), :Wᵢ])

p2 = plot(spp_data[indx, :δ], spp_data[indx, :eqN_norm], seriestype = :scatter,
     color = :black, markersize = 5, frame = :box, grid = :false, legend = :none)

eq_plot = plot(p1, p2, layout = (2,1), size = (450, 800))

savefig("../figures/figure6/eq_pop_norm.pdf")


ρ = 0.02
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 15.0])

ρ = 0.03
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 25.0])

ρ = 0.04
calc_eqN(spp_data, 1 / P, ρ, F, E, θ_fc, μ, false, false, uf, 3.0)
plot(spp_data.δ, spp_data.eqN, seriestype = :scatter, ylim = [0.0, 500.0])





##---------------------------------------------------------------
## h* as a function of rainfall
##---------------------------------------------------------------

Random.seed!(1)
spp_data = generate_spp_data(100, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 1.5, 30.0, 0.9, 1.1, true)
spp_data[1,:δ] = 1.5

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (4.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)


rho_list = [0.001:0.001:0.12;]
hs_list = Vector{Float64}(undef, length(rho_list))
for i in 1:length(rho_list)
    hs_list[i] = calc_hstar(spp_data, rho_list[i], F, μ, 1 / 4, uf, E)
end

p1 = plot(rho_list, hs_list, linewidth = 5, color = "#54278f",
      frame = :box, grid = :false, widen = :false, legend = :false, ylim = [0.0, 0.9])

rho_list = [0.001:0.001:0.12;]
hs_list = Vector{Float64}(undef, length(rho_list))
for i in 1:length(rho_list)
    hs_list[i] = calc_hstar(spp_data, rho_list[i], F, μ, 1 / 10, uf, E)
end

p1 = plot!(rho_list, hs_list, linewidth = 5, color = "#756bb1",
      frame = :box, grid = :false, widen = :false, legend = :false, ylim = [0.0, 0.9])

rho_list = [0.001:0.001:0.12;]
hs_list = Vector{Float64}(undef, length(rho_list))
for i in 1:length(rho_list)
    hs_list[i] = calc_hstar(spp_data, rho_list[i], F, μ, 1 / 20, uf, E)
end

p1 = plot!(rho_list, hs_list, linewidth = 5, color = "#9e9ac8",
      frame = :box, grid = :false, widen = :false, legend = :false, ylim = [0.0, 0.9])

rho_list = [0.001:0.001:0.12;]
hs_list = Vector{Float64}(undef, length(rho_list))
for i in 1:length(rho_list)
    hs_list[i] = calc_hstar(spp_data, rho_list[i], F, μ, 1 / 40, uf, E)
end

p1 = plot!(rho_list, hs_list, linewidth = 5, color = "#cbc9e2",
      frame = :box, grid = :false, widen = :false, legend = :false, ylim = [0.0, 0.9])


og_alpha = spp_data.α[1]
alpha_list = [0.2:0.2:15.0;]
hs_list_alpha = Vector{Float64}(undef, length(alpha_list))
for i in 1:length(alpha_list)
    spp_data.α[1] = alpha_list[i]
    hs_list_alpha[i] = calc_hstar(spp_data, 0.04, F, μ, 1 / 10, uf, E)
end
spp_data.α[1] = og_alpha
hs_list_P[hs_list_P .== 0.1] .= 0.0

p2 = plot(alpha_list, hs_list_alpha, linewidth = 5, color = "#016c59",
      frame = :box, grid = :false, widen = :false, legend = :false, ylim = [0.0, 0.9])


og_alpha = spp_data.α[1]
alpha_list = [0.2:0.2:15.0;]
hs_list_alpha = Vector{Float64}(undef, length(alpha_list))
for i in 1:length(alpha_list)
    spp_data.α[1] = alpha_list[i]
    hs_list_alpha[i] = calc_hstar(spp_data, 0.05, F, μ, 1 / 10, uf, E)
end
spp_data.α[1] = og_alpha
hs_list_P[hs_list_P .== 0.1] .= 0.0

p2 = plot!(alpha_list, hs_list_alpha, linewidth = 5, color = "#1c9099",
      frame = :box, grid = :false, widen = :false, legend = :false, ylim = [0.0, 0.9])


og_alpha = spp_data.α[1]
alpha_list = [0.2:0.2:15.0;]
hs_list_alpha = Vector{Float64}(undef, length(alpha_list))
for i in 1:length(alpha_list)
    spp_data.α[1] = alpha_list[i]
    hs_list_alpha[i] = calc_hstar(spp_data, 0.06, F, μ, 1 / 10, uf, E)
end
spp_data.α[1] = og_alpha
hs_list_P[hs_list_P .== 0.1] .= 0.0

p2 = plot!(alpha_list, hs_list_alpha, linewidth = 5, color = "#67a9cf",
      frame = :box, grid = :false, widen = :false, legend = :false, ylim = [0.0, 0.9])
p2 = plot!(ylim = [0.0, 1.4])

hs_plot = plot(p1, p2, layout = (2,1), size = (450, 800))


savefig("../figures/figure6/hstar.pdf")


##---------------------------------------------------------------
## break-even time as function of alpha
##---------------------------------------------------------------

alpha_list = [0.5:0.2:15.0;]
tau_list = Vector{Float64}(undef, length(alpha_list))
for i in 1:length(alpha_list)
    tau_list[i] = calc_τ(spp_data.Aₘ[1], spp_data.r[1], alpha_list[i], spp_data.δ[1],F,
                         μ, 0.1, 3.0)
end
p1 = plot(alpha_list, tau_list ./ 0.1, linewidth = 5, color = "#980043",
      frame = :box, grid = :false, widen = :false, legend = :false)


hs = 0.3
alpha_list = [0.5:0.2:15.0;]
tau_list = Vector{Float64}(undef, length(alpha_list))
for i in 1:length(alpha_list)
    tau_list[i] = calc_τ_cc(hs, 0.1, spp_data.Aₘ[1], spp_data.r[1], alpha_list[i], spp_data.δ[1],
                            μ, uf)
end
p1 = plot!(alpha_list, tau_list ./ 0.1, linewidth = 5, color = "#dd1c77",
      frame = :box, grid = :false, widen = :false, legend = :false)

hs = 0.5
alpha_list = [0.5:0.2:15.0;]
tau_list = Vector{Float64}(undef, length(alpha_list))
for i in 1:length(alpha_list)
    tau_list[i] = calc_τ_cc(hs, 0.1, spp_data.Aₘ[1], spp_data.r[1], alpha_list[i], spp_data.δ[1],
                            μ, uf)
end
p1 = plot!(alpha_list, tau_list ./ 0.1, linewidth = 5, color = "#df65b0",
      frame = :box, grid = :false, widen = :false, legend = :false)

hs = 0.8
alpha_list = [0.5:0.2:15.0;]
tau_list = Vector{Float64}(undef, length(alpha_list))
for i in 1:length(alpha_list)
    tau_list[i] = calc_τ_cc(hs, 0.1, spp_data.Aₘ[1], spp_data.r[1], alpha_list[i], spp_data.δ[1],
                            μ, uf)
end
p1 = plot!(alpha_list, tau_list ./ 0.1, linewidth = 5, color = "#d7b5d8",
      frame = :box, grid = :false, widen = :false, legend = :false)


##---------------------------------------------------------------
## optimum height investment
##---------------------------------------------------------------

using JuMP

using Ipopt
using Plots

## constant definition
hs = 0.17
b = 3.0
δ = 1.0
F = 100.0
μ = 0.11
uf = 0.1


## species with minimum feasible investment in stem:

μ = 0.09
Random.seed!(1)
spp_data = generate_spp_data(100, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 1.5, 30.0, 0.9, 1.1, true)
spp_data[1,:δ] = 1.5

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (4.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)



function opt_alpha(δ, hs = 0.17, F = 100.0, μ = 0.11, uf = 0.1, b = 3.0)

    ## model definition
    model = Model(Ipopt.Optimizer)

    @variable(model, a >= 0)
    @variable(model, 0.1 >= g >= 0)

    @objective(model, Min, g)

    @NLconstraint(model,
                (F * g^(b-1) * gamma(b)) / (μ^b * a * δ) *
                    exp(-hs^2 * δ / a / g * (μ * (1 - uf) / uf)) == 1)

    ## model execution
    optimize!(model);
    value(g)
    return value(a)

end

opt_alpha(1.0, 0.08)


hs = calc_hstar(spp_data, 0.04, F, μ, 1 / P, uf, E)

dlist = [1.5:0.2:10.0;]
opt_a = Vector{Float64}(undef, length(dlist))

for i in 1:length(dlist)
    opt_a[i] = opt_alpha(dlist[i], hs)
end

p2 = plot(dlist[opt_a .< 100.0], opt_a[opt_a .< 100.0], linewidth = 5, color = "#006d2c")


hs = calc_hstar(spp_data, 0.05, F, μ, 1 / P, uf, E)
dlist = [1.5:0.2:10.0;]
opt_a = Vector{Float64}(undef, length(dlist))

for i in 1:length(dlist)
    opt_a[i] = opt_alpha(dlist[i], hs, F)
end

p2 = plot!(dlist[opt_a .< 100.0], opt_a[opt_a .< 100.0], linewidth = 5, color = "#2ca25f")


hs = calc_hstar(spp_data, 0.06, F, μ, 1 / P, uf, E)
dlist = [1.5:0.2:10.0;]
opt_a = Vector{Float64}(undef, length(dlist))

for i in 1:length(dlist)
    opt_a[i] = opt_alpha(dlist[i], hs, F)
end

p2 = plot!(dlist[opt_a .< 100.0], opt_a[opt_a .< 100.0], linewidth = 5, color = "#66c2a4")
p2 = plot!(frame = :box, grid = :false, widen = :false, legend = :false)

alpha_plot = plot(p1, p2, layout = (2,1), size = (450, 800))

plot(eq_plot, hs_plot, alpha_plot, layout = (1,3), size = (1200, 700))
savefig("../figures/figure6/full_plot.pdf")
savefig("../figures/figure6/full_plot.svg")
