##---------------------------------------------------------------
## Figure 7
## Equilibrium stand structure plots
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

function plot_stand_struct(results::Vector{Any}, ntree_fac::Float64 = 1.0)

    final_pops = collect(results[2][nrow(results[2]), 2:ncol(results[2])])
    pos_spp = findall(final_pops .> 0.1)

    p = plot()
    for s in pos_spp


        hs = results[9][length(results[9])]
        canopy_rows = findall(results[7][:,s+1] .> hs)


        ntree = Int(round(final_pops[s] * ntree_fac))
        pop = DataFrame(cohort = Vector{Int64}(undef, ntree), height = Vector{Float64}(undef, ntree))
        pop.cohort .= sample(collect(canopy_rows), FrequencyWeights(results[3][canopy_rows,s+1]), ntree)

        for i in 1:ntree
            pop[i,:height] = results[7][pop[i,:cohort],s+1]
        end


        if s == pos_spp[1]
            x = rand(Uniform(0,1), 1)
            p = plot([x[1], x[1]], [0.0, pop[1,:height]], color = cg[s],
                     linewidth = 5, legend = :none)
            for i in 2:nrow(pop)
                x = rand(Uniform(0,1), 1)
                p = plot!([x[1], x[1]], [0.0, pop[i,:height]], color = cg[s],
                          linewidth = 5)
            end
        else
            for i in 1:nrow(pop)
                x = rand(Uniform(0,1), 1)
                p = plot!([x[1], x[1]], [0.0, pop[i,:height]], color = cg[s],
                          linewidth = 5)
            end
        end
    end

    return p
end



F = 10.0
μ = 0.09

Random.seed!(1)
spp_data = generate_spp_data(1, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 1.0, 30.0)
#spp_data = generate_spp_data(2, 0.7, 1, 1.0 / P, F, μ,
#                             3.0, 0.4, 0.0, 0.06, 0.03, 1.5, 1.5, 0.5, 4.0)

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (3.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)




## canopy open species 2 wins
cg = cgrad(["#0868ac", "#7bccc4","#bae4bc"], 2, rev = false)

ρ = 0.001
out_1_low = sim_water_ppa(spp_data, 200, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
out_1_low[9][2000] > 0

Random.seed!(6)
p_1_low = plot_stand_struct(out_1_low, 2.0)
p_1_low = plot!(ylim = [0.0, 2.1], xlim = [0.0, 1.0])


## canopy open, coexistence
ρ = 0.0125
out_1_lomid = sim_water_ppa(spp_data, 200, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
out_1_lomid[9][2000] > 0
Random.seed!(6)
p_1_lomid = plot_stand_struct(out_1_lomid)
p_1_lomid = plot!(ylim = [0.0, 2.1], xlim = [0.0, 1.0])


## canopy closed, coexistence
ρ = 0.018
out_1_himid = sim_water_ppa(spp_data, 200, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
Random.seed!(6)
p_1_himid = plot_stand_struct(out_1_himid)
p_1_himid = plot!(ylim = [0.0, 2.1], xlim = [0.0, 1.0])


## canopy closed, species 1 wins
ρ = 0.035
out_1_high = sim_water_ppa(spp_data, 200, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

Random.seed!(6)
p_1_high = plot_stand_struct(out_1_high)
p_1_high = plot!(ylim = [0.0, 2.1])

plot(p_1_low, p_1_lomid, p_1_himid, p_1_high)


## 2 spp

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
cg = cgrad(["#0868ac", "#7bccc4","#bae4bc"], 2, rev = false)

ρ = 0.001
out_2_low = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
out_2_low[9][4000] > 0

Random.seed!(6)
p_2_low = plot_stand_struct(out_2_low, 2.0)
p_2_low = plot!(ylim = [0.0, 2.1], xlim = [0.0, 1.0])


## canopy open, coexistence
ρ = 0.0125
out_2_lomid = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
Random.seed!(6)
p_2_lomid = plot_stand_struct(out_2_lomid)
p_2_lomid = plot!(ylim = [0.0, 2.1], xlim = [0.0, 1.0])


## canopy closed, coexistence
ρ = 0.018
out_2_himid = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
Random.seed!(6)
p_2_himid = plot_stand_struct(out_2_himid)
p_2_himid = plot!(ylim = [0.0, 2.1], xlim = [0.0, 1.0])


## canopy closed, species 1 wins
ρ = 0.035
out_2_high = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

Random.seed!(6)
p_2_high = plot_stand_struct(out_2_high)
p_2_high = plot!(ylim = [0.0, 2.1])

plot(p_2_low, p_2_lomid, p_2_himid, p_2_high)




## 10 spp

F = 10.0
μ = 0.09

cg = cgrad(["#0868ac", "#7bccc4","#bae4bc"], 10, rev = true)
Random.seed!(1)
spp_data = generate_spp_data(10, 0.15, 1, 1.0 / P, F, μ,
                             3.0, 0.15, 0.0, 0.2, 0.03, 1.5, 0.5, 1.0, 30.0)
#spp_data = generate_spp_data(2, 0.7, 1, 1.0 / P, F, μ,
#                             3.0, 0.4, 0.0, 0.06, 0.03, 1.5, 1.5, 0.5, 4.0)

function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

psi = -10.0 .+ (3.5 ./ spp_data.δ) .^ 1
spp_data.Wᵢ .= w_psi.(psi)


ρ = 0.001
out_10_low = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

Random.seed!(6)
p_10_low = plot_stand_struct(out_10_low, 2.0)
p_10_low = plot!(ylim = [0.0, 2.1], xlim = [0.0, 1.0])


ρ = 0.0125
out_10_lomid = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
Random.seed!(6)
p_10_lomid = plot_stand_struct(out_10_lomid)
p_10_lomid = plot!(ylim = [0.0, 2.1], xlim = [0.0, 1.0])


ρ = 0.018
out_10_himid = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)
Random.seed!(6)
p_10_himid = plot_stand_struct(out_10_himid)
p_10_himid = plot!(ylim = [0.0, 2.1], xlim = [0.0, 1.0])


ρ = 0.035
out_10_high = sim_water_ppa(spp_data, 400, nrow(spp_data),
                       Ninit, μ, F, P, ρ*P, θ_fc, zeros(1,1), true, 0.4, 3.0, uf, false)

Random.seed!(6)
p_10_high = plot_stand_struct(out_10_high)
p_10_high = plot!(ylim = [0.0, 2.1])

plot(p_10_low, p_10_lomid, p_10_himid, p_10_high)


plot!(p_1_low, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_1_lomid, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_1_himid, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_1_high, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_2_low, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_2_lomid, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_2_himid, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_2_high, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_10_low, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_10_lomid, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_10_himid, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)
plot!(p_10_high, frame = :box, grid = :none, xtickfontcolor = :white, ytickfontcolor = :white)



plot(p_1_low, p_1_lomid, p_1_himid, p_1_high)
savefig("../figures/figure7/figure7a.pdf")
plot(p_2_low, p_2_lomid, p_2_himid, p_2_high)
savefig("../figures/figure7/figure7b.pdf")
plot(p_10_low, p_10_lomid, p_10_himid, p_10_high)
savefig("../figures/figure7/figure7c.pdf")
