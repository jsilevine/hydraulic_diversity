

## convert from volumetric water content to soil water potential
function w_psi(psi, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    wmin + (wmax - wmin) * (psimin / psi) ^ lambda
end

function psi_w(w, wmin = 0.01, wmax = 0.05, psimin = -10, lambda = 1.5)
    psimin / (((w - wmin)/(wmax - wmin))^(1/lambda))
end

function calc_phi(m = 5.6)
    m / ((Cₐ - big_gamma) * (1 + (Dₛ / D₀)))
end

function calc_Aₘ()
    Cᵢ = Cₐ - (1 / calc_phi())
    ((Vmax * Cᵢ) / (ω + Cᵢ))
end

function calc_Wᵢ_first(psiₛ, phi)
    w_psi(((phi * Aₘ * Dₛ) / (γ * k)) + psiₛ)
end


function calc_g(t::Float64, Aₘ::Float64, r::Float64, α::Float64, δ::Float64, T::Float64)
    (t * Aₘ - T * r) / (T * α * δ)
end;

## calculate time at which species first become water limited
function calc_t_first!(t_first::Vector{Float64}, biomass_data::DataFrame,
                 n_data::DataFrame, v::Vector{Float64},
                 Wᵢ_first::Vector{Float64}, W₀::Float64,
                 T::Float64, ht_data::DataFrame, zstar::Float64)

    ## determine the proportion of the canopy occupied by each species
    x = canopy_proportion(ht_data, n_data, biomass_data, zstar)

    ## only perform calculations for spp with Wᵢ < W₀
    g0 = findall(W₀ .> Wᵢ_first)
    t_first[1:g0[1]-1] .= 0.0
    t_first[g0] .= T
    if length(g0) != 0
        for s in g0[1]:length(x)
            if s == g0[1] ## special case for earliest species
                t_first[s] = (W₀ - Wᵢ_first[s]) / (E * sum(x[g0]))
            else
                t_first[s] = ((Wᵢ_first[s-1] - Wᵢ_first[s]) / (E * sum(x[s:length(v)]))) + t_first[s-1]
            end
            if t_first[s] > T
                t_first[s] = T
                break
            end
        end
    end
    nothing
end

function calc_a_nostep(psi, psiₛ)
    function fn(aa)
        ((aa .* (Vmax .- aa)) ./ (((γ .* k) ./ Dₛ) .* (Cₐ .* (Vmax .- aa) .- ω .* aa))) .+ psiₛ .- psi
    end

    nlsolve(fn, [0.01]).zero[1]

end

function calc_g_nostep(t::Vector{Float64}, t_first::Vector{Float64}, Aₘ::Float64, r::Float64,
                       α::Float64, δ::Float64, T::Float64,
                       biomass_data::DataFrame, n_data::DataFrame, v::Vector{Float64},
                       Wᵢ_first::Vector{Float64}, Wᵢ::Vector{Float64}, W₀::Float64,
                       ht_data::DataFrame, zstar::Float64)

    ## determine the proportion of the canopy occupied by each species
    x = canopy_proportion(ht_data, n_data, biomass_data, zstar)
    g = Vector{Float64}(undef, length(x))

    fullin = findall(t_first .< T)
    partin = findall(t_first .< T .&& t .> T)
    g = (Aₘ .- r) ./ (α .* δ)
    if length(fullin) > 0
        for s in 1:fullin[length(fullin)]
            if s in partout
                a = t_first[s] * Aₘ
                w_fin = Wᵢ_first[s] - E*sum(x[s:length(x)])*(T-t_first[s])
                step = (Wᵢ_first[s] - w_fin) / 40
                wlist = [Wᵢ_first[s]:step:w_fin]
                for i in wlist
                    a += calc_a_nostep(psi_w(wlist[i]), psi_w(Wᵢ[s]))
                end
            else
                a = t_first[s] * Aₘ
                step = (Wᵢ_first[s] - Wᵢ[s]) / 40
                tstep = (t[s] - t_first[s]) / 40
                wlist = collect(range(Wᵢ_first[s], Wᵢ[s], 40))
                for i in 1:length(wlist)
                    a += tstep * calc_a_nostep(psi_w(wlist[i]), psi_w(Wᵢ[s]))
                    println(calc_a_nostep(psi_w(wlist[i]), psi_w(Wᵢ[s])))
                end
            end

            g[s] = (a - T * r) / (T * α[s] * δ[s])
            a = 0.0
        end
    end

    return g

end;

##---------------------------------------------------------------
## water and light competition simulator, relaxed step-function approx.
##---------------------------------------------------------------

"""
    iterate_water_ppa!(yr::Int64, biomass_data::DataFrame, biomass_dynamics::DataFrame,
                       n_data::DataFrame, n_dynamics::DataFrame, r_data::DataFrame,
                       w_data::DataFrame, height_data::DataFrame, canopy_dynamics::DataFrame,
                       g::Vector{Float64}, t::Vector{Float64}, v::Vector{Float64}, Nspp::Int64,
                       μ::Float64, F::Float64, mt::Matrix{Float64}, W₀vec::Vector{Float64},
                       Tvec::Vector{Float64}, understory_factor::Float64, θ_fc::Float64, b::Float64)

Iterates through years of a simulation and returns results. Called by simulation wrappers.

"""
function iterate_water_ppa_nostep(Nyr::Int64, spp_data::DataFrame,
                           biomass_data::DataFrame, biomass_dynamics::DataFrame,
                           n_data::DataFrame, n_dynamics::DataFrame,
                           r_data::DataFrame, w_data::Vector{Float64},
                           height_data::DataFrame, canopy_dynamics::DataFrame,
                           g::Vector{Float64}, t::Vector{Float64}, t_first::Vector{Float64},
                           v::Vector{Float64}, Nspp::Int64,
                           μ::Float64 = 0.1, F::Float64 = 10.0, mt::Matrix{Float64} = zeros(1,1),
                           W₀vec::Vector{Float64} = repeat([0.6], Nyr),
                           Tvec::Vector{Float64} = repeat([40.0], inner = Nyr),
                           understory_factor::Float64 = 0.1,
                           θ_fc::Float64 = 0.4, pb::Bool = true, w_init::Float64 = 0.4, b::Float64 = 2.5,
                           perturb::Bool = false, perturb_iter::Int64 = 200, perturb_fac::Float64 = 0.9)

    ## initialize soil water content
    w = minimum([θ_fc, w_init])
    w_in_data = copy(w_data)
    zstar_data = copy(w_data)

    t_data = DataFrame(k = [1:1:Nyr;])
    for i in 1:nrow(spp_data)
        t_data[:,string(spp_data[i,:spp])] .= Vector{Float64}(undef, Nyr)
    end

    ## begin iterating
    if pb
        prog = ProgressBar(total = Nyr)
    end

    zstar = 0.0

    for yr in 1:Nyr

        ## calculate canopy closure height
        zstar = calc_zstar(biomass_data, height_data, n_data)
        zstar_data[yr] = zstar

        ## calculate initial water content
        w = minimum([w + W₀vec[yr], θ_fc])
        w_in_data[yr] = w

        ## calculate growth time for each species
        if any(spp_data.Wᵢ .< w)
            calc_t!(t, biomass_data, n_data, v, spp_data.Wᵢ,
                    w, Tvec[yr], height_data, zstar)
            calc_t_first!(t, biomass_data, n_data, v, spp_data.Wᵢ_first,
                          w, Tvec[yr], height_data, zstar)
        else
            t = repeat([0.0], inner = Nspp)
        end

        t_data[yr, 2:Nspp+1] .= t

        ## recalculate soil water content
        w = maximum([calc_w(w, t, Tvec[yr], biomass_data, n_data, v, spp_data.Wᵢ, height_data, zstar), 0.0])
        w_data[yr] = w

        ## calculate average growth rate for each species
        g = calc_g.(t, spp_data[:, :Aₘ], spp_data[:, :r], spp_data[:,:α], spp_data[:,:δ], Tvec[yr])

        for s in spp_data.spp
            biomass_data[1:yr, s+1] =
                grow(biomass_data[1:yr, s+1], height_data[1:yr, s+1], zstar, b, g[s],
                     g[s] * understory_factor, Tvec[yr])

        end

        die!(n_data, r_data, yr, mt)
        height_data = calc_ht(biomass_data, height_data, spp_data)

        v = ΣB(biomass_data, n_data, v)
        biomass_dynamics[yr, [2:1:Nspp+1;]] = v

        canopy_dynamics[yr,2:Nspp+1] = canopy_proportion(height_data, n_data, biomass_data, zstar)

        ## generate new cohort and record population dynamics
        n_data = birth(n_data, canopy_dynamics, yr, F, Tvec[yr])
        for i in 2:ncol(r_data)
            r_data[yr+1, i] = n_data[yr+1, i]
        end
        Σn!(n_data, n_dynamics, yr)

        ## extinction cutoff
        r_data[:,Matrix(n_dynamics)[yr,:] .< 1e-3] .= 0.0

        if pb
            update(prog)
        end

    end

    return [biomass_data, n_dynamics, n_data, r_data, w_data,
            w_in_data, height_data, canopy_dynamics, zstar_data, g, t, biomass_dynamics, Tvec, t_data]

end;

"""
    sim_water_ppa(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                       Ninit::Any, μ::Float64 = 0.15, F::Float64 = 10.0,
                       P::Float64 = 40.0, mean_p::Float64 = 16.0, θ_fc = 0.4, mt::Matrix{Float64} = zeros(1,1),
                       pb::Bool = true, w_init::Float64 = 0.4, b::Float64 = 2.5, understory_factor::Float64 = 0.1,
                       perturb::Bool = false, perturb_frac::Float64 = 0.5, perturb_factor::Float64 = 0.9,
                       perturb_water::Bool = false, perturb_factor_water::Float64 = 0.7,
                       perturb_water_return::Bool = false, perturb_water_return_frac::Float64 = 0.7)

Perform simulation of a community of species limited by both water and (potentially) light. The community
of species must be described by `spp_data`, a data frame generated by the `generate_spp_data()` function.
"""
function sim_water_ppa_nostep(spp_data::DataFrame, Nyr::Int64, Nspp::Int64,
                       Ninit::Any, μ::Float64 = 0.15, F::Float64 = 10.0,
                       P::Float64 = 40.0, mean_p::Float64 = 16.0, θ_fc = 0.4, mt::Matrix{Float64} = zeros(1,1),
                       pb::Bool = true, w_init::Float64 = 0.4, b::Float64 = 2.5, understory_factor::Float64 = 0.1,
                       perturb::Bool = false, perturb_frac::Float64 = 0.5, perturb_factor::Float64 = 0.9,
                       perturb_water::Bool = false, perturb_factor_water::Float64 = 0.7,
                       perturb_water_return::Bool = false, perturb_water_return_frac::Float64 = 0.7)

    Nyr = Nyr * Int(round(P))

    ## setup for simulations
    ## generate biomass and population data frames
    biomass_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                             spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                             B = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    n_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    r_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                       spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                       N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))
    height_data = DataFrame(rowkey = repeat([1:1:Nyr+1;], inner = Nspp),
                            spp = repeat(string.(spp_data.spp), outer = Nyr+1),
                            N = repeat(repeat(Float64[0], inner = Nspp), inner = Nyr+1))

    biomass_data = unstack(biomass_data, :spp, :B)
    biomass_data[!,[2:ncol(biomass_data);]] .= convert.(Float64,biomass_data[!,[2:ncol(biomass_data);]])
    n_data = unstack(n_data, :spp, :N)
    n_data[!,[2:ncol(n_data);]] .= convert.(Float64,n_data[!,[2:ncol(n_data);]])
    r_data = unstack(r_data, :spp, :N)
    r_data[!,[2:ncol(r_data);]] .= convert.(Float64,r_data[!,[2:ncol(r_data);]])
    height_data = unstack(height_data, :spp, :N)
    height_data[!,[2:ncol(height_data);]] .= convert.(Float64,height_data[!,[2:ncol(height_data);]])

    ## dynamics data frames
    biomass_dynamics = copy(biomass_data[1:Nyr,:])
    n_dynamics = copy(n_data[1:Nyr,:])
    canopy_dynamics = copy(n_data[1:Nyr,:])

    ## inital population:
    if typeof(Ninit) == Float64
        n_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)
        r_data[1, 2:Nspp+1] = repeat([Ninit], inner = Nspp)
    elseif typeof(Ninit) == Vector{Float64}
        n_data[1, 2:Nspp+1] = Ninit
        r_data[1, 2:Nspp+1] = Ninit
    else
        println("Please supply Ninit value that is a Float64 or Vector{Float64}")
    end

    g = Vector{Float64}(undef, Nspp)
    t = Vector{Float64}(undef, Nspp)
    t_first = Vector{Float64}(undef, Nspp)
    v = Vector{Float64}(undef, Nspp)

    w_data = Vector{Float64}(undef, Nyr)

    T = 1.0 / P
    W₀ = mean_p / P

    if mt == zeros(1,1)
        mt = mortality_table(Nyr, μ, repeat([T], inner = Nyr))
    end

    ρ_list = repeat([W₀], inner = Nyr)
    if perturb_water
        ρ_list[Int(round(perturb_frac * Nyr)):Nyr] .= perturb_factor_water * W₀
        if perturb_water_return
            ρ_list[Int(round(perturb_water_return_frac * Nyr)):Nyr] .= W₀
        end
    end

    iterate_water_ppa_nostep(Nyr, spp_data,
                      biomass_data, biomass_dynamics,
                      n_data, n_dynamics,
                      r_data, w_data,
                      height_data, canopy_dynamics,
                      g, t, t_first, v, Nspp, μ, F, mt,
                      ρ_list,
                      repeat([T], inner = Nyr), understory_factor, θ_fc,
                      pb, w_init, b, perturb, Int(round(perturb_frac * Nyr)),
                      perturb_factor)

end;



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
