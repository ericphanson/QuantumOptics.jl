using ..ode_dopri, ..metrics, ..operators

import OrdinaryDiffEq, DiffEqCallbacks, StochasticDiffEq

function recast! end

"""
    integrate(tspan::Vector{Float64}, df::Function, x0::Vector{Complex128},
            state::T, dstate::T, fout::Function; kwargs...)

Integrate using OrdinaryDiffEq
"""
function integrate(tspan::Vector{Float64}, df::Function, x0::Vector{Complex128},
            state::T, dstate::T, fout::Function;
            alg::OrdinaryDiffEq.OrdinaryDiffEqAlgorithm = OrdinaryDiffEq.DP5(),
            steady_state = false, tol = 1e-3, save_everystep = false,
            callback = nothing, kwargs...) where T

    function df_(dx::Vector{Complex128}, x::Vector{Complex128}, p, t)
        recast!(x, state)
        recast!(dx, dstate)
        df(t, state, dstate)
        recast!(dstate, dx)
    end
    function fout_(x::Vector{Complex128}, t::Float64, integrator)
        recast!(x, state)
        fout(t, state)
    end

    out_type = pure_inference(fout, Tuple{eltype(tspan),typeof(state)})

    out = DiffEqCallbacks.SavedValues(Float64,out_type)

    scb = DiffEqCallbacks.SavingCallback(fout_,out,saveat=tspan,
                                         save_everystep=save_everystep,
                                         save_start = false)

    prob = OrdinaryDiffEq.ODEProblem{true}(df_, x0,(tspan[1],tspan[end]))

    if steady_state
        affect! = function (integrator)
            !save_everystep && scb.affect!(integrator,true)
            OrdinaryDiffEq.terminate!(integrator)
        end
        _cb = OrdinaryDiffEq.DiscreteCallback(
                                SteadyStateCondtion(copy(state),tol,state),
                                affect!;
                                save_positions = (false,false))
        cb = OrdinaryDiffEq.CallbackSet(_cb,scb)
    else
        cb = scb
    end

    full_cb = OrdinaryDiffEq.CallbackSet(callback,cb)

    sol = OrdinaryDiffEq.solve(
                prob,
                alg;
                reltol = 1.0e-6,
                abstol = 1.0e-8,
                save_everystep = false, save_start = false,
                save_end = false,
                callback=full_cb, kwargs...)
    out.t,out.saveval
end

function integrate(tspan::Vector{Float64}, df::Function, x0::Vector{Complex128},
            state::T, dstate::T, ::Void; kwargs...) where T
    function fout(t::Float64, state::T)
        copy(state)
    end
    integrate(tspan, df, x0, state, dstate, fout; kwargs...)
end

struct SteadyStateCondtion{T,T2,T3}
    rho0::T
    tol::T2
    state::T3
end
function (c::SteadyStateCondtion)(rho,t,integrator)
    timeevolution.recast!(rho,c.state)
    dt = integrator.dt
    drho = metrics.tracedistance(c.rho0, c.state)
    c.rho0.data[:] = c.state.data
    drho/dt < c.tol
end



"""
    integrate_stoch(tspan::Vector{Float64}, df::Function, dg::Vector{Function}, x0::Vector{Complex128},
            state::T, dstate::T, fout::Function; kwargs...)

Integrate scalar noise problem using StochasticDiffEq
"""
function integrate_stoch(tspan::Vector{Float64}, df::Function, dg::Function, x0::Vector{Complex128},
            state::T, dstate::Union{T, Array{T}}, fout::Function;
            alg::StochasticDiffEq.StochasticDiffEqAlgorithm = StochasticDiffEq.RKMil(interpretation=:Stratonovich),
            save_everystep = false,
            callback = nothing, kwargs...) where T

    function df_(dx::Vector{Complex128}, x::Vector{Complex128}, p, t)
        recast!(x, state)
        recast!(dx, dstate)
        df(t, state, dstate)
        recast!(dstate, dx)
    end

    function dg_(dx::Vector{Complex128}, x::Vector{Complex128}, p, t)
        recast!(x, state)
        recast!(dx, dstate)
        dg(t, state, dstate)
        recast!(dstate, dx)
    end

    function fout_(x::Vector{Complex128}, t::Float64, integrator)
        recast!(x, state)
        fout(t, state)
    end

    out_type = pure_inference(fout, Tuple{eltype(tspan),typeof(state)})

    out = DiffEqCallbacks.SavedValues(Float64,out_type)

    scb = DiffEqCallbacks.SavingCallback(fout_,out,saveat=tspan,
                                         save_everystep=save_everystep,
                                         save_start = false)

    prob = StochasticDiffEq.SDEProblem{true}(df_, dg_, x0,(tspan[1],tspan[end]))

    sol = StochasticDiffEq.solve(
                prob,
                alg;
                reltol = 1.0e-6,
                abstol = 1.0e-8,
                save_everystep = false, save_start = false,
                save_end = false,
                callback=scb, kwargs...)

    out.t,out.saveval
end

"""
    integrate_stoch(tspan::Vector{Float64}, df::Function, dg::Vector{Function}, x0::Vector{Complex128},
            state::T, dstate::T, fout::Function, n::Int; kwargs...)

Integrate non-diagonal noise problem using StochasticDiffEq
"""
function integrate_stoch(tspan::Vector{Float64}, df::Function, dg::Function, x0::Vector{Complex128},
            state::T, dstate::Union{T, Array{T}}, fout::Function, n::Int;
            alg::StochasticDiffEq.StochasticDiffEqAlgorithm = StochasticDiffEq.EulerHeun(),
            save_everystep = false,
            callback = nothing, kwargs...) where T

    function df_(dx::Vector{Complex128}, x::Vector{Complex128}, p, t)
        recast!(x, state)
        @inbounds for i=1:size(dstate)[2]
            recast!(dx, dstate[1, i])
            df(t, state, dstate[1, i])
        end
        recast!(dstate, dx)
    end

    function dg_(dx::Array{Complex128, 2}, x::Vector{Complex128}, p, t)
        recast!(x, state)
        @inbounds for i=1:size(dx)[2]
            recast!(dx[:, i], dstate[1, i])
            dg(t, state, dstate[1, i])
        end
        recast!(dstate, dx)
    end

    function fout_(x::Vector{Complex128}, t::Float64, integrator)
        recast!(x, state)
        fout(t, state)
    end

    out_type = pure_inference(fout, Tuple{eltype(tspan),typeof(state)})

    out = DiffEqCallbacks.SavedValues(Float64,out_type)

    scb = DiffEqCallbacks.SavingCallback(fout_,out,saveat=tspan,
                                         save_everystep=save_everystep,
                                         save_start = false)

    prob = StochasticDiffEq.SDEProblem{true}(df_, dg_, x0,(tspan[1],tspan[end]), noise_rate_prototype=Array{Complex128}(length(state), n))

    sol = StochasticDiffEq.solve(
                prob,
                alg;
                reltol = 1.0e-6,
                abstol = 1.0e-8,
                save_everystep = false, save_start = false,
                save_end = false,
                callback=scb, kwargs...)

    out.t,out.saveval
end

"""
    integrate_stoch

Define fout if it was omitted.
"""
function integrate_stoch(tspan::Vector{Float64}, df::Function, dg::Function, x0::Vector{Complex128},
    state::T, dstate::Union{T, Array{T}}, ::Void, n::Int; kwargs...) where T
    function fout(t::Float64, state::T)
        copy(state)
    end

    if n == 1
        integrate_stoch(tspan, df, dg, x0, state, dstate, fout; kwargs...)
    else
        integrate_stoch(tspan, df, dg, x0, state, dstate, fout, n; kwargs...)
    end
end


Base.@pure pure_inference(fout,T) = Core.Inference.return_type(fout, T)
