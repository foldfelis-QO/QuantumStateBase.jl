export
    q_pdf,
    q_pdf!

real_tr_mul(𝐚, 𝐛) = sum(real(𝐚[i, :]' * 𝐛[:, i]) for i in 1:size(𝐚, 1))

"""
    q_pdf(state::AbstractState, θ::Real, x::Real; T=Float64)

Quadrature prabability at point (θ, x)
"""
function q_pdf(state::AbstractState, θ::Real, x::Real; T=Float64)
    𝛑̂_res = Matrix{complex(T)}(undef, state.dim, state.dim)

    return q_pdf!(𝛑̂_res, state, θ, x)
end

function q_pdf!(𝛑̂_res::Matrix{Complex{T}}, state::StateMatrix, θ::Real, x::Real) where {T}
    return real_tr_mul(𝛑̂!(𝛑̂_res, T(θ), T(x), dim=state.dim), state.𝛒)
end

function q_pdf!(𝛑̂_res::Matrix{Complex{T}}, state::StateVector, args...; kwargs...) where {T}
    return q_pdf!(𝛑̂_res, StateMatrix(state), args...; kwargs...)
end

"""
    q_pdf(state::AbstractState, θ::Real, x::Real; T=Float64)

Quadrature prabability at points (θs, xs)
"""
function q_pdf(state::AbstractState, θs, xs; T=Float64)
    𝛑̂_res_vec = [Matrix{complex(T)}(undef, state.dim, state.dim) for _ in 1:Threads.nthreads()]
    𝐩 = Matrix{T}(undef, length(θs), length(xs))

    return q_pdf!(𝛑̂_res_vec, 𝐩, state, θs, xs)
end

function q_pdf!(𝛑̂_res_vec::Vector{Matrix{Complex{T}}}, 𝐩::Matrix{T}, state::StateMatrix, θs, xs) where {T}
    @sync for (j, x) in enumerate(xs)
        for (i, θ) in enumerate(θs)
            Threads.@spawn 𝐩[i, j] = q_pdf!(𝛑̂_res_vec[Threads.threadid()], state, θ, x)
        end
    end

    return 𝐩
end

function q_pdf!(𝛑̂_res_vec::Vector{Matrix{Complex{T}}}, 𝐩::Matrix{T}, state::StateVector, args...; kwargs...) where {T}
    q_pdf!(𝛑̂_res_vec, 𝐩, StateMatrix(state), args...; kwargs...)
end
