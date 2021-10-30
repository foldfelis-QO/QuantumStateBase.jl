export
    q_pdf,
    q_pdf!

real_tr_mul(𝐚, 𝐛) = sum(real(𝐚[i, :]' * 𝐛[:, i]) for i in 1:size(𝐚, 1))

"""
    q_pdf(state::AbstractState, θ::Real, x::Real; T=Float64)

Quadrature prabability at point (θ, x)
"""
function q_pdf(T::Type{<:Real}, state::AbstractState, θ::Real, x::Real)
    𝛑̂_res = Matrix{Complex{T}}(undef, state.dim, state.dim)

    return q_pdf!(𝛑̂_res, state, θ, x)
end

q_pdf(state::AbstractState, θ::Real, x::Real) = q_pdf(Float64, state, θ, x)

function q_pdf!(𝛑̂_res::AbstractMatrix, state::StateMatrix, θ::Real, x::Real)
    return real_tr_mul(𝛑̂!(𝛑̂_res, θ, x, dim=state.dim), state.𝛒)
end

function q_pdf!(𝛑̂_res::AbstractMatrix, state::StateVector, args...; kwargs...)
    return q_pdf!(𝛑̂_res, StateMatrix(state), args...; kwargs...)
end

"""
    q_pdf(state::AbstractState, θ::Real, x::Real; T=Float64)

Quadrature prabability at points (θs, xs)
"""
function q_pdf(T::Type{<:Real}, state::AbstractState, θs, xs)
    𝛑̂_res_vec = [Matrix{Complex{T}}(undef, state.dim, state.dim) for _ in 1:Threads.nthreads()]
    𝐩 = Matrix{T}(undef, length(θs), length(xs))

    return q_pdf!(𝛑̂_res_vec, 𝐩, state, θs, xs)
end

q_pdf(state::AbstractState, θs, xs) = q_pdf(Float64, state, θs, xs)

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
