export
    gaussian_state_sampler,
    gaussian_state_sampler!

###########################
# gaussian data generator #
###########################

"""
    gaussian_state_sampler(state::AbstractState, n::Integer; θ_offset=0.)

Random points sampled from quadrature probability density function of Gaussian `state`.

* `state`: Quantum state.
* `n`: N points.
* `θ_offset`: The offset of the θ coordinate
"""
function gaussian_state_sampler(T::Type{<:Real}, state::AbstractState, n::Integer; θ_offset=zero(T))
    points = Matrix{T}(undef, 2, n)

    return gaussian_state_sampler!(points, state, θ_offset)[1]
end

gaussian_state_sampler(state::AbstractState, n::Integer; θ_offset=0.) = gaussian_state_sampler(Float64, state, n, θ_offset=θ_offset)

function gaussian_state_sampler!(
    points::AbstractMatrix{T},
    state::StateMatrix, θ_offset::T
) where {T}
    n = size(points, 2)

    # θs
    points[1, :] .= sort!(2π*rand(n) .+ θ_offset)

    # μ and σ given θ
    μ = π̂ₓ_μ(points[1, :], state)
    σ = real(sqrt.(π̂ₓ²_μ(points[1, :], state) - μ.^2))

    # xs
    points[2, :] .= real(μ) + σ .* randn(n)

    return points, real(μ), σ
end

function gaussian_state_sampler!(
    points::AbstractMatrix{T},
    state::StateVector, θ_offset::T
) where {T}
    return gaussian_state_sampler!(points, StateMatrix(state), θ_offset)
end
