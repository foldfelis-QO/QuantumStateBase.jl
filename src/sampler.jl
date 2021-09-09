using KernelDensity

export
    gaussian_state_sampler,
    gaussian_state_sampler!,
    state_sampler,
    state_sampler!,
    IsGaussian

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
function gaussian_state_sampler(state::AbstractState, n::Integer; θ_offset=0.)
    points = Matrix{Float64}(undef, 2, n)

    return gaussian_state_sampler!(points, state, θ_offset)
end

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

##############################
# nongaussian data generator #
##############################

function ranged_rand(n, range::Tuple{T, T}) where {T <: Number}
    return range[1] .+ (range[2]-range[1]) * rand(T, n)
end

function ranged_rand(range::Tuple{T, T}) where {T <: Number}
    return range[1] + (range[2]-range[1]) * rand(T)
end

"""
    state_sampler(
        state::AbstractState, n::Integer;
        c=1, θ_range=(0., 2π), x_range=(-10., 10.),
    )

Random points sampled from quadrature probability density function of Gaussian `state`.

* `state`: Quantum state.
* `n`: N points.
* `θ_range`: Sampling range of θ.
* `x_range`: Sampling range of x.
"""
function state_sampler(
    state::AbstractState, n::Integer;
    c=1, θ_range=(0., 2π), x_range=(-10., 10.),
)
    sampled_points = Matrix{Float64}(undef, 2, n)
    𝛑̂_res_vec = [
        Matrix{complex(Float64)}(undef, state.dim, state.dim)
        for _ in 1:Threads.nthreads()
    ]

    return state_sampler!(sampled_points, 𝛑̂_res_vec, state, c, θ_range, x_range)
end

# rejection method: \frac{f(x)}{c * g(x)} ≥ u
# here, say g(x) is a uniform distribution, and c = 1
function state_sampler!(
    sampled_points::AbstractMatrix{T}, 𝛑̂_res_vec::Vector{Matrix{Complex{T}}},
    state::StateMatrix, c::Real, θ_range, x_range
) where {T}
    n = size(sampled_points, 2)
    sampled_points[1, :] .= sort!(ranged_rand(n, θ_range))

    Threads.@threads for i in 1:n
        p = (thread_id, x) -> q_pdf!(𝛑̂_res_vec[thread_id], state, sampled_points[1, i], x)

        sampled_points[2, i] = ranged_rand(x_range)
        while (p(Threads.threadid(), sampled_points[2, i]) < c*rand())
            sampled_points[2, i] = ranged_rand(x_range)
        end
    end

    return sampled_points
end

function state_sampler!(
    sampled_points::AbstractMatrix{T}, 𝛑̂_res_vec::Vector{Matrix{Complex{T}}},
    state::StateVector,
    args...
) where {T}
    return state_sampler!(sampled_points, 𝛑̂_res_vec, StateMatrix(state), args...)
end

############
# wrapping #
############

struct IsGaussian end

"""
    Base.rand(state::AbstractState, n::Integer, ::Type{IsGaussian}; kwargs...)

Random points sampled from quadrature probability density function of Gaussian `state`.

* `state`: Quantum state.
* `n`: n points.
* `IsGaussian`: To declare `state` is a Gaussian state.
* `kwargs...`: see `gaussian_state_sampler`
"""
function Base.rand(state::StateMatrix, n::Integer, ::Type{IsGaussian}; kwargs...)
    return gaussian_state_sampler(state, n; kwargs...)
end

"""
    Base.rand(state::AbstractState, ::Type{IsGaussian}; kwargs...)

One random point sampled from quadrature probability density function of Gaussian `state`.

* `state`: Quantum state.
* `IsGaussian`: To declare `state` is a Gaussian state.
* `kwargs...`: see `gaussian_state_sampler`
"""
function Base.rand(state::StateMatrix, ::Type{IsGaussian}; kwargs...)
    return rand(state, 1, IsGaussian; kwargs...)
end

function Base.rand(state::StateVector, n::Integer, ::Type{IsGaussian}; kwargs...)
    return rand(StateMatrix(state), n, IsGaussian; kwargs...)
end

function Base.rand(state::StateVector, ::Type{IsGaussian}; kwargs...)
    return rand(state, 1, IsGaussian; kwargs...)
end

"""
    Base.rand(state::AbstractState, n::Integer; kwargs...)

Random points sampled from quadrature probability density function of `state`.

* `state`: Quantum state.
* `n`: n points.
* `kwargs...`: see `state_sampler`
"""
function Base.rand(state::StateMatrix, n::Integer; kwargs...)
    return state_sampler(state, n; kwargs...)
end

"""
    Base.rand(state::AbstractState, n::Integer; kwargs...)

One random point sampled from quadrature probability density function of `state`.

* `state`: Quantum state.
* `n`: n points.
* `kwargs...`: see `state_sampler`
"""
function Base.rand(state::StateMatrix; kwargs...)
    return rand(state, 1; kwargs...)
end

function Base.rand(state::StateVector, n::Integer; kwargs...)
    return rand(StateMatrix(state), n; kwargs...)
end

function Base.rand(state::StateVector; kwargs...)
    return rand(state, 1; kwargs...)
end
