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
    gaussian_state_sampler(state::AbstractState, n::Integer; bias_phase=0.)

Random points sampled from quadrature probability density function of Gaussian `state`.

* `state`: Quantum state.
* `n`: N points.
* `bias_phase`: The offset of the θ coordinate
"""
function gaussian_state_sampler(state::AbstractState, n::Integer; bias_phase=0.)
    points = Matrix{Float64}(undef, 2, n)

    return gaussian_state_sampler!(points, state, bias_phase)
end

function gaussian_state_sampler!(
    points::Matrix{T},
    state::StateMatrix, bias_phase::T
) where {T}
    n = size(points, 2)

    # θs
    points[1, :] .= sort!(2π*rand(n) .+ bias_phase)

    # μ and σ given θ
    μ = π̂ₓ_μ(points[1, :], state)
    σ = real(sqrt.(π̂ₓ²_μ(points[1, :], state) - μ.^2))

    # xs
    points[2, :] .= real(μ) + σ .* randn(n)

    return points
end

function gaussian_state_sampler!(
    points::Matrix{T},
    state::StateVector, bias_phase::T
) where {T}
    return gaussian_state_sampler!(points, StateMatrix(state), bias_phase)
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

function reject!(new_points, gen_point, p, g, c)
    Threads.@threads for i in 1:size(new_points, 2)
        new_points[:, i] .= gen_point()
        while p(Threads.threadid(), new_points[:, i]...) / g(new_points[:, i]...) < c
            new_points[:, i] .= gen_point()
        end
    end

    return new_points
end

"""
    state_sampler(
        state::AbstractState, n::Integer;
        warm_up_n=128, batch_size=64, c=0.9, θ_range=(0., 2π), x_range=(-10., 10.),
        show_log=false
    )

Random points sampled from quadrature probability density function of Gaussian `state`.

* `state`: Quantum state.
* `n`: N points.
* `warm_up_n`: N points sampled from uniform random and accepted by rejection method.
* `batch_size`: Adapt `g` for every `batch_size` points.
* `θ_range`: Sampling range of θ.
* `x_range`: Sampling range of x.
"""
function state_sampler(
    state::AbstractState, n::Integer;
    warm_up_n=128, batch_size=64, c=0.9, θ_range=(0., 2π), x_range=(-10., 10.),
    show_log=false
)
    sampled_points = Matrix{Float64}(undef, 2, n)
    𝛑̂_res_vec = [Matrix{complex(Float64)}(undef, state.dim, state.dim) for _ in 1:Threads.nthreads()]

    return state_sampler!(
        sampled_points, 𝛑̂_res_vec,
        state,
        warm_up_n, batch_size, c, θ_range, x_range,
        show_log
    )
end

function state_sampler!(
    sampled_points::Matrix{T}, 𝛑̂_res_vec::Vector{Matrix{Complex{T}}},
    state::StateMatrix,
    warm_up_n::Integer, batch_size::Integer, c::Real, θ_range, x_range,
    show_log::Bool
) where {T}
    n = size(sampled_points, 2)
    p = (thread_id, θ, x) -> q_pdf!(𝛑̂_res_vec[thread_id], state, θ, x)

    show_log && @info "Warm up"
    warm_up_n = n < warm_up_n ? n : warm_up_n
    warm_up_points = view(sampled_points, :, 1:warm_up_n)
    gen_rand_point = () -> [ranged_rand(θ_range), ranged_rand(x_range)]
    kde_result = kde((ranged_rand(n, θ_range), ranged_rand(n, x_range)))
    g = (θ, x) -> pdf(kde_result, θ, x)
    reject!(warm_up_points, gen_rand_point, p, g, c)

    show_log && @info "Start to generate data"
    batch = div(n-warm_up_n, batch_size)
    for b in 1:batch
        ref_range = 1:(warm_up_n+(b-1)*batch_size)
        ref_points = view(sampled_points, :, ref_range)
        new_range = (warm_up_n+(b-1)*batch_size+1):(warm_up_n+b*batch_size)
        new_points = view(sampled_points, :, new_range)

        h = KernelDensity.default_bandwidth((ref_points[1, :], ref_points[2, :]))
        gen_point_from_g = () -> ref_points[:, rand(ref_range)] + randn(2)./h
        kde_result = kde((ref_points[1, :], ref_points[2, :]), bandwidth=h)
        g = (θ, x) -> pdf(kde_result, θ, x)
        reject!(new_points, gen_point_from_g, p, g, c)

        show_log && @info "progress: $b/$(batch+1)"
    end
    rem_n = rem(n-warm_up_n, batch_size)
    if  rem_n > 0
        ref_range = 1:(n-rem_n)
        ref_points = view(sampled_points, :, ref_range)
        new_range = (n-rem_n+1):n
        new_points = view(sampled_points, :, new_range)

        h = KernelDensity.default_bandwidth((ref_points[1, :], ref_points[2, :]))
        gen_point_from_g = () -> ref_points[:, rand(ref_range)] + randn(2)./h
        kde_result = kde((ref_points[1, :], ref_points[2, :]), bandwidth=h)
        g = (θ, x) -> pdf(kde_result, θ, x)
        reject!(new_points, gen_point_from_g, p, g, c)
    end
    show_log && @info "progress: $(batch+1)/$(batch+1)"

    sampled_points .= sampled_points[:, sortperm(sampled_points[1, :])]

    return sampled_points
end

function state_sampler!(
    sampled_points::Matrix{T}, 𝛑̂_res_vec::Vector{Matrix{Complex{T}}},
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
