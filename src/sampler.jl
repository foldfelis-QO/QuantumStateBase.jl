using KernelDensity

export
    gaussian_state_sampler,
    gaussian_state_sampler!,
    nongaussian_state_sampler,
    nongaussian_state_sampler!

###########################
# gaussian data generator #
###########################

function gaussian_state_sampler(state::StateMatrix, n::Integer; bias_phase=0.)
    points = Vector{Float64}(undef, n)

    return gaussian_state_sampler!(points, state, bias_phase)
end

function gaussian_state_sampler!(
    points::AbstractVector{Float64},
    state::StateMatrix, bias_phase::Float64
)
    n = length(points)

    # θs
    view(points, :) .= sort!(2π*rand(n) .+ bias_phase)

    # μ and σ given θ
    μ = π̂ₓ_μ(view(points, :), state)
    σ = real(sqrt.(π̂ₓ²_μ(view(points, :), state) - μ.^2))

    # xs
    view(points, :) .= real(μ) + σ .* randn(n)

    return points
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

function nongaussian_state_sampler(
    state;
    n=4096, warm_up_n=128, batch_size=64, c=0.9, θ_range=(0., 2π), x_range=(-10., 10.),
    show_log=true
)
    sampled_points = Matrix{Float64}(undef, 2, n)
    𝛑̂_res_vec = [Matrix{complex(Float64)}(undef, state.dim, state.dim) for _ in 1:Threads.nthreads()]

    return nongaussian_state_sampler!(
        sampled_points, 𝛑̂_res_vec,
        state,
        warm_up_n, batch_size, c, θ_range, x_range,
        show_log
    )
end

function nongaussian_state_sampler!(
    sampled_points::Matrix{T}, 𝛑̂_res_vec::Vector{Matrix{Complex{T}}},
    state::StateMatrix,
    warm_up_n::Integer, batch_size::Integer, c::Real, θ_range, x_range,
    show_log::Bool
) where {T}
    n = size(sampled_points, 2)
    p = (thread_id, θ, x) -> q_pdf!(𝛑̂_res_vec[thread_id], state, θ, x)

    show_log && @info "Warm up"
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

        show_log && @info "progress: $b/$batch"
    end

    sampled_points .= sampled_points[:, sortperm(sampled_points[1, :])]

    return sampled_points
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
