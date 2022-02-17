export
    WignerFunction,
    WignerSurface,
    wigner

#=
    Wigner function by Laguerre Polynominal
=#

function calc_wigner(
    x_range::AbstractRange,
    p_range::AbstractRange,
    dim::Integer
)
    𝐰 = Array{ComplexF64,4}(undef, dim, dim, length(x_range), length(p_range))
    xs = collect(x_range)
    ps = collect(p_range)
    @sync for m in 1:dim
        Threads.@spawn for n in 1:dim
            # 𝐰[m, n, :, :] = wigner(m, n, collect(x_range), collect(p_range))
            𝐰[m, n, :, :] .=
                gaussian_function(xs, ps) .*
                coefficient_of_wave_function(m, n) .*
                z_to_power(m, n, xs, ps) .*
                laguerre(m, n, xs, ps)
        end
    end

    return 𝐰
end

function load_wigner(
    bin_path::String,
    x_range::AbstractRange,
    p_range::AbstractRange,
    dim::Integer
)
    return load_𝐰(bin_path, x_range, p_range, dim)
end

function create_wigner(
    x_range::AbstractRange,
    p_range::AbstractRange,
    dim::Integer
)
    bin_path = gen_wigner_bin_path(x_range, p_range, dim)

    if isfile(bin_path)
        return load_wigner(bin_path, x_range, p_range, dim)
    end

    𝐰 = calc_wigner(x_range, p_range, dim)
    save_𝐰(bin_path, 𝐰)

    return 𝐰
end

mutable struct WignerFunction{T<:Integer, U<:AbstractRange}
    𝐰::Array{ComplexF64,4}
    x_range::U
    p_range::U
    dim::T

    function WignerFunction(
        x_range::U,
        p_range::U,
        dim::T
    ) where {T<:Integer, U<:AbstractRange}
        # !check_argv(m_dim, n_dim, x_range, p_range) && throw(ArgumentError)
        𝐰 = create_wigner(x_range, p_range, dim)

        return new{T, U}(𝐰, x_range, p_range, dim)
    end
end

function WignerFunction(x_range::AbstractRange, p_range::AbstractRange, dim)
    return WignerFunction(x_range, p_range, dim)
end

struct WignerSurface{T<:AbstractRange}
    x_range::T
    p_range::T
    𝐰_surface::Matrix{Float64}

    function WignerSurface(
        x_range::T,
        p_range::T,
        𝐰_surface::Matrix{Float64}
    ) where {T<:AbstractRange}
        return new{T}(x_range, p_range, 𝐰_surface)
    end
end

function (wf::WignerFunction)(ρ::AbstractMatrix{T}) where {T}
    ρ = collect(ρ) # this is due to LazyArrays
    𝐰_surface = Matrix{Float64}(undef, length(wf.x_range), length(wf.p_range))
    @sync for i in 1:length(wf.x_range)
        Threads.@spawn for j in 1:length(wf.p_range)
            𝐰_surface[i, j] = real(sum(ρ .* wf.𝐰[:, :, i, j]))
        end
    end

    return WignerSurface(wf.x_range, wf.p_range, 𝐰_surface)
end

# TODO: (wf::WignerFunction)(state::StateVector) = wf(ρ(state))

function wigner(ρ::AbstractMatrix{T}, x_range::AbstractRange, p_range::AbstractRange) where {T}
    dim = size(ρ, 1)
    wf = WignerFunction(x_range, p_range, dim)

    return wf(ρ)
end

# TODO: wigner(v::AbstractVector{T}, x_range::AbstractRange, p_range::AbstractRange) where {T}
