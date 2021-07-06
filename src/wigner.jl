using Mmap

export
    wigner,
    WignerFunction,
    WignerSurface

#=
    Wigner function by Laguerre Polynominal
=#

abstract type CreateWignerMethod end
struct Load𝐖 <: CreateWignerMethod end
struct Calc𝐖 <: CreateWignerMethod end

function wigner(m::Integer, n::Integer, x::Vector{<:Real}, p::Vector{<:Real})
    w = gaussian_function(x, p) .*
        coefficient_of_wave_function(m, n) .*
        z_to_power(m, n, x, p) .*
        laguerre(m, n, x, p)

    return w
end

function create_wigner(
    m_dim::Integer,
    n_dim::Integer,
    x_range::AbstractRange,
    p_range::AbstractRange,
    ::Type{Calc𝐖}
)
    𝐰 = Array{ComplexF64,4}(undef, m_dim, n_dim, length(x_range), length(p_range))
    @sync for m in 1:m_dim
        Threads.@spawn for n in 1:n_dim
            𝐰[m, n, :, :] = wigner(m, n, collect(x_range), collect(p_range))
        end
    end

    return 𝐰
end

function create_wigner(
    m_dim::Integer,
    n_dim::Integer,
    x_range::AbstractRange,
    p_range::AbstractRange,
    bin_path::String,
    ::Type{Load𝐖}
)
    return load_𝐰(m_dim, n_dim, x_range, p_range, bin_path)
end

function create_wigner(
    m_dim::Integer,
    n_dim::Integer,
    x_range::AbstractRange,
    p_range::AbstractRange,
)
    bin_path = gen_wigner_bin_path(m_dim, n_dim, x_range, p_range)

    if isfile(bin_path)
        return create_wigner(m_dim, n_dim, x_range, p_range, bin_path, Load𝐖)
    end

    𝐰 = create_wigner(m_dim, n_dim, x_range, p_range, Calc𝐖)
    save_𝐰(bin_path, 𝐰)

    return 𝐰
end

mutable struct WignerFunction{T<:Integer, U<:AbstractRange}
    m_dim::T
    n_dim::T
    x_range::U
    p_range::U
    𝐰::Array{ComplexF64,4}

    function WignerFunction(
        m_dim::T,
        n_dim::T,
        x_range::U,
        p_range::U
    ) where {T<:Integer, U<:AbstractRange}
        !check_argv(m_dim, n_dim, x_range, p_range) && throw(ArgumentError)

        𝐰 = create_wigner(m_dim, n_dim, x_range, p_range)
        return new{T, U}(m_dim, n_dim, x_range, p_range, 𝐰)
    end
end

function WignerFunction(x_range::AbstractRange, p_range::AbstractRange; dim=DIM)
    return WignerFunction(dim, dim, x_range, p_range)
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

function (wf::WignerFunction)(state::StateMatrix)
    𝛒 = state.𝛒

    𝐰_surface = Matrix{Float64}(undef, length(wf.x_range), length(wf.p_range))
    @sync for i in 1:length(wf.x_range)
        Threads.@spawn for j in 1:length(wf.p_range)
            𝐰_surface[i, j] = real(sum(𝛒 .* wf.𝐰[:, :, i, j]))
        end
    end

    return WignerSurface(wf.x_range, wf.p_range, 𝐰_surface)
end

(wf::WignerFunction)(state::StateVector) = wf(StateMatrix(state))
