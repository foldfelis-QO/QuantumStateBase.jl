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
    bin_name = gen_wigner_bin_name(x_range, p_range, dim)
    hash = my_artifact_hash(bin_name, my_artifacts[])

    if !isnothing(hash) && my_artifact_exists(hash)
        return load_wigner(my_artifact_path(hash), x_range, p_range, dim)
    end

    𝐰 = calc_wigner(x_range, p_range, dim)
    save_𝐰(bin_name, 𝐰)

    return 𝐰
end

"""
    WignerFunction

`WignerFunction` to calculate Wigner function for a quantum state in gevin `x` and `p` range.

## Arguments

* `x_range`: First physical domain range in phase space.
* `p_range`: Second physical domain range in phase space.
* `dim`: Photon number truncation of the quantum state.

## Example

```jldoctest
julia> xs = -1:0.1:1; ps = -1:0.1:1; dim=35;

julia> wf = WignerFunction(xs, ps, dim=dim);

julia> wf(VacuumState(Matrix, dim=dim))
```
"""
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
        𝐰 = create_wigner(x_range, p_range, dim)

        return new{T, U}(𝐰, x_range, p_range, dim)
    end
end

function WignerFunction(x_range::AbstractRange, p_range::AbstractRange; dim)
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
    𝐰_surface = Matrix{Float64}(undef, length(wf.x_range), length(wf.p_range))
    @sync for i in 1:length(wf.x_range)
        Threads.@spawn for j in 1:length(wf.p_range)
            𝐰_surface[i, j] = real(sum(ρ .* wf.𝐰[:, :, i, j]))
        end
    end

    return WignerSurface(wf.x_range, wf.p_range, 𝐰_surface)
end

# TODO: (wf::WignerFunction)(state::StateVector) = wf(ρ(state))

"""
    wigner(ρ::AbstractMatrix{T}, x_range::AbstractRange, p_range::AbstractRange)

Calculate the Wigner fuction of the quantum state in the gevin `x` and `p` range.

## Arguments

* `ρ`: Density matrix of a quantum state.
* `x_range`: First physical domain range in phase space.
* `p_range`: Second physical domain range in phase space.

## Example

```jldoctest
julia> xs = -1:0.1:1; ps = -1:0.1:1;

julia> wigner(VacuumState(Matrix, dim=35), xs, ps);
```
"""
function wigner(ρ::AbstractMatrix{T}, x_range::AbstractRange, p_range::AbstractRange) where {T}
    dim = size(ρ, 1)
    wf = WignerFunction(x_range, p_range, dim=dim)

    return wf(ρ)
end

# TODO: wigner(v::AbstractVector{T}, x_range::AbstractRange, p_range::AbstractRange) where {T}
