export
    create,
    annihilate,
    displace,
    squeeze

############
# a† and a #
############
"""
    QuantumStateBase.Creation(T; dim::Integer)

Construct a creation operator in matrix form.

``\\hat{a}^{\\dagger}``
"""
Creation(T; dim::Integer) = .√(diagm(-1 => T.(Base.OneTo(dim-1))))

"""
    create(s::AbstractArray{T})

Apply creation operator on a quantum state.

## Arguments

* `s`: `s` can be a quantum state vector or a density matrix of a quantum state.

## Examples

```jldoctest
julia> state = VacuumState(dim=5)
5-element Vector{Float64}:
 1.0
 0.0
 0.0
 0.0
 0.0

julia> create(state)
5-element Vector{Float64}:
 0.0
 1.0
 0.0
 0.0
 0.0
```

```jldoctest
julia> state = VacuumState(Matrix, dim=5)
5×5 Matrix{Float64}:
 1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> create(state)
5×5 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
```
"""
function create(v::AbstractVector{T}) where {T}
    dim = length(v)

    return Creation(T, dim=dim) * v
end

function create(ρ::AbstractMatrix{T}) where {T}
    dim = size(ρ, 1)
    c = Creation(T, dim=dim)

    return c * ρ * c'
end

"""
    QuantumStateBase.Annihilation(T; dim::Integer)

Construct an annihilation operator in matrix form.

``\\hat{a}``
"""
Annihilation(T; dim::Integer) = .√(diagm(1 => T.(Base.OneTo(dim-1))))

"""
    annihilate(s::AbstractArray{T})

Apply annihilation operator on a quantum state.

## Arguments

* `s`: `s` can be a quantum state vector or a density matrix of a quantum state.

## Examples

```jldoctest
julia> state = SinglePhotonState(dim=5)
5-element Vector{Float64}:
 0.0
 1.0
 0.0
 0.0
 0.0

julia> annihilate(state)
5-element Vector{Float64}:
 1.0
 0.0
 0.0
 0.0
 0.0
```

```jldoctest
julia> state = SinglePhotonState(Matrix, dim=5)
5×5 Matrix{Float64}:
 0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> annihilate(state)
5×5 Matrix{Float64}:
 1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
```
"""
function annihilate(v::AbstractVector{T}) where {T}
    dim = length(v)

    return Annihilation(T, dim=dim) * v
end

function annihilate(ρ::AbstractMatrix{T}) where {T}
    dim = size(ρ, 1)
    a = Annihilation(T, dim=dim)

    return a * ρ * a'
end

###########
# α and ξ #
###########

struct ComplexVec{T<:Real}
    r::T
    θ::T
end

"""
    QuantumStateBase.ComplexVec(r::T, θ::T)

Construct a vector in polar coordinate for complex plane.

``v = re^{-i\\theta}``

## Arguments

* `r`: Radius of the complex number.
* `θ`: Argument of the complex number.
"""
ComplexVec(r::T, θ::T) where {T} = ComplexVec{T}(r, θ)

function Base.show(io::IO, complexvec::ComplexVec{T}) where {T}
    print(io, "ComplexVec{$T}($(complexvec.r)exp(-$(complexvec.θ)im))")
end

"""
    QuantumStateBase.z(complexvec::ComplexVec)

Evaluate the complex vector.
"""
z(complexvec::ComplexVec) = complexvec.r * exp(-im * complexvec.θ)

################
# displacement #
################
"""
    QuantumStateBase.Displacement(T::Type{<:Complex}, r::Real, θ::Real; dim::Integer)

Construct a displacement operator in matrix form.

``\\hat{D}(\\alpha) = exp(\\alpha \\hat{a}^{\\dagger} - \\alpha^{*} \\hat{a})``
"""
function Displacement(T::Type{<:Complex}, r::Real, θ::Real; dim)
    U = real(T)
    α = ComplexVec(U(r), U(θ))

    return exp(
        z(α) * Creation(U, dim=dim) - z(α)' * Annihilation(U, dim=dim)
    )
end

"""
    displace(s::AbstractArray{T}, r, θ)

Apply displacement operator on a quantum state.

## Arguments

* `s`: `s` can be a quantum state vector or a density matrix of a quantum state.
* `r`: Radius of displacement.
* `θ`: Phase of displacement.

## Examples

```jldoctest
julia> state = VacuumState(dim=5)
5-element Vector{Float64}:
 1.0
 0.0
 0.0
 0.0
 0.0

julia> displace(state, 2, π/4)
5-element Vector{ComplexF64}:
   0.14864196494651452 - 5.551115123125783e-17im
   0.15309939720606983 - 0.15309939720606955im
 8.326672684688674e-17 - 0.5201876236174409im
   -0.1359447277638306 - 0.13594472776383068im
   -0.7896009204837534 - 3.0531133177191805e-16im
```
"""
function displace(v::AbstractVector{T}, r, θ) where {T}
    dim = length(v)

    return Displacement(complex(T), r, θ, dim=dim) * v
end

function displace(ρ::AbstractMatrix{T}, r, θ) where {T}
    dim = size(ρ, 1)
    d = Displacement(complex(T), r, θ, dim=dim)

    return d * ρ * d'
end

#############
# squeezing #
#############

"""
    QuantumStateBase.Squeezing(T::Type{<:Complex}, r::Real, θ::Real; dim::Integer)

Squeezing operator in matrix form.

## Arguments

* `r`: Squeezing level.
* `θ`: Squeezing phase.

``\\hat{S}(\\xi) = exp(\\frac{1}{2}(\\xi^{*} \\hat{a}^{2} - \\xi \\hat{a}^{\\dagger 2}))``
"""
function Squeezing(T::Type{<:Complex}, r::Real, θ::Real; dim)
    U = real(T)
    ξ = ComplexVec(U(r), U(θ))

    return exp(
        (z(ξ)' * Annihilation(U, dim=dim)^2)/2 - (z(ξ) * Creation(U, dim=dim)^2)/2
    )
end

"""
    squeeze(s::AbstractArray{T})

Apply squeezing operator on a quantum state.

## Arguments

* `s`: `s` can be a quantum state vector or a density matrix of a quantum state.
* `r`: Squeezing level.
* `θ`: Squeezing phase.

## Examples

```jldoctest
julia> state = VacuumState(dim=5);

julia> new_state = squeeze(state, 0.5, π/4);

julia> new_state == SqueezedState(0.5, π/4, dim=5)
true
```

```jldoctest
julia> state = VacuumState(Matrix, dim=5)
5×5 Matrix{Float64}:
 1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> squeeze(state, 0.5, π/4)
5×5 Matrix{ComplexF64}:
    0.887225+0.0im       0.0+0.0im   -0.202612-0.202612im   0.0+0.0im  2.84731e-17+0.13399im
         0.0+0.0im       0.0+0.0im         0.0+0.0im        0.0+0.0im          0.0+0.0im
   -0.202612+0.202612im  0.0+0.0im   0.0925394+0.0im        0.0+0.0im   -0.0305987-0.0305987im
         0.0+0.0im       0.0+0.0im         0.0+0.0im        0.0+0.0im          0.0+0.0im
 2.84731e-17-0.13399im   0.0+0.0im  -0.0305987+0.0305987im  0.0+0.0im    0.0202352+0.0im
```
"""
function squeeze(v::AbstractVector{T}, r, θ) where {T}
    dim = length(v)

    return Squeezing(complex(T), r, θ, dim=dim) * v
end

function squeeze(ρ::AbstractMatrix{T}, r, θ) where {T}
    dim = size(ρ, 1)
    s = Squeezing(complex(T), r, θ, dim=dim)

    return s * ρ * s'
end
