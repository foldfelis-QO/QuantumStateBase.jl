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
[...]
5-element Vector{ComplexF64}:
    0.14864 - 5.5511e-17im
     0.1531 - 0.1531im
 8.3267e-17 - 0.52019im
   -0.13594 - 0.13594im
    -0.7896 - 3.0531e-16im
```

```jldoctest
julia> state = VacuumState(Matrix, dim=5)
5×5 Matrix{Float64}:
 1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> displace(state, 2, π/4)
[...]
5×5 Matrix{ComplexF64}:
   0.0221+0.0im          0.0228+0.0228im  …    -0.117+8.92e-17im
   0.0228-0.0228im       0.0469+0.0im          -0.121+0.121im
 4.13e-17-0.0773im       0.0796-0.0796im     9.31e-17+0.411im
  -0.0202-0.0202im    -2.43e-17-0.0416im        0.107+0.107im
   -0.117-8.92e-17im     -0.121-0.121im         0.623+0.0im
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

Construct a squeezing operator in matrix form.

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
julia> state = VacuumState(dim=5)
5-element Vector{Float64}:
 1.0
 0.0
 0.0
 0.0
 0.0

julia> squeeze(state, 0.5, π/4)
[...]
5-element Vector{ComplexF64}:
    0.94193 + 0.0im
        0.0 + 0.0im
    -0.2151 + 0.2151im
        0.0 + 0.0im
 3.0229e-17 - 0.14225im
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
[...]
5×5 Matrix{ComplexF64}:
    0.887+0.0im    0.0+0.0im   -0.203-0.203im   0.0+0.0im  2.85e-17+0.134im
      0.0+0.0im    0.0+0.0im      0.0+0.0im     0.0+0.0im       0.0+0.0im
   -0.203+0.203im  0.0+0.0im   0.0925+0.0im     0.0+0.0im   -0.0306-0.0306im
      0.0+0.0im    0.0+0.0im      0.0+0.0im     0.0+0.0im       0.0+0.0im
 2.85e-17-0.134im  0.0+0.0im  -0.0306+0.0306im  0.0+0.0im    0.0202+0.0im
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
