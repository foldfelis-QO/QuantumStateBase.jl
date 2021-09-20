export
    Creation,
    create!,
    create,
    Annihilation,
    annihilate!,
    annihilate,

    ComplexVec,
    α,
    ξ,

    Displacement,
    displace!,

    Squeezing,
    squeeze!

############
# a† and a #
############

"""
    Creation(; dim=DIM)

Creation operator in matrix representation

``\\hat{a}^{\\dagger}``
"""
Creation(T::Type{<:Number}; dim=DIM) = diagm(-1 => sqrt.(T.(1:dim-1)))
Creation(; dim=DIM) = Creation(ComplexF64, dim=dim)

"""
    create!(state::AbstractState)

Apply creation operator on the quantum state.

# Examples
```jldoctest
julia> state = VacuumState();

julia> create!(state);

julia> vec(state) == vec(SinglePhotonState())
true
```
"""
function create!(state::StateVector{T}) where {T}
    dim = state.dim
    state.v = Creation(T, dim=dim) * state.v

    return state
end

function create!(state::StateMatrix{T}) where {T}
    dim = state.dim
    𝐜 = Creation(T, dim=dim)
    state.𝛒 = 𝐜 * state.𝛒 * 𝐜'

    return state
end

"""
    create(state::AbstractState)

Apply creation operator on the new instance of the quantum state.

# Examples
```jldoctest
julia> state = VacuumState();

julia> new_state = create(state);

julia> vec(state) == vec(VacuumState())
true

julia> vec(new_state) == vec(SinglePhotonState())
true
```
"""
create(state::AbstractState) = create!(copy(state))

"""
    Annihilation(; dim=DIM)

Annihilation operator in matrix representation

``\\hat{a}``
"""
Annihilation(T::Type{<:Number}; dim=DIM) = diagm(1 => sqrt.(T.(1:dim-1)))
Annihilation(; dim=DIM) = Annihilation(ComplexF64, dim=dim)

"""
    annihilate!(state::AbstractState)

Apply annihilation operator on the quantum state.

# Examples
```jldoctest
julia> state = SinglePhotonState();

julia> annihilate!(state);

julia> vec(state) == vec(VacuumState())
true
```
"""
function annihilate!(state::StateVector{T}) where {T}
    dim = state.dim
    state.v = Annihilation(T, dim=dim) * state.v

    return state
end

function annihilate!(state::StateMatrix{T}) where {T}
    dim = state.dim
    𝐚 = Annihilation(T, dim=dim)
    state.𝛒 = 𝐚 * state.𝛒 * 𝐚'

    return state
end

"""
    annihilate!(state::AbstractState)

Apply annihilation operator on the new instance of the quantum state.

# Examples
```jldoctest
julia> state = SinglePhotonState();

julia> new_state = annihilate(state);

julia> vec(state) == vec(SinglePhotonState())
true

julia> vec(new_state) == vec(VacuumState())
true
```
"""
annihilate(state::AbstractState) = annihilate!(copy(state))

###########
# α and ξ #
###########

"""
    ComplexVec{T <: Real}(r::T, θ::T)

Vector in polar coordinate for complex plane.

``v = r e^{-i\\theta}``
"""
struct ComplexVec{T<:Real}
    r::T
    θ::T
end

function Base.show(io::IO, complexvec::ComplexVec{T}) where {T}
    print(io, "ComplexVec{$T}($(complexvec.r)exp(-$(complexvec.θ)im))")
end

z(complexvec::ComplexVec) = complexvec.r * exp(-im * complexvec.θ)

"""
    α(r::Real, θ::Real)

Eigenvalue of annihilation operator.

``\\hat{a} | \\alpha \\rangle = \\alpha | \\alpha \\rangle``

# Examples
```jldoctest
julia> α(1.5, π/4)
ComplexVec{Float64}(1.5exp(-0.7853981633974483im))
```
"""
α(T::Type{<:Real}, r::Real, θ::Real) = ComplexVec{T}(T(r), T(θ))
α(r, θ) = α(Float64, r, θ)

"""
    ξ(r::Real, θ::Real)

# Examples
```jldoctest
julia> ξ(1.5, π/4)
ComplexVec{Float64}(1.5exp(-0.7853981633974483im))
```
"""
const ξ = α

################
# displacement #
################

"""
    Displacement(α::ComplexVec{<:Real}; dim=DIM)

Displacement operator in matrix representation

``\\hat{D}(\\alpha) = exp(\\alpha \\hat{a}^{\\dagger} - \\alpha^{*} \\hat{a})``
"""
function Displacement(T::Type{<:Number}, α::ComplexVec; dim=DIM)
    return exp(z(α) * Creation(T, dim=dim) - z(α)' * Annihilation(T, dim=dim))
end

Displacement(α::ComplexVec; dim=DIM) = Displacement(ComplexF64, α, dim=dim)

"""
    displace!(state::AbstractState)

Apply displacement operator on the quantum state.

# Examples
```jldoctest
julia> state = VacuumState();

julia> displace!(state, α(5., π/4));

julia> vec(state) == vec(CoherentState(α(5., π/4)))
true
```
"""
function displace!(state::StateVector{T}, α::ComplexVec) where {T}
    dim = state.dim
    state.v = Displacement(T, α, dim=dim) * state.v

    return state
end

function displace!(state::StateMatrix{T}, α::ComplexVec) where {T}
    dim = state.dim
    𝐝 = Displacement(T, α, dim=dim)
    state.𝛒 = 𝐝 * state.𝛒 * 𝐝'

    return state
end

#############
# squeezing #
#############

"""
    Squeezing(ξ::ComplexVec{<:Real}; dim=DIM)

Squeezing operator in matrix representation

``\\hat{S}(\\xi) = exp(\\frac{1}{2} (\\xi^{*} \\hat{a}^{2} - \\xi \\hat{a}^{\\dagger 2}))``
"""
function Squeezing(T::Type{<:Number}, ξ::ComplexVec; dim=DIM)
    return exp(0.5 * z(ξ)' * Annihilation(T, dim=dim)^2 - 0.5 * z(ξ) * Creation(T, dim=dim)^2)
end

Squeezing(ξ::ComplexVec; dim=DIM) = Squeezing(ComplexF64, ξ, dim=dim)

"""
    squeeze!(state::AbstractState)

Apply squeezing operator on the quantum state.

# Examples
```jldoctest
julia> state = VacuumState();

julia> squeeze!(state, α(0.5, π/4));

julia> vec(state) == vec(SqueezedState(ξ(0.5, π/4)))
true
```
"""
function squeeze!(state::StateVector{T}, ξ::ComplexVec) where {T}
    dim = state.dim
    state.v = Squeezing(T, ξ, dim=dim) * state.v

    return state
end

function squeeze!(state::StateMatrix{T}, ξ::ComplexVec) where {T}
    dim = state.dim
    𝐬 = Squeezing(T, ξ, dim=dim)
    state.𝛒 = 𝐬 * state.𝛒 * 𝐬'

    return state
end

# ###############
# # measurement #
# ###############

# # ##### for arb. statein θ-x quadrature coordinate #####

# # |θ, x⟩ = ∑ₙ |n⟩ ⟨n|θ, x⟩ = ∑ₙ ψₙ(θ, x) |n⟩
# # ⟨n|θ, x⟩ = ψₙ(θ, x) = exp(im n θ) (2/π)^(1/4) exp(-x^2) Hₙ(√2 x)/√(2^n n!)
# # coeff_ψₙ = (2/π)^(1/4)/√(2^n n!)
# # ψₙ = coeff_ψₙ(n) exp(im n θ) exp(-x^2) Hₙ(√2 x)
# calc_coeff_ψₙ(n::BigInt) = (2/π)^(1/4) / sqrt(2^n * factorial(n))
# COEFF_ψₙ = [calc_coeff_ψₙ(big(n)) for n in 0:(DIM-1)]

# function extend_coeff_ψₙ!(n::Integer)
#     while length(COEFF_ψₙ)-1 < n
#         push!(COEFF_ψₙ, calc_coeff_ψₙ(big(length(COEFF_ψₙ))))
#     end
# end

# function coeff_ψₙ(n::Integer)
#     (n < length(COEFF_ψₙ)) && (return COEFF_ψₙ[n+1])

#     return calc_coeff_ψₙ(big(n))
# end

# function ψₙ(n::Integer, θ::Real, x::Real)
#     return coeff_ψₙ(n) * exp(im * n * θ - x^2) * hermiteh(n, sqrt(2)x)
# end

# function 𝛑̂!(result::Matrix{<:Complex}, θ::Real, x::Real; dim=DIM)
#     view(result, :, 1) .= ψₙ.(0:dim-1, θ, x)
#     result .= view(result, :, 1) * view(result, :, 1)'

#     return result
# end

# function 𝛑̂(θ::Real, x::Real; dim=DIM, T=ComplexF64)
#     result = Matrix{T}(undef, dim, dim)
#     U = T.parameters[1]

#     return 𝛑̂!(result, U(θ), U(x), dim=dim)
# end

# # ##### for Gaussian state in θ-x quadrature coordinate #####

# # π̂ₓ = (â exp(-im θ) + â† exp(im θ)) / 2

# tr_mul(𝐚, 𝐛) = sum(𝐚[i, :]' * 𝐛[:, i] for i in 1:size(𝐚, 1))
# create_μ(state::StateMatrix) = tr_mul(Creation(dim=state.dim), state.𝛒)
# create²_μ(state::StateMatrix) = tr_mul(Creation(dim=state.dim)^2, state.𝛒)
# annihilate_μ(state::StateMatrix) = tr_mul(Annihilation(dim=state.dim), state.𝛒)
# annihilate²_μ(state::StateMatrix) = tr_mul(Annihilation(dim=state.dim)^2, state.𝛒)
# create_annihilate_μ(state::StateMatrix) = tr_mul(
#     Creation(dim=state.dim) * Annihilation(dim=state.dim),
#     state.𝛒
# )

# # ⟨π̂ₓ²⟩ = ⟨ââ exp(-2im θ) + â†â† exp(2im θ) + ââ† + â†â⟩ / 4
# # ⟨π̂ₓ²⟩ = (exp(-2im θ)⟨â²⟩ + exp(2im θ)⟨â†²⟩ + 1 + 2⟨ââ†⟩) / 4
# # here, ⟨ââ† + â†â⟩ = 1 + 2⟨ââ†⟩ due to the commutation relation
# function π̂ₓ²_μ(θs::AbstractVector{<:Number}, state::StateMatrix)
#     return (
#         exp.(-2im*θs) .* annihilate²_μ(state) .+
#         exp.(2im*θs) .* create²_μ(state) .+
#         1 .+ 2create_annihilate_μ(state)
#     ) ./ 4
# end

# # ⟨π̂ₓ⟩ = ⟨â exp(-im θ) + â† exp(im θ)⟩ / 2
# # ⟨π̂ₓ⟩ = (exp(-im θ)⟨â⟩ + exp(im θ)⟨â†⟩) / 2
# function π̂ₓ_μ(θs::AbstractVector{<:Number}, state::StateMatrix)
#     return (
#         exp.(-im*θs) .* annihilate_μ(state) .+
#         exp.(im*θs) .* create_μ(state)
#     ) ./ 2
# end
