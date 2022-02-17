export
    create,
    annihilate,
    displace,
    squeeze

function check_finite_size(s::AbstractArray)
    ℵ₀ in size(s) && throw(BoundsError("infinite size"))

    return true
end

############
# a† and a #
############

Creation(T::Type{<:Number}) = .√(Diagonal(T.(0:∞))[:, 2:end])
Creation(T; dim::Integer) = view(Creation(T), 1:dim, 1:dim)

function create(v::AbstractVector{T}) where {T}
    dim = length(v)

    check_finite_size(v) && (return Creation(T, dim=dim) * v)
end

function create(ρ::AbstractMatrix{T}) where {T}
    dim = size(ρ, 1)
    c = Creation(T, dim=dim)

    check_finite_size(ρ) && (return c * ρ * c')
end

Annihilation(T::Type{<:Number}) = .√(Diagonal(T.(0:∞))[2:end, :])
Annihilation(T; dim::Integer) = view(Annihilation(T), 1:dim, 1:dim)

function annihilate(v::AbstractVector{T}) where {T}
    dim = length(v)

    check_finite_size(v) && (return Annihilation(T, dim=dim) * v)
end

function annihilate(ρ::AbstractMatrix{T}) where {T}
    dim = size(ρ, 1)
    a = Annihilation(T, dim=dim)

    check_finite_size(ρ) && (return a * ρ * a')
end

###########
# α and ξ #
###########

struct ComplexVec{T<:Real}
    r::T
    θ::T
end

ComplexVec(r::T, θ::T) where {T} = ComplexVec{T}(r, θ)

function Base.show(io::IO, complexvec::ComplexVec{T}) where {T}
    print(io, "ComplexVec{$T}($(complexvec.r)exp(-$(complexvec.θ)im))")
end

z(complexvec::ComplexVec) = complexvec.r * exp(-im * complexvec.θ)

################
# displacement #
################

function Displacement(T::Type{<:Complex}, r::Real, θ::Real; dim)
    U = real(T)
    α = ComplexVec(U(r), U(θ))

    return exp(collect( # `collect` applied due to limitation of `exp(A::AbstractMatrix)`
        z(α) * Creation(U, dim=dim) - z(α)' * Annihilation(U, dim=dim)
    ))
end

function displace(v::AbstractVector{T}, r, θ) where {T}
    dim = length(v)

    check_finite_size(v) && (return Displacement(complex(T), r, θ, dim=dim) * v)
end

function displace(ρ::AbstractMatrix{T}, r, θ) where {T}
    dim = size(ρ, 1)
    d = Displacement(complex(T), r, θ, dim=dim)

    check_finite_size(ρ) && (return d * ρ * d')
end

#############
# squeezing #
#############

function Squeezing(T::Type{<:Complex}, r::Real, θ::Real; dim)
    U = real(T)
    ξ = ComplexVec(U(r), U(θ))

    return exp(collect( # `collect` applied due to limitation of `exp(A::AbstractMatrix)`
        (z(ξ)' * Annihilation(U, dim=dim)^2)/2 - (z(ξ) * Creation(U, dim=dim)^2)/2
    ))
end

function squeeze(v::AbstractVector{T}, r, θ) where {T}
    dim = length(v)

    check_finite_size(v) && (return Squeezing(complex(T), r, θ, dim=dim) * v)
end

function squeeze(ρ::AbstractMatrix{T}, r, θ) where {T}
    dim = size(ρ, 1)
    s = Squeezing(complex(T), r, θ, dim=dim)

    check_finite_size(s) && (return s * ρ * s')
end
