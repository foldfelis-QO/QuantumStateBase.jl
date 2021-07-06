using QuantumStateBase
using Test
using LinearAlgebra
using ClassicalOrthogonalPolynomials

Base.:(==)(s1::StateVector, s2::StateVector) = (s1.v == s2.v) && (s1.dim == s2.dim)
Base.:(==)(s1::StateMatrix, s2::StateMatrix) = (s1.𝛒 == s2.𝛒) && (s1.dim == s2.dim)
Base.:(≈)(s1::StateVector, s2::StateVector) = (s1.v ≈ s2.v) && (s1.dim == s2.dim)
Base.:(≈)(s1::StateMatrix, s2::StateMatrix) = (s1.𝛒 ≈ s2.𝛒) && (s1.dim == s2.dim)

@testset "QuantumStateBase.jl" begin
    include("polynomial.jl")
    include("representation.jl")
    include("basis.jl")
    include("operator.jl")
    include("state.jl")
end
