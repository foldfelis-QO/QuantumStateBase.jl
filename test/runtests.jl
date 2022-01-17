using QuantumStateBase
using Test
using LinearAlgebra
using ClassicalOrthogonalPolynomials

const QSB = QuantumStateBase
# const DIM = QSB.DIM

Base.:(==)(s1::StateVector, s2::StateVector) = (s1.v == s2.v) && (s1.dim == s2.dim)
Base.:(==)(s1::StateMatrix, s2::StateMatrix) = (s1.𝛒 == s2.𝛒) && (s1.dim == s2.dim)
Base.:(≈)(s1::StateVector, s2::StateVector) = (s1.v ≈ s2.v) && (s1.dim == s2.dim)
Base.:(≈)(s1::StateMatrix, s2::StateMatrix) = (s1.𝛒 ≈ s2.𝛒) && (s1.dim == s2.dim)

@testset "QuantumStateBase.jl" begin
    # polynomial
    include("polynomial.jl")

    # state
    # include("representation.jl")
    # include("basis.jl")
    # include("operator.jl")
    # include("state.jl")

    # wigner
    # include("wigner_util.jl")
    # include("wigner.jl")

    # pdf of quadrature
    # include("quadrature_pdf.jl")
    # include("sampler.jl")
end
