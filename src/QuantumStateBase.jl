module QuantumStateBase
    using LinearAlgebra
    using ClassicalOrthogonalPolynomials
    using OhMyArtifacts

    const my_artifacts = Ref{String}()

    function __init__()
        my_artifacts[] = @my_artifacts_toml!()
        return
    end

    # polynomial
    include("polynomial.jl")

    # state
    include("operator.jl")
    include("state.jl")

    # wigner
    include("wigner_util.jl")
    include("wigner.jl")
end
