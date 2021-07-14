@testset "pdf of quadrature" begin

    gaussian(x, μ, σ) = exp(-0.5((x-μ)/σ)^2) / (σ*sqrt(2π))

    vacuum_state = VacuumState(rep=StateMatrix)
    v_μ = QSB.π̂ₓ_μ([0], vacuum_state)[1]
    v_σ = real(sqrt.(QSB.π̂ₓ²_μ([0], vacuum_state) .- v_μ.^2))[1]

    v(x) = gaussian(x, v_μ, v_σ)

    @test all(
        q_pdf(
            VacuumState(rep=StateMatrix, dim=100),
            LinRange(0, 2π, 2),
            LinRange(-10, 10, 100),
        )[i, :] ≈ v.(LinRange(-10, 10, 100))
        for i in 1:2
    )
end
