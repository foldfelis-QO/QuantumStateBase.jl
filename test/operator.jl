@testset "a† and a" begin
    dim = 70

    @test create!(VacuumState(dim=dim)) ≈ SinglePhotonState(dim=dim)
    @test annihilate!(SinglePhotonState(dim=dim)) ≈ VacuumState(dim=dim)
    @test create!(VacuumState(dim=dim, rep=StateMatrix)) ≈ SinglePhotonState(dim=dim, rep=StateMatrix)
    @test annihilate!(SinglePhotonState(dim=dim, rep=StateMatrix)) ≈ VacuumState(dim=dim, rep=StateMatrix)

    @test create(VacuumState(dim=dim)) ≈ SinglePhotonState(dim=dim)
    @test annihilate(SinglePhotonState(dim=dim)) ≈ VacuumState(dim=dim)
    @test create(VacuumState(dim=dim, rep=StateMatrix)) ≈ SinglePhotonState(dim=dim, rep=StateMatrix)
    @test annihilate(SinglePhotonState(dim=dim, rep=StateMatrix)) ≈ VacuumState(dim=dim, rep=StateMatrix)
end

@testset "α and ξ" begin
    @test repr(Arg(2., π/4)) == "Arg{Float64}(2.0exp($(π/4)im))"
    @test QSB.z(α(2., π/4)) ≈ 2 * exp(-im * π/4)
    @test QSB.z(ξ(2., π/4)) ≈ 2 * exp(-im * π/4)
end

@testset "Displacement" begin
    dim = 70
    r = 2.
    θ = π/4

    @test displace!(VacuumState(dim=dim), α(r, θ)).v ≈ exp(
        QSB.z(α(r, θ)) * Creation(dim=dim) -
        QSB.z(α(r, θ))' * Annihilation(dim=dim)
    ) * VacuumState().v
    @test displace!(VacuumState(dim=dim, rep=StateMatrix), α(r, θ)).𝛒 ≈ exp(
        QSB.z(α(r, θ)) * Creation(dim=dim) -
        QSB.z(α(r, θ))' * Annihilation(dim=dim)
    ) * VacuumState(rep=StateMatrix).𝛒 * exp(
        QSB.z(α(r, θ)) * Creation(dim=dim) -
        QSB.z(α(r, θ))' * Annihilation(dim=dim)
    )'
end

@testset "squeezing" begin
    dim = 70
    r = 2.
    θ = π/4

    @test squeeze!(VacuumState(dim=dim), α(r, θ)).v ≈ exp(
        0.5 * QSB.z(ξ(r, θ))' * Annihilation(dim=dim)^2 -
        0.5 * QSB.z(ξ(r, θ)) * Creation(dim=dim)^2
    ) * VacuumState().v
    @test squeeze!(VacuumState(dim=dim, rep=StateMatrix), α(r, θ)).𝛒 ≈ exp(
        0.5 * QSB.z(ξ(r, θ))' * Annihilation(dim=dim)^2 -
        0.5 * QSB.z(ξ(r, θ)) * Creation(dim=dim)^2
    ) * VacuumState(rep=StateMatrix).𝛒 * exp(
        0.5 * QSB.z(ξ(r, θ))' * Annihilation(dim=dim)^2 -
        0.5 * QSB.z(ξ(r, θ)) * Creation(dim=dim)^2
    )'
end

@testset "measurement" begin
    @test QSB.COEFF_ψₙ == [QSB.calc_coeff_ψₙ(big(n)) for n in 0:(QSB.DIM-1)]
    @test QSB.coeff_ψₙ(100) == QSB.calc_coeff_ψₙ(big(100))
    @test QSB.COEFF_ψₙ == [QSB.calc_coeff_ψₙ(big(n)) for n in 0:100]
    @test QSB.coeff_ψₙ(500) == QSB.calc_coeff_ψₙ(big(500))
    @test QSB.COEFF_ψₙ == [QSB.calc_coeff_ψₙ(big(n)) for n in 0:500]

    ψₙs = QSB.ψₙ.(0:QSB.DIM-1, 2., 3.)
    @test QSB.𝛑̂(2, 3) ≈ ψₙs * ψₙs'
end

@testset "Gaussian state" begin
    𝐚 = rand(10, 10)
    𝐛 = rand(10, 10)

    @test QSB.tr_mul(𝐚, 𝐛) ≈ tr(𝐚 * 𝐛)

    state = SqueezedThermalState(ξ(1., π/4), 0.5)

    @test QSB.create_μ(state) ≈ tr(Creation(dim=state.dim) * state.𝛒)
    @test QSB.create²_μ(state) ≈ tr(Creation(dim=state.dim)^2 * state.𝛒)
    @test QSB.annihilate_μ(state) ≈ tr(Annihilation(dim=state.dim) * state.𝛒)
    @test QSB.annihilate²_μ(state) ≈ tr(Annihilation(dim=state.dim)^2 * state.𝛒)
    @test QSB.create_annihilate_μ(state) ≈ tr(Creation(dim=state.dim) * Annihilation(dim=state.dim) * state.𝛒)
end
