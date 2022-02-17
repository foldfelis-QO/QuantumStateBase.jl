@testset "a† and a" begin
    v3 = zeros(Float64, 5); v3[4] = 1
    v4 = zeros(Float64, 5); v4[5] = 1
    @test create(v3) ≈ √(3+1) * v4
    @test annihilate(v4) ≈ √4 * v3
    @test eltype(create(v3)) === Float64
    @test eltype(annihilate(v4)) === Float64

    ρ3 = zeros(Float64, 5, 5); ρ3[4, 4] = 1
    ρ4 = zeros(Float64, 5, 5); ρ4[5, 5] = 1
    @test create(ρ3) ≈ (3+1) * ρ4
    @test annihilate(ρ4) ≈ 4 * ρ3
    @test eltype(create(ρ3)) === Float64
    @test eltype(annihilate(ρ4)) === Float64

    v3 = zeros(Float32, 5); v3[4] = 1
    v4 = zeros(Float32, 5); v4[5] = 1
    @test create(v3) ≈ √(3+1) * v4
    @test annihilate(v4) ≈ √4 * v3
    @test eltype(create(v3)) === Float32
    @test eltype(annihilate(v4)) === Float32

    ρ3 = zeros(Float32, 5, 5); ρ3[4, 4] = 1
    ρ4 = zeros(Float32, 5, 5); ρ4[5, 5] = 1
    @test create(ρ3) ≈ (3+1) * ρ4
    @test annihilate(ρ4) ≈ 4 * ρ3
    @test eltype(create(ρ3)) === Float32
    @test eltype(annihilate(ρ4)) === Float32
end

@testset "α and ξ" begin
    @test repr(QSB.ComplexVec(2., π/4)) == "ComplexVec{Float64}(2.0exp(-$(π/4)im))"
    @test repr(QSB.ComplexVec(2, 3)) == "ComplexVec{Int64}(2exp(-3im))"
    @test QSB.z(QSB.ComplexVec(2., π/4)) ≈ 2 * exp(-im * π/4)
end

@testset "Displacement" begin
    T = Float64
    v0 = zeros(T, DIM); v0[1] = 1
    α = QSB.ComplexVec(2., π/4)
    @test displace(v0, 2, π/4) ≈
        exp(collect(QSB.z(α) * QSB.Creation(T, dim=DIM) - QSB.z(α)' * QSB.Annihilation(T, dim=DIM))) * v0
    @test eltype(displace(v0, 2, π/4)) == complex(T)

    T = Float32
    v0 = zeros(T, DIM); v0[1] = 1
    α = QSB.ComplexVec(2., π/4)
    @test displace(v0, 2, π/4) ≈
        exp(collect(QSB.z(α) * QSB.Creation(T, dim=DIM) - QSB.z(α)' * QSB.Annihilation(T, dim=DIM))) * v0
    @test eltype(displace(v0, 2, π/4)) == complex(T)
end

@testset "Squeezing" begin
    T = Float64
    v0 = zeros(T, DIM); v0[1] = 1
    ξ = QSB.ComplexVec(1., π/4)
    @test squeeze(v0, 1, π/4) ≈
        exp(collect(0.5QSB.z(ξ)' * QSB.Annihilation(T, dim=DIM)^2 - 0.5QSB.z(ξ) * QSB.Creation(T, dim=DIM)^2)) * v0
    @test eltype(squeeze(v0, 1, π/4)) == complex(T)

    T = Float32
    v0 = zeros(T, DIM); v0[1] = 1
    ξ = QSB.ComplexVec(1., π/4)
    @test squeeze(v0, 1, π/4) ≈
        exp(collect(0.5QSB.z(ξ)' * QSB.Annihilation(T, dim=DIM)^2 - 0.5QSB.z(ξ) * QSB.Creation(T, dim=DIM)^2)) * v0
    @test eltype(squeeze(v0, 1, π/4)) == complex(T)
end
