@testset "Fock basis" begin
    @test FockState(Float64, 3, Vector, dim=5) ≈ [0, 0, 0, 1, 0]
    @test FockState(Float32, 3, Vector, dim=5) ≈ [0, 0, 0, 1, 0]
    @test eltype(FockState(Float64, 3, Vector, dim=5)) === Float64
    @test eltype(FockState(Float32, 3, Vector, dim=5)) === Float32

    @test FockState(Float64, 3, Matrix, dim=5) ≈ [
        0 0 0 0 0;
        0 0 0 0 0;
        0 0 0 0 0;
        0 0 0 1 0;
        0 0 0 0 0;
    ]
    @test FockState(Float32, 3, Matrix, dim=5) ≈ [
        0 0 0 0 0;
        0 0 0 0 0;
        0 0 0 0 0;
        0 0 0 1 0;
        0 0 0 0 0;
    ]
    @test eltype(FockState(Float64, 3, Matrix, dim=5)) === Float64
    @test eltype(FockState(Float32, 3, Matrix, dim=5)) === Float32

    @test FockState(0, Vector, dim=DIM) ≈ FockState(Float64, 0, Vector, dim=DIM)
    @test FockState(0, Matrix, dim=DIM) ≈ FockState(Float64, 0, Matrix, dim=DIM)
    @test FockState(0, dim=DIM) ≈ FockState(Float64, 0, Vector, dim=DIM)

    @test NumberState === FockState

    @test VacuumState(Float64, Vector, dim=DIM) == FockState(Float64, 0, Vector, dim=DIM)
    @test VacuumState(Float64, Matrix, dim=DIM) == FockState(Float64, 0, Matrix, dim=DIM)
    @test VacuumState(Vector, dim=DIM) == FockState(Float64, 0, Vector, dim=DIM)
    @test VacuumState(Matrix, dim=DIM) == FockState(Float64, 0, Matrix, dim=DIM)
    @test VacuumState(dim=DIM) == FockState(0, Vector, dim=DIM)

    @test SinglePhotonState(Float64, Vector, dim=DIM) == FockState(Float64, 1, Vector, dim=DIM)
    @test SinglePhotonState(Float64, Matrix, dim=DIM) == FockState(Float64, 1, Matrix, dim=DIM)
    @test SinglePhotonState(Vector, dim=DIM) == FockState(Float64, 1, Vector, dim=DIM)
    @test SinglePhotonState(Matrix, dim=DIM) == FockState(Float64, 1, Matrix, dim=DIM)
    @test SinglePhotonState(dim=DIM) == FockState(1, Vector, dim=DIM)
end

@testset "pure state" begin
    v = displace(VacuumState(Float64, Vector, dim=DIM), 2, π/4)
    @test CoherentState(ComplexF64, 2, π/4, Vector, dim=DIM) == v
    @test CoherentState(2, π/4, Vector, dim=DIM) == v
    @test CoherentState(2, π/4, dim=DIM) == v
    v = displace(VacuumState(Float32, Vector, dim=DIM), 2, π/4)
    @test CoherentState(ComplexF32, 2, π/4, Vector, dim=DIM) == v


    ρ = displace(VacuumState(Float64, Matrix, dim=DIM), 2, π/4)
    @test CoherentState(ComplexF64, 2, π/4, Matrix, dim=DIM) == ρ
    @test CoherentState(2, π/4, Matrix, dim=DIM) == ρ
    ρ = displace(VacuumState(Float32, Matrix, dim=DIM), 2, π/4)
    @test CoherentState(ComplexF32, 2, π/4, Matrix, dim=DIM) == ρ

    v = squeeze(VacuumState(Float64, Vector, dim=DIM), 2, π/4)
    @test SqueezedState(ComplexF64, 2, π/4, Vector, dim=DIM) == v
    @test SqueezedState(2, π/4, Vector, dim=DIM) == v
    @test SqueezedState(2, π/4, dim=DIM) == v
    v = squeeze(VacuumState(Float32, Vector, dim=DIM), 2, π/4)
    @test SqueezedState(ComplexF32, 2, π/4, Vector, dim=DIM) == v


    ρ = squeeze(VacuumState(Float64, Matrix, dim=DIM), 2, π/4)
    @test SqueezedState(ComplexF64, 2, π/4, Matrix, dim=DIM) == ρ
    @test SqueezedState(2, π/4, Matrix, dim=DIM) == ρ
    ρ = squeeze(VacuumState(Float32, Matrix, dim=DIM), 2, π/4)
    @test SqueezedState(ComplexF32, 2, π/4, Matrix, dim=DIM) == ρ

end

@testset "mixed state" begin
    ρ = diagm(QSB.bose_einstein(0.3).(Float64.(0:DIM)))
    @test ThermalState(Float64, 0.3, dim=DIM) == ρ
    @test ThermalState(0.3, dim=DIM) == ρ

    ρ = squeeze(ThermalState(Float64, 0.3, dim=DIM), 1, π/4)
    @test SqueezedThermalState(ComplexF64, 1, π/4, 0.3, dim=DIM) == ρ
    @test SqueezedThermalState(1, π/4, 0.3, dim=DIM) == ρ
end
