@testset "Fock basis" begin
    @test FockState(Float64, 3, Vector)[1:5] ≈ [0, 0, 0, 1, 0]
    @test FockState(Float32, 3, Vector)[1:5] ≈ [0, 0, 0, 1, 0]
    @test eltype(FockState(Float64, 3, Vector)) === Float64
    @test eltype(FockState(Float32, 3, Vector)) === Float32

    @test FockState(Float64, 3, Matrix)[1:5, 1:5] ≈ [
        0 0 0 0 0;
        0 0 0 0 0;
        0 0 0 0 0;
        0 0 0 1 0;
        0 0 0 0 0;
    ]
    @test FockState(Float32, 3, Matrix)[1:5, 1:5] ≈ [
        0 0 0 0 0;
        0 0 0 0 0;
        0 0 0 0 0;
        0 0 0 1 0;
        0 0 0 0 0;
    ]
    @test eltype(FockState(Float64, 3, Matrix)) === Float64
    @test eltype(FockState(Float32, 3, Matrix)) === Float32

    @test FockState(0, Vector, DIM) ≈ FockState(Float64, 0, Vector, DIM)
    @test FockState(0, Matrix, DIM) ≈ FockState(Float64, 0, Matrix, DIM)
    @test FockState(0, DIM) ≈ FockState(Float64, 0, Vector, DIM)

    @test NumberState === FockState

    @test VacuumState(Float64, Vector, DIM) == FockState(Float64, 0, Vector, DIM)
    @test VacuumState(Float64, Matrix, DIM) == FockState(Float64, 0, Matrix, DIM)
    @test VacuumState(Vector, DIM) == FockState(Float64, 0, Vector, DIM)
    @test VacuumState(Matrix, DIM) == FockState(Float64, 0, Matrix, DIM)
    @test VacuumState(DIM) == FockState(0, Vector, DIM)
    
    @test SinglePhotonState(Float64, Vector, DIM) == FockState(Float64, 1, Vector, DIM)
    @test SinglePhotonState(Float64, Matrix, DIM) == FockState(Float64, 1, Matrix, DIM)
    @test SinglePhotonState(Vector, DIM) == FockState(Float64, 1, Vector, DIM)
    @test SinglePhotonState(Matrix, DIM) == FockState(Float64, 1, Matrix, DIM)
    @test SinglePhotonState(DIM) == FockState(1, Vector, DIM)
end

@testset "pure state" begin
end

@testset "mixed state" begin
end
