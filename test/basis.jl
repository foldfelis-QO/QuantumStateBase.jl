@testset "FockState in StateVector" begin
    dim = DIM
    T = ComplexF64

    vacuum_state_v = zeros(T, dim)
    vacuum_state_v[0 + 1] = 1

    vacuum_state = StateVector{T}(vacuum_state_v, dim)

    @test FockState(0, dim=dim) == vacuum_state
    @test NumberState(0, dim=dim) == vacuum_state
    @test VacuumState(dim=dim) == vacuum_state
    @test SinglePhotonState(dim=dim) == FockState(1, dim=dim)

    @test FockState(T, 0, dim=dim) == vacuum_state
    @test NumberState(T, 0, dim=dim) == vacuum_state
    @test VacuumState(T, dim=dim) == vacuum_state
    @test SinglePhotonState(T, dim=dim) == FockState(1, dim=dim)
end

@testset "FockState in StateMatrix" begin
    dim = DIM
    T = ComplexF64

    vacuum_state_𝛒 = zeros(T, dim, dim)
    vacuum_state_𝛒[0+1, 0+1] = 1

    vacuum_state = StateMatrix{T}(vacuum_state_𝛒, dim)

    @test FockState(0, dim=dim, rep=StateMatrix) == vacuum_state
    @test NumberState(0, dim=dim, rep=StateMatrix) == vacuum_state
    @test VacuumState(dim=dim, rep=StateMatrix) == vacuum_state
    @test SinglePhotonState(dim=dim, rep=StateMatrix) == FockState(1, dim=dim, rep=StateMatrix)

    @test FockState(T, 0, dim=dim, rep=StateMatrix) == vacuum_state
    @test NumberState(T, 0, dim=dim, rep=StateMatrix) == vacuum_state
    @test VacuumState(T, dim=dim, rep=StateMatrix) == vacuum_state
    @test SinglePhotonState(T, dim=dim, rep=StateMatrix) == FockState(1, dim=dim, rep=StateMatrix)
end
