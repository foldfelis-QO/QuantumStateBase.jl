@testset "FockState in StateVector" begin
    dim = 70
    T = ComplexF64

    vacuum_state_v = zeros(T, dim)
    vacuum_state_v[0 + 1] = 1

    vacuum_state = StateVector{T}(vacuum_state_v, dim)

    @test FockState(0, T=T, dim=dim) == vacuum_state
    @test NumberState(0, T=T, dim=dim) == vacuum_state
    @test VacuumState(T=T, dim=dim) == vacuum_state
    @test SinglePhotonState(T=T, dim=dim) == FockState(1, dim=dim)
end

@testset "FockState in StateMatrix" begin
    dim = 70
    T = ComplexF64

    vacuum_state_𝛒 = zeros(T, dim, dim)
    vacuum_state_𝛒[0+1, 0+1] = 1

    vacuum_state = StateMatrix{T}(vacuum_state_𝛒, dim)

    @test FockState(0, T=T, dim=dim, rep=StateMatrix) == vacuum_state
    @test NumberState(0, T=T, dim=dim, rep=StateMatrix) == vacuum_state
    @test VacuumState(T=T, dim=dim, rep=StateMatrix) == vacuum_state
    @test SinglePhotonState(T=T, dim=dim, rep=StateMatrix) == FockState(1, T=T, dim=dim, rep=StateMatrix)
end
