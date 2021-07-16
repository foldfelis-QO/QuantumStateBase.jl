@testset "StateVector" begin
    dim = 70
    T = ComplexF64

    vacuum_state_v = zeros(T, dim)
    vacuum_state_v[0 + 1] = 1

    state = StateVector{T}(vacuum_state_v, dim)
    @test vec(state) == vacuum_state_v
    @test state.dim == dim
    @test 𝛒(state) == vacuum_state_v * vacuum_state_v'
    @test purity(state) ≈ 1

    @test repr(VacuumState(dim=5)) == "StateVector{ComplexF64}(dim=5, vec=[\n" *
        " 1.0 + 0.0im\n 0.0 + 0.0im\n 0.0 + 0.0im\n 0.0 + 0.0im\n 0.0 + 0.0im\n])"

    # new_state = create!(copy(state))
    # @test new_state == SinglePhotonState(dim=dim)
    # @test state == VacuumState(dim=dim)
end

@testset "StateMatrix" begin
    dim = 70
    T = ComplexF64

    vacuum_state_𝛒 = zeros(T, dim, dim)
    vacuum_state_𝛒[0+1, 0+1] = 1

    state = StateMatrix{T}(vacuum_state_𝛒, dim)
    @test state.dim == dim
    @test 𝛒(state) == vacuum_state_𝛒
    @test purity(state) ≈ 1

    @test repr(VacuumState(dim=5, rep=StateMatrix)) == "StateMatrix{ComplexF64}(dim=5, 𝛒=[\n" *
        " 1.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im\n" *
        " 0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im\n" *
        " 0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im\n" *
        " 0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im\n" *
        " 0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im  0.0 + 0.0im\n])"

    # new_state = create!(copy(state))
    # @test new_state == SinglePhotonState(rep=StateMatrix, dim=dim)
    # @test state == VacuumState(rep=StateMatrix, dim=dim)

    # constructor for state vector
    vacuum_state_v = zeros(T, dim)
    vacuum_state_v[0 + 1] = 1
    state_vector = StateVector{T}(vacuum_state_v, dim)
    @test state == StateMatrix(state_vector)
end
