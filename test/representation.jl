@testset "StateVector" begin
    dim = 70
    T = ComplexF64

    vacuum_state_v = zeros(T, dim)
    vacuum_state_v[0 + 1] = 1

    state = StateVector{T}(vacuum_state_v, dim)
    # @test repr(state) == "StateVector{ComplexF64}( " *
    #     "\e[38;2;255;102;102m\u2587" *
    #     "\e[38;2;178;178;178m\u2587"^(dim-1) *
    #     "\e[0m )"
    @test vec(state) == vacuum_state_v
    @test state.dim == dim
    @test 𝛒(state) == vacuum_state_v * vacuum_state_v'
    @test purity(state) ≈ 1

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
    # @test repr(state) == "StateMatrix{ComplexF64}(\n" *
    #     "\e[38;2;255;102;102m\u2587" * "\e[38;2;178;178;178m\u2587"^(dim-1) * "\n" *
    #     ("\e[38;2;178;178;178m\u2587"^dim * "\n")^(dim-1) *
    #     "\e[0m)"
    @test state.dim == dim
    @test 𝛒(state) == vacuum_state_𝛒
    @test purity(state) ≈ 1

    # new_state = create!(copy(state))
    # @test new_state == SinglePhotonState(rep=StateMatrix, dim=dim)
    # @test state == VacuumState(rep=StateMatrix, dim=dim)

    # constructor for state vector
    vacuum_state_v = zeros(T, dim)
    vacuum_state_v[0 + 1] = 1
    state_vector = StateVector{T}(vacuum_state_v, dim)
    @test state == StateMatrix(state_vector)
end
