using KernelDensity

@testset "pdf and Gaussian state data generator" begin
    # StateMatrix
    state = SqueezedThermalState(ξ(1., π/4), 0.5)
    θs = LinRange(0, 2π, 10)
    xs = LinRange(-10, 10, 10)

    ground_truth_pdf = q_pdf(state, θs, xs)

    single_point_pdf = (θ, x) -> q_pdf(state, θ, x)
    @test single_point_pdf.(θs, xs') ≈ ground_truth_pdf

    n = 100000
    data = gaussian_state_sampler(state, n)
    sampled_pdf = pdf(kde((data[1, :], data[2, :])), θs, xs)
    @test sum(abs.(sampled_pdf .- ground_truth_pdf)) / n < 5e-5

    # StateVector
    state = VacuumState()
    θs = LinRange(0, 2π, 10)
    xs = LinRange(-10, 10, 10)

    ground_truth_pdf = q_pdf(state, θs, xs)

    single_point_pdf = (θ, x) -> q_pdf(state, θ, x)
    @test single_point_pdf.(θs, xs') ≈ ground_truth_pdf

    n = 100000
    data = gaussian_state_sampler(state, n)
    sampled_pdf = pdf(kde((data[1, :], data[2, :])), θs, xs)
    @test sum(abs.(sampled_pdf .- ground_truth_pdf)) / n < 5e-5
end
