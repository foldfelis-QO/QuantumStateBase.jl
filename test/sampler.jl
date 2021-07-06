using KernelDensity

@testset "pdf and Gaussian state data generator" begin
    state = SqueezedThermalState(ξ(1., π/4), 0.5)
    θs = LinRange(0, 2π, 10)
    xs = LinRange(-10, 10, 10)

    ground_truth_pdf = q_pdf(state, θs, xs)

    single_point_pdf = (θ, x) -> q_pdf(state, θ, x)
    @test single_point_pdf.(θs, xs') ≈ ground_truth_pdf

    n = 100000
    data = gen_gaussian_training_data(state, n)
    sampled_pdf = pdf(kde((LinRange(0, 2π, n), data)), θs, xs)

    @test sum(abs.(sampled_pdf .- ground_truth_pdf)) / n < 5e-5
end
