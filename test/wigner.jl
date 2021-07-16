@testset "wigner" begin
    m, n = 40, 3
    x = p = collect(-2:0.1:2)
    @test QSB.wigner(m, n, x, p) ==
        QSB.gaussian_function(x, p) .*
        QSB.coefficient_of_wave_function(m, n) .*
        QSB.z_to_power(m, n, x, p) .*
        laguerre(m, n, x, p)
end

@testset "WignerFunction" begin
    m_dim = 10
    n_dim = 10
    xs = -1:0.1:1
    ps = -1:0.1:1

    # w/o mmap
    for _ in 1:2
        wf = WignerFunction(xs, ps)
        state = VacuumState(rep=StateMatrix)

        w = wf(state)
        ans = real(sum(state.𝛒 .* wf.𝐰, dims=(1, 2)))
        is_correct = []
        for (i, e) in enumerate(w.𝐰_surface)
            push!(is_correct, e==ans[i])
        end
        @test all(is_correct)

        wf = WignerFunction(xs, ps, dim=m_dim)
        @test size(wf.𝐰) == (m_dim, n_dim, length(xs), length(ps))
    end
end

@testset "WignerSurface" begin
    _range = -10:0.1:10
    len = length(_range)
    @test WignerSurface(_range, _range, ones(Float64, len, len)) isa
        WignerSurface{typeof(_range)}
end
