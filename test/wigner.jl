@testset "WignerFunction" begin
    dim1 = 10
    dim2 = 11
    xs = -1:0.1:1
    ps = -1:0.1:1

    # w/o mmap
    for _ in 1:2
        wf = WignerFunction(xs, ps, dim=dim1)
        state = VacuumState(Matrix, dim=dim1)

        w = wf(state)
        ans = real(sum(state .* wf.𝐰, dims=(1, 2)))
        @test all(e == ans[i] for (i, e) in enumerate(w.𝐰_surface))

        @test size(wigner(VacuumState(Matrix, dim=dim2), xs, ps).𝐰_surface) == (length(xs), length(ps))
    end
end

@testset "WignerSurface" begin
    _range = -10:0.1:10
    len = length(_range)
    @test WignerSurface(_range, _range, ones(Float64, len, len)) isa
        WignerSurface{typeof(_range)}
end
