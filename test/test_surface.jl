@testset "surface_mesh" begin
    dist = SurfaceDistribution(1e-3, 1e-4, .025, .05)
    af = NACA"0012"s
    x0, y0 = gen_airfoil(af)
    points = surface_mesh(dist, x0, y0)
    xs = getindex.(points, 1)
    ys = getindex.(points, 2)
    im = argmin(xs)
    x1 = xs[2:im-1]
    y1 = ys[2:im-1]
    dy = abs.(y_lower(af, x1)-y1)
    @test all(dy .< 1e-4)
    @test ys[2] < ys[end-1]
    @test points[1] === points[end]
    af = NACA"4412"
    x0, y0 = gen_airfoil(af)
    points = surface_mesh(dist, x0, y0)
    @test points[2][2] < points[end-1][2]
    @test points[1] === points[end]
end