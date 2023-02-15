@testset "surface_mesh" begin
    dist = SurfaceDistribution(1e-3, 1e-4, 0.025, 0.05)
    x0 = copy(xaf)
    y0 = copy(yaf)
    points = surface_mesh(dist, x0, y0)
    xs = getindex.(points, 1)
    ys = getindex.(points, 2)
    im = argmin(xs)
    xgrd = xs[im+1:end-1]
    ygrd = ys[im+1:end-1]
    yext = map(x -> yt(x, t), xgrd)
    dy = abs.(ygrd .- yext)
    @test all(dy .< 1e-4)
    @test ys[2] < ys[end-1]
    @test points[1] === points[end]
end