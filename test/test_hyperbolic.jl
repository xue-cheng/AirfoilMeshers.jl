@testset "hyperbolic_mesh" begin
    @testset "L" begin
        np1 = 33
        np = 2np1 - 1
        sp = AirfoilMeshers.TanhSpacing(0.01, -1, np1)
        s = AirfoilMeshers.nodes(sp)
        x = vcat(fill(0.0, np - np1), s)
        y = vcat(reverse(s), fill(0.0, np - np1))
        p0 = [MVector{2,Float64}(x[i], y[i]) for i in 1:np]
        mesher = HyperbolicMesher(
            :H,
            ExpStepSize(1e-3, 1.1),
            SymmetryBoundary(),
            SymmetryBoundary();
            f_blend=1,
            n_smooth=25
        )
        mesh = gen_mesh(mesher, p0, 0.3; mg_level=2)
        save_mesh("output/L_inside.dat", mesh, :TECPLOT)
        reverse!(p0)
        mesher = HyperbolicMesher(
            :H,
            ExpStepSize(1e-3, 1.1),
            SymmetryBoundary(),
            SymmetryBoundary();
            f_blend=1,
            n_smooth=25
        )
        mesh = gen_mesh(mesher, p0, 0.3; mg_level=2)
        save_mesh("output/L_outside.dat", mesh, :TECPLOT)
        save_mesh("output/L_out.cfl.dat", mesh, :PLOT3D; cfl3d=true)
        save_mesh("output/L_out.p3d.dat", mesh, :PLOT3D)
    end
    @testset "C0012" begin
        dist = SurfaceDistribution(5e-4, 1e-5, 0.01, 0.025)
        x0 = copy(xaf)
        y0 = copy(yaf)
        points = surface_mesh(dist, x0, y0; mg_level=2)
        mesher = HyperbolicMesher(
            :C,
            MinStepSize(ExpStepSize(5e-6, 1.14), 2),
            SymmetryBoundary(),
            SymmetryBoundary();
            f_blend=1e-5,
            n_smooth=10,
            wake_len=25,
            wake_maxl=1,
            wake_ratio=1.14,
            wake_aoa=10
        )

        mesh = gen_mesh(mesher, points, 20; mg_level=2)
        save_mesh("output/C0012.dat", mesh, :TECPLOT)

        save_mesh("output/C0012.cas", mesh, :FLUENT)
    end
    @testset "Surface Jet" begin
        x0 = copy(xaf)
        y0 = copy(yaf)
        xjet = 0.90
        wjet = 0.01
        xlb = xjet - wjet / 2
        ylb = yt(xlb, t)
        xrb = xjet + wjet / 2
        yrb = yt(xrb, t)
        npjet = 33
        dsjet = sqrt((yrb - ylb)^2 + (xrb - xlb)^2) / (npjet - 1)
        ds_le = 5e-4
        ds_te = 1e-5
        ds_up_max = 0.01
        ds_lo_max = 0.025
        rs_max = 1.14
        snode = [
            TrailingNode(:L),
            LeadingNode(),
            SurfaceNode(:U, xlb),
            SurfaceNode(:U, xrb),
            TrailingNode(:U)
        ]
        sdist = [
            TanhSpacing(ds_te, ds_le; dsmax=ds_lo_max, rsmax=rs_max),
            TanhSpacing(ds_le, dsjet; dsmax=ds_up_max, rsmax=rs_max),
            EqualSpacing(npjet),
            TanhSpacing(dsjet, ds_te; dsmax=ds_up_max, rsmax=rs_max)
        ]
        dist = SurfaceDistribution(snode, sdist)
        points = surface_mesh(dist, x0, y0; mg_level=2)
        mesher = HyperbolicMesher(
            :C,
            MinStepSize(ExpStepSize(5e-6, 1.14), 2),
            SymmetryBoundary(),
            SymmetryBoundary();
            f_blend=1e-5,
            n_smooth=10,
            wake_len=25,
            wake_maxl=1,
            wake_ratio=1.14,
            wake_aoa=5
        )

        mesh = gen_mesh(mesher, points, 20; mg_level=2)
        save_mesh("output/surfjet.dat", mesh, :TECPLOT)
    end
end
