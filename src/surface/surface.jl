include("distributions.jl")
include("control.jl")

function surface_mesh(
    dist::SurfaceDistribution, x::AbstractVector, y::AbstractVector; mg_level::Int=0
)::Vector{MVector{2,Float64}}
    # check inputs
    @assert length(x) == length(y)
    @assert mg_level in 0:6
    if y[end - 1] < y[2]
        reverse!(x)
        reverse!(y)
    end
    x_te = (x[end] + x[1]) / 2
    y_te = (y[end] + y[1]) / 2
    t_te = sqrt((y[end] - y[1])^2 + (x[end] - x[1])^2)
    closed_te = t_te < 1e-3
    @assert abs(x_te - 1) < 0.01
    if closed_te
        y[end] = y[1] = y_te
        x[end] = x[1] = x_te
        np_te = 0
    end
    i_le = argmin(x)
    @assert issorted(x[i_le:end]) && issorted(x[1:i_le]; rev=true)

    s = Vector{Float64}(undef, length(x))
    s[1] = 0
    @inbounds for i in 2:length(s)
        s[i] = s[i - 1] + sqrt((x[i] - x[i - 1])^2 + (y[i] - y[i - 1])^2)
    end

    smax = s[end]
    s_le = s[i_le]
    sx = Spline1D(s, x)
    sy = Spline1D(s, y)
    xus = Spline1D(x[i_le:end], s[i_le:end])
    xls = Spline1D(x[i_le:-1:1], s[i_le:-1:1])

    s_nodes = map(n -> to_s(n, s_le, smax, xus, xls), dist.node)
    #
    sm = nodes(dist, s_nodes; mg_level=mg_level)
    xsurf, ysurf = sx(sm), sy(sm)
    if !closed_te
        ds1 = sqrt((xsurf[1]-xsurf[2])^2 + (ysurf[1]-ysurf[2])^2)
        ds2 = sqrt((xsurf[end]-xsurf[end-1])^2 + (ysurf[end]-ysurf[end-1])^2)
        ds_te = sqrt(ds1*ds2)
        ns_te = round(Int, t_te/ds_te)
        mg_te = 2^max(mg_level,1)
        np_te = ceil(Int, ns_te//mg_te)*mg_te + 1
        halfte = (np_te + 1) ÷ 2
        xsurf = vcat(
            LinRange(x_te, x[1], halfte),
            xsurf[2:(end - 1)],
            LinRange(x[end], x_te, halfte),
        )
        ysurf = vcat(
            LinRange(y_te, y[1], halfte),
            ysurf[2:(end - 1)],
            LinRange(y[end], y_te, halfte),
        )
    else
        xsurf[1] = xsurf[end] = x_te
        ysurf[1] = ysurf[end] = y_te
    end
    points = map(zip(xsurf, ysurf)) do (x, y)
        MVector{2,Float64}(x, y)
    end
    points[end] = points[1]
    return points
end
