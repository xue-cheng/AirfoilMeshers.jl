const I2 = SMatrix{2,2,Float64}(1, 0, 0, 1)

function check_points(topo::Symbol, p0::AbstractVector{MVector{2,Float64}})
    if topo == :O || topo == :C
        @assert p0[1] === p0[end] "the first and last node must be the same object for `:O` and `:C` type grids"
    end
    np = length(p0)
    for i in 2:np
        @assert !isapprox(p0[i],p0[i - 1],atol=1e-10)  "duplicated nodes $(i-1) and $(i) (distance < 1e-10)"
    end
    return nothing
end

function gen_steps(ss::StepSize, height::Float64, mg_level::Int)::Vector{Float64}
    i::Int = 0
    h::Float64 = 0.0
    s = Vector{Float64}()
    while h < height
        i += 1
        push!(s, ss(i))
        h += s[end]
    end
    if mg_level > 0
        mg_step = 2^mg_level
        ns1 = ceil(Int, length(s)//mg_step) * mg_step
        dnp = ns1 - length(s)
        append!(s, fill(s[end], dnp))
    end
    return s
end

@generated function diag_add(a::A, val::T) where {N, T, A<:MMatrix{N,N,T}}
    quote
        $(map(i->:(@inbounds a[$i,$i] += val), 1:N)...)
        :(return a)
    end
end

@generated function diag_sub(a::A, val::T) where {N, T, A<:MMatrix{N,N,T}}
    quote
        $(map(i->:(@inbounds a[$i,$i] -= val), 1:N)...)
        :(return a)
    end
end


function create_wake(
    w::WakeInfo, surf::V; mg_level::Int=0
) where {V<:AbstractVector{MVector{2,Float64}}}
    ds0 = norm(surf[1] - surf[2])
    ds1 = norm(surf[end] - surf[end - 1])
    @assert isapprox(ds0, ds1, rtol=0.5)
    ds = sqrt(ds0 * ds1)
    ss = ExpStepSize(ds, w.ratio)
    if w.dsmax > 0
        ss = MinStepSize(ss, w.dsmax)
    end
    s = gen_steps(ss, w.len, mg_level)
    dir = SVector(cosd(w.aoa), sind(w.aoa))
    ns = length(s)
    points = Vector{MVector{2,Float64}}(undef, ns)
    points[1] = surf[end] + s[1]*dir
    for i in 2:ns
        points[i] = points[i-1] + s[i]*dir
    end

    return points
end

join_wake(surf, wake) = vcat(reverse(wake), surf, wake)
