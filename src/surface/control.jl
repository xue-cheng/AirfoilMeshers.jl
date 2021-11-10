
abstract type AbstractSurfaceNode end

struct LeadingNode <: AbstractSurfaceNode end

to_c(::LeadingNode)::Float64 = 0
to_s(::LeadingNode, sle::Float64, ::Float64, ::Spline1D, ::Spline1D)::Float64 = sle

function Base.isless(a::A, b::B) where {A<:AbstractSurfaceNode,B<:AbstractSurfaceNode}
    return isless(to_c(a), to_c(b))
end

struct TrailingNode{S} <: AbstractSurfaceNode
    TrailingNode(side::Symbol) = (@assert side in (:L, :U); new{side}())
end
to_c(::TrailingNode{:U})::Float64 = 1
to_c(::TrailingNode{:L})::Float64 = -1
to_s(::TrailingNode{:U}, ::Float64, smax::Float64, ::Spline1D, ::Spline1D)::Float64 = smax
to_s(::TrailingNode{:L}, ::Float64, smax::Float64, ::Spline1D, ::Spline1D)::Float64 = 0

struct SurfaceNode{S} <: AbstractSurfaceNode
    x::Float64
    SurfaceNode(side::Symbol, x::Real) = (@assert side in (:L, :U); new{side}(x))
end
to_c(n::SurfaceNode{:U})::Float64 = n.x
to_c(n::SurfaceNode{:L})::Float64 = -n.x
to_s(n::SurfaceNode{:U}, ::Float64, ::Float64, su::Spline1D, ::Spline1D)::Float64 = su(n.x)
to_s(n::SurfaceNode{:L}, ::Float64, ::Float64, ::Spline1D, sl::Spline1D)::Float64 = sl(n.x)

struct SurfaceDistribution
    node::Vector{AbstractSurfaceNode}
    dist::Vector{AbstractDist}
    function SurfaceDistribution(
        nodes::N, dist::D
    )where{NN<:AbstractSurfaceNode, DD<:AbstractDist, N<:AbstractVector{NN}, D<:AbstractVector{DD}}
        @assert (length(nodes)-1) == length(dist)
        @assert length(dist) > 0 
        return new(nodes, dist)
    end
    function SurfaceDistribution(ds_le, ds_te, ds_up_max, ds_lo_max, rs_max=1.2)
        return new(
            [TrailingNode(:L), LeadingNode(), TrailingNode(:U)],
            [
                TanhSpacing(ds_te, ds_le; dsmax=ds_lo_max, rsmax=rs_max),
                TanhSpacing(ds_le, ds_te; dsmax=ds_up_max, rsmax=rs_max),
            ],
        )
    end
end

function nodes(
    d::SurfaceDistribution, s::S; mg_level::Int=0
) where {T<:Real,S<:AbstractVector{T}}
    nseg = length(d.dist)
    @assert nseg > 1
    @assert length(s) == length(d.node) == nseg + 1
    @assert issorted(d.node)
    @assert issorted(s) && s[1] == 0
    ssall = []
    for i in 1:nseg
        s0 = s[i]
        s1 = s[i + 1]
        len = s1 - s0
        ss = nodes(d.dist[i]; len=len) .+ s0
        @assert length(ss) > 1
        push!(ssall, ss)
    end
    if mg_level > 0
        all(np_fixed.(d.dist)) && error(
            "Cannot ensure multi-grid level = $mg_level, all segments have fixed `np`"
        )
        mg_step = 2^mg_level
        np_tot = sum(length, ssall) - nseg + 1 # nseg has (nseg-1) common nodes
        n, m = divrem(np_tot - 1, mg_step)
        if m != 0
            for p in 1:(mg_step - m)
                ratios = map(1:nseg) do i
                    if np_fixed(d.dist[i])
                        return 0
                    end
                    ss = ssall[i]
                    ds = ss[2:end] - ss[1:(end - 1)]
                    mr = 1.0
                    for j in 2:length(ds)
                        r = ds[j] / ds[j - 1]
                        if r < 1
                            r = 1 / r
                        end
                        if r > mr
                            mr = r
                        end
                    end
                    return mr
                end
                i = argmax(ratios)
                np0 = length(ssall[i])
                s0 = s[i]
                slen = s[i + 1] - s0
                ssall[i] = nodes(d.dist[i], np0 + 1; len=slen) .+ s0
            end
        end
    end
    ret = ssall[1]
    for i in 2:nseg
        append!(ret, ssall[i][2:end])
    end
    return ret
end
