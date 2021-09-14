
abstract type AbstractDist{E} end

np_float(::AbstractDist{E}) where {E} = E
np_fixed(::AbstractDist{E}) where {E} = !E

include("dist_tanh.jl")

struct EqualSpacing{E} <: AbstractDist{E}
    np::Int
    dsmax::Float64
end

function EqualSpacing(np::Int=0; dsmax::Real=-1)
    E = np < 2
    if E
        dsmax <= 0 && "`dsmax`(> 0) must be specified if `np` < 2 (auto determin)"
    else
        dsmax <= 0 || @warn "`dsmax` is ignored since `np` has been specified"
    end
    return EqualSpacing{E}(ds0, ds1, dsmax, rsmax, np)
end

isvalid(d::EqualSpacing{true}) = d.dsmax > 0
isvalid(d::EqualSpacing{false}) = d.np > 1

function nodes(d::EqualSpacing{false}; len::L=1.0) where {L<:Real}
    return collect(LinRange(0, len, d.np))
end

function nodes(::EqualSpacing, np::Int; len::L=1.0) where {L<:Real}
    return collect(LinRange(0, len, np))
end

function nodes(d::EqualSpacing{true}; len::L=1.0) where {L<:Real}
    np = ceil(Int, len / d.dsmax) + 1
    return nodes(d, np; len=len)
end
