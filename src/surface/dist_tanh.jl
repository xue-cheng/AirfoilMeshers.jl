const _tanh_iter = 100
const _tanh_btol = 1e-8
const _tanh_dmin = 1e-4

struct TanhSpacing{E} <: AbstractDist{E}
    ds0::Float64
    ds1::Float64
    dsmax::Float64
    rsmax::Float64
    np::Int
end

function TanhSpacing(ds0::Real, ds1::Real, np::Int=0; dsmax::Real=-1, rsmax::Real=-1.0)
    E = np < 2
    if E
        if dsmax <= 0 && rsmax < 1.01
            error("`dsmax`(> 0) or `rsmax`(>= 1.01) must be specified if `np` < 2")
        end
        if dsmax > 0 && (dsmax < ds1 || dsmax < ds0)
            error("`dsmax` must be greater than `ds0` and `ds1`")
        end
    else
        dsmax <= 0 || @warn "`dsmax` is ignored since `np` has been specified"
        rsmax < 1.01 || @warn "`rsmax` is ignored since `np` has been specified"
    end
    return TanhSpacing{E}(ds0, ds1, dsmax, rsmax, np)
end

isvalid(d::TanhSpacing{true}) = d.dsmax > max(d.ds0, d.ds1, 0)
isvalid(d::TanhSpacing{false}) = d.np > 1

function nodes(d::TanhSpacing{false}; len::L=1.0) where {L<:Real}
    return tanh_dist(d.np, d.ds0 / len, d.ds1 / len) * len
end
function nodes(d::TanhSpacing, np::Int; len::L=1.0) where {L<:Real}
    return tanh_dist(np, d.ds0 / len, d.ds1 / len) * len
end

ratio_max(ds) =
    maximum(2:length(ds)) do i
        r = ds[i] / ds[i - 1]
        r < 1 ? 1 / r : r
    end

function nodes(d::TanhSpacing{true}; len::L=1.0) where {L<:Real}
    np = max(round(Int, len / d.dsmax) + 1, 3)
    s = nodes(d, np; len=len)
    ds = s[2:end] - s[1:(end - 1)]
    if (d.dsmax > 0)
        dsm = maximum(ds)
        while dsm > d.dsmax
            np = max(round(Int, np * sqrt(dsm / d.dsmax)), np + 1)
            s = nodes(d, np; len=len)
            ds = s[2:end] - s[1:(end - 1)]
            dsm = maximum(ds)
        end
    end
    if d.rsmax > 1.01
        rsm = ratio_max(ds)
        if (rsm > d.rsmax)
            np = max(round(Int, np * (rsm / d.rsmax)^0.45), np + 1)
            s = nodes(d, np; len=len)
            ds = s[2:end] - s[1:(end - 1)]
            rsm = ratio_max(ds)
            while rsm > d.rsmax
                np = max(round(Int, np * (rsm / d.rsmax)^0.45), np + 1)
                s = nodes(d, np; len=len)
                ds = s[2:end] - s[1:(end - 1)]
                rsm1 = ratio_max(ds)
                rsm1 > rsm && error("cannot determin `np`, check inputs.")
                rsm = rsm1
            end
        end
    end
    return s
end

function tanh_dist(np::Int, s0::Real, s1::Real)
    @assert np >= 2
    if np == 2
        return Float64[0, 1]
    end
    fix0 = s0 > 0
    fix1 = s1 > 0
    if fix0 && fix1
        s = tanh_2(np, convert(Float64, s0), convert(Float64, s1))
    elseif fix0
        s = tanh_1(np, convert(Float64, s0))
    elseif fix1
        s = tanh_1(np, convert(Float64, s1))
        reverse!(s)
        @. s = 1 - s
    else
        s = collect(LinRange(0, 1, np))
    end
    return s
end

function tanh_1(np::Int, ds0::Float64)::Vector{Float64}
    np1 = np - 1
    ds = 1 / np1
    b = ds / ds0
    if b < 1
        ?? = tanh_??2(b)
        ds1 = ds / (sinh(??) / ??)
        return tanh_2(np, ds0, ds1)
    else
        s = Vector{Float64}(undef, np)
        ?? = tanh_??(b)
        ??_2 = ?? / 2
        tnh2 = tanh(??_2)
        s[1] = 0
        s[end] = 1
        for i in 2:np1
            s[i] = (1+tanh(??_2*((i - 1)*ds - 1))/tnh2)
        end
        return s
    end
end

function tanh_2(np::Int, ds0::Float64, ds1::Float64)::Vector{Float64}
    np1 = np - 1
    ds = 1 / np1
    a = sqrt(ds1 / ds0)
    b = ds / sqrt(ds0 * ds1)
    s = Vector{Float64}(undef, np)
    s[1] = 0
    s[end] = 1
    if b < 1
        ds00 = ds0
        b0 = ds / ds00
        if b0 < 1
            ??0 = tanh_??2(b0)
            ds10 = ds * ??0 / sinh(??0)
        else
            ??0 = tanh_??(b0)
            ??0_2 = ??0 / 2
            ds10 = ds * ??0_2 / tanh(??0_2)
        end
        a0 = sqrt(ds10 / ds00)
        b0 = ds / sqrt(ds00 * ds10)
        ??0 = tanh_??(b0)
        ??0_2 = ??0 / 2
        tnh20 = tanh(??0_2)
        ds11 = ds1
        b1 = ds / ds11

        if b1 < 1.0
            ??1 = tanh_??2(b1)
            ds01 = ds * ??1 / (sinh(??1))
        else
            ??1 = tanh_??(b1)
            ??1_2 = ??1 / 2
            ds01 = ds * ??1_2 / tanh(??1_2)
        end
        a1 = sqrt(ds11 / ds01)
        b1 = ds / sqrt(ds01 * ds11)
        ??1 = tanh_??(b1)
        ??1_2 = ??1 / 2
        tnh21 = tanh(??1_2)
        for i in 2:np1
            u0 = (1 + tanh(??0_2 * (2 * (i - 1) * ds - 1.0)) / tnh20) / 2
            u1 = (1 + tanh(??1_2 * (2 * (i - 1) * ds - 1.0)) / tnh21) / 2
            d0 = u0 / (a0 + (1 - a0) * u0)
            d1 = u1 / (a1 + (1 - a1) * u1)
            s[i] = ((np1 - (i - 1)) * d0 + (i - 1) * d1) * ds
        end
    else
        ?? = tanh_??(b)
        ??_2 = ?? / 2
        tnh2 = tanh(??_2)
        for i in 2:np1
            u = 0.5 * (1 + tanh(??_2 * (2 * (i - 1) * ds - 1.0)) / tnh2)
            s[i] = (u / (a + (1 - a) * u))
        end
    end
    return s
end

"""
    Compute a value for `??` given `b` (b>=-1)
    b = sinh(??)/??
"""
function tanh_??(b::Float64)::Float64
    isapprox(b, 1; atol=_tanh_btol) && return _tanh_dmin
    @assert b > 1

    ?? = max(min(sqrt(6 * (b - 1)), 30), _tanh_dmin)

    b0 = sinh(??) / ??
    resid = (b - b0) / b
    @debug "[tanh_??] initial guess: ??=$??, b=$b0, resisual=$resid"
    for i in 1:_tanh_iter
        @assert ?? > 0
        db = (cosh(??) - b0) / ??
        @assert db != 0
        ?? += (b - b0) / db
        b0 = sinh(??) / ??
        resid = (b - b0) / b
        isapprox(resid, 0; atol=_tanh_btol) && break
        @debug "[tanh_??] iter $i: ??=$??, b=$b0, resisual=$resid"
    end
    @assert isapprox(resid, 0; atol=_tanh_btol) "Newton iteration converged poorly."
    return ??
end
"""
    Compute a value for `??` given `b` at the opposite end. (b<=1)
    b = tanh(??/2)/(??/2)
"""
function tanh_??2(b::Float64)::Float64
    isapprox(b, 1; atol=_tanh_btol) && return _tanh_dmin
    @assert b < 1
    hmin = _tanh_dmin / 2

    ??_2 = max(sqrt(6 * (1 - b)), hmin)
    b0 = tanh(??_2) / ??_2
    resid = (b - b0) / b
    @debug "[tanh_??2] initial guess: ??/2=$??_2, b=$b0, resisual=$resid"
    for i in 1:_tanh_iter
        @assert ??_2 > 0
        db = (cosh(??_2)^-2 - b0) / ??_2
        @assert db != 0
        ??_2 += (b - b0) / db
        b0 = tanh(??_2) / ??_2
        resid = (b - b0) / b
        isapprox(resid, 0; atol=_tanh_btol) && break
        @debug "[tanh_??2] iter $i: ??/2=$??_2, b=$b0, resisual=$resid"
    end
    @assert isapprox(resid, 0; atol=_tanh_btol) "Newton iteration converged poorly."
    return ??_2 * 2
end
