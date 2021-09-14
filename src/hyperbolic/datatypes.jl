
abstract type StepSize end

struct UserStepSize{F} <: StepSize
    f::F
    function UserStepSize(f::F) where {F}
        @assert isa(f(1), Float64) "Invalid return type, must be Float64"
        return new{F}(f)
    end
end
(s::UserStepSize)(η::Int) = s.f(max(η, 1))

struct SpecifiedStepSize <: StepSize
    s::Vector{Float64}
    function SpecifiedStepSize(a::A) where {T<:Real,A<:AbstractVector{T}}
        return new(convert(Vector{Float64}, a))
    end
end
(s::SpecifiedStepSize)(η::Int) = s.s[min(max(η, 1), length(s.s))]

struct ExpStepSize <: StepSize
    ds0::Float64
    ratio::Float64
end
(s::ExpStepSize)(η::Int) = s.ds0 * s.ratio^(max(η, 1) - 1)

struct ConstStepSize <: StepSize
    ds::Float64
end
(s::ConstStepSize)(::Int) = s.ds

struct FilterStepSize <: StepSize
    η::Int
    sa::Float64
    sb::Float64
end
(s::FilterStepSize)(η::Int) = η <= s.η ? s.sa : s.sb

struct MinStepSize <: StepSize
    children::Vector{StepSize}
    MinStepSize(a...) =  new(StepSize[StepSize(f) for f in a])
end
(s::MinStepSize)(η::Int) = minimum(map(ss->ss(η), s.children))



StepSize(a::Real) = ConstStepSize(a)
StepSize(a::A) where {A<:StepSize} = a
StepSize(a::Function) = UserStepSize(a)
StepSize(a::AbstractVector) = SpecifiedStepSize(a)

abstract type BoundaryCondition end

struct FreeBoundary <: BoundaryCondition end

struct SymmetryBoundary <: BoundaryCondition end

struct LineBoundary <: BoundaryCondition
    dir::SVector{2,Float64}
    LineBoundary(x, y) = new(normalize(SVector{2,Float64}(x, y)))
end

ConstantX() = LineBoundary(0, 1)
ConstantY() = LineBoundary(1, 0)
ConstantAngle(α) = LineBoundary(cosd(α), sind(α))

one_step(::FreeBoundary, h::Float64, d0::MVector{2,Float64}) = normalize(d0) * h
function one_step(b::LineBoundary, h::Float64, d0::MVector{2,Float64})
    return dot(b.dir, d0) < 0 ? -h * b.dir : h * b.dir
end

struct WakeInfo
    len::Float64
    ratio::Float64
    dsmax::Float64
    aoa::Float64
end
