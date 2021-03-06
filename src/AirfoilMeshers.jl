module AirfoilMeshers

using Statistics, Printf, LinearAlgebra
using StaticArrays
using Dierckx
using TimerOutputs
using FortranFiles
struct AirfoilMesh
    topo::Symbol
    grid::Matrix{MVector{2,Float64}}
    ncon::Int
end

const to = TimerOutput()

include("tridiag/solver.jl")
include("surface/surface.jl")
include("hyperbolic/hyperbolic.jl")
include("output/output.jl")

export SurfaceDistribution, TanhSpacing, EqualSpacing
export LeadingNode, TrailingNode, SurfaceNode
export surface_mesh
export WakeInfo, create_wake, join_wake
export StepSize, UserStepSize, SpecifiedStepSize, ExpStepSize, ConstStepSize, MinStepSize, FilterStepSize
export ConstantX,ConstantY,ConstantAngle,FreeBoundary,SymmetryBoundary,LineBoundary
export HyperbolicMesher
export gen_mesh, save_mesh


end # module
