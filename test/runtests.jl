using Test
using AirfoilMeshers
using LinearAlgebra
using Statistics
using StaticArrays

mkpath(joinpath(@__DIR__, "output"))

# NACA0012 airfoil
yt(x, t) = 5t * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x^2 + 0.2843 * x^3 - 0.1036x^4)
t = 0.12
xx = map(a -> (1 - cos(a)) / 2, LinRange(0, pi, 65))
yy = map(x -> yt(x, t), xx)
xaf = vcat(reverse(xx), xx[2:end])
yaf = vcat(-reverse(yy), yy[2:end])

include("test_tridiag.jl")
include("test_surface.jl")
include("test_hyperbolic.jl")