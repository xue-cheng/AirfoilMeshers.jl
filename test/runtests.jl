using Test
using AirfoilMeshers
using LinearAlgebra
using Statistics
using ParametricAirfoils
using StaticArrays

mkpath(joinpath(@__DIR__, "output"))

include("test_tridiag.jl")
include("test_surface.jl")
include("test_hyperbolic.jl")