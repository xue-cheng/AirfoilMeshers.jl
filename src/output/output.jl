
function save_mesh(file::AbstractString, m::AirfoilMesh, format::Symbol; kwargs...)
    open(file, "w") do io
        save_mesh(io, m, format; kwargs...)
    end
end

const _valid_formats = (:TECPLOT, :FLUENT, :PLOT3D)

function save_mesh(io::IO, m::AirfoilMesh, format::Symbol; kwargs...) 
    @assert format in _valid_formats "invalid format :$format, must be one of $(_valid_formats)"
    save_mesh(Val(format), io, m; kwargs...)
end

include("tecplot.jl")
include("fluent.jl")
include("plot3d.jl")