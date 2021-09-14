
function save_mesh(file::AbstractString, m::AirfoilMesh, format::Symbol)
    open(file, "w") do io
        save_mesh(io, m, format)
    end
end

const _valid_formats = (:TECPLOT, :FLUENT)

function save_mesh(io::IO, m::AirfoilMesh, format::Symbol) 
    @assert format in _valid_formats "invalid format :$format, must be one of $(_valid_formats)"
    save_mesh(Val(format), io, m)
end

include("tecplot.jl")
include("fluent.jl")