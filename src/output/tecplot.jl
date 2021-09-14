function save_mesh(::Val{:TECPLOT}, io::IO, m::AirfoilMesh)
    NI, NJ = size(m.grid)
    println(io, """T="Exported from AirfoilMeshers.jl"
    VARIABLES = "X" "Y"
    ZONE T="FLOW" I=$NI J=$NJ DATAPACKING=POINT""")
    for p in m.grid
        @printf io "%24.16E %24.16E\n" p[1] p[2]
    end
    return nothing
end