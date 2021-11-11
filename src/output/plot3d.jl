function save_mesh(::Val{:PLOT3D}, io::IO, m::AirfoilMesh; cfl3d::Bool=false)
    f = FortranFile(io)
    write(f, Int32(1))
    NI, NJ = size(m.grid)
    if cfl3d
        write(f, Int32[2, NI, NJ])
        blk = Array{Float64, 4}(undef, 2, NI, NJ, 3)
        for k in 1:NJ
            for j in 1:NI
                blk[:,j,k,1] .= m.grid[j,k][1] # x
                blk[:,j,k,3] .= m.grid[j,k][2] # z
            end
        end
        blk[1,:,:,2] .= 1 # y 
        blk[2,:,:,2] .= 0
        write(f, blk)
    else
        write(f, Int32[NI, NJ])
        blk = Array{Float64, 3}(undef, NI, NJ, 2)
        for j in 1:NJ
            for i in 1:NI
                blk[i,j,1] = m.grid[i,j][1]
                blk[i,j,2] = m.grid[i,j][2]
            end
        end
        write(f, blk)
    end
    return nothing
end