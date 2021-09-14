function save_mesh(::Val{:FLUENT}, io::IO, m::AirfoilMesh)
    NI, NJ = size(m.grid)
    NV = NI * NJ - m.ncon
    NF = NI * (NJ - 1) + NJ * (NI - 1) - max(m.ncon - 1, 0)
    NC = (NI - 1) * (NJ - 1)

    NLI = Matrix(LinearIndices((NI, NJ)))
    NCI = Dict{Int,CartesianIndex{2}}()
    # try do remove duplicated
    if m.ncon > 0
        if m.topo == :O
            @assert m.ncon == NJ
            NLI[:, NJ] .= NLI[:, 1]
        elseif m.topo == :C
            NLI[NI:-1:(NI - m.ncon + 1), 1] .= NLI[1:(m.ncon), 1]
        end
        old_indicies = sort(collect(Set(NLI)))
        @assert length(old_indicies) == NV
        o2n = Dict{Int,Int}()
        for i in eachindex(old_indicies)
            o2n[old_indicies[i]] = i
        end
        for i in eachindex(NLI)
            NLI[i] = o2n[NLI[i]]
        end
    end
    for j in 1:NJ
        for i in 1:NI
            li = NLI[i, j]
            if !haskey(NCI, li)
                NCI[li] = CartesianIndex(i, j)
            end
        end
    end

    CLI = LinearIndices((NI - 1, NJ - 1))
    f_int = Vector{NTuple{4,Int}}()
    f_wall = Vector{NTuple{4,Int}}()
    f_far = Vector{NTuple{4,Int}}()
    f_min = Vector{NTuple{4,Int}}()
    f_max = Vector{NTuple{4,Int}}()
    # INTERIOR
    for j in 1:(NJ - 1)
        for i in 2:(NI - 1)
            p0 = NLI[i, j]
            p1 = NLI[i, j + 1]
            c0 = CLI[i - 1, j]
            c1 = CLI[i, j]
            push!(f_int, (p0, p1, c0, c1))
        end
    end
    for j in 2:(NJ - 1)
        for i in 1:(NI - 1)
            p0 = NLI[i, j]
            p1 = NLI[i + 1, j]
            c0 = CLI[i, j]
            c1 = CLI[i, j - 1]
            push!(f_int, (p0, p1, c0, c1))
        end
    end
    # J = 1
    if m.topo == :C
        for i in 1:(m.ncon - 1)
            p0 = NLI[i, 1]
            p1 = NLI[i + 1, 1]
            c0 = CLI[i, 1]
            c1 = CLI[NI - i, 1]
            push!(f_int, (p0, p1, c0, c1))
        end
        for i in (m.ncon):(NI - m.ncon)
            p0 = NLI[i, 1]
            p1 = NLI[i + 1, 1]
            c0 = CLI[i, 1]
            push!(f_wall, (p0, p1, c0, 0))
        end
    else
        for i in 1:(NI - 1)
            p0 = NLI[i, 1]
            p1 = NLI[i + 1, 1]
            c0 = CLI[i, 1]
            push!(f_wall, (p0, p1, c0, 0))
        end
    end
    # J = NJ
    for i in 1:(NI - 1)
        p0 = NLI[i + 1, NJ] # reversed at J=NJ
        p1 = NLI[i, NJ]
        c0 = CLI[i, NJ - 1]
        c1 = 0
        push!(f_far, (p0, p1, c0, c1))
    end
    # I = 1
    for j in 1:(NJ - 1)
        p0 = NLI[1, j + 1] # reversed
        p1 = NLI[1, j]
        c0 = CLI[1, j]
        if m.topo == :O
            c1 = CLI[NI - 1, j]
            push!(f_int, (p0, p1, c0, c1))
        else
            push!(f_min, (p0, p1, c0, 0))
        end
    end
    # I = NI
    if m.topo != :O
        for j in 1:(NJ - 1)
            p0 = NLI[NI, j]
            p1 = NLI[NI, j + 1]
            c0 = CLI[NI - 1, j]
            push!(f_max, (p0, p1, c0, 0))
        end
    end
    NFINT = length(f_int)
    I_INT = (1, NFINT)
    NFWAL = length(f_wall)
    I_WAL = (I_INT[2] + 1, I_INT[2] + NFWAL)
    NFFAR = length(f_far)
    I_FAR = (I_WAL[2] + 1, I_WAL[2] + NFFAR)
    NFMIN = length(f_min)
    I_MIN = (I_FAR[2] + 1, I_FAR[2] + NFMIN)
    NFAMX = length(f_max)
    I_MAX = (I_MIN[2] + 1, I_MIN[2] + NFAMX)
    @assert I_MAX[2] == NF
    # informations
    _fluent_comment(io, "Exported from AirfoilMeshers.jl")
    _fluent_comment(io, "Dim = 2")
    println(io, "(2 2)")
    _fluent_comment(io, "# of nodes: $(NV)")
    println(io, "(10 (0 1 $(string(NV, base=16)) 0 2))")
    _fluent_comment(io, "# of cells: $(NC)")
    println(io, "(12 (0 1 $(string(NC, base=16)) 0))")
    _fluent_comment(io, "# of faces: $(NF)")
    println(io, "(13 (0 1 $(string(NF, base=16)) 0))")
    # Zone 1: Nodes
    izone = 1
    _fluent_comment(io, "Zone $izone: nodes [$NV]")
    println(io, "(10 ($izone 1 $(string(NV, base=16)) 1 2)(")
    for li in 1:NV
        p = m.grid[CartesianIndex(NCI[li])]
        @printf io "%24.16E %24.16E\n" p[1] p[2]
    end
    println(io, "))")
    # Zone 2: Cells
    izone += 1
    _fluent_comment(io, "Zone $izone: cells [$NC]")
    println(io, "(12 ($izone 1 $(string(NC, base=16)) 1 3))")
    println(io, "(45 ($izone fluid flowfield)())")
    # Zone 3: Interior
    izone += 1
    _fluent_comment(io, "Zone $izone: interior faces $(I_INT[1]):$(I_INT[2])")
    println(
        io, "(13 ($izone $(string(I_INT[1], base=16)) $(string(I_INT[2], base=16)) 2 2)("
    )
    for f in f_int
        println(io, join(string.(f, base=16), ' '))
    end
    println(io, "))")
    println(io, "(45 ($izone interior interior-flow)())")
    # Zone 4: Wall
    izone += 1
    _fluent_comment(io, "Zone $izone: wall faces $(I_WAL[1]):$(I_WAL[2])")
    println(
        io, "(13 ($izone $(string(I_WAL[1], base=16)) $(string(I_WAL[2], base=16)) 3 2)("
    )
    for f in f_wall
        println(io, join(string.(f, base=16), ' '))
    end
    println(io, "))")
    println(io, "(45 ($izone wall wall)())")
    # Zone 5: Far
    izone += 1
    _fluent_comment(io, "Zone $izone: far faces $(I_FAR[1]):$(I_FAR[2])")
    println(
        io, "(13 ($izone $(string(I_FAR[1], base=16)) $(string(I_FAR[2], base=16)) a 2)("
    )
    for f in f_far
        println(io, join(string.(f, base=16), ' '))
    end
    println(io, "))")
    println(io, "(45 ($izone velocity-inlet far)())")
    if NFMIN > 0
        # Zone 6: IMin
        izone += 1
        _fluent_comment(io, "Zone $izone: imin faces $(I_MIN[1]):$(I_MIN[2])")
        println(
            io,
            "(13 ($izone $(string(I_MIN[1], base=16)) $(string(I_MIN[2], base=16)) 5 2)(",
        )
        for f in f_min
            println(io, join(string.(f, base=16), ' '))
        end
        println(io, "))")
        println(io, "(45 ($izone pressure-outlet imin)())")
        # Zone 7: IMax
        izone += 1
        _fluent_comment(io, "Zone $izone: imax faces $(I_MAX[1]):$(I_MAX[2])")
        println(
            io,
            "(13 ($izone $(string(I_MAX[1], base=16)) $(string(I_MAX[2], base=16)) 5 2)(",
        )
        for f in f_max
            println(io, join(string.(f, base=16), ' '))
        end
        println(io, "))")
        println(io, "(45 ($izone pressure-outlet imax)())")
    end

    return nothing
end

function _fluent_comment(io::IO, msg::AbstractString, tag::Int=0)
    return println(io, "($tag \"$msg\")")
end
