
struct HyperMeshData
    α::Vector{Float64}
    β::Vector{Float64}
    γ::Vector{Float64}
    λ::Vector{Float64}
    λi::Vector{Float64}
    Area0::Vector{Float64}
    Area1::Vector{Float64}
    rξ::Vector{MVector{2,Float64}}
    rξr::Vector{MVector{2,Float64}}
    rη::Vector{MVector{2,Float64}}
    S::Vector{MVector{2,Float64}}
    A::Vector{MMatrix{2,2,Float64,4}}
    B::Vector{MMatrix{2,2,Float64,4}}
    C::Vector{MMatrix{2,2,Float64,4}}
    function HyperMeshData(NI::Int)
        rξ = [MVector{2,Float64}(undef) for _ in 1:NI]
        rξr = [MVector{2,Float64}(undef) for _ in 1:NI]
        rη = [MVector{2,Float64}(undef) for _ in 1:NI]
        α = Vector{Float64}(undef, NI)
        β = Vector{Float64}(undef, NI)
        γ = Vector{Float64}(undef, NI)
        λ = Vector{Float64}(undef, NI)
        λi = Vector{Float64}(undef, NI)
        Area0 = Vector{Float64}(undef, NI)
        Area1 = Vector{Float64}(undef, NI)
        A = [MMatrix{2,2,Float64}(undef) for _ in 1:NI]
        B = [MMatrix{2,2,Float64}(undef) for _ in 1:NI]
        C = [MMatrix{2,2,Float64}(undef) for _ in 1:NI]
        S = [MVector{2,Float64}(undef) for _ in 1:NI]
        return new(α, β, γ, λ, λi, Area0, Area1, rξ, rξr, rη, S, A, B, C)
    end
end

function gen_mesh(
    mshr::HyperbolicMesher,
    p0::AbstractVector{MVector{2,Float64}},
    height::Real;
    mg_level::Int=0,
    timer::Bool=true,
)::AirfoilMesh
    @timeit to "gen_mesh" begin
        @assert height > 0
        @assert mg_level >= 0
        check_points(mshr.topo, p0)
        np_wake = 0
        if mshr.topo == :C && mshr.wake.len > 0
            wake = create_wake(mshr.wake, p0; mg_level=mg_level)
            np_wake = length(wake) + 1
            p0 = join_wake(p0, wake)
        end
        steps = gen_steps(mshr.step, convert(Float64, height), mg_level)
        @assert length(steps) % 2^mg_level == 0
        NJ = length(steps) + 1
        NI = length(p0)
        grid = Matrix{MVector{2,Float64}}(undef, NI, NJ)
        data = HyperMeshData(NI)
        grid[:, 1] .= p0

        if mshr.topo == :O # BCs are ignored in :O mesh
            bc0 = mshr.bc0
            bc1 = mshr.bc1
        else
            bc0 = if isa(mshr.bc0, SymmetryBoundary)
                dx, dy = p0[2] - p0[1]
                LineBoundary(-dy, dx)
            else
                mshr.bc0
            end
            bc1 = if isa(mshr.bc1, SymmetryBoundary)
                dx, dy = p0[end] - p0[end - 1]
                LineBoundary(-dy, dx)
            else
                mshr.bc1
            end
        end
        for j in 2:NJ
            for i in 1:NI
                grid[i, j] = MVector{2,Float64}(undef)
            end
            h0 = mshr.step(j - 2)
            h1 = mshr.step(j - 1)
            @timeit to "extrusion" hyp_run!(
                mshr,
                data,
                mshr.f_dissi,
                (1 + (1 - 0.4mshr.f_blend)^(j - 2)) / 2,
                h0,
                h1,
                bc0,
                bc1,
                view(grid, :, j - 1),
                view(grid, :, j),
            )
            if j == 2
                for i in 1:NI
                    rη = normalize(data.rη[i])
                    dη = grid[i, j] - grid[i, j - 1]
                    hc = norm(dη)
                    dc = dot(dη, rη)
                    if dc < 0 # negative direction
                        grid[i, j] .= rη * 0.5 * h1 + grid[i, j - 1]
                    elseif hc < h1
                        grid[i, j] .= (h1 + hc) / (2hc) * dη + grid[i, j - 1]
                    end
                end
            elseif (mshr.n_smooth > 0 && mshr.f_smooth > 0)
                α0 = mshr.f_smooth / 2
                α = min(α0, α0 * (j - 2) / NJ)
                if α > 0
                    β = 1 - α
                    for _ in 1:(mshr.n_smooth)
                        @timeit to "smooth sweept" for i in 2:(NI - 1)
                            dη = grid[i, j] - grid[i, j - 1]
                            hc = norm(dη)
                            if hc < 0.5h1
                                continue
                            end
                            p₋ = grid[i - 1, j]
                            p₊ = grid[i + 1, j]
                            p = grid[i, j]
                            lm = norm(p - p₋)
                            lp = norm(p₊ - p)
                            lmp = lm + lp
                            lm /= lmp
                            lp /= lmp
                            dss = lm * p₊ + lp * p₋
                            grid[i, j] .= β * p + α * dss
                        end
                    end
                end
            end
        end
    end
    if timer
        @info "$NI×$NJ Grid Extruded"
        show(to)
        println()
    end
    reset_timer!(to)
    nconn = if mshr.topo == :C
        np_wake
    elseif mshr.topo == :O
        NJ
    else
        0
    end
    return AirfoilMesh(mshr.topo, grid, nconn)
end

function hyp_run!(
    mshr::HyperbolicMesher,
    dat::HyperMeshData,
    fdis::Float64,
    fbld::Float64,
    h0::Float64,
    h1::Float64,
    bc0::BoundaryCondition,
    bc1::BoundaryCondition,
    r0::L0,
    r1::L1,
) where {L0<:AbstractVector{MVector{2,Float64}},L1<:AbstractVector{MVector{2,Float64}}}
    nξ = length(r0)
    # rξ
    @inbounds for ξ in 2:(nξ - 1)
        dat.rξ[ξ] .= (r0[ξ + 1] - r0[ξ - 1]) / 2
    end
    if mshr.topo == :O
        dat.rξ[1] .= (r0[2] - r0[nξ - 1]) / 2
        dat.rξ[nξ] .= dat.rξ[1]
    else
        dat.rξ[1] .= r0[2] - r0[1]
        dat.rξ[nξ] .= r0[nξ] - r0[nξ - 1]
    end
    # rξr, γ, Area0, Area1
    for ξ in 1:nξ
        xξ, yξ = dat.rξ[ξ]
        dat.rξr[ξ] .= (-yξ, xξ)
        dat.γ[ξ] = γ = xξ^2 + yξ^2
        γ_sqr = sqrt(γ)
        dat.Area0[ξ] = A0 = γ_sqr * h0
        dat.Area1[ξ] = γ_sqr * h1
        dat.rη[ξ] .= dat.rξr[ξ] * (A0 / γ)
    end
    # extrude ends
    if mshr.topo == :O
        r1[end] = r1[1] = r0[1] + dat.rξr[1] / sqrt(dat.γ[1]) * h1
    else
        r1[1] = r0[1] + one_step(bc0, h1, dat.rξr[1])
        r1[nξ] = r0[nξ] + one_step(bc1, h1, dat.rξr[nξ])
    end
    # vol_blending
    if fbld != 1
        AM1 = mean(dat.Area1)
        @. dat.Area1 = fbld * dat.Area1 + (1 - fbld) * AM1
    end
    # α, β, λ
    for ξ in 1:nξ
        dxξ0, dyξ0 = dat.rξ[ξ]
        dxη0, dyη0 = dat.rη[ξ]
        α = dxξ0 .* dxη0 - dyξ0 .* dyη0
        β = dxξ0 .* dyη0 + dxη0 .* dyξ0
        γ = dat.γ[ξ]
        dat.λ[ξ] = sqrt(α^2 + β^2) / γ * fdis
        dat.C[ξ][1, 1] = α
        dat.C[ξ][2, 1] = β
        dat.C[ξ][1, 2] = β
        dat.C[ξ][2, 2] = -α
        dat.C[ξ] ./= 2γ
        dat.A[ξ] .= -dat.C[ξ]
        dat.B[ξ] .= I2
        AS = (dat.Area0[ξ] + dat.Area1[ξ]) / γ
        dat.S[ξ] .= dat.rξr[ξ] * AS + r0[ξ]
    end
    # update dissipation at interface
    if fdis > 0
        for ξ in 1:(nξ - 1)
            dat.λi[ξ] = sqrt(dat.λ[ξ] * dat.λ[ξ + 1])
        end
        # apply dissipation
        for ξ in 2:(nξ - 1)
            λ₋ = 0.5(dat.λi[ξ - 1])
            λ₊ = 0.5(dat.λi[ξ])
            diag_sub(dat.A[ξ], λ₋)
            diag_sub(dat.C[ξ], λ₊)
            diag_add(dat.B[ξ], λ₋ + λ₊)
        end
    end
    # solve
    dat.S[2] -= dat.A[2] * r1[1]
    dat.S[nξ - 1] -= dat.C[nξ - 1] * r1[nξ]
    solve_block_tridiag!(
        view(dat.A, 2:(nξ - 1)),
        view(dat.B, 2:(nξ - 1)),
        view(dat.C, 2:(nξ - 1)),
        view(dat.S, 2:(nξ - 1)),
    )
    for ξ in 2:(nξ - 1)
        r1[ξ] .= dat.S[ξ]
    end
    return nothing
end
