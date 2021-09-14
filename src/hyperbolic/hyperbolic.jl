include("datatypes.jl")
include("utils.jl")

struct HyperbolicMesher
    topo::Symbol
    step::StepSize
    bc0::BoundaryCondition
    bc1::BoundaryCondition
    f_dissi::Float64
    f_blend::Float64
    n_smooth::Float64
    f_smooth::Float64
    wake::WakeInfo
    function HyperbolicMesher(
        topo::Symbol,
        step,
        bc0::BoundaryCondition=FreeBoundary(),
        bc1::BoundaryCondition=FreeBoundary();
        f_dissi::Real=1,
        f_blend::Real=0,
        n_smooth::Int=0,
        f_smooth::Real=1,
        wake_len::Real=0,
        wake_aoa::Real=0,
        wake_ratio::Real=1.1,
        wake_maxl::Real=Inf,
    )
        @assert topo in (:C, :O, :H)
        @assert 0 <= f_dissi <= 1
        @assert 0 <= f_blend <= 1
        @assert 0 <= f_smooth <= 1
        if topo == :C && wake_len > 0
            @assert wake_ratio >= 1.01
            @assert wake_maxl > 0
        else
            @assert wake_len <= 0 "'wake_*' is only supported in :C type grids"
        end

        return new(
            topo,
            StepSize(step),
            bc0,
            bc1,
            f_dissi,
            f_blend,
            n_smooth,
            f_smooth,
            WakeInfo(wake_len, wake_ratio, wake_maxl, wake_aoa),
        )
    end
end
include("mesh.jl")
