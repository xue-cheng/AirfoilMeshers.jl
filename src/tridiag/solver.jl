
function solve_block_tridiag!(
    A::AA, B::BB, C::CC, D::DD
) where {
    AA<:AbstractVector,
    BB<:AbstractVector,
    CC<:AbstractVector,
    DD<:AbstractVector,
}
    dim = length(A)
    @assert length(A) == length(B) == length(C) == length(D)
    B[1] .= inv(B[1])
    C[1] .= B[1] * C[1]
    D[1] .= B[1] * D[1]
    for i in 2:dim
        B[i] -= A[i] * C[i - 1]
        D[i] -= A[i] * D[i - 1]
        B[i] .= inv(B[i])
        C[i] .= B[i] * C[i]
        D[i] .= B[i] * D[i]
    end
    for i in (dim - 1):-1:1
        D[i] -= C[i] * D[i + 1]
    end
end
