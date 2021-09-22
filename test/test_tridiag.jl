@testset "Block Tridiagonal" begin
    solve = AirfoilMeshers.solve_block_tridiag!
    DIM = 5
    A = [rand(2, 2) for i in 1:DIM]
    B = [I + rand(2, 2) for i in 1:DIM]
    C = [rand(2, 2) for i in 1:DIM]
    D = [rand(2) for i in 1:DIM]
    MM = zeros(2DIM, 2DIM)
    blk(MM, i, j) = @view MM[(2i - 1):(2i), (2j - 1):(2j)]
    for i in 1:DIM
        blk(MM, i, i) .= B[i]
        if i > 1
            blk(MM, i, i - 1) .= A[i]
        end
        if i < DIM
            blk(MM, i, i + 1) .= C[i]
        end
    end
    DD = vcat(D...)
    X0 = MM \ DD
    solve(A,B,C,D)
    X1 = vcat(D...)
    @test X0 ≈ X1
end
