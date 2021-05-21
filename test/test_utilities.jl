@testset "utilities" begin



# Test _schurstructure
R = [-4.4][:,:]
@test ([1], [1:1], 1) == ControlMatrixEquations._schurstructure(R, I(1), Val(:U))
@test ([1], [1:1], 1) == ControlMatrixEquations._schurstructure(R, I(1), Val(:L))

@test ([2], [1:2], 1) == ControlMatrixEquations._schurstructure(ones(2,2), I(2), Val(:U))
@test ([2], [1:2], 1) == ControlMatrixEquations._schurstructure(ones(2,2), I(2), Val(:L))
@test ([2], [1:2], 1) == ControlMatrixEquations._schurstructure(I(2), ones(2,2), Val(:U))
@test ([2], [1:2], 1) == ControlMatrixEquations._schurstructure(I(2), ones(2,2), Val(:L))

R = [1 2; 0 3]
@test ([1,1], [1:1, 2:2], 2) == ControlMatrixEquations._schurstructure(R, I(2), Val(:U))
@test ([2], [1:2], 1) == ControlMatrixEquations._schurstructure(R, I(2), Val(:L))

R = diagm(0 => [1, 1, 1, 1, 1, 1, 1, 1], -1 => [-3, 0, 0, 0, -2, 0, 0])
S = diagm(0 => [1, 1, 1, 1, 1, 1, 1, 1], -1 => [0, 0, 4, 0, 0, 0, 0])

d_R = [2, 1, 1, 2, 1, 1]
b_R = [1:2, 3:3, 4:4, 5:6, 7:7, 8:8]

d_S = [1, 1, 2, 1, 1, 1, 1]
b_S = [1:1, 2:2, 3:4, 5:5, 6:6, 7:7, 8:8]

d_RS = [2, 2, 2, 1, 1]
b_RS = [1:2, 3:4, 5:6, 7:7, 8:8]

@test (d_R, b_R, 6) == ControlMatrixEquations._schurstructure(R, I(8), Val(:U))
@test (d_R, b_R, 6) == ControlMatrixEquations._schurstructure(R', I(8), Val(:L))
@test (d_S, b_S, 7) == ControlMatrixEquations._schurstructure(I(8), S, Val(:U))
@test (d_S, b_S, 7) == ControlMatrixEquations._schurstructure(I(8), S', Val(:L))
@test (d_RS, b_RS, 5) == ControlMatrixEquations._schurstructure(R, S, Val(:U))
@test (d_RS, b_RS, 5) == ControlMatrixEquations._schurstructure(R', S', Val(:L))


Random.seed!(1)
T = diagm(0 => ones(15), -1 => [-2, 0, -1, 0, 0, -3, 0, -4, 0, 0, 2.3, 0, 0, 2.0] )
@test ControlMatrixEquations._schurstructure(T, Val(:U)) == ControlMatrixEquations._schurstructure(T, I(15), Val(:U))
@test ControlMatrixEquations._schurstructure(T', Val(:L)) == ControlMatrixEquations._schurstructure(T', I(15), Val(:L))

end


