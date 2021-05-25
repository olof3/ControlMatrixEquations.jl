@testset "utilities" begin

# Test _schurstructure
S = [-4.4][:,:]
@test ([1:1], 1) == ControlMatrixEquations._schurstructure(S, I(1), Val(:U))
@test ([1:1], 1) == ControlMatrixEquations._schurstructure(S, I(1), Val(:L))

@test ([1:2], 1) == ControlMatrixEquations._schurstructure(ones(2,2), I(2), Val(:U))
@test ([1:2], 1) == ControlMatrixEquations._schurstructure(ones(2,2), I(2), Val(:L))
@test ([1:2], 1) == ControlMatrixEquations._schurstructure(I(2), ones(2,2), Val(:U))
@test ([1:2], 1) == ControlMatrixEquations._schurstructure(I(2), ones(2,2), Val(:L))

S = [1 2; 0 3]
@test ([1:1, 2:2], 2) == ControlMatrixEquations._schurstructure(S, I(2), Val(:U))
@test ([1:2], 1) == ControlMatrixEquations._schurstructure(S, I(2), Val(:L))

S = diagm(0 => [1, 1, 1, 1, 1, 1, 1, 1], -1 => [-3, 0, 0, 0, -2, 0, 0])
T = diagm(0 => [1, 1, 1, 1, 1, 1, 1, 1], -1 => [0, 0, 4, 0, 0, 0, 0])

b_S = [1:2, 3:3, 4:4, 5:6, 7:7, 8:8]
b_T = [1:1, 2:2, 3:4, 5:5, 6:6, 7:7, 8:8]
b_ST = [1:2, 3:4, 5:6, 7:7, 8:8]

@test (b_S, 6) == ControlMatrixEquations._schurstructure(S, I(8), Val(:U))
@test (b_S, 6) == ControlMatrixEquations._schurstructure(S', I(8), Val(:L))
@test (b_T, 7) == ControlMatrixEquations._schurstructure(I(8), T, Val(:U))
@test (b_T, 7) == ControlMatrixEquations._schurstructure(I(8), T', Val(:L))
@test (b_ST, 5) == ControlMatrixEquations._schurstructure(S, T, Val(:U))
@test (b_ST, 5) == ControlMatrixEquations._schurstructure(S', T', Val(:L))


S = diagm(0 => ones(15), -1 => [-2, 0, -1, 0, 0, -3, 0, -4, 0, 0, 2.3, 0, 0, 2.0] )
@test ControlMatrixEquations._schurstructure(S, Val(:U)) == ControlMatrixEquations._schurstructure(S, I(15), Val(:U))
@test ControlMatrixEquations._schurstructure(S', Val(:L)) == ControlMatrixEquations._schurstructure(S', I(15), Val(:L))

end


