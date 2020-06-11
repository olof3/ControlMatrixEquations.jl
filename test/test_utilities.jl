@testset "utilities" begin

# Test _schurstructure
R = diagm(0 => [1, 1, 1, 1, 1, 1, 1, 1], -1 => [1, 0, 0, 0, 1, 0, 0])
d0 = [2, 1, 1, 2, 1, 1]
b0 = [1:2, 3:3, 4:4, 5:6, 7:7, 8:8]
@test (d0, b0, 6) == SylvesterEquations._schurstructure(R, Val(:U))
@test (d0, b0, 6) == SylvesterEquations._schurstructure(R', Val(:L))


end
