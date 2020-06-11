@testset "small_eqs_naive_explicit" begin

A = reshape(1:9, 3, 3) + I # Without the I there are numerical problems..
Ac = reshape((1:9) .+ im*(9:-1:1), 3, 3) + I
B = reshape(1:4, 2, 2)
X = reshape(1:6, 3, 2)
Xs = Matrix(Symmetric(reshape(1:9, 3, 3)))


@test sylvc(A, B, sylvc_rhs(A, B, X), Val(:naive)) ≈ X
@test sylvc(Ac, B, sylvc_rhs(Ac, B, X), Val(:naive)) ≈ X
@test lyapc(A, lyapc_rhs(A, Xs), Val(:naive)) ≈ Xs
@test lyapc(Ac, lyapc_rhs(Ac, Xs), Val(:naive)) ≈ Xs

@test sylvd(A, B, sylvd_rhs(A, B, X), Val(:naive)) ≈ X
@test sylvd(Ac, B, sylvd_rhs(Ac, B, X), Val(:naive)) ≈ X
@test lyapc(A, lyapc_rhs(A, Xs), Val(:naive)) ≈ Xs
@test lyapc(Ac, lyapc_rhs(Ac, Xs), Val(:naive)) ≈ Xs

end
