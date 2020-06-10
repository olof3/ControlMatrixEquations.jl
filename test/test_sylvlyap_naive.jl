using SymPy


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



@vars a1 a2 real=true
@vars b1 b2 real=true
@vars x1 x2 x3 x4 real=true

# Need the definition tranpose(a) = a to make things work for non-real numbers, possibly a problem with SymPy
Base.transpose(f::Sym) = f

A = [a1 a2; 2 a2+1]
B = [b1 1; b2 2]
X = [x1 x2; x3 x4]
Xs = [x1 x2; x2 x3]


@test sylvc(A, B, sylvc_rhs(A, B, X)) == X
@test lyapc(A, lyapc_rhs(A, Xs)) == Xs
@test sylvd(A, B, sylvd_rhs(A, B, X)) == X
