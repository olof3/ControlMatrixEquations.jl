A = reshape(1:9, 3, 3) + 2I # Without the I there are numerical problems..
Ac = reshape((1:9) .+ im*(9:-1:1), 3, 3) + 2I
B = reshape(1:4, 2, 2)
Bc = reshape(1:4, 2, 2) .+ 2im
C = reshape(1:6, 3, 2)
Cc = reshape(1:6, 3, 2) .- 1im

Q = Matrix(Symmetric(reshape(1:9, 3, 3)))

E = reshape(9:-1:1, 3, 3) + 2I # Without the I there are numerical problems..
F = reshape(2:5, 2, 2)

Qr = Symmetric(reshape(1:9, 3, 3))
Qc = Hermitian(reshape(1:9, 3, 3) .+ im)



@testset "sylv*_naive" begin

X = sylvc(A, B, C, Val(:naive))
@test norm(sylvc_residual(X, A, B, C)) < 5e-15

X = sylvc(Ac, B, C, Val(:naive))
@test norm(sylvc_residual(X, Ac, B, C)) < 5e-15

X = sylvc(Ac, Bc, Cc, Val(:naive))
@test norm(sylvc_residual(X, Ac, Bc, Cc)) < 5e-15


X = sylvd(A, B, C, Val(:naive))
@test norm(sylvd_residual(X, A, B, C)) < 5e-15

X = sylvd(Ac, B, C, Val(:naive))
@test norm(sylvd_residual(X, Ac, B, C)) < 5e-15

X = sylvd(Ac, Bc, Cc, Val(:naive))
@test norm(sylvd_residual(X, Ac, Bc, Cc)) < 1e-14


X = sylvg(A, B, C, E, F, Val(:naive))
@test norm(sylvg_residual(X, A, B, C, E, F)) < 1e-14

X = sylvg(Ac, B, C, E, F, Val(:naive))
@test norm(sylvg_residual(X, Ac, B, C, E, F)) < 1e-14

X = sylvg(Ac, Bc, Cc, E, F, Val(:naive))
@test norm(sylvg_residual(X, Ac, Bc, Cc, E, F)) < 1e-14

X = sylvg(Ac, im*B, C, E, F, Val(:naive))
@test norm(sylvg_residual(X, Ac, im*B, C, E, F)) < 2e-14

end




X = lyapc(A, Qr, Val(:naive))
@test norm(lyapc_residual(X, A, Qr)) < 1e-13

X = lyapc(Ac, Qc, Val(:naive))
@test norm(lyapc_residual(X, Ac, Qc)) < 1e-13

X = lyapc(A, Qr, E, Val(:naive))
@test norm(lyapcg_residual(X, A, Qr, E)) < 1e-13

X = lyapc(Ac, Qc, E, Val(:naive))
@test norm(lyapcg_residual(X, Ac, Qc, E)) < 1e-13


X = lyapd(A, Qr, Val(:naive))
@test norm(lyapd_residual(X, A, Qr)) < 2e-13

X = lyapd(Ac, Qc, Val(:naive))
@test norm(lyapd_residual(X, Ac, Qc)) < 1e-13

# Modified some elements to get meaningful solutions
X = lyapd(2A, Qr, E, Val(:naive))
@test norm(lyapdg_residual(X, 2A, Qr, E)) < 1e-13

X = lyapd(Ac, Qc, 2E, Val(:naive))
@test norm(lyapdg_residual(X, Ac, Qc, 2E)) < 1e-13



#@edit lyapd(Ac, Qc, E, Val(:naive))
X = sylvg(Ac, Ac', -Qc, -E, E', Val(:naive))
X = sylvg(Ac, Ac', -Qc, -E, E', Val(:bartstew))

MatrixEquations.gsylv(Ac, Ac', -E, E', -Qc)
norm(sylvg_residual(X, Ac, Ac', -Qc, -E, E'))


@testset "lyap*_naive" begin


X = sylvc(Ac, B, C, Val(:naive))
@test norm(sylvc_residual(X, Ac, B, C)) < 5e-15

X = sylvd(A, B, C, Val(:naive))
@test norm(sylvd_residual(X, A, B, C)) < 5e-15

X = sylvd(Ac, B, C, Val(:naive))
@test norm(sylvd_residual(X, Ac, B, C)) < 5e-15

X = sylvg(A, B, C, E, F, Val(:naive))
@test norm(sylvg_residual(X, A, B, C, E, F)) < 1e-14

X = sylvg(Ac, B, C, E, F, Val(:naive))
@test norm(sylvg_residual(X, Ac, B, C, E, F)) < 1e-14
    
end