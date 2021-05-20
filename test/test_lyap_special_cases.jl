@testset "lyap*, Special cases" begin

n = 7

AI = (2.5 + im)I
Q = Hermitian(randn(n,n))
A = randn(n, n)
EI = (1-2im)I
E = randn(n,n)


X = lyapc(AI, Q)
@test norm(lyapc_residual(X, AI, Q)) < 1e-12

X = lyapd(AI, Q)
@test norm(lyapd_residual(X, AI, Q)) < 1e-12

X = lyapc(A, Q, EI)
@test norm(lyapcg_residual(X, A, Q, EI)) < 1e-12

X = lyapd(A, Q, EI)
@test norm(lyapdg_residual(X, A, Q, EI)) < 1e-12

X = lyapc(AI, Q, E)
@test norm(lyapcg_residual(X, AI, Q, E)) < 1e-12

X = lyapd(AI, Q, E)
@test norm(lyapdg_residual(X, AI, Q, E)) < 1e-12

X = lyapc(AI, Q, EI)
@test norm(lyapcg_residual(X, AI, Q, EI)) < 1e-12

X = lyapd(AI, Q, EI)
@test norm(lyapdg_residual(X, AI, Q, EI)) < 1e-12



end