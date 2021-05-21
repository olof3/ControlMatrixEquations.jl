Random.seed!(2000)

n = 7

AI = (2.5 + im)I
Q = Hermitian(randn(n,n))
A = randn(n, n)
EI = (1-2im)I
E = randn(n,n)




@testset "lyap*, UniformScaling inputs" begin

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


@testset "lyap*, UniformScaling inputs" begin

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


@testset "lyap*, scalar inputs" begin

ac = 3 + im
q = 1.5
ec = 2 - im

x = lyapc(ac, q)
@test abs(lyapc_residual(x, ac, q)) < 1e-15

x = lyapc(ac, q, ec)
@test abs(lyapcg_residual(x, ac, q, ec)) < 1e-15

x = lyapd(ac, q)
@test abs(lyapd_residual(x, ac, q)) < 1e-15

x = lyapd(ac, q, ec)
@test abs(lyapdg_residual(x, ac, q, ec)) < 1e-15

end