using ControlMatrixEquations, Random
using Test

sylvc_residual(X, A, B, C) = A*X + X*B - C
sylvd_residual(X, A, B, C) = A*X*B - X - C
sylvg_residual(X, A, B, C, E, F) = A*X*B + E*X*F - C

lyapc_residual(X, A, Q) = A*X + X*A' + Q
lyapd_residual(X, A, Q) = A*X*A' - X + Q
lyapcg_residual(X, A, Q, E) = A*X*E' + E*X*A' + Q
lyapdg_residual(X, A, Q, E) = A*X*A' - E*X*E' + Q


##

Random.seed!(5000)
n = 30
Ar = randn(n, n)
Br = randn(n, n)
Cr = randn(n, n)
Er = randn(n, n)
Fr = randn(n, n)

Ac = randn(ComplexF64, n, n)
Bc = randn(ComplexF64, n, n)
Cc = randn(ComplexF64, n, n)
Ec = randn(ComplexF64, n, n)
Fc = randn(ComplexF64, n, n)

Qr = Symmetric(randn(n,n))
Qc = Hermitian(randn(ComplexF64,n,n))

##


function test_sylvc(A, B, C; tol=1e-9)
    X = sylvc(A, B, C)
    @test norm(sylvc_residual(X, A, B, C)) <= tol
end

@testset "sylvc" begin

test_sylvc(Ar, Br, Cr)
test_sylvc(Ac, Br, Cr)
test_sylvc(Ar, Bc, Cr)
test_sylvc(Ar, Br, Cc)
test_sylvc(Ac, Bc, Cr)
test_sylvc(Ac, Br, Cc)
test_sylvc(Ar, Bc, Cc)
test_sylvc(Ac, Bc, Cc)

end

##


function test_sylvd(A, B, C; tol=1e-9)
    X = sylvd(A, B, C)
    @test norm(sylvd_residual(X, A, B, C)) <= tol
end

@testset "sylvd" begin

test_sylvd(Ar, Br, Cr)
test_sylvd(Ac, Br, Cr)
test_sylvd(Ar, Bc, Cr)
test_sylvd(Ar, Br, Cc)
test_sylvd(Ac, Bc, Cr)
test_sylvd(Ac, Br, Cc)
test_sylvd(Ar, Bc, Cc)
test_sylvd(Ac, Bc, Cc)

end

##

function test_sylvg(A, B, C, E, F; tol=1e-9)
    X = sylvg(A, B, C, E, F)
    @test norm(sylvg_residual(X, A, B, C, E, F)) <= tol
end

@testset "sylvg" begin

# A subset of all the combinations

test_sylvg(Ar, Br, Cr, Er, Fr)
test_sylvg(Ar, Br, Cc, Er, Fr)
test_sylvg(Ar, Br, Cr, Ec, Fr)
test_sylvg(Ar, Br, Cc, Ec, Fr)

test_sylvg(Ac, Br, Cc, Er, Fc)
test_sylvg(Ac, Bc, Cc, Ec, Fr)
test_sylvg(Ac, Bc, Cr, Ec, Fc)
test_sylvg(Ac, Bc, Cc, Ec, Fc)

end

##


function test_lyapc(A, Q; tol=1e-9)
    X = lyapc(A, Q)
    @test norm(lyapc_residual(X, A, Q)) <= tol
end

@testset "lyapc" begin

test_lyapc(Ar, Qr)
test_lyapc(Ar, Qc)
test_lyapc(Ac, Qr)
test_lyapc(Ac, Qc)

end

##

function test_lyapcg(A, Q, E; tol=1e-9)
    X = lyapc(A, Q, E)
    @test norm(lyapcg_residual(X, A, Q, E)) <= tol
end

@testset "lyapcg" begin

test_lyapcg(Ar, Qr, Er)
test_lyapcg(Ar, Qc, Er)
test_lyapcg(Ac, Qr, Er)
test_lyapcg(Ac, Qc, Er)

test_lyapcg(Ar, Qr, Ec)
test_lyapcg(Ar, Qc, Ec)
test_lyapcg(Ac, Qr, Ec)
test_lyapcg(Ac, Qc, Ec)

end

##

function test_lyapd(A, Q; tol=1e-10)
    X = lyapd(A, Q)
    @test norm(lyapd_residual(X, A, Q)) <= tol
end

@testset "lyapd" begin

test_lyapd(Ar, Qr)
test_lyapd(Ar, Qc, tol=2e-10)
test_lyapd(Ac, Qr)
test_lyapd(Ac, Qc)

end

##

function test_lyapdg(A, Q, E; tol=1e-9)
    X = lyapd(A, Q, E)
    @test norm(lyapdg_residual(X, A, Q, E)) <= tol
end

@testset "lyapdg" begin

test_lyapdg(Ar, Qr, Er)
test_lyapdg(Ar, Qc, Er)
test_lyapdg(Ac, Qr, Er)
test_lyapdg(Ac, Qc, Er)

test_lyapdg(Ar, Qr, Ec)
test_lyapdg(Ar, Qc, Ec)
test_lyapdg(Ac, Qr, Ec)
test_lyapdg(Ac, Qc, Ec)

end

