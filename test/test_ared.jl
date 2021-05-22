function test_ared(A, B, Q, R, S=nothing; tol=1e-11)
    @testset "ared(::$(typeof(A)), ::$(typeof(B)), ::$(typeof(Q)), ::$(typeof(R)), ::$(typeof(S)))" begin
        X, cl_eigvals = ared(A, B, Q, R, S)
        residual = ared_residual(X, A, B, Q, R, S)    
        @test norm(residual) <= tol
        @test all(abs.(cl_eigvals) .< 1)
        @test isposdef(X)

        R isa UniformScaling && return
        Λ = isnothing(S) ? eigvals(A - B * ((B'*X*B .+ R) \ (B'*X*A))) : eigvals(A - B * ((B'*X*B .+ R) \ (B'*X *A + S')))
        @test sort(real(Λ)) ≈ sort(real(cl_eigvals)) && sort(imag(Λ)) ≈ sort(imag(cl_eigvals)) # Sloppy but convenient        
    end
end

function test_aredg(E, A, B, Q, R, S=nothing; tol=1e-11)
    @testset "aredg(::$(typeof(A)), ::$(typeof(A)), ::$(typeof(B)), ::$(typeof(Q)), ::$(typeof(R)), ::$(typeof(S)))" begin
        X, cl_eigvals = aredg(E, A, B, Q, R, S)
        residual = aredg_residual(X, E, A, B, Q, R, S)
        @test norm(residual) <= tol
        @test all(abs.(cl_eigvals) .< 1)
        @test isposdef(X)
    end
end

##

Random.seed!(900)

n = 30
m = 5
Ar = 0.1*randn(n, n)
Ac = 0.1*randn(ComplexF64, n, n)

br = randn(n)
bc = randn(ComplexF64, n)
Br = randn(n, m)
Bc = randn(ComplexF64, n, m)

Q = I(n)

Qr = Symmetric(randn(n,n)) + 20I
Qc = Hermitian(randn(n,n)) + 20I

r = 1.0
R = I(m)
Rr = Symmetric(randn(n,n)) + 20I
Rc = Hermitian(randn(n,n)) + 20I

E = randn(n, n) + 4I

sr = randn(size(br))
sc = randn(ComplexF64, size(br))
Sr = 0.1randn(size(Bc))
Sc = 0.1randn(ComplexF64, size(Bc))

##

@testset "ared" begin

test_ared(Ar, br, Q, r, tol=1e-9)
test_ared(Ar, Br, Q, R, tol=1e-9)
test_ared(Ar, bc, Q, r, tol=1e-9)
test_ared(Ar, Bc, Q, R, tol=1e-9)

test_ared(Ac, br, Q, r, tol=1e-9)
test_ared(Ac, Br, Q, R, tol=1e-9)
test_ared(Ac, bc, Q, r, tol=1e-9)
test_ared(Ac, Bc, Q, R, tol=1e-9)

test_ared(Ac, br, Qc, r, tol=1e-9)
test_ared(Ac, Br, Qc, R, tol=1e-9)
test_ared(Ac, bc, Qc, r, tol=1e-9)
test_ared(Ac, Bc, Qc, R, tol=1e-9)

test_ared(Ac, br, Qr, r, sr, tol=1e-9)
test_ared(Ac, Br, Qr, R, Sc, tol=1e-9)
test_ared(Ac, bc, Qr, r, sr, tol=1e-9)
test_ared(Ac, Bc, Qr, R, Sc, tol=1e-9)

# With uniform scaling
test_ared(Ar, Bc, Q, I, tol=1e-9)
test_ared(Ar, Bc, I, R, tol=1e-9)

end

##

# Apparently difficult to get good accuracy when there is only one input signal
@testset "aredg" begin

test_aredg(E, Ar, br, Q, r, tol=1e-9)
test_aredg(E, Ar, Br, Q, R, tol=1e-9)
test_aredg(E, Ar, bc, Q, r, tol=1e-9)
test_aredg(E, Ar, Bc, Q, R, tol=1e-9)

test_aredg(E, Ac, br, Q, r, tol=1e-9)
test_aredg(E, Ac, Br, Q, R, tol=1e-9)
test_aredg(E, Ac, bc, Q, r, tol=1e-9)
test_aredg(E, Ac, Bc, Q, R, tol=1e-9)

test_aredg(E, Ac, br, Qc, r, tol=1e-9)
test_aredg(E, Ac, Br, Qc, R, tol=1e-9)
test_aredg(E, Ac, bc, Qc, r, tol=1e-9)
test_aredg(E, Ac, Bc, Qc, R, tol=1e-9)

test_aredg(E, Ac, br, Qr, r, tol=1e-9)
test_aredg(E, Ac, Br, Qr, R, tol=1e-9)
test_aredg(E, Ac, bc, Qr, r, tol=1e-9)
test_aredg(E, Ac, Bc, Qr, R, tol=1e-9)

test_aredg(E, Ac, Bc, I, R, tol=1e-9)
test_aredg(I, Ac, Bc, Qr, R, tol=1e-9)
test_aredg(I, Ac, Bc, I, I, Sr, tol=1e-9)

end
