function test_arec(A, B, Q, R, S=nothing; tol=1e-11)
    @testset "arec(::$(typeof(A)), ::$(typeof(B)), ::$(typeof(Q)), ::$(typeof(R)), ::$(typeof(S)))" begin
        X, cl_eigvals = arec(A, B, Q, R, S)
        residual = arec_residual(X, A, B, Q, R, S)    
        @test norm(residual) <= tol
        @test all(real.(cl_eigvals) .< 0)
        
        Λ = isnothing(S) ? eigvals(A - B*inv(R)*B'*X) : eigvals(A - B*inv(R)*(S' + B'*X))
        @test sort(real(Λ)) ≈ sort(real(cl_eigvals)) && sort(imag(Λ)) ≈ sort(imag(cl_eigvals)) # Sloppy but convenient
    end
end

function test_arecg(E, A, B, Q, R, S=nothing; tol=1e-11)
    @testset "arecg(::$(typeof(A)), ::$(typeof(A)), ::$(typeof(B)), ::$(typeof(Q)), ::$(typeof(R)), ::$(typeof(S)))" begin
        X, cl_eigvals = arecg(E, A, B, Q, R, S)
        residual = arecg_residual(X, E, A, B, Q, R, S)
        @test norm(residual) <= tol
        @test all(real.(cl_eigvals) .< 0)
        @test isposdef(X)
    end
end

##

Random.seed!(1000)

n = 30
m = 5
Ar = randn(n, n) - 10I
Ac = randn(ComplexF64, n, n) - 10I

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

@testset "arec" begin

test_arec(Ar, br, Q, r)
test_arec(Ar, Br, Q, R)
test_arec(Ar, bc, Q, r)
test_arec(Ar, Bc, Q, R)

test_arec(Ac, br, Q, r)
test_arec(Ac, Br, Q, R)
test_arec(Ac, bc, Q, r)
test_arec(Ac, Bc, Q, R)

test_arec(Ac, br, Qc, r)
test_arec(Ac, Br, Qc, R)
test_arec(Ac, bc, Qc, r)
test_arec(Ac, Bc, Qc, R)

test_arec(Ac, br, Qr, r, sr)
test_arec(Ac, Br, Qr, R, Sc)
test_arec(Ac, bc, Qr, r, sr)
test_arec(Ac, Bc, Qr, R, Sc)

# With uniform scaling
test_arec(Ar, Bc, Q, I)
test_arec(Ar, Bc, I, R)

end

##

# Apparently difficult to get good accuracy when there is only one input signal
@testset "arecg" begin

test_arecg(E, Ar, br, Q, r, tol=1e-5)
test_arecg(E, Ar, Br, Q, R, tol=1e-10)
test_arecg(E, Ar, bc, Q, r, tol=1e-5)
test_arecg(E, Ar, Bc, Q, R, tol=1e-10)

test_arecg(E, Ac, br, Q, r, tol=1e-5)
test_arecg(E, Ac, Br, Q, R, tol=1e-10)
test_arecg(E, Ac, bc, Q, r, tol=1e-5)
test_arecg(E, Ac, Bc, Q, R, tol=1e-10)

test_arecg(E, Ac, br, Qc, r, tol=1e-5)
test_arecg(E, Ac, Br, Qc, R, tol=1e-10)
test_arecg(E, Ac, bc, Qc, r, tol=1e-5)
test_arecg(E, Ac, Bc, Qc, R, tol=1e-10)

test_arecg(E, Ac, br, Qr, r, tol=1e-5)
test_arecg(E, Ac, Br, Qr, R, tol=1e-10)
test_arecg(E, Ac, bc, Qr, r, tol=1e-5)
test_arecg(E, Ac, Bc, Qr, R, tol=1e-10)

test_arecg(E, Ac, Bc, I, R, tol=1e-10)
test_arecg(I, Ac, Bc, Qr, R, tol=1e-12)
test_arecg(I, Ac, Bc, I, I, Sr, tol=1e-10)

end
