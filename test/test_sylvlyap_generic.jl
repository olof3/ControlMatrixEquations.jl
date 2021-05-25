Random.seed!(200)

A, B, C, E, F = [big.(randn(3,3)) for k=1:5]
Q = Symmetric(randn(3,3))


@testset "Test fallback to naive for unsupported types" begin

X = sylvc(A, B, C)
@test norm(sylvc_residual(X, A, B, C)) < 1e-50
@test eltype(X) <: BigFloat

X = sylvd(A, B, C)
@test norm(sylvd_residual(X, A, B, C)) < 1e-50
@test eltype(X) <: BigFloat

X = sylvg(A, B, C, E, F)
@test norm(sylvg_residual(X, A, B, C, E, F)) < 1e-50
@test eltype(X) <: BigFloat

X = lyapc(A, Q)
@test norm(lyapc_residual(X, A, Q)) < 1e-50
@test eltype(X) <: BigFloat

X = lyapd(A, Q)
@test norm(lyapd_residual(X, A, Q)) < 1e-50
@test eltype(X) <: BigFloat

end


# FIXME: Add testcases with generic schur