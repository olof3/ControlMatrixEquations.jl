using ControlMatrixEquations
include("framework.jl")


@testset "Small equations (max 2 x 2) handled with static arrays" begin

for m=1:2
    for n=1:2
        A = reshape(1:m*m, m, m)
        B = reshape(2 .+ (1:n*n), n, n)
        E = reshape(-4 .+ (1:m*m), m, m)
        F = reshape(-2 .+ (1:n*n), n, n)

        C = reshape(10.0 .+ (1:m*n), m, n)

        X = ControlMatrixEquations._sylvc!(A, B, copy(C))
        @test norm(sylvc_residual(X, A, B, C)) < 1e-14

        X = ControlMatrixEquations._sylvd!(A, B, copy(C))
        @test norm(sylvd_residual(X, A, B, C)) < 1e-14

        X = ControlMatrixEquations._sylvg!(A, B, copy(C), E, F)
        @test norm(sylvg_residual(X, A, B, C, E, F)) < 1e-14
    end
end


end