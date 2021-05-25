@testset "errors" begin

    @test_throws ErrorException ControlMatrixEquations.lyapc_schur!(ones(2,2), [0 1; 0 1])
    @test_throws ErrorException ControlMatrixEquations.lyapd_schur!(ones(2,2), [0 1; 0 1])

    
    A, B, C, E, F = [randn(3,3) for k=1:5]
    # The following functions will only handle 2x2 matrices and smaller
    @test_throws ErrorException ControlMatrixEquations._sylvc!(A, B, C)
    @test_throws ErrorException ControlMatrixEquations._sylvd!(A, B, C)
    @test_throws ErrorException ControlMatrixEquations._sylvg!(A, B, C, E, F)
end
    