@testset "errors" begin

    @test_throws ErrorException ControlMatrixEquations.lyapc_schur!(ones(2,2), [0 1; 0 1])
    @test_throws ErrorException ControlMatrixEquations.lyapd_schur!(ones(2,2), [0 1; 0 1])
    
end
    