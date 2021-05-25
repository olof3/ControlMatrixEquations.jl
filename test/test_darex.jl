@testset "DAREX tests" begin

darex_tests = [("1.2", Dict(), 1e-11),
               ("1.3", Dict(), 1e-14),
               ("1.4", Dict(), 1e-15),
               ("1.5", Dict(), 1e-11),
               ("1.6", Dict(), 2e-13),
               ("2.1", Dict(:δ=>1), 1e-10),
               ("2.1", Dict(:δ=>1e6), 1e-7),
               ("4.1", Dict(:n=>10), 2e-13)
               ]

for (id,kwargs,atol) in darex_tests
    @testset "DAREX $id $kwargs" begin
        A, B, Q, R, S, Xtrue = darex(id; kwargs...)
        X = ared(A, B, Q, R, S)[1]
        @test norm(ared_residual(X, A, B, Q, R, S)) ≈ 0 atol=atol # Doesn't check stabilizing solution...
    end

end

end