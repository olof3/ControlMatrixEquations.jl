@testset "arec" begin

# Scalars/complex scalars
@test arec(1, 2, 4, 2)[1] ≈ fill(2.0, (1,1))
@test arec(1, 2, 1, 5, 1)[1] ≈ fill(2.0, (1,1))
@test arec(1+im, 2, 1, 5, 1)[1] ≈ fill(2.0, (1,1))
@test arec(1, 2-im, 1, 5, 1+2im)[1] ≈ fill(2.0, (1,1))


@testset "CAREX" begin

carex_tests = [("1.4", Dict(), 1e-12),
               ("1.5", Dict(), 1e-12),
               ("2.1", Dict(:ϵ=>1e-2), 1e-7), # A little bit high tolerance..
               ("2.2", Dict(), 1e-7),
               ("2.3", Dict(), 1e-12),
               ("4.1", Dict(:n=>10), 1e-8) # Also not too impressive
               ]

for (id,kwargs,atol) in carex_tests
    @testset "CAREX $id $kwargs" begin
        A, B, Q, R, S, Xtrue = carex(id; kwargs...)
        X = arec(A, B, Q, R, S)[1]
        @test norm(arec_residual(X, A, B, Q, R, S)) ≈ 0 atol=atol # Doesn't check stabilizing solution...
    end
end

end # testset CAREX

end # testset
