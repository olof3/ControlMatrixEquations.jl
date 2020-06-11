#
# Tests with small trinagular matrices of various types
#
A3 = [2 0 0; 3 4 0; 5 6 7]
B3 = [2 3 4; 0 1 1; 0 0 2]
X3 = reshape(1:9, 3, 3)

for (A, B, X, tsname) in [(fill(2.0, 1, 1), fill(3.0, 1, 1), fill(1.0, 1, 1), "1x1"),
                  (diagm([2.0,3]), diagm([4.0,5]), [1.0 2; 3 4], "2x2 diagonal"),
                  ([-2 0; 2 3], [3.0 4; 0 5], [1.0 2; 3 4], "2x2 tringular"),
                  ([-2.0 0; 2 3], fill(5.0, 1, 1), [1.0; 2][:,:], "2x1"),
                  ([-2 0; 2 3], fill(5, 1, 1), [1; 2][:,:], "2x1"),
                  (fill(2.0, 1, 1), [3.0 4; 0 5], [1.0 2], "1x2"),
                  (float(A3), float(B3), float(X3), "3x3 (Float64)"),
                  (float(A3), float(B3), X3, "3x3 (Float64)"),
                  (A3, B3, X3, "3x3"),
                  ]

    # Generate complex version, while keeping the structure
    Ac = A .* (1 + im)
    Bc = B .* (1 + im)
    Xc = X .* (1 + im)

    @testset "sylv(c/d)_schur!, $(tsname) $(eltype.((A, B, X)))" begin
        @test SylvesterEquations.sylvc_schur!(A, B, sylvc_rhs(A, B, X)) ≈ X
        @test SylvesterEquations.sylvc_schur!(Ac, B, sylvc_rhs(Ac, B, X)) ≈ X
        @test SylvesterEquations.sylvc_schur!(A, B, sylvc_rhs(A, B, Xc)) ≈ Xc
        @test SylvesterEquations.sylvc_schur!(Ac, Bc, sylvc_rhs(Ac, Bc, Xc)) ≈ Xc

        @test SylvesterEquations.sylvd_schur!(A, B, sylvd_rhs(A, B, X)) ≈ X
        @test SylvesterEquations.sylvd_schur!(Ac, B, sylvd_rhs(Ac, B, X)) ≈ X
        @test SylvesterEquations.sylvd_schur!(A, B, sylvd_rhs(A, B, Xc)) ≈ Xc
        @test SylvesterEquations.sylvd_schur!(Ac, Bc, sylvd_rhs(Ac, Bc, Xc)) ≈ Xc

    end

    if size(X,1) == size(X,2)
        Xh = X + X'
        Xch = Xc + Xc'

        @testset "lyap(c/d)_schur!, $(tsname) $(eltype.((A, Xh)))" begin
            @test SylvesterEquations.lyapc_schur!(A, lyapc_rhs(A, Xh)) ≈ Xh
            @test SylvesterEquations.lyapc_schur!(Ac, lyapc_rhs(Ac, Xh)) ≈ Xh
            @test SylvesterEquations.lyapc_schur!(A, lyapc_rhs(A, Xch)) ≈ Xch
            @test SylvesterEquations.lyapc_schur!(Ac, lyapc_rhs(Ac, Xch)) ≈ Xch

            @test SylvesterEquations.lyapd_schur!(A, lyapd_rhs(A, Xh)) ≈ Xh
            @test SylvesterEquations.lyapd_schur!(Ac, lyapd_rhs(Ac, Xh)) ≈ Xh
            @test SylvesterEquations.lyapd_schur!(A, lyapd_rhs(A, Xch)) ≈ Xch
            @test SylvesterEquations.lyapd_schur!(Ac, lyapd_rhs(Ac, Xch)) ≈ Xch
        end
    end
end
