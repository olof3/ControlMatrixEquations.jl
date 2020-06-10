# Tests of the top levels functions

A = 0.01*reshape(1:81, 9, 9) - 0.2I
B = 0.01*reshape(1:81, 9, 9) - 0.3I
X = diagm(-1 => 2ones(8), 0 => ones(9), 1 => 0.5ones(8))
Xs = (1:9)*(1:9)'

Ac = A + diagm(1 => (1:8)*im)
Bc = B + diagm(2 => (1:7)*im)
Xc = X + diagm(2 => (1:7)*im)
Xcs = Xs + diagm(1 => -(1:8)*im, 1 => (1:8)*im)

@test sylvc(A, B, sylvc_rhs(A, B, X)) ≈ X atol=1e-10
@test lyapc(A, lyapc_rhs(A, Xs)) ≈ Xs atol=1e-10
@test sylvc(Ac, Bc, sylvc_rhs(Ac, Bc, Xc)) ≈ Xc atol=1e-10
@test lyapc(Ac, lyapc_rhs(Ac, Xcs)) ≈ Xcs atol=1e-10

@test sylvd(A, B, sylvd_rhs(A, B, X)) ≈ X
@test lyapd(A, lyapd_rhs(A, Xs)) ≈ Xs
@test sylvd(Ac, Bc, sylvd_rhs(Ac, Bc, Xc)) ≈ Xc
@test lyapd(Ac, lyapd_rhs(Ac, Xcs)) ≈ Xcs



# Special cases
@test lyapd(0.5I, lyapd_rhs(0.5I, Xs)) ≈ Xs
@test lyapc(0.5I, lyapc_rhs(0.5I, Xs)) ≈ Xs
