A, B, Q, R, S, Xtrue = carex(5, ϵ=1e-10)
@test arec(A, B[:,:], Q, R, S)[1] ≈ Xtrue
