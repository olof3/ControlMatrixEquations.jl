import Random




function carex(id; ϵ=nothing,n::Int=-1,p=1,seed=0)
    Random.seed!(seed)
    if id == 0
        A = randn(n,n)
        B = randn(n,p)
        C = randn(p,n)
        Q = C'*C
        R = fill(1, (1, 1))
        S = randn(n,p)
        return A, B, Q, R, S, nothing
    elseif id == 1 # CAREX 2.1
        if ϵ === nothing; ϵ = 1e-3; end
        A = [1 0; 0 -2]
        B = [ϵ; 0]
        C = [1 1]
        Q = C'*C
        R = fill(1, (1, 1))
        return A, B, Q, R, nothing, nothing
    elseif id == 2 # CAREX 2.2
        if ϵ === nothing; ϵ = 1e-3; end
        A = [-0.01 0; 0 -0.02]
        B = [0.1 0; 0.001 0.01]
        C = [10 100]
        Q = C'*C
        R = [1+ϵ 1; 1 1]
        return A, B, Q, R, nothing, nothing
    elseif id == 3 # CAREX 2.3
        if ϵ === nothing; ϵ = 1e-3; end
        A = [0 0; ϵ 0]
        B = [0; 1]
        Q = Matrix{Float64}(I, 2, 2)
        R = fill(1, (1, 1))
        return A, B, Q, R, nothing, nothing
    elseif id == 5 # CAREX 2.5
        if ϵ === nothing; ϵ = 1e-3; end
        A = [3-ϵ 1; 4 2-ϵ]
        B = [1; 1]
        Q = [4ϵ - 11   2ϵ - 5; 2ϵ - 5 2ϵ - 2]
        R = fill(1, (1, 1))
        X = [2 1; 1 1]
        return A, B, Q, R, nothing, X
    elseif id == 6 # CAREX 4.1
        if n == -1; n = 20; end
        r = 0.001
        A = diagm(1 => ones(n-1))
        B = [zeros(n-1, 1); 1]
        Q = Matrix(I, n, n)
        R = fill(r, 1, 1)
        return A, B, Q, R, nothing, nothing
    elseif id == 7 # CAREX 1.4 (binary distillation column)
        A = [-0.991   0.529   0.0     0.0     0.0     0.0     0.0     0.0
              0.522  -1.051   0.596   0.0     0.0     0.0     0.0     0.0
              0.0     0.522  -1.118   0.596   0.0     0.0     0.0     0.0
              0.0     0.0     0.522  -1.548   0.718   0.0     0.0     0.0
              0.0     0.0     0.0     0.922  -1.64    0.799   0.0     0.0
              0.0     0.0     0.0     0.0     0.922  -1.721   0.901   0.0
              0.0     0.0     0.0     0.0     0.0     0.922  -1.823   1.021
              0.0     0.0     0.0     0.0     0.0     0.0     0.922  -1.943]

        B = [ 3.84   4.0   37.6   3.08   2.36   2.88   3.08   3.0
             -2.88  -3.04  -2.8  -2.32  -3.32  -3.82  -4.12  -3.96]'

        Q = [1.0  0.0  0.0  0.0  0.5  0.0  0.0  0.1
             0.0  1.0  0.0  0.0  0.1  0.0  0.0  0.0
             0.0  0.0  1.0  0.0  0.0  0.5  0.0  0.0
             0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0
             0.5  0.1  0.0  0.0  0.1  0.0  0.0  0.0
             0.0  0.0  0.5  0.0  0.0  0.1  0.0  0.0
             0.0  0.0  0.0  0.0  0.0  0.0  0.1  0.0
             0.1  0.0  0.0  0.0  0.0  0.0  0.0  0.1]

        R = I(2)
        return A, B, Q, R, nothing, nothing
    elseif id == 8 # CAREX 1.5 (tubular ammonia reactor)
        A = [ -4.019   5.12    0.0     0.0    -2.082     0.0     0.0    0.0    0.87
              -0.346   0.986   0.0     0.0    -2.34      0.0     0.0    0.0    0.97
              -7.909  15.407  -4.069   0.0    -6.45      0.0     0.0    0.0    2.68
             -21.816  35.606  -0.339  -3.87  -17.8       0.0     0.0    0.0    7.39
             -60.196  98.188  -7.907   0.34  -53.008     0.0     0.0    0.0   20.4
               0.0     0.0     0.0     0.0    94.0    -147.2     0.0   53.2    0.0
               0.0     0.0     0.0     0.0     0.0      94.0  -147.2    0.0    0.0
               0.0     0.0     0.0     0.0     0.0      12.8     0.0  -31.6    0.0
               0.0     0.0     0.0     0.0    12.8       0.0     0.0   18.8  -31.6]


        B = [0.01    0.003   0.009   0.024   0.068  0.0  0.0  0.0  0.0
            -0.011  -0.021  -0.059  -0.162  -0.445  0.0  0.0  0.0  0.0
            -0.151   0.0     0.0     0.0     0.0    0.0  0.0  0.0  0.0]'

        Q = I(9)
        R = I(3)
        return A, B, Q, R, nothing, nothing
    else
        error("Unknown test case")        
    end
end






function darex(id)
    if id == 1
        A = [0  1; 0 -1]
        B = [1 0; 2 1]
        Q = 1/11*[-4 -4; -4 7]
        R = [9 3; 3 1]
        S = [3 1; -1 7]
        return A, B, Q, R, S, nothing
    elseif id == 2
        A = [0  0; 0 1]
        B = [0; 1][:,:]
        Q = [1 2; 2 4]
        R = fill(1, 1, 1)
        return A, B, Q, R, nothing, nothing
    elseif id == 3
        A = [0  0.1 0; 0 0 0.1; 0 0 0]
        B = [1 0; 0 0; 0 1]
        Q = diagm([1e5, 1e3, -10])
        R = diagm([0, 1])
        return A, B, Q, R, nothing, nothing
    elseif id == 4
        A = [0.998   0.067   0     0
             -0.067  0.098   0     0
                0      0   0.998  0.153
                0      0  -0.153  0.998]
        B = [0.0033    0.02
             0.1     -7e-4
             0.4       0.0073
             -0.0028    0.1]
        Q = [1.87     0      0     -0.244
             0        0.744  0.205   0
             0        0.205  0.589   0
             -0.244  0      0      1.048]
        R = diagm([1, 1])
        return A, B, Q, R, nothing, nothing
    elseif id == 5
        A = [4 3; -9/2 -7/2]
        B = [1; -1][:,:]
        Q = [9 6; 6 4]
        δ = 1e-15
        R = fill(δ, 1, 1)
        return A, B, Q, R, nothing, nothing
    elseif id == 44
        A = 1e-3 * [
            984.75 -79.903 0.9054 -1.0765
            41.588  998.99 -35.855 12.684
            546.62  44.916 329.91    193.18
            2662.4 -100.45 -924.55 -263.25]
        B = 1e-3 * [
            37.112    7.361
            -870.51   0.093411
            -11984.0  -4.1378
            -31927.0  9.2535]
        R = Matrix(I, 2, 2)
        Q = 0.01*Matrix(I, 4, 4)
        return A, B, Q, R, nothing, nothing
    elseif id == 6 # CAREX 4.1, Hmm, this is essentially a time-delay, not integrators
        n = 50
        r = 0.001
        A = diagm(1 => ones(n-1))
        B = [zeros(n-1, 1); 1]
        Q = Matrix(I, n, n)
        R = fill(r, 1, 1)
        return A, B, Q, R, nothing, nothing
    elseif id == 8
        A, B, C, D = ControlSystems.ssdata(c2d(DemoSystems.fotd(), 0.1))
        R = fill(1e-7, (1, 1))
        Q = C'*C
        S = nothing
        return A, B, Q, R, nothing, nothing
    elseif id == 9
        n = 20
        m = 4
        A = randn(n, n)
        B = randn(n, m)
        Q = Symmetric(200I + randn(n, n))
        R = Symmetric(15I + randn(m, m))
        S = randn(n, m)
        return A, B, Q, R, S, nothing
    end
end
