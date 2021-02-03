using ControlMatrixEquations


n = 50
A = triu(randn(n,n)) - 50I
b = randn(n,3)


#C = lyapc(A, b, Val(:hammarling))

C = ControlMatrixEquations.lyapham(A, b)

X = C*C'


maximum(eigvals(A))

residual = A*X + X*A' + b*b'

norm(residual)

@btime C = lyapc(A, copy(b), Val(:hammarling))

Juno.@profiler lyapc(A, copy(b), Val(:hammarling))

C = lyapc(A, copy(b), Val(:hammarling))
X = C*C'
@time X2 = lyapc(A, b*b')



#display(X - X2)
norm(X - X2)
