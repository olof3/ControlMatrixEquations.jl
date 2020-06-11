using SymPy

@vars a1 a2 real=true
@vars b1 b2 real=true
@vars c1 c2 c3 c4 real=true
@vars q1 q2 q3 real=true

# Need the definition tranpose(a) = a to make things work for complex (non-real) numbers, possibly a problem with SymPy
Base.transpose(f::Sym) = f

A = [a1 a2; 2 a2+1]
B = [b1 1; b2 2]
C = [c1 c2; c3 c4]
Q = [q1 q2; q2 q3]

X = sylvc(A, B, C)
simplify.(A*X + X*B) == C

X = lyapc(A, Q)
simplify.(A*X + X*A') == -Q

X = sylvd(A, B, C)
simplify.(A*X*B - X) == C

X = lyapd(A, Q)
simplify.(A*X*A' - X) == -Q
