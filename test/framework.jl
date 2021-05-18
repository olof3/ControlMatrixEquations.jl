
arecg_residual(X, E, A, B, Q, R, S) = A'*X*E + E'*X*A - (E'*X*B + S)/ControlMatrixEquations.to_matrix(eltype(R), R)*(E'*X*B + S)' + Q
arec_residual(X, A, B, Q, R, S) = A'*X + X*A - (X*B + S)/ControlMatrixEquations.to_matrix(eltype(R), R)*(X*B + S)' + Q
arec_residual(X, A, B, Q, R, S::Nothing=nothing) = arec_residual(X, A, B, Q, R, zeros(size(B)))

ared_residual(X, A, B, Q, R, S) = A'*X*A - X - (A'*X*B + S)/(B'*X*B + R)*(A'*X*B + S)' + Q
ared_residual(X, A, B, Q, R, S::Nothing=nothing) = ared_residual(X, A, B, Q, R, zeros(size(B)))
