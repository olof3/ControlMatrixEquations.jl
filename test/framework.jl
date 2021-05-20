sylvc_residual(X, A, B, C) = A*X + X*B - C
sylvd_residual(X, A, B, C) = A*X*B - X - C
sylvg_residual(X, A, B, C, E, F) = A*X*B + E*X*F - C

lyapc_residual(X, A, Q) = A*X + X*A' + Q
lyapd_residual(X, A, Q) = A*X*A' - X + Q
lyapcg_residual(X, A, Q, E) = A*X*E' + E*X*A' + Q
lyapdg_residual(X, A, Q, E) = A*X*A' - E*X*E' + Q



arec_residual(X, A, B, Q, R, S) = A'*X + X*A - (X*B[:,:] + S[:,:])/ControlMatrixEquations.to_matrix(eltype(R), R)*(X*B[:,:] + S[:,:])' + Q
arec_residual(X, A, B, Q, R, S::Nothing=nothing) = arec_residual(X, A, B, Q, R, zeros(size(B)))
arecg_residual(X, E, A, B, Q, R, S) = A'*X*E + E'*X*A - (E'*X*B[:,:] + S[:,:])/ControlMatrixEquations.to_matrix(eltype(R), R)*(E'*X*B[:,:] + S[:,:])' + Q
arecg_residual(X, E, A, B, Q, R, S::Nothing=nothing) = arecg_residual(X, E, A, B, Q, R, zeros(size(B)))

ared_residual(X, A, B, Q, R, S) = A'*X*A - X - (A'*X*B[:,:] + S[:,:])/(B'*X*B + R)*(A'*X*B[:,:] + S[:,:])' + Q
ared_residual(X, A, B, Q, R, S::Nothing=nothing) = ared_residual(X, A, B, Q, R, zeros(size(B)))
