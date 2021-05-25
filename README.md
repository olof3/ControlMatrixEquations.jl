# ControlMatrixEquations.jl
Solvers for Sylvester, Lyapunov and Riccati Equations


[![Build Status](https://github.com/olof3/ControlMatrixEquations.jl/workflows/CI/badge.svg)](https://github.com/olof3/ControlMatrixEquations.jl/actions?query=workflow%3ACI)

[![codecov](https://codecov.io/gh/olof3/ControlMatrixEquations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/olof3/ControlMatrixEquations.jl)

Numerous forms of the matrix equations below (in terms of symbols, signs, and transposes) occur in the literature and other software packages. The ones used for this package were chosen to be sensible in a control setting, but is WIP.


### Linear matrix equations

The packages provides direct solvers based a vanilla implementation of [Bartels–Stewart's algorithm](https://en.wikipedia.org/wiki/Bartels%E2%80%93Stewart_algorithm).
If there is no method `schur` for the `A` or the `B` matrix, there is a fallback to the "naive" (Kronecker product) method, this is useful for, e.g., symbolic equations.

#### Sylvester Equations
* `sylvc(A, B, C)` solves `AX + BX = C`
* `sylvd(A, B, C)` solves `AXB - X = C`
* `sylvg(A, B, C, E, F)` solves `AXB + EXF = C`

#### Lyapunov Equations
* `lyapc(A, Q)` solves `AX + XA' = -Q`
* `lyapd(A, Q)` solves `AXA' - X = -Q`
* `lyapc(A, Q, E)` solves `AXE' + EXA' = -Q` (will possibly be changed to to `lyapcg(E, A, Q)`)
* `lyapd(A, Q, E)` solves `AXA' - EXE' = -Q`

### Riccati Equations

Schur-factorization based Riccati solvers, including extended pencil versions that handle singular and near singular `R` matrices.

* `arec(A, B, Q, R, S)` solves the equation `A'X + XA - (XB + S)/R(XB + S)' + Q = 0`
* `ared(A, B, Q, R, S)` solves the equation `A'XA - X - (A'XB + S)/(B'XB + R)(A'XB + S)' + Q = 0`
* `arecg(E, A, B, Q, R, S)` solves the equation `A'XE + E'XA - (E'XB + S)/R(E'XB + S)' + Q = 0`
* `aredg(E, A, B, Q, R, S)` solves the equation `A'XA - E'XE - (A'XB + S)/(B'XB + R)(A'XB + S)' + Q = 0`

See Arnold & Laub (1984) "Generalized eigenproblem algorithms and software for algebraic Riccati equations."
