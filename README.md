# ControlMatrixEquations.jl
Solvers for Sylvester, Lyapunov and Riccati Equations

[![codecov](https://codecov.io/gh/olof3/ControlMatrixEquations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/olof3/ControlMatrixEquations.jl)



`sylvc(A, B, C)` solves `AX + BX = C`  
`sylvd(A, B, C)` solves `AXB - X = C`  
`lyapc(A, Q)` solves `AX + BA' = -Q`  
`lyapd(A, Q)` solves `AXA' - X = -Q`  

The packages provides direct solvers based on straight-forward implementations of [Bartels-Stewart's algorithm](https://en.wikipedia.org/wiki/Bartels%E2%80%93Stewart_algorithm).
If there is no method `schur` for the `A` or the `B` matrix, there is a fallback to the "naive" method, this is useful for, e.g., symbolic equations.
