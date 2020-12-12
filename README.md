# SylvesterEquations.jl
Solvers for Sylvester and Lyapunov Equations

[![Build Status](https://travis-ci.org/olof3/SylvesterEquations.jl.svg?branch=master)](https://travis-ci.org/olof3/SylvesterEquations.jl)
[![codecov](https://codecov.io/gh/olof3/SylvesterEquations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/olof3/SylvesterEquations.jl)



`sylvc(A, B, C)` solves `AX + BX = C`  
`sylvd(A, B, C)` solves `AXB - X = C`  
`lyapc(A, Q)` solves `AX + BA' = -Q`  
`lyapd(A, Q)` solves `AXA' - X = -Q`  

The packages provides direct solvers based on straight-forward implementations of [Bartels-Stewart's algorithm](https://en.wikipedia.org/wiki/Bartels%E2%80%93Stewart_algorithm).
If there is no method `schur` for the `A` or the `B` matrix, there is a fallback to the "naive" method, this is useful for, e.g., symbolic equations.
