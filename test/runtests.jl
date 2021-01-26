using ControlMatrixEquations
using Test, LinearAlgebra, Random


include("framework.jl")

my_tests = [
            "test_utilities",
            "test_sylvlyap_schur",
            "test_sylvlyap_naive",
            "test_sylvlyap",
            ]


# Compute C/Q from X to allow convenient check of solution
sylvc_rhs = (A, B, X) -> (A*X + X*B)
sylvd_rhs = (A, B, X) -> (A*X*B - X)
lyapc_rhs = (A, X) -> -Matrix(Hermitian(A*X + X*A'))
lyapd_rhs = (A, X) -> -Matrix(Hermitian(A*X*A' - X))

for test in my_tests
    println("In $test.jl:")
    include("$(test).jl")
end
