module ControlMatrixEquations

using LinearAlgebra
using StaticArrays

export sylvc, sylvd
export sylvg

export lyapc, lyapd
export lyapcg, lyapdg

export arec, ared
export arecg, aredg

const AbstractNumOrArray = Union{Number, AbstractVecOrMat}

include("utilities.jl")

include("sylv_2x2.jl")
include("sylvlyap_bartstew.jl")
include("sylvlyap_naive.jl")

include("riccati.jl")


# Infer the algorithm based on the type of the problem data
function _infer_sylvalg(A, B, C)
    T = promote_type(eltype(A), eltype(B), Float64)
    if hasmethod(schur!, (Matrix{T},))
        return Val(:bartstew)
    else
        return Val(:naive)
    end
end

sylvc(A, B, C) = sylvc(A, B, C, _infer_sylvalg(A,B,C))
sylvd(A, B, C) = sylvd(A, B, C, _infer_sylvalg(A,B,C))
lyapc(A, Q) = lyapc(A, Q, _infer_sylvalg(A,A,Q))
lyapd(A, Q) = lyapd(A, Q, _infer_sylvalg(A,A,Q))

sylvg(A, B, C, E, F) = sylvg(A, B, C, E, F, Val(:bartstew))
lyapc(A, Q, E) = lyapc(A, Q, E, Val(:bartstew))
lyapd(A, Q, E) = lyapd(A, Q, E, Val(:bartstew))


## Special cases
lyapc(a::Number, q::Number) = sylvc(a, a', -q)
lyapd(a::Number, q::Number) = sylvd(a, a', -q)
lyapc(a::Number, q::Number, e::Number) = sylvg(a, e', -q, e, a')
lyapd(a::Number, q::Number, e::Number) = sylvg(a, a', -q, -e, e')


lyapc(A::UniformScaling, Q) = Q ./ (-2*real(A.λ))
lyapd(A::UniformScaling, Q) = Q ./ (1 - abs2(A.λ))
lyapc(A, Q, E::UniformScaling) = lyapc(conj(E.λ) * A, Q)
lyapd(A, Q, E::UniformScaling) = lyapd(A / E.λ, Q/abs2(E.λ))
lyapc(A::UniformScaling, Q, E) = lyapc(conj(A.λ) * E, Q)
lyapd(A::UniformScaling, Q, E) = lyapd(E / A.λ, -Q/abs2(A.λ))
lyapc(A::UniformScaling, Q, E::UniformScaling) = -Q / (2real(A.λ * conj(E.λ)))
lyapd(A::UniformScaling, Q, E::UniformScaling) = -Q / (abs2(A.λ) - abs2(E.λ))

# Could probably add a few more of these, i.e., for sylv, FIXME: add further tests




# The following should preferably be fixed in LinearAlgebra, there is an issue posted...
LinearAlgebra.schur(A::AbstractMatrix{T}) where T = schur!(LinearAlgebra.copy_oftype(A, LinearAlgebra.eigtype(T)))
function LinearAlgebra.schur(A::AbstractMatrix{T1}, B::AbstractMatrix{T2}) where {T1, T2}
    T = promote_type(LinearAlgebra.eigtype(T1), LinearAlgebra.eigtype(T2))
    schur!(LinearAlgebra.copy_oftype(A, T), LinearAlgebra.copy_oftype(B, T))
end


end # module
