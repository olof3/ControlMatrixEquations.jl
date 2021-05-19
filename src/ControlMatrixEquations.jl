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

sylvc(A, B, C) = sylvc(A, B, C, Val(:infer))
sylvd(A, B, C) = sylvd(A, B, C, Val(:infer))
lyapc(A, Q) = lyapc(A, Q, Val(:infer))
lyapd(A, Q) = lyapd(A, Q, Val(:infer))

lyapc(A, Q, E) = lyapc(A, Q, E, Val(:bartstew))
lyapd(A, Q, E) = lyapd(A, Q, E, Val(:bartstew))

lyapc(A::UniformScaling, Q) = Q ./ (-2*real(A.λ))
lyapd(A::UniformScaling, Q) = Q ./ (1 - abs2(A.λ))


sylvc(A, B, C, ::Val{:infer}; kwargs...) = sylvc(A, B, C, _infer_sylvalg(A,B,C); kwargs...)
sylvd(A, B, C, ::Val{:infer}; kwargs...) = sylvd(A, B, C, _infer_sylvalg(A,B,C); kwargs...)
lyapc(A, Q, ::Val{:infer}; kwargs...) = lyapc(A, Q, _infer_sylvalg(A,A,Q); kwargs...)
lyapd(A, Q, ::Val{:infer}; kwargs...) = lyapd(A, Q, _infer_sylvalg(A,A,Q); kwargs...)


## Special cases
lyapc(a::Number, q::Number) = sylvc(a, a', -q)
lyapd(a::Number, q::Number) = sylvd(a, a', -q)
lyapc(a::Number, q::Number, e::Number) = sylvg(a, e', -q, e, a')
lyapd(a::Number, q::Number, e::Number) = sylvg(a, a', -q, -e, e')



lyapc(A::UniformScaling, Q::AbstractMatrix) = -Q / (2real(A.λ))
lyapd(A::UniformScaling, Q::AbstractMatrix) = -Q / (abs2(A.λ) - 1)
lyapc(A, Q, E::UniformScaling) = lyapc(conj(E.λ) * A, Q)
lyapd(A, Q, E::UniformScaling) = lyapd(A / E.λ, Q/abs2(E.λ))
lyapc(A::UniformScaling, Q, E) = lyapc(conj(A.λ) * E, Q)
lyapd(A::UniformScaling, Q, E) = lyapd(E / A.λ, -Q/abs2(A.λ))
lyapc(A::UniformScaling, Q, E::UniformScaling) = -Q / (2real(A.λ * conj(E.λ)))
lyapd(A::UniformScaling, Q, E::UniformScaling) = -Q / (abs2(A.λ) - abs2(E.λ))

# Could probably add a few more of these, FIXME: add further tests






# The following should preferably be fixed in LinearAlgebra, there is an issue posted...
LinearAlgebra.schur(A::LinearAlgebra.AdjOrTrans{T}) where T = schur!(LinearAlgebra.copy_oftype(A, LinearAlgebra.eigtype(T)))
LinearAlgebra.schur(A::AbstractMatrix{T}) where T = schur!(LinearAlgebra.copy_oftype(A, LinearAlgebra.eigtype(T)))
LinearAlgebra.schur(A::AbstractMatrix{T1}, B::AbstractMatrix{T2}) where {T1, T2} = schur!(LinearAlgebra.copy_oftype(A, LinearAlgebra.eigtype(T1)), LinearAlgebra.copy_oftype(B, LinearAlgebra.eigtype(T2)))

# Note: not using inplace
LinearAlgebra.schur(A::LinearAlgebra.AdjOrTrans{T1}, B::LinearAlgebra.AdjOrTrans{T2}) where {T1, T2} = schur(LinearAlgebra.copy_oftype(A, LinearAlgebra.eigtype(T1)), LinearAlgebra.copy_oftype(B, LinearAlgebra.eigtype(T2)))


end # module
