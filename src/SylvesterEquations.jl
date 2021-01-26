module SylvesterEquations

using LinearAlgebra
using StaticArrays

export sylvc, sylvd
export lyapc, lyapd
export arec, ared

export HermitianLowRank

const AbstractNumOrArray = Union{Number, AbstractVecOrMat}

include("HermitianLowRank.jl")
include("utilities.jl")

include("sylv_2x2.jl")
include("sylvlyap_bartstew.jl")
include("sylvlyap_naive.jl")

include("riccati.jl")


# Infer the algorithm based on the type of the problem data
function _infer_sylvalg(A, B, C)
    if hasmethod(schur!, (typeof(A),)) && hasmethod(schur!, (typeof(B),))
        return Val(:bartstew)
    else
        return Val(:naive)
    end
end

sylvc(A, B, C) = sylvc(A, B, C, Val(:infer))
sylvd(A, B, C) = sylvd(A, B, C, Val(:infer))
lyapc(A, Q) = lyapc(A, Q, Val(:infer))
lyapd(A, Q) = lyapd(A, Q, Val(:infer))


lyapc(A::UniformScaling, Q) = Q ./ (-2*real(A.λ))
lyapd(A::UniformScaling, Q) = Q ./ (1 - abs2(A.λ))


sylvc(A, B, C, ::Val{:infer}; kwargs...) = sylvc(A, B, C, _infer_sylvalg(A,B,C); kwargs...)
sylvd(A, B, C, ::Val{:infer}; kwargs...) = sylvd(A, B, C, _infer_sylvalg(A,B,C); kwargs...)
lyapc(A, Q, ::Val{:infer}; kwargs...) = lyapc(A, Q, _infer_sylvalg(A,A,Q); kwargs...)
lyapd(A, Q, ::Val{:infer}; kwargs...) = lyapd(A, Q, _infer_sylvalg(A,A,Q); kwargs...)


# The following should preferably be fixed in LinearAlgebra, there is an issue posted...
LinearAlgebra.schur(A::AbstractMatrix{T}) where T = schur!(LinearAlgebra.copy_oftype(A, LinearAlgebra.eigtype(T)))

end # module
