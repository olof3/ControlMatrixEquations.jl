# There is a similar type in LowRankApprox.jl, perhaps we should use that instead of this.
"""
    MatrixLowRank{T, MT}

Represents a symmetric positive semideifnite matrix of the form

P = factor1 * factor2'

"""
struct MatrixLowRank{T <: Number, MT <: AbstractVecOrMat{T}} <: AbstractMatrix{T}
    factor1::MT
    factor2::MT

    function MatrixLowRank{T,MT}(factor1, factor2) where {T,MT}
        size(factor1)[1] != size(factor2)[1] &&  error("Incompatible dimensions of factor1 and factor2")
        new{T,MT}(factor1, factor2)
    end
end
function MatrixLowRank(factor1::AbstractVecOrMat, factor2::AbstractVecOrMat)
    MT = promote_type(typeof(factor1), typeof(factor2))
    MatrixLowRank{eltype(MT),MT}(MT(factor1), MT(factor2))
end

"""
    HermitianLowRank{T, MT}

Represents a symmetric positive semideifnite matrix of the form

P = factor * factor'

"""
struct HermitianLowRank{T <: Number, MT <: AbstractMatrix{T}} <: AbstractMatrix{T}
    factor::MT
end
Base.convert(Matrix, Z::HermitianLowRank) = Z.factor * Z.factor'
Base.adjoint(Z::HermitianLowRank) = Z

Base.size(Z::HermitianLowRank) = (size(Z.factor,1), size(Z.factor,1))
Base.getindex(Z::HermitianLowRank, i, j) = conj(dot(Z.factor[i, :], Z.factor[j, :]))
