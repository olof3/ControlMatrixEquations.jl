# There is a similar type in LowRankApprox.jl, perhaps we should use that instead of this.
"""
    MatrixLowRank{T, MT}

Represents a symmetric positive semideifnite matrix of the form

P = factor * factor'

"""
struct MatrixLowRank{T <: Number, MT <: AbstractMatrix{T}} <: AbstractMatrix{T}
    factor1::MT
    factor2::MT
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
