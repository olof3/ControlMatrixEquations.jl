"""
    _similarity_transform!(A, S::Diagonal)

Compute `S\\A*S` in place. (currently only supports Diagonal S)
"""
function _similarity_transform!(A, S::Diagonal)
    rmul!(A, S)
    ldiv!(S, A)
    return A
end


to_matrix(T, A::AbstractVector) = Matrix{T}(reshape(A, length(A), 1))
to_matrix(T, A::AbstractMatrix) = convert(Matrix{T}, A)
to_matrix(T, A::Number) = fill(T(A), 1, 1)
to_matrix(T, A::Adjoint{R, MT}) where {R<:Number, MT<:AbstractMatrix} = to_matrix(T, MT(A))



"""
    `_schurstructure(S::AbstractMatrix, ul=Union{Val{:U}, Val{:L}}) -> (block_indicies, nblocks)`
    `_schurstructure(S::AbstractMatrix, T::AbstractMatrix, ul=Union{Val{:U}, Val{:L}}) -> (block_indicies, nblocks)`

Return the indicies of the diagonal blocks (1x1 or 2x2) of a quasi-traingular Schur matrix `S`.
If `T` is provided, the joint structure of `S` and `T` is returned.

`ul` indicates if `S` (and `T`) are upper (`ul=Val(:U)`) or lower (`ul=Val(:L)`) triangular.

`block_indicies` contains the indices of the blocks

`nblocks` is the number of blocks

"""
_schurstructure(S::AbstractMatrix, ul::Union{Val{:U}, Val{:L}}) = _schurstructure(S, nothing, ul)
function _schurstructure(S::AbstractMatrix, T::Union{AbstractMatrix,Nothing}, ul::Union{Val{:U}, Val{:L}})
    n = size(S,1)

    block_indicies = Vector{UnitRange{Int}}(undef, n) # indices of each block (do not use UnitRange{Int32}, julia issue #40931)

    # Create function for checking if the offdiagonal element below (:U) / to the right (:L) of element (j,j) is zero
    is_offdiag_zero =
        if ul === Val(:U)
            T === nothing ? j -> iszero(S[j+1, j]) : j -> iszero(S[j+1, j]) && iszero(T[j+1, j])
        else
            T === nothing ? j -> iszero(S[j, j+1]) : j -> iszero(S[j, j+1]) && iszero(T[j, j+1])
        end

    j = 1 # current column (if :U) or row (if :L) index
    k = 0 # block number
    while j <= n
        j >= 2 && !is_offdiag_zero(j-1) && error("Matrix is not on Schur form / Incompatible Schur structures")

        k = k + 1
        # Check if the blocksize is 1 or 2
        if j == n || is_offdiag_zero(j)
            block_indicies[k] = j:j
            j += 1
        else
            block_indicies[k] = j:j+1
            j += 2
        end
    end

    # k now contains the total number of blocks
    return resize!(block_indicies, k), k
end


# FIXME: better handling of uniform scaling?!
issquare(A::Number) = true
issquare(A::AbstractMatrix) = size(A,1) == size(A,2)
function _check_lyap_inputs(A, Q, E=nothing)
    if !issquare(A); error("The A matrix must be square"); end
    if !ishermitian(Q); error("The Q matrix must be Hermitian"); end
    if size(Q, 1) != size(A, 1); error("The A and Q matrices must have the same dimensions"); end

    if E === nothing; return; end
    if size(A) != size(E); error("A and E matrices must have the same size"); end
end

function _check_sylv_inputs(A, B, C, E=nothing, F=nothing)
    if !issquare(A); error("The A matrix must be square"); end
    if !issquare(B); error("The B matrix must be square"); end
    if size(C, 1) != size(A, 1); error("The A and C matrices have inconsistent dimensions"); end
    if size(C, 2) != size(B, 2); error("The B and C matrices have inconsistent dimensions"); end

    if E === nothing; return; end
    if size(A) != size(E); error("A and E matrices must have the same size"); end
    if size(B) != size(F); error("B and F matrices must have the same size"); end
end


sylvcsoltype(A, B, C) = Base.promote_op((a,b,c) -> c / (a + b), eltype(A), eltype(B), eltype(C))
#sylvcsoltype(A, B, C) = Base.promote_op(sylvc, eltype(A), eltype(B), eltype(C))
sylvdsoltype(A, B, C) = Base.promote_op(sylvd, eltype(A), eltype(B), eltype(C))
