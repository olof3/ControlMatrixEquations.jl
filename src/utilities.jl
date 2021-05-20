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
# Handle Adjoint Matrices
to_matrix(T, A::Adjoint{R, MT}) where {R<:Number, MT<:AbstractMatrix} = to_matrix(T, MT(A))
to_matrix(T, ::Nothing) = nothing # For matrices that are not supplied


"""
    `_schurstructure(R::AbstractMatrix, ul=Union{Val{:U}, Val{:L}}) -> (b, d, nblocks)`

Return the block strucutre of an upper quasi-traingular Schur matrix `R`.
`ul` indicates if R is upper (`ul=Val(:U)`) or lower (`ul=Val(:L)`) triangular.


`d` contains the block sizees of each diagonal block (`1` or `2`)

`b` contains the indices of the blocks

`nblocks` is the number of blocks

"""
_schurstructure(R::AbstractMatrix, ul::Union{Val{:U}, Val{:L}}) = _schurstructure(R, nothing, ul)
function _schurstructure(R::AbstractMatrix, S::Union{AbstractMatrix,Nothing}, ul::Union{Val{:U}, Val{:L}})
    n = size(R,1)

    d = Vector{Int32}(undef, n) # block sizes
    b = Vector{UnitRange{Int32}}(undef, n) # block indices

    # Create function for checking if the offdiagonal element below (:U) / to the right (:L) of element (j,j) is zero
    is_offdiag_zero =
        if ul === Val(:U)
            S === nothing ? j -> iszero(R[j+1, j]) : j -> iszero(R[j+1, j]) && iszero(S[j+1, j])
        else
            S === nothing ? j -> iszero(R[1, j+1]) : j -> iszero(R[j, j+1]) && iszero(S[j, j+1])
        end

    j = 1 # current column (if :U) or row (if :L) index
    k = 0 # block number
    while j <= n
        j >= 2 && !is_offdiag_zero(j-1) && error("Matrix is not on Schur form")
        k += 1
        if j == n
            d[k] = 1
        else
            d[k] = is_offdiag_zero(j) ? 1 : 2
        end
        b[k] = j:j+d[k]-1
        j += d[k]
    end
    resize!(d, k)
    resize!(b, k)
    # k now contains the total number of blocks
    return d, b, k
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
