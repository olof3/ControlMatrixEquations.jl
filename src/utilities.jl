"""
    `_schurstructure(R::AbstractMatrix, ul=Union{Val{:U}, Val{:L}}) -> (b, d, nblocks)`

Return the block strucutre of an upper quasi-traingular Schur matrix `R`.
`ul` indicates if R is upper (`ul=Val(:U)`) or lower (`ul=Val(:L)`) triangular.


`d` contains the block sizees of each diagonal block (`1` or `2`)

`b` contains the indices of the blocks

`nblocks` is the number of blocks

"""
function _schurstructure(R::AbstractMatrix, ul=Val(:U)::Union{Val{:U}, Val{:L}})
    n = size(R,1)

    d = Vector{Int}(undef, n) # block sizes
    b = Vector{UnitRange{Int64}}(undef, n) # block indices

    j = 1 # column if ul=:U, row if ul=:L
    k = 0 # block number
    while j <= n
        k += 1
        if j == n
            d[k] = 1
        else
            if ul === Val(:U)
                d[k] = iszero(R[j+1, j]) ? 1 : 2
            else
                d[k] = iszero(R[j, j+1]) ? 1 : 2
            end
        end
        b[k] = j:j+d[k]-1
        j += d[k]
    end
    resize!(d, k)
    resize!(b, k)
    return d, b, k
end

# FIXME: better handling of uniform scaling?!
issquare(A::Number) = true
issquare(A::AbstractMatrix) = size(A,1) == size(A,2)
function _check_lyap_inputs(A, Q)
    if !issquare(A); error("The A matrix must be square"); end
    if !ishermitian(Q); error("The Q matrix must be Hermitian"); end
    if size(Q, 1) != size(A, 1); error("The A and Q matrices must have the same dimensions"); end
end

function _check_sylv_inputs(A, B, C)
    if !issquare(A); error("The A matrix must be square"); end
    if !issquare(B); error("The B matrix must be square"); end
    if size(C, 1) != size(A, 1); error("The A and C matrices have inconsistent dimensions"); end
    if size(C, 2) != size(B, 2); error("The B and C matrices have inconsistent dimensions"); end
end


sylvcsoltype(A, B, C) = Base.promote_op((a,b,c) -> c / (a + b), eltype(A), eltype(B), eltype(C))
#sylvcsoltype(A, B, C) = Base.promote_op(sylvc, eltype(A), eltype(B), eltype(C))
sylvdsoltype(A, B, C) = Base.promote_op(sylvd, eltype(A), eltype(B), eltype(C))
