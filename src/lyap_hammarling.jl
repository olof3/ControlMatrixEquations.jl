# Based on Sorensen, Algorithm 3.1
"""
    lyapc(A, B, ::Val{:hammarling}) -> C

returns an upper-trinagular Cholesky factor `C` such that
`X = C'*C` solves `A*X + X*B + B'*B = 0`.

The matrix `A` must be Hurwitz.

Based on

[1] **Sorensen, D. C., & Zhou, Y.** (2003). Direct methods for matrix Sylvester and
Lyapunov equations. Journal of Applied Mathematics, 2003(6), 277-303.
"""

# FIXME, double check transposes/conjugation etc

# FIXME: Only handles the case that A is Upper triangular
function lyapc(A, B0, ::Val{:hammarling})

    B = copy(B0)

    n = size(A,1)

    U = zeros(size(A))
    bh = Vector{eltype(B)}(undef, size(B,2))

    for j=n:-1:2
        bh .= B[j, :]

        μ = norm(bh)
        μ1 = sqrt(-2*real(A[j,j]))
        τ = μ / μ1

        if μ == 0; continue; end

        u = view(U, 1:j-1, j)

        # store -btmp in the space allocated for u
        @views u .= -τ .* A[1:j-1, j]
        @views mul!(u, B[1:j-1, :], conj!(bh), -1/τ, 1)

        # replace btmp so that U[1:j-1,j] contains the solution
        @views _ldiv!(UpperTriangular(A[1:j-1, 1:j-1]), u, shift=conj(A[j,j]))

        U[j, j] = τ

        # rank-1 update of B
        @views mul!(B[1:j-1, :], u, transpose(B[j,:]), -1/τ, 1)
    end

    U[1,1] = norm(B[1, :]) / sqrt(-2*real(A[1,1]))

    return Cholesky(U, :U, 0)
end


#U[1:j-1,j] .= (A[1:j-1, 1:j-1] + A[j,j]*I) \ ( -(B[1:j-1, :] * B[j, :] * (1/τ) .+ τ * A[1:j-1, j]) )


function _ldiv!(A::LowerTriangular, b::AbstractVector; shift=0)
    b[1] /= (A[1,1] + shift)
    @inbounds for k=2:size(A,1)
        for l=1:k-1
            b[k] -= A[k,l] * b[l]
        end
        b[k] /= (A[k,k] + shift)
    end
    b
end

function _ldiv!(A::UpperTriangular, b::AbstractVector; shift=0)
    b[end] /= (A[end,end] + shift)
    @inbounds for k=size(A,1)-1:-1:1
        for l=k+1:size(A,2)
            b[k] -= A[k,l] * b[l]
        end
        b[k] /= (A[k,k] + shift)
    end
    b
end
