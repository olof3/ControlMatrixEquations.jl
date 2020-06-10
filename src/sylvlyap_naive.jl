# Methods for the naive method where the solution is simply vectorized
# This is necessary for matrices with element types that don't allow a Schur factorization,
# e.g., symbolic types

function sylvc(A, B, C, ::Val{:naive})
    Xv = ( kron(I(size(B,1)), A) + kron(transpose(B), I(size(A,1))) ) \ C[:]
    return reshape(Xv, size(C))
end
function sylvd(A, B, C, ::Val{:naive})
    # FIXME: just use I when bug in SymPy.jl has been fixed
    Xv = ( kron(transpose(B), A) - I(prod(size(C))) ) \ C[:] # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)
    return reshape(Xv, size(C))
end

# Mapping from Cartesian index into vector represenentation of Symmetric matrix
@inline sub2triidx(i,j) = (j >= i ? (j*(j-1))>>>1 + i : (i*(i-1))>>>1 + j)

function lyapc(A, Q, ::Val{:naive}) # Only works for real matrices A
    # Sets up and solves a system of equations for the upper triangular part of X
    # and solves that equation. This gives an n(n+1)/2 system instead of n^2
    # ONLY WORKS FOR REAL A!
    # Should be able to to base the discrete time version on this as well
    _check_lyap_inputs(A, Q)

    if !isreal(A); return sylvc(A, A', -Q); end # Should call sylvc

    n = size(Q, 1)
    nR = (n*(n+1)) >>> 1 # Ssize of the system to be solved

    R = zeros(eltype(A), nR, nR)
    minusQv = Vector{eltype(Q)}(undef, nR)

    for j=1:n
        for i=1:j
            m = sub2triidx(i,j) # Set up equation for X[i,j], m corresponds to index (i,j)
            minusQv[m] = -Q[i,j]
            # (AX + XA')[i,j] = sum(A[i,k]*X[k,j]) + sum(X[i,k]*A'[k,j])
            # the r = trinum[i]+r gives the kth element in the upper traingle
            # which correpsonds to X[i,j]
            for k=1:n
                R[m, sub2triidx(k,j)] += A[i,k]
                R[m, sub2triidx(i,k)] += A[j,k]
            end
        end
    end

    Xv = R \ minusQv

    # Fill the upper traingle of X from Xv
    X = [Xv[sub2triidx(i,j)] for i=1:n, j=1:n]

    return X
end
lyapd(A, Q, ::Val{:naive}) = sylvd(A, A', lmul!(-1,-Q), Val(:naive)) # No specilized method yet
