# (A, B/R*B) must be stabilizable
"""
    arec(A, B, Q, R, S=nothing; force_esp=false)

Find the stabilizing solution `X` to the continuous-time Riccati equation
`A'X + XA - (XB + S)/R(XB + S)' + Q`
"""
function arec(A::AbstractNumOrArray, B::AbstractNumOrArray, Q::AbstractNumOrArray, R::AbstractNumOrArray, S=nothing; force_esp=false)
    A, B, Q, R, S = _check_ARE_inputs(A, B, Q, R, S)

    if cond(R) <= 1e4 && !force_esp# FIXME: Think this condition through
        if S == nothing
            return arec_noinv(A, B/R*B', Q)
        else
            return arec_noinv(A - B/R*S', B/R*B', Hermitian(Q - S/R*S'))
        end
    else
        return _are_extended_pencil(Val(:c), I, A, B, Q, R, S)
    end
end

"""
    ared(A, B, Q, R, S=nothing; force_esp=false)

Find the stabilizing solution `X` to the discrete-time Riccati equation
`A'XA - X - (A'XB + S)/(B'XB + R)(A'XB + S)' + Q`
"""
function ared(A::AbstractNumOrArray, B::AbstractNumOrArray, Q::AbstractNumOrArray, R::AbstractNumOrArray, S=nothing; force_esp=false)
    A, B, Q, R, S = _check_ARE_inputs(A, B, Q, R, S)

    if cond(R) <= 1e4 && !force_esp
        if S == nothing
            return _ared(A, B, Q, R)
        else
            return _ared(A - B/R*S', B, Q - S/R*S', R)
        end
    else
        return _are_extended_pencil(Val(:d), I, A, B, Q, R, S)
    end
end


function _check_ARE_inputs(A, B, Q, R, S)
    # Convert all inputs to Matrix{T}
    T = promote_type(eltype(A), eltype(B), eltype(Q), isnothing(S) ? Union{} : eltype(S))

    A, B, Q, R, S = to_matrix(T, A), to_matrix(T, B), to_matrix(T, Q), to_matrix(T, R), to_matrix(T, S)

    n, m = size(B)

    size(A) != (n, n) && error("A and B have mismatched sizes $(size(A)) vs $(size(B))")
    size(Q) != (n, n) && error("A and Q have mismatched sizes $(size(A)) vs $(size(Q))")
    size(R) != (m, m) && error("B and R have mismatched sizes $(size(B)) vs $(size(R))")
    !isnothing(S) && size(S) != (n,m) && error("A and S have mismatched sizes $(size(A)) vs $(size(S))")

    !ishermitian(Q) && error("Q must be Hermitian")
    !ishermitian(R) && error("R must be Hermitian")

    # return all inputs converted to matrices
    return A, B, Q, R, S
end
#_check_ARE_inputs(A, B, Q, R) = _check_ARE_inputs(A, B, Q, R, nothing)

# Solve the dicrete-time Riccati (with S=0) by setting
# up the basic matrix pencil
function _ared(A::Matrix{T}, B::Matrix{T}, Q::Matrix{T}, R::Matrix{T}; stable_solution::Bool=true) where {T <: Number}
    n, m = size(B)

    L = [Matrix{T}(I, n, n) B*(R\B');
         zeros(n,n) A']
    M = [A zeros(n,n);
         -Q Matrix{T}(I, n, n)]

    return _sovle_ARE_pencil(L, M, Val(:c), stable_solution)
end

"""
    arec_noinv(A, G, Q)

Find the solution `X` to the Riccati equation
A'X + XA - XGX + Q = 0
"""
function arec_noinv(A::Matrix{T}, G::Matrix{T}, Q::Matrix{T}) where {T <: Number}
    L = [A  -G;
        -Q  -A']

    return _sovle_ARE_pencil(L, I, Val(:c))
end


# This is a the same balancing used in scipy.linalg,
# but modified (hopefully correctly) for better readability
function _balance_extended_pencil!(H::Matrix, J::Matrix, n::Int, m::Int)
    W = abs.(H) + abs.(J)
    for k=1:size(H,1); W[k,k] = 0; end # Put diagonal entries to zero

    _, _, scaling0 = LAPACK.gebal!('S', W) # Only need the scaling vector

    # Make  scaling maintains symplectic structure and are powers of 2
    scaling0_log2 = [frexp.(el)[2] for el in scaling0]
    s = (scaling0_log2[1:n] - scaling0_log2[n+1:2n]) .>> 1 # Div 2
    scaling_pow2rounded = [exp2.(s); exp2.(-s); ones(eltype(s), m)]

    # Apply the scaling
    De = Diagonal(scaling_pow2rounded)
    _similarity_transform!(H, De)
    _similarity_transform!(J, De)

    D = Diagonal(scaling_pow2rounded[1:n]) # Transformation for recovering the original solution
    return D
end

# This method works for poorly conditioned R matrices
function _are_extended_pencil(timetype::Union{Val{:c},Val{:d}}, E, A::Matrix{T}, B::Matrix{T}, Q::Matrix{T}, R::Matrix{T}, S=nothing; balance=true) where {T <: Number}
    n, m = size(B)

    (isnothing(E) || E == I) && (E = I(n))
    isnothing(S) && (S = zeros(n, m))

    if timetype === Val(:c)
        J = zeros(2n+m, 2n+m)
        J[1:n, 1:n] .= E
        J[n+1:2n, n+1:2n] .= E

        H = [A zeros(n,n) B
            -Q  -A' -S
            S'   B' R]
    else
        H = [A zeros(n,n)
        -Q E'
        S' zeros(m,n)]

        J = [E zeros(n,n)
        zeros(n,n) A'
        zeros(m,n) -B']
    end

    if balance
        D = _balance_extended_pencil!(H, J, n, m)
    end

    # Compute P_compress that keeps everything except
    # the deflating subspace correpsonding to eigenvalues at infinity
    F_compress = qr(Matrix(H[:,2n+1:end]))
    P_compress = [zeros(2n, m) Matrix(I, 2n, 2n)] * F_compress.Q'

    L = P_compress * H[:, 1:2n]
    M = P_compress * J[:, 1:2n]

    X, cl_eigvals = _sovle_ARE_pencil(L, M, timetype, stable_solution=true)

    if balance # X -> D*X*D
        ldiv!(D, X)
        rdiv!(X, D)
    end

    return X, cl_eigvals
end


# Find the c/d - stabilizing/antistabilizing subspace to a Riccati matrix pencil
# `L - λM`
function _sovle_ARE_pencil(L, M, timetype::Union{Val{:c},Val{:d}}; stable_solution::Bool=true)
    n = Int(size(L,1) / 2)

    schurfact = (M == I) ? schur(L) : schur(L, M)

    # Find the basis vectors corresponding to the chosen subspace (c-stable/antistable, d-stable/antistable)
    # and reorder the Schur/QZ factorization so that these come first
    if timetype === Val(:c) && stable_solution
        select = real(schurfact.values) .< 0 # Cirrect? Best for readability to loo at the eigenvalues?
        #select = schurfact.β .> 0
    elseif timetype === Val(:d) && stable_solution
        select = abs2.(schurfact.values) .< 1
        #select = abs2.(schurfact.α) .< abs2.(schurfact.β)
    else
        error("Unknown/unsupported subspace to select when solving ARE pencil")
    end

    count(select) != n && error("Unequal numbers of stable and anti-stable eigenvalues of ARE pencil")
    ordschur!(schurfact, select)

    X = schurfact.Z[n+1:end, 1:n] / schurfact.Z[1:n, 1:n]

    return X, schurfact.values[1:n]
end
