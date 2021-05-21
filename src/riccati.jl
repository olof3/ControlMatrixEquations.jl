"""
    arec(A, B, Q, R, S=nothing; force_esp=false)

Find the stabilizing solution `X` to the continuous-time Riccati equation
`A'X + XA - (XB + S)/R(XB + S)' + Q = 0`

Arnold, W. F., & Laub, A. J. (1984). Generalized eigenproblem algorithms and software
for algebraic Riccati equations. Proceedings of the IEEE, 72(12), 1746-1754.

Bini, D. A., Iannazzo, B., & Meini, B. (2011). Numerical solution of algebraic Riccati equations.
"""
function arec(A::AbstractNumOrArray, B::AbstractNumOrArray, Q::Union{AbstractNumOrArray,UniformScaling}, R::Union{AbstractNumOrArray,UniformScaling}, S=nothing; force_extended_pencil=false, stabsol::Bool=true)
    _, A, B, Q, R, S = _check_ARE_inputs(I, A, B, Q, R, S)

    if cond(R) <= 1e4 && !force_extended_pencil # FIXME: What condition to use?
        if S == nothing
            return arec_noinv(A, B/R*B', Q, stabsol=stabsol)
        else
            return arec_noinv(A - B/R*S', B/R*B', Q - Hermitian(S/R*S'), stabsol=stabsol)
        end
    else
        return _ARE_extended_pencil(Val(:c), I, A, B, Q, R, S, stabsol=stabsol)
    end
end

"""
    arecg(E, A, B, Q, R, S=nothing; force_esp=false) -> X

Find the stabilizing solution `X` to the continuous-time generalized Riccati equation
`A'XE + E'XA - (E'XB + S)/R(E'XB + S)' + Q = 0`
"""
function arecg(E, A, B, Q, R, S=nothing; balance=false, stabsol=true)
    E, A, B, Q, R, S = _check_ARE_inputs(E, A, B, Q, R, S)
    _ARE_extended_pencil(Val(:c), E, A, B, Q, R, S, stabsol=stabsol, balance=balance)
end

"""
    arec_noinv(A, G, Q)

Find the solution `X` to the Riccati equation
`A'X + XA - XGX + Q = 0`
"""
function arec_noinv(A::Matrix{T}, G::Matrix{T}, Q::Matrix{T}; stabsol=true) where {T <: Number}
    M = [A  -G;
        -Q  -A']

    return _sovle_ARE_pencil(M, I, Val(:c), stabsol=stabsol)
end


"""
    ared(A, B, Q, R, S=nothing; force_esp=false) -> X

Find the stabilizing solution `X` to the discrete-time Riccati equation
`A'XA - X - (A'XB + S)/(B'XB + R)(A'XB + S)' + Q = 0`

Arnold, W. F., & Laub, A. J. (1984). Generalized eigenproblem algorithms and software
for algebraic Riccati equations. Proceedings of the IEEE, 72(12), 1746-1754.

Bini, D. A., Iannazzo, B., & Meini, B. (2011). Numerical solution of algebraic Riccati equations.
"""
function ared(A::AbstractNumOrArray, B::AbstractNumOrArray, Q::Union{AbstractNumOrArray,UniformScaling}, R::Union{AbstractNumOrArray,UniformScaling}, S=nothing; force_extended_pencil=false, stabsol::Bool=true)
    _, A, B, Q, R, S = _check_ARE_inputs(I, A, B, Q, R, S)

    if cond(R) <= 1e4 && !force_extended_pencil
        if S == nothing
            return _ared(A, B, Q, R, stabsol=stabsol)
        else
            return _ared(A - B/R*S', B, Q - S/R*S', R, stabsol=stabsol)
        end
    else
        return _ARE_extended_pencil(Val(:d), I, A, B, Q, R, S, stabsol=stabsol)
    end
end
function aredg(E, A, B, Q, R, S=nothing; balance=false, stabsol=true)
    E, A, B, Q, R, S = _check_ARE_inputs(E, A, B, Q, R, S)
    _ARE_extended_pencil(Val(:d), E, A, B, Q, R, S, stabsol=stabsol, balance=balance)
end



# Solve the dicrete-time Riccati (with S=0) by setting
# up the standard matrix pencil
function _ared(A::Matrix{T}, B::Matrix{T}, Q::Matrix{T}, R::Matrix{T}; stabsol::Bool=true) where {T <: Number}
    n, m = size(B)

    M = [Matrix{T}(I, n, n) B*(R\B');
         zeros(n,n) A']
    L = [A zeros(n,n);
         -Q Matrix{T}(I, n, n)]

    return _sovle_ARE_pencil(L, M, Val(:d), stabsol=stabsol)
end


"""
    _check_ARE_inputs(E, A, B, Q, R, S)

Check that the arguments are valid input to a Riccati equation and convert all of them
to `Matrix{T}`. `S` may equal `nothing`

"""
function _check_ARE_inputs(E, A, B, Q, R, S)
    n, m = size(A,1), size(B,2)
    T = float(promote_type(eltype(E), eltype(A), eltype(B), eltype(Q), isnothing(S) ? Union{} : eltype(S)))

    # Convert all inputs to Matrix{T}
    A, B = to_matrix(T, A), to_matrix(T, B)
    Q = (Q isa UniformScaling) ? Matrix{T}(Q,n,n) : to_matrix(T, Q)
    R = (R isa UniformScaling) ? Matrix{T}(R,m,m) : to_matrix(T, R) # There would be small improvements for cases where this conversion could be avoided
    !isnothing(S) && (S = to_matrix(T, S))
    if !(E isa UniformScaling); E = to_matrix(T, E); end

    # Check matrix sizes
    size(A) != (n, n) && error("A and B have mismatched sizes $(size(A)) vs $(size(B))")
    !(E isa UniformScaling) && size(E) != (n, n) && error("A and E have mismatched sizes $(size(A)) vs $(size(E))")
    size(Q) != (n, n) && error("A and Q have mismatched sizes $(size(A)) vs $(size(Q))")
    size(R) != (m, m) && error("B and R have mismatched sizes $(size(B)) vs $(size(R))")
    !isnothing(S) && size(S) != (n,m) && error("B and S have mismatched sizes $(size(B)) vs $(size(S))")

    # Check structure of R and Q
    !ishermitian(Q) && error("Q must be Hermitian")
    !ishermitian(R) && error("R must be Hermitian")

    # return all inputs converted to matrices
    return E, A, B, Q, R, S
end
#_check_ARE_inputs(A, B, Q, R) = _check_ARE_inputs(A, B, Q, R, nothing)


# This is the same balancing used in scipy.linalg,
# but modified (hopefully correctly) for better readability
function _balance_extended_pencil!(H::Matrix{T}, J::Matrix{T}, n::Int, m::Int) where {T <: Number}
    W = abs.(H) + abs.(J)
    for k=1:size(H,1); W[k,k] = 0; end # Put diagonal entries to zero

    _, _, scaling0 = LAPACK.gebal!('S', W) # Only need the vector of scaling factors

    # Make sure the scaling factors maintain the symplectic structure and are powers of 2
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



# The extended pencil method can handle poorly conditioned R matrices
function _ARE_extended_pencil(timetype::Union{Val{:c},Val{:d}}, E, A::Matrix{T}, B::Matrix{T}, Q::Matrix{T}, R::Matrix{T}, S=nothing; balance=true, stabsol=stabsol) where {T <: Number}
    n, m = size(B)

    (isnothing(E) || E == I) && (E = Matrix{T}(I, n, n))
    isnothing(S) && (S = zeros(n, m))

    if timetype === Val(:c)
        H = [A zeros(T, (n,n)) B
            -Q  -A' -S
            S'   B' R]

        J = zeros(T, (2n+m, 2n+m))
        J[1:n, 1:n] .= E
        J[n+1:2n, n+1:2n] .= E'
    else
        H = [A zeros(n,n) B
        -Q E' -S
        S' zeros(m,n) R]

        J = [E zeros(n,n) zeros(n,m)
        zeros(n,n) A' zeros(n,m)
        zeros(m,n) -B' zeros(m,m)]
    end


    if balance
        D = _balance_extended_pencil!(H, J, n, m)
    end


    # Compute P_compress that keeps everything except
    # the deflating subspace correpsonding to eigenvalues at infinity
    F_compress = qr(Matrix(H[:,2n+1:end]))
    P_compress = [zeros(2n, m) Matrix(I, 2n, 2n)] * F_compress.Q'

    M = P_compress * H[:, 1:2n]
    L = P_compress * J[:, 1:2n]

    X, cl_eigvals = _sovle_ARE_pencil(M, L, timetype, E, stabsol=true)

    if balance # X -> D*X*D
        ldiv!(D, X)
        #rdiv!(X, D) # Method missing in LinearAlgebra
        rmul!(X, inv(D))
    end

    return X, cl_eigvals
end


"""
    _sovle_ARE_pencil(M, L, timetype::Union{Val{:c},Val{:d}}; stabsol::Bool=true)

    Find the c/d - stabilizing/antistabilizing subspace to a Riccati matrix pencil `L - λM`
"""
function _sovle_ARE_pencil(M, L, timetype::Union{Val{:c},Val{:d}}, E = I; stabsol::Bool=true)
    n = Int(size(M,1) / 2)

    schurfact = (L == I) ? schur(M) : schur(M, L)

    # Find the basis vectors corresponding to the chosen subspace (c-stable/antistable, d-stable/antistable)
    # so that the Schur/QZ factorization can be reordered to have these first
    if timetype === Val(:c) && stabsol
        select = real(schurfact.values) .< 0 # Correct? Best for readability to look at the eigenvalues rather than α and β?
    elseif timetype === Val(:c) && !stabsol
        select = real(schurfact.values) .> 0
    elseif timetype === Val(:d) && stabsol
        select = abs2.(schurfact.values) .< 1
    elseif timetype === Val(:d) && !stabsol
        select = abs2.(schurfact.values) .> 1
    else
        error("Unknown/unsupported subspace to select when solving ARE pencil")
    end

    count(select) != n && error("Unequal numbers of stable and anti-stable eigenvalues of ARE pencil")
    ordschur!(schurfact, select)

    if E == I
        Z11 = schurfact.Z[1:n, 1:n]
        Z21 = schurfact.Z[n+1:end, 1:n]
    else
        Z0 = [E * schurfact.Z[1:n, 1:n];
             schurfact.Z[n+1:end, 1:n]]
        Z = qr(Z0).Q * [Matrix(I, n, n); zeros(n, n)]
        Z11 = Z[1:n, :]
        Z21 = Z[n+1:2n, :]
    end

    X = Z21 / Z11

    return X, schurfact.values[1:n] # FIXME: Consier symmetrizing X
end
