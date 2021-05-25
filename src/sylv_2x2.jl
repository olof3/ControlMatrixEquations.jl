# These funcitons are used by the Bertels-Stewart algoirthm for complex arithmetic,
# they might also have some value in their own right
@inline function sylvd(a::Number, b::Number, c::Number)
    x = c / (a * b - 1)

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    x
end
@inline function sylvc(a::Number, b::Number, c::Number)
    x = c / (a + b)

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    x
end
@inline function sylvg(a::Number, b::Number, c::Number, e::Number, f::Number)
    x = c / (a*b + e*f)

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvg / ?lyapcg / ?lyapdg"); end

    x
end


"""
    _sylvc!(A, B, C) -> X

Find the solution `X` to the continuous-time Sylvester equation

`AX + XB = C`

for small matrices (1x1, 1x2, 2x1, 2x2),
overwriting the input `C` with the solution `X`.
"""
@inline function _sylvc!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    M, N = size(C)
    if M == 2 && N == 2
        C[:] .= _sylvc(SMatrix{2,2}(A), SMatrix{2,2}(B), SMatrix{2,2}(C))
    elseif M == 2 && N == 1
        C[:] .= _sylvc(SMatrix{2,2}(A), SMatrix{1,1}(B), SMatrix{2,1}(C))
    elseif M == 1 && N == 2
        C[:] .= _sylvc(SMatrix{1,1}(A), SMatrix{2,2}(B), SMatrix{1,2}(C))
    elseif M == 1 && N == 1
        C[:] .= _sylvc(SMatrix{1,1}(A), SMatrix{1,1}(B), SMatrix{1,1}(C))
    else
        error("Problem size must not be greater than (2,2)")
    end
    return C
end
@inline function _sylvc(A::SMatrix, B::SMatrix, C::SMatrix)
    M, N = size(C)

    Xv = lu(kron(SMatrix{N,N}(I), A) + kron(transpose(B), SMatrix{M,M}(I))) \ C[:] # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X) (with A = I or B = I)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    return Xv
end


"""
    _sylvd!(A, B, C) -> X

Find the solution `X` to the discrete-time Sylvester equation

`AXB - X = C`

for small matrices (1x1, 1x2, 2x1, 2x2),
overwriting the input `C` with the solution `X`.
"""
function _sylvd!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)    
    M, N = size(C)
    if M == 2 && N == 2
        C[:] .= _sylvd(SMatrix{2,2}(A), SMatrix{2,2}(B), SMatrix{2,2}(C))
    elseif M == 2 && N == 1
        C[:] .= _sylvd(SMatrix{2,2}(A), SMatrix{1,1}(B), SMatrix{2,1}(C))
    elseif M == 1 && N == 2
        C[:] .= _sylvd(SMatrix{1,1}(A), SMatrix{2,2}(B), SMatrix{1,2}(C))
    elseif M == 1 && N == 1
        C[:] .= _sylvd(SMatrix{1,1}(A), SMatrix{1,1}(B), SMatrix{1,1}(C))
    else
        error("Problem size must not be greater than (2,2)")
    end
    return C
end
function _sylvd(A::SMatrix, B::SMatrix, C::SMatrix)
    M, N = size(C)

    Xv = (kron(transpose(B), A) - SMatrix{M*N,M*N}(I)) \ C[:] # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    return Xv
end



"""
    _sylvd!(A, B, C, E, F) -> X

Find the solution `X` to the discrete-time Sylvester equation

`AXB + EXF = C`

for small matrices (1x1, 1x2, 2x1, 2x2)
overwriting the input `C` with the solution `X`.
"""
function _sylvg!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, E::AbstractMatrix, F::AbstractMatrix)
    M, N = size(C)
    if M == 2 && N == 2
        C[:] .= _sylvg(SMatrix{2,2}(A), SMatrix{2,2}(B), SMatrix{2,2}(C), SMatrix{2,2}(E), SMatrix{2,2}(F))
    elseif M == 2 && N == 1
        C[:] .= _sylvg(SMatrix{2,2}(A), SMatrix{1,1}(B), SMatrix{2,1}(C), SMatrix{2,2}(E), SMatrix{1,1}(F))
    elseif M == 1 && N == 2
        C[:] .= _sylvg(SMatrix{1,1}(A), SMatrix{2,2}(B), SMatrix{1,2}(C), SMatrix{1,1}(E), SMatrix{2,2}(F))
    elseif M == 1 && N == 1
        C[:] .= _sylvg(SMatrix{1,1}(A), SMatrix{1,1}(B), SMatrix{1,1}(C), SMatrix{1,1}(E), SMatrix{1,1}(F))
    else
        error("Problem size must not be greater than (2,2)")
    end
    return C
end
function _sylvg(A::SMatrix, B::SMatrix, C::SMatrix, E::SMatrix, F::SMatrix)

    Xv = (kron(transpose(B), A) + kron(transpose(F), E)) \ C[:] # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvg / ?lyapcg / ?lyapdg"); end

    return Xv
end