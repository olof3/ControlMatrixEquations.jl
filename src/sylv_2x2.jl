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

    if !isfinite(x); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    x
end


"""
    _sylvc!(A, B, C) -> X

Find the solution `X` to the continuous-time Sylvester equation

`AX + XB = C`

for small matrices (1x1, 1x2, 2x1, 2x2), overwriting the input `C`.
"""
@inline function _sylvc!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    M, N = size(C)
    if M == 2 && N == 2
        _sylvc!(A, B, C, Val(2), Val(2))
    elseif M == 2 && N == 1
        _sylvc!(A, B, C, Val(2), Val(1))
    elseif M == 1 && N == 2
        _sylvc!(A, B, C, Val(1), Val(2))
    elseif M == 1 && N == 1
        _sylvc!(A, B, C, Val(1), Val(1))
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
    return C
end
@inline function _sylvc!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, ::Val{M}, ::Val{N}) where {T <: Number, M, N}
    As = SMatrix{M,M}(A)
    Bs = SMatrix{N,N}(B)
    Cvs = SMatrix{M,N}(C)[:] # vectorization of C

    Xv = lu(kron(SMatrix{N,N}(I), As) + kron(transpose(Bs), SMatrix{M,M}(I))) \ Cvs # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X) (with A = I or B = I)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvc or ?lyapc"); end

    C .= reshape(Xv, M, N)
end


"""
    _sylvd!(A, B, C) -> X

Find the solution `X` to the discrete-time Sylvester equation

`AXB - X = C`

for small matrices (1x1, 1x2, 2x1, 2x2), overwriting the input `C`.
"""
function _sylvd!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix)
    M, N = size(C)
    if M == 2 && N == 2
        _sylvd!(A, B, C, Val(2), Val(2))
    elseif M == 2 && N == 1
        _sylvd!(A, B, C, Val(2), Val(1))
    elseif M == 1 && N == 2
        _sylvd!(A, B, C, Val(1), Val(2))
    elseif M == 1 && N == 1
        _sylvd!(A, B, C, Val(1), Val(1))
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
    return C
end
function _sylvd!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, ::Val{M}, ::Val{N}) where {T <: Number, M, N}
    As = SMatrix{M,M}(A)
    Bs = SMatrix{N,N}(B)
    Cvs = SMatrix{M,N}(C)[:] # vectorization of C

    Xv = (kron(transpose(Bs), As) - SMatrix{M*N,M*N}(I)) \ Cvs # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    C .= reshape(Xv, M, N)
end
