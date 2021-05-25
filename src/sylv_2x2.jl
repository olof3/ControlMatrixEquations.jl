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
    @inbounds if M == 2 && N == 2
        Xv = _sylvd2(SMatrix{2,2}(A), SMatrix{2,2}(B), SMatrix{2,2}(C))
        C[:] .= Xv
    elseif M == 2 && N == 1
        Xv = _sylvd2(SMatrix{2,2}(A), SMatrix{1,1}(B), SMatrix{2,1}(C))
        C[:] .= Xv
    elseif M == 1 && N == 2
        Xv = _sylvd2(SMatrix{1,1}(A), SMatrix{2,2}(B), SMatrix{1,2}(C))
        C[:] .= Xv
    elseif M == 1 && N == 1
        Xv = _sylvd2(SMatrix{1,1}(A), SMatrix{1,1}(B), SMatrix{1,1}(C))
        C[:] .= Xv
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
end
function _sylvd2(A::SMatrix, B::SMatrix, C::SMatrix)

    M, N = size(C)
    Xv = (kron(transpose(B), A) - SMatrix{M*N,M*N}(I)) \ C[:] # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvd or ?lyapd"); end

    return Xv
end



"""
    _sylvd!(A, B, C, E, F) -> X

Find the solution `X` to the discrete-time Sylvester equation

`AXB + EXF = C`

for small matrices (1x1, 1x2, 2x1, 2x2), overwriting the input `C`.
"""
function _sylvg!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, E::AbstractMatrix, F::AbstractMatrix)
    M, N = size(C)
    if M == 2 && N == 2
        _sylvg!(A, B, C, E, F, Val(2), Val(2))
    elseif M == 2 && N == 1
        _sylvg!(A, B, C, E, F, Val(2), Val(1))
    elseif M == 1 && N == 2
        _sylvg!(A, B, C, E, F, Val(1), Val(2))
    elseif M == 1 && N == 1
        _sylvg!(A, B, C, E, F, Val(1), Val(1))
    else
        error("Matrix dimensionsins should not be greater than 2")
    end
    return C
end
function _sylvg!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, E::AbstractMatrix, F::AbstractMatrix, ::Val{M}, ::Val{N}) where {T <: Number, M, N}
    As = SMatrix{M,M}(A)
    Bs = SMatrix{N,N}(B)
    Es = SMatrix{M,M}(E)
    Fs = SMatrix{N,N}(F)
    Cvs = SMatrix{M,N}(C)[:] # vectorization of C

    Xv = (kron(transpose(Bs), As) + kron(transpose(Fs), Es)) \ Cvs # using the vectorization identity vec(AXB) = kron(B'*A)*vec(X)

    if any(!isfinite, Xv); error("Matrix equation has no solution, see ?sylvg / ?lyapcg / ?lyapdg"); end

    C .= reshape(Xv, M, N)
end