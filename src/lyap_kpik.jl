
#using IteratveSolvers

"""
`kpik(A,B,E=1;<keyword arguments>)`

Julia code for K-PIK (Krylov-plus-inverted-Krylov)
Based on kpik.m avalible from V. Simoncini's website
(http://www.dm.unibo.it/~simoncin/software.html)
and includes many of the same comments and much of the same description.
Approximately solve

AX + XA' + BB' = 0

by means of the extended Krylov subspace method
ARGUMENTS
`A` : coeff matrix, A < 0
`B` : factor of rhs,   nxk matrix with k << n

keyword arguments
`m` : max dimension of the Krylov space
`tol` : stopping tolerance based on the backwards error,
with stopping criterion
 ‖AX + XA' - BB'‖
-----------------------     < tol
‖BB'‖ + ‖A‖ ‖X‖
computed in a cheap manner. Defult: 1e-9
Output
`Z` : solution factor s.t.  X = ZZ'
`er2` : history of scaled residual, as above
Comments
1. The projected solution is computed at each iteration
As an alternative, a periodic computation could be considered.
2. This code performs a factorization of A. As an alternative,
iterative solves could be considered.

    **V. Simoncini** A new iterative method for solving large-scale Lyapunov matrix equations,
    SIAM J.  Scient. Computing, v.29, n.3 (2007), pp. 1268-1288.
"""

# Perhaps use r or q instead of s


eye(m,n) = I(max(m,n))[1:m,1:n]
#function kpik(A, B, E=I; LE=I,m=100,tol=1e-9,tolY=1e-12,infoV=false)
#function lyapc(A::AbstractMatrix, B, ::Val{:kpik}; m=100,tol=1e-9,tolY=1e-12,infoV=false)
function lyapc(A::AbstractMatrix, Q::HermitianLowRank, ::Val{:kpik}; m=100,tol=1e-9,tolY=1e-12,infoV=false)
    B = Q.factor

    normB = norm(B)^2
    normA = norm(A)
    er2 = zeros(m,1)

    n, sh = size(B)

    Y = []
    odds = []
    er2 = [] # why the 2? Error squared, square root of error, or 2-norm?


    if norm(A-A',1)<1e-14
        if issparse(A)
            cfA = cholfact(-A,perm=1:size(A,1))
            LA = sparse(cfA[:L])
            UA = -LA'
        else
            UA = chol(-A)
            LA = -UA'
        end
        infoV && println("A sym. Completed Chol factorization\n")
        k_max = 2
    else
        #A_lu = lu(A,Val(false))
        A_lu = lu(A)
        k_max = m
    end

    s = 2*sh

    # Orthogonalize [B,A^{-1}B] with an economy-size QR
    # Julia's qr decomposition is always "economy size" (no true any longer)
    U, beta = qr([ B A_lu\B ])
    U = U[1:n,1:s]

    ibeta = inv(beta[1:s,1:s])
    beta  =  beta[1:sh,1:sh]
    beta2 = beta*beta'

    Up = zeros(n,s)
    H = zeros((m+1)*s,m*s)
    T = zeros((m+1)*s,m*s)
    L = zeros((m+1)*s,m*s)
    local js, j, rho


    for j=1:m
        jms=(j-1)*s+1
        j1s=(j+1)*s
        js=j*s
        js1=js+1
        jsh=(j-1)*s+sh

        # Expand the basis
        # multiply by A

        # Fill the matrix Up with the new columns from [A^j*B  inv(A)^j*B]
        # (should be doen in-place)
        Up[:,1:sh] = A*U[:,jms:jsh]
        Up[:,sh+1:s] = A_lu\U[:,jsh+1:js] # Doesn't U have n rows?

        # Orthogonalize the new basis block wrt all the previous ones by modified gram
        for l=1:2
            k_min = max(1,j-k_max)
            for kk = k_min:j
                k1 = (kk-1)*s+1
                k2 = kk*s
                coef =  U[1:n,k1:k2]'*Up
                H[k1:k2,jms:js] = H[k1:k2,jms:js]+ coef
                Up = Up - U[:,k1:k2]*coef
            end
        end

        if j <= m
            Up_qr = qr(Up)

            H[js1:j1s,jms:js] = Up_qr.R
            Hinv = inv(H[js1:j1s,jms:js])

            Up = Matrix(Up_qr.Q)
        end


        ###############################################################
        # Recover the columns of T=U'*A*U (projection of A onto the space) from
        # the colums of H.
        # REMARK: we need T as coefficient matrix of the projected problem.
        Iden=I(js+s)

        if j == 1
            L[1:j*s+sh,(j-1)*sh+1:j*sh] = [H[1:s+sh,1:sh]/ibeta[1:sh,1:sh] eye(s+sh,sh)/ibeta[1:sh,1:sh]]*ibeta[1:s,sh+1:s]
        else
            L[1:j*s+s,(j-1)*sh+1:j*sh] = L[1:j*s+s,(j-1)*sh+1:j*sh] + H[1:j*s+s,jms:jms-1+sh]*rho
        end

        odds = [odds; jms:(jms-1+sh)]   # store the odd block columns
        evens = 1:js
        flag = trues(size(evens))
        flag[odds] .= false
        evens = evens[flag]
        T[1:js+s,odds]=H[1:js+s,odds]   #odd columns

        T[1:js+sh,evens]=L[1:js+sh,1:j*sh]   #even columns
        L[1:j*s+s,j*sh+1:(j+1)*sh] = ( Iden[1:j*s+s,(js-sh+1):js]- T[1:js+s,1:js]*H[1:js,js-sh+1:js])*Hinv[sh+1:s,sh+1:s]
        rho = Hinv[1:sh,1:sh] \ Hinv[1:sh,sh+1:s]

        #################################################################

        # Solve the projected problem by Bartels-Stewart
        Y = lyap(T[1:js,1:js], eye(j*s,sh)*beta2*eye(j*s,sh)')

        # safeguard to preserve symmetry
        Y = (Y+Y')/2

        # Compute the residual norm. See the article by Valeria

        cc = [H[js1:j1s,js-s+1:js-sh] L[js1:j1s,(j-1)*sh+1:j*sh]]

        normY = norm(Y)

        er2=[er2;sqrt(2)*norm(cc*Y[js-s+1:js,:])/(normB+normA*normY)]

        infoV && println("KPIK It: $j -- Current Backwards Error: $(er2[j])")

        if (er2[j]<tol)
            break
        end

        U = [U Up]
    end

    # Done
    # reduce solution rank if needed


    # We do not reduce the solution for trouble-shooting reasons..
    #return U*Y*U'
    return HermitianLowRank(U*cholesky(Y).L)

    sY,uY=eigen(Y)

    id=sortperm(sY) # What sorting is this supposed to do?
    sY=sort(sY)
    sY=flipdim(sY,1)


    uY .= uY[:,id[end:-1:1]]
    is = 0
    for ii in 1:size(sY)[1]
        if abs(sY[ii])>tolY
            is = is+1
        end
    end

    Y0 = uY[:,1:is]*diagm(sqrt(sY[1:is])) # This is a Cholesky factor of Y, should have another name

    Z = U[1:n,1:js]*Y0 # X = Z*Z'
    er2=er2[1:j]

    #if infoV
    #    println("its  Back. Error            space dim. CPU Time")
    #    println("$j    $(er2[j])  $js          $(toq())")
    #end

    return HermitianLowRank(Z), er2

end
