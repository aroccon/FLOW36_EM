
! !***********************************************************************
! !************ Matrix Diagonalization/Inversion Subroutine **************
! !***********************************************************************
subroutine Diagonalization(Eigenv,Pmat,Matrice,NN)
    !Perform Matrix Diagonalization: Matrice=Pmat*Lambda*PmatInv
    !Inputs: Matrice, NN(dimension of Matrice)
    !Outputs are: Eigenv,Pmat,PmatInv
    !Eigenv is the diagonal of Lambda
    implicit none 
    
    real(kind=8)::Eigenv(1:NN)
    real(kind=8)::Pmat(1:NN,1:NN)
    real(kind=8)::Matrice(1:NN,1:NN)
    real(kind=8),allocatable,dimension(:,:)::AA
    integer(kind=8)::NN
    integer::LDA, LDVL, LDVR
    INTEGER::LWMAX
    !PARAMETER        ( LWMAX = 4*NN )	
    !....Local Scalars ..	
    INTEGER          INFO, LWORK
    !....Local Arrays ..
    DOUBLE PRECISION:: VL(NN,NN), VR(NN,NN), WR(NN), WI(NN)
    double precision,allocatable,dimension(:)::WORK
    !....External Subroutines ..
    EXTERNAL         DGEEV
    !....Intrinsic Functions ..
    INTRINSIC        INT, MIN
    !....Query the optimal workspace.
         double precision :: temp
         double precision :: temp_v(NN)
         integer :: j
  
#define flow_stream flow_dir 
        
    LWMAX = int(4*NN) 
        
    allocate(AA(1:NN,1:NN))
    allocate(WORK(LWMAX))
    AA=Matrice
        
    LDA=int(NN)
    LDVL=int(NN)
    LDVR=int(NN)
        
    LWORK = -1
    CALL DGEEV( 'Vectors', 'Vectors', NN, AA, LDA, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )
    LWORK = MIN( LWMAX, INT(WORK( 1 )))
        
    !....Solve eigenproblem.
    
    CALL DGEEV( 'Vectors', 'Vectors', NN, AA, LDA, WR, WI, VL, LDVL,VR, LDVR, WORK, LWORK, INFO )
    
    !....Check for convergence.
    
    IF( INFO.GT.0 ) THEN
    WRITE(*,*)'The algorithm failed to compute eigenvalues.'
    STOP
    END IF
    
    !....Right-Eigenvector matrix
    Pmat=VR
    
    !!....Left-Eigenvector matrix
    !call InvertiMatrice(PmatInv,Pmat,NN)
    
    !....Eigenvalues
    Eigenv=WR
  
#if flow_stream == 0 
  
    temp = Eigenv(2)
    Eigenv(2) = Eigenv(3)
    Eigenv(3) = temp
  
    do j = 1, NN
        temp_v(j) = Pmat(j, 2)   
        Pmat(j, 2) = Pmat(j, 3)      
        Pmat(j, 3) = temp_v(j)    
    enddo
  
#elif flow_stream == 1
    temp = Eigenv(2)
    Eigenv(2) = Eigenv(3)
    Eigenv(3) = temp
  
    do J = 1, 3
        temp_v = Pmat(J, 2)
        Pmat(J, 2) = Pmat(J, 3)
        Pmat(J, 3) = temp_v
    enddo
#endif
  
    deallocate(AA)
    deallocate(work)
    
    return
    end subroutine
  
  !!***********************************************************************
  !!***********************************************************************
  
  ! module Eigen3x3
  !     implicit none
  !     double precision, parameter :: EPS = 2.2204460492503131D-16
  !     double precision, parameter :: SQRT3 = 1.73205080756887729352744634151d0
  !   contains
  
  !   subroutine Diagonalization(A, Q, W)   !(CTOT,Eigenvectors,Eigenvalues)
  !     double precision, intent(in) :: A(3,3)
  !     double precision, intent(out) :: Q(3,3), W(3)
  !     double precision :: T, U, ERROR, NORM
  !     integer :: J,k
  !     double precision :: E(3)  
  !     double precision :: temp
  !     double precision :: temp_v(3)
  
  ! #define flow_stream flow_dir
  
  !     CALL Cardano_Eig(A, W)
  
  !     T = max(abs(W(1)), abs(W(2)), abs(W(3)))
  !     U = max(T, T**2.d0)
  !     ERROR = 256.0d0 * EPS * U**2.d0
  
  !     Q(1, 2) = A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2)
  !     Q(2, 2) = A(1, 3) * A(1, 2) - A(2, 3) * A(1, 1)
  !     Q(3, 2) = A(1, 2)**2.d0
  
  !     Q(1, 1) = Q(1, 2) + A(1, 3) * W(1)
  !     Q(2, 1) = Q(2, 2) + A(2, 3) * W(1)
  !     Q(3, 1) = (A(1,1) - W(1)) * (A(2,2) - W(1)) - Q(3,2)
  !     NORM = Q(1, 1)**2.d0 + Q(2, 1)**2.d0 + Q(3, 1)**2.d0
  
  !     if(NORM<=ERROR)then
  !         call QL_Eig(A, Q, W, E)  ! Pass E as an argument
  !         return
  !     else
  !         NORM = SQRT(1.0d0 / NORM)
  !         do J = 1, 3
  !             Q(J, 1) = Q(J, 1) * NORM
  !         enddo
  !     endif
  
  !     Q(1, 2) = Q(1, 2) + A(1, 3) * W(2)
  !     Q(2, 2) = Q(2, 2) + A(2, 3) * W(2)
  !     Q(3, 2) = (A(1,1) - W(2)) * (A(2,2) - W(2)) - Q(3,2)
  !     NORM = Q(1, 2)**2.d0 + Q(2, 2)**2.d0 + Q(3, 2)**2.d0
  
  !     if(NORM<=ERROR)then
  !         call QL_Eig(A, Q, W, E)  
  !         return
  !     else
  !         NORM = SQRT(1.0d0 / NORM)
  !         do J = 1, 3
  !             Q(J, 2) = Q(J, 2) * NORM
  !         enddo
  !     endif
  
  !     Q(1, 3) = Q(2, 1) * Q(3, 2) - Q(3, 1) * Q(2, 2)
  !     Q(2, 3) = Q(3, 1) * Q(1, 2) - Q(1, 1) * Q(3, 2)
  !     Q(3, 3) = Q(1, 1) * Q(2, 2) - Q(2, 1) * Q(1, 2)
  
  !     ! After computing the eigenvectors in QL_Eig
  !     do J = 1, 3
  !       ! Ensure the first non-zero component is positive for each eigenvector
  !       if (Q(1, J) < 0.0d0) then
  !         Q(1, J) = -Q(1, J)
  !         Q(2, J) = -Q(2, J)
  !         Q(3, J) = -Q(3, J)
  !       endif
  !     enddo
  
  
  !   ! sorting based on the flow
  !     !if X direction flow
  ! !#if flow_stream == 0 
  ! !   temp = W(2)
  ! !   W(2) = W(3)
  ! !   W(3) = temp
  ! !
  ! !   do J = 1, 3
  ! !       temp = Q(J, 2)
  ! !       Q(J, 2) = Q(J, 3)
  ! !       Q(J, 3) = temp
  ! !   enddo
  ! !
  ! !   !if Y direction flow
  ! !#elif flow_stream == 1
  ! !   temp = W(1)
  ! !   W(1) = W(3)
  ! !   W(3) = W(2)
  ! !   W(2) = temp
  ! !
  ! !   do J = 1, 3
  ! !       temp_v(J) = Q(J, 1)   
  ! !       Q(J, 1) = Q(J, 3)      
  ! !       Q(J, 3) = Q(J, 2)     
  ! !       Q(J, 2) = temp_v(J)    
  ! !   enddo
  ! !#endif
      
  !   end subroutine Diagonalization
  
    
  !   subroutine Cardano_Eig(A, W)
  !     double precision, intent(in) :: A(3,3)
  !     double precision, intent(out) ::W(3)
  !     double precision :: M, C1, C0, DE, DD, EE, FF, P, SQRTP, Q, C, S, PHI
  
  !     DE = A(1,2) * A(2,3)
  !     DD = A(1,2)**2.d0
  !     EE = A(2,3)**2.d0
  !     FF = A(1,3)**2.d0
  !     M = A(1,1) + A(2,2) + A(3,3)
  !     C1 = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) - (DD + EE + FF)
  !     C0 = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) - 2.0d0 * A(1,3)*DE
  !     P = M**2 - 3.0d0 * C1
  !     Q = M*(P - (3.0d0/2.0d0)*C1) - (27.0d0/2.0d0)*C0
  !     SQRTP = sqrt(abs(P))
  !     PHI = atan2(sqrt(abs(27.0d0 * ( 0.25d0 * C1**2.d0 * (P - C1) + C0 * (Q + (27.0d0/4.0d0)*C0) ))), Q) / 3.0d0
  
  !     C = SQRTP * cos(PHI)
  !     S = (1.0d0 / SQRT3) * SQRTP * sin(PHI)
  
  !     W(2) = (1.0d0/3.0d0) * (M - C)
  !     W(3) = W(2) + S
  !     W(1) = W(2) + C
  !     W(2) = W(2) - S
        
  !   end subroutine Cardano_Eig
  
  !     subroutine QL_Eig(A, Q, W, E)
  !       double precision, intent(in) :: A(3,3)
  !       double precision, intent(out) :: Q(3,3), W(3)
  !       double precision :: E(3), G, R, P, F, B, S, C, T
  !       integer :: L, M, I, J, K, NITER
    
  !       call Householder(A, Q, W, E) 
    
  !       do L = 1, 2
  !         NITER = 0
  !         do I = 1, 50
  !           do M = L, 2
  !             G = abs(W(M)) + abs(W(M+1))
  !             if(abs(E(M)) + G == G) exit
  !           enddo
  !           if(M == L) exit
  !           NITER = NITER + 1
  !           if(NITER >= 30) then
  !             print *, 'QL_Eig: No convergence.'
  !             return
  !           endif
    
  !           G = (W(L+1) - W(L)) / (2.0d0 * E(L))
  !           R = SQRT(1.0d0 + G**2.d0)
  !           G = W(M) - W(L) + E(L)/(G + R * MERGE(1.0d0, -1.0d0, G >= 0.0d0))
  !           S = 1.0d0
  !           C = 1.0d0
  !           P = 0.0d0
    
  !           do J = M - 1, L, -1
  !             F = S * E(J)
  !             B = C * E(J)
  !             if (ABS(F) > ABS(G)) then
  !               C = G / F
  !               R = SQRT(1.0d0 + C**2.d0)
  !               E(J+1) = F * R
  !               S = 1.0d0 / R
  !               C = C * S
  !             else
  !               S = F / G
  !               R = SQRT(1.0d0 + S**2.d0)
  !               E(J+1) = G * R
  !               C = 1.0d0 / R
  !               S = S * C
  !             endif
    
  !             G = W(J+1) - P
  !             R = (W(J) - G) * S + 2.0d0 * C * B
  !             P = S * R
  !             W(J+1) = G + P
  !             G = C * R - B
    
  !             do K = 1, 3
  !               T = Q(K, J+1)
  !               Q(K, J+1) = S * Q(K, J) + C * T
  !               Q(K, J) = C * Q(K, J) - S * T
  !             enddo
  !           enddo
  !           W(L) = W(L) - P
  !           E(L) = G
  !           E(M) = 0.0d0
  !         enddo
  !       enddo
  !     end subroutine QL_Eig
    
  !     subroutine Householder(A, Q, D, E)
  !       double precision, intent(in) :: A(3,3)
  !       double precision, intent(out) :: Q(3,3), D(3), E(2)
  !       double precision :: U(3), P(3), OMEGA, F, G, K, H
  !       integer :: I, J
    
  !       ! Initialize Q as the identity matrix
  !       Q = 0.0d0
  !       do I = 1, 3
  !           Q(I, I) = 1.0d0
  !       end do
    
  !       ! Step 1: Create the first Householder reflection for the (1,2) and (1,3) elements
  !       H = A(1,2)**2.d0 + A(1,3)**2.d0
  !       if (H > 0.0d0) then
  !           G = -sign(sqrt(H), A(1,2))  ! Set G to the negative norm of the vector
  !           E(1) = G                    ! E(1) will store the first off-diagonal element
  !           F = G * A(1,2)
  !           U(2) = A(1,2) - G           ! Construct Householder vector U
  !           U(3) = A(1,3)
            
  !           OMEGA = H - F
  !           if (OMEGA > 0.0d0) then
  !               OMEGA = 1.0d0 / OMEGA   ! Normalize the Householder vector
    
  !               ! Calculate P = A * U / (U^T * A * U)
  !               do I = 2, 3
  !                   F = A(2,I)*U(2) + A(I,3)*U(3)
  !                   P(I) = OMEGA * F
  !               end do
    
  !               ! Calculate K for adjusting P
  !               K = 0.0d0
  !               do I = 2, 3
  !                   K = K + U(I) * P(I)
  !               end do
  !               K = 0.5d0 * K * OMEGA
    
  !               do I = 2, 3
  !                   P(I) = P(I) - K * U(I)
  !               end do
    
  !               ! Update A to be tridiagonal in the first column and row
  !               D(1) = A(1,1)
  !               D(2) = A(2,2) - 2.0d0 * P(2) * U(2)
  !               D(3) = A(3,3) - 2.0d0 * P(3) * U(3)
                
  !               ! Set the off-diagonal element in E(2)
  !               E(2) = A(2,3) - P(2) * U(3) - U(2) * P(3)
                
  !               ! Apply the transformation to update Q
  !               do I = 1, 3
  !                   F = OMEGA * (Q(I,2) * U(2) + Q(I,3) * U(3))
  !                   Q(I,2) = Q(I,2) - F * U(2)
  !                   Q(I,3) = Q(I,3) - F * U(3)
  !               end do
  !           else
  !               D(1) = A(1,1)
  !               D(2) = A(2,2)
  !               D(3) = A(3,3)
  !               E(2) = 0.0d0
  !           end if
  !       else
  !           D(1) = A(1,1)
  !           D(2) = A(2,2)
  !           D(3) = A(3,3)
  !           E(1) = 0.0d0
  !           E(2) = 0.0d0
  !       end if
  !   end subroutine Householder
    
    
  ! end module Eigen3x3
    
  
  