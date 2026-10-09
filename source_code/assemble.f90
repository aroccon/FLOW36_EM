module assemble

use commondata
use par_size
use polymer
use sim_par
use derivatives
!use Eigen3x3

implicit none

double precision, allocatable, dimension(:,:,:) :: diag11,diag22,diag33
double precision, allocatable, dimension(:,:,:) :: erre11,erre12,erre13
double precision, allocatable, dimension(:,:,:) :: erre21,erre22,erre23
double precision, allocatable, dimension(:,:,:) :: erre31,erre32,erre33
double precision, allocatable, dimension(:,:,:) :: emme11,emme12,emme13
double precision, allocatable, dimension(:,:,:) :: emme21,emme22,emme23
double precision, allocatable, dimension(:,:,:) :: emme31,emme32,emme33
double precision, allocatable, dimension(:,:,:) :: omega1,omega2,omega3
double precision :: CTOT(3,3),Eigenvectors(3,3),Eigenvalues(3)

private 
public diag11,diag22,diag33,assemble_eigen,assemble_matrices
public erre11,erre12,erre13,erre21,erre22,erre23,erre31,erre32,erre33
public omega1,omega2,omega3,emme11,emme22,emme33


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assemble_matrices

  integer :: i,j,k
  integer :: ntime
  double precision :: epsnum,mask
  
  ! matrix M = R gradU R^T  
  
  allocate(emme11(nx,fpz,fpy))
  allocate(emme12(nx,fpz,fpy))
  allocate(emme13(nx,fpz,fpy))
  allocate(emme21(nx,fpz,fpy))
  allocate(emme22(nx,fpz,fpy))
  allocate(emme23(nx,fpz,fpy))
  allocate(emme31(nx,fpz,fpy))
  allocate(emme32(nx,fpz,fpy))
  allocate(emme33(nx,fpz,fpy))
  
  !$acc parallel loop collapse(3)
  do j=1,fpy
   do k=1,fpz
     do i=1,nx
      emme11(i,k,j)=erre11(i,k,j)*(dudx(i,k,j)*erre11(i,k,j)+dudy(i,k,j)*erre21(i,k,j)+dudz(i,k,j)*erre31(i,k,j))+ &
              &     erre21(i,k,j)*(dvdx(i,k,j)*erre11(i,k,j)+dvdy(i,k,j)*erre21(i,k,j)+dvdz(i,k,j)*erre31(i,k,j))+ &
              &     erre31(i,k,j)*(dwdx(i,k,j)*erre11(i,k,j)+dwdy(i,k,j)*erre21(i,k,j)+dwdz(i,k,j)*erre31(i,k,j)) 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      emme12(i,k,j)=erre11(i,k,j)*(dudx(i,k,j)*erre12(i,k,j)+dudy(i,k,j)*erre22(i,k,j)+dudz(i,k,j)*erre32(i,k,j))+ &
              &     erre21(i,k,j)*(dvdx(i,k,j)*erre12(i,k,j)+dvdy(i,k,j)*erre22(i,k,j)+dvdz(i,k,j)*erre32(i,k,j))+ &
              &     erre31(i,k,j)*(dwdx(i,k,j)*erre12(i,k,j)+dwdy(i,k,j)*erre22(i,k,j)+dwdz(i,k,j)*erre32(i,k,j)) 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      emme13(i,k,j)=erre11(i,k,j)*(dudx(i,k,j)*erre13(i,k,j)+dudy(i,k,j)*erre23(i,k,j)+dudz(i,k,j)*erre33(i,k,j))+ &
              &     erre21(i,k,j)*(dvdx(i,k,j)*erre13(i,k,j)+dvdy(i,k,j)*erre23(i,k,j)+dvdz(i,k,j)*erre33(i,k,j))+ &
              &     erre31(i,k,j)*(dwdx(i,k,j)*erre13(i,k,j)+dwdy(i,k,j)*erre23(i,k,j)+dwdz(i,k,j)*erre33(i,k,j)) 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      emme21(i,k,j)=erre12(i,k,j)*(dudx(i,k,j)*erre11(i,k,j)+dudy(i,k,j)*erre21(i,k,j)+dudz(i,k,j)*erre31(i,k,j))+ &
              &     erre22(i,k,j)*(dvdx(i,k,j)*erre11(i,k,j)+dvdy(i,k,j)*erre21(i,k,j)+dvdz(i,k,j)*erre31(i,k,j))+ &
              &     erre32(i,k,j)*(dwdx(i,k,j)*erre11(i,k,j)+dwdy(i,k,j)*erre21(i,k,j)+dwdz(i,k,j)*erre31(i,k,j)) 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      emme22(i,k,j)=erre12(i,k,j)*(dudx(i,k,j)*erre12(i,k,j)+dudy(i,k,j)*erre22(i,k,j)+dudz(i,k,j)*erre32(i,k,j))+ &
              &     erre22(i,k,j)*(dvdx(i,k,j)*erre12(i,k,j)+dvdy(i,k,j)*erre22(i,k,j)+dvdz(i,k,j)*erre32(i,k,j))+ &
              &     erre32(i,k,j)*(dwdx(i,k,j)*erre12(i,k,j)+dwdy(i,k,j)*erre22(i,k,j)+dwdz(i,k,j)*erre32(i,k,j)) 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      emme23(i,k,j)=erre12(i,k,j)*(dudx(i,k,j)*erre13(i,k,j)+dudy(i,k,j)*erre23(i,k,j)+dudz(i,k,j)*erre33(i,k,j))+ &
              &     erre22(i,k,j)*(dvdx(i,k,j)*erre13(i,k,j)+dvdy(i,k,j)*erre23(i,k,j)+dvdz(i,k,j)*erre33(i,k,j))+ &
              &     erre32(i,k,j)*(dwdx(i,k,j)*erre13(i,k,j)+dwdy(i,k,j)*erre23(i,k,j)+dwdz(i,k,j)*erre33(i,k,j)) 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      emme31(i,k,j)=erre13(i,k,j)*(dudx(i,k,j)*erre11(i,k,j)+dudy(i,k,j)*erre21(i,k,j)+dudz(i,k,j)*erre31(i,k,j))+ &
              &     erre23(i,k,j)*(dvdx(i,k,j)*erre11(i,k,j)+dvdy(i,k,j)*erre21(i,k,j)+dvdz(i,k,j)*erre31(i,k,j))+ &
              &     erre33(i,k,j)*(dwdx(i,k,j)*erre11(i,k,j)+dwdy(i,k,j)*erre21(i,k,j)+dwdz(i,k,j)*erre31(i,k,j)) 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      emme32(i,k,j)=erre13(i,k,j)*(dudx(i,k,j)*erre12(i,k,j)+dudy(i,k,j)*erre22(i,k,j)+dudz(i,k,j)*erre32(i,k,j))+ &
              &     erre23(i,k,j)*(dvdx(i,k,j)*erre12(i,k,j)+dvdy(i,k,j)*erre22(i,k,j)+dvdz(i,k,j)*erre32(i,k,j))+ &
              &     erre33(i,k,j)*(dwdx(i,k,j)*erre12(i,k,j)+dwdy(i,k,j)*erre22(i,k,j)+dwdz(i,k,j)*erre32(i,k,j)) 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      emme33(i,k,j)=erre13(i,k,j)*(dudx(i,k,j)*erre13(i,k,j)+dudy(i,k,j)*erre23(i,k,j)+dudz(i,k,j)*erre33(i,k,j))+ &
              &     erre23(i,k,j)*(dvdx(i,k,j)*erre13(i,k,j)+dvdy(i,k,j)*erre23(i,k,j)+dvdz(i,k,j)*erre33(i,k,j))+ &
              &     erre33(i,k,j)*(dwdx(i,k,j)*erre13(i,k,j)+dwdy(i,k,j)*erre23(i,k,j)+dwdz(i,k,j)*erre33(i,k,j)) 
     enddo
   enddo
  enddo
  
  
  
  ! omega antisymmetric matrix 
  
  allocate(omega1(nx,fpz,fpy))
  allocate(omega2(nx,fpz,fpy))
  allocate(omega3(nx,fpz,fpy))
  
  epsnum=0.000000001d0 
  
  !$acc parallel loop collapse(3)
  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        mask=merge(1.0d0, 0.0d0, &
                  (dabs(diag22(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                  (dabs(diag33(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                  (dabs(diag33(i,k,j)-diag22(i,k,j))<=epsnum))
        !!!!!!!!!!!!!!!!
        omega1(i,k,j)=(1.d0-mask)*((diag22(i,k,j)*emme12(i,k,j)+diag11(i,k,j)*emme21(i,k,j))/ &
                                   (diag22(i,k,j)-diag11(i,k,j)+epsnum))
        omega2(i,k,j)=(1.d0-mask)*((diag33(i,k,j)*emme13(i,k,j)+diag11(i,k,j)*emme31(i,k,j))/ &
                                   (diag33(i,k,j)-diag11(i,k,j)+epsnum))
        omega3(i,k,j)=(1.d0-mask)*((diag33(i,k,j)*emme23(i,k,j)+diag22(i,k,j)*emme32(i,k,j))/ &
                                   (diag33(i,k,j)-diag22(i,k,j)+epsnum))
      enddo
    enddo
  enddo
  
  deallocate(emme12)
  deallocate(emme13)
  deallocate(emme21)
  deallocate(emme23)
  deallocate(emme31)
  deallocate(emme32)
    
  
  !$acc parallel loop collapse(3)
  do j=1,fpy
   do k=1,fpz
     do i=1,nx
        omega1(i,k,j)=omega1(i,k,j)*(erre22(i,k,j)*erre11(i,k,j)-erre21(i,k,j)*erre12(i,k,j))+&
                    & omega2(i,k,j)*(erre23(i,k,j)*erre11(i,k,j)-erre21(i,k,j)*erre13(i,k,j))+&
                    & omega3(i,k,j)*(erre23(i,k,j)*erre12(i,k,j)-erre22(i,k,j)*erre13(i,k,j))
        omega2(i,k,j)=omega1(i,k,j)*(erre32(i,k,j)*erre11(i,k,j)-erre31(i,k,j)*erre12(i,k,j))+&
                    & omega2(i,k,j)*(erre33(i,k,j)*erre11(i,k,j)-erre31(i,k,j)*erre13(i,k,j))+&
                    & omega3(i,k,j)*(erre33(i,k,j)*erre12(i,k,j)-erre32(i,k,j)*erre13(i,k,j))
        omega3(i,k,j)=omega1(i,k,j)*(erre32(i,k,j)*erre21(i,k,j)-erre31(i,k,j)*erre22(i,k,j))+&
                    & omega2(i,k,j)*(erre33(i,k,j)*erre21(i,k,j)-erre31(i,k,j)*erre23(i,k,j))+&
                    & omega3(i,k,j)*(erre33(i,k,j)*erre22(i,k,j)-erre32(i,k,j)*erre23(i,k,j))
     enddo
   enddo
  enddo
  
  ! get (logC)^n=psi^n
  
  !$acc kernels
  diag11=log(diag11)
  diag22=log(diag22)
  diag33=log(diag33)
  !$acc end kernels
  
  !$acc parallel loop collapse(3)
  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        C_xxf(i,k,j)=diag11(i,k,j)*erre11(i,k,j)**2.d0+ diag22(i,k,j)*erre12(i,k,j)**2.d0+ diag33(i,k,j)*erre13(i,k,j)**2.d0 
        !!!!!!!!!!!!!!!!!!!!!!!!
        C_xyf(i,k,j)=diag11(i,k,j)*erre11(i,k,j)*erre21(i,k,j) + &
              &      diag22(i,k,j)*erre12(i,k,j)*erre22(i,k,j)+diag33(i,k,j)*erre13(i,k,j)*erre23(i,k,j)
        !!!!!!!!!!!!!!!!!!!!!!!!
        C_xzf(i,k,j)=diag11(i,k,j)*erre11(i,k,j)*erre31(i,k,j) + &
              &      diag22(i,k,j)*erre12(i,k,j)*erre32(i,k,j)+diag33(i,k,j)*erre13(i,k,j)*erre33(i,k,j)
        !!!!!!!!!!!!!!!!!!!!!!!!
        C_yyf(i,k,j)=diag11(i,k,j)*erre21(i,k,j)**2.d0+ diag22(i,k,j)*erre22(i,k,j)**2.d0+ diag33(i,k,j)*erre23(i,k,j)**2.d0
        !!!!!!!!!!!!!!!!!!!!!!!!
        C_yzf(i,k,j)=diag11(i,k,j)*erre21(i,k,j)*erre31(i,k,j) + &
              &      diag22(i,k,j)*erre22(i,k,j)*erre32(i,k,j)+diag33(i,k,j)*erre23(i,k,j)*erre33(i,k,j)
        !!!!!!!!!!!!!!!!!!!!!!!!
        C_zzf(i,k,j)=diag11(i,k,j)*erre31(i,k,j)**2.d0+ diag22(i,k,j)*erre32(i,k,j)**2.d0+ diag33(i,k,j)*erre33(i,k,j)**2.d0
      enddo
    enddo
   enddo
  
  
  !$acc kernels
  diag11=exp(diag11)
  diag22=exp(diag22)
  diag33=exp(diag33)
  !$acc end kernels  
  
   call phys_to_spectral(C_xxf,C_xxs,0)
   call phys_to_spectral(C_xyf,C_xys,0)
   call phys_to_spectral(C_xzf,C_xzs,0)
   call phys_to_spectral(C_yyf,C_yys,0)
   call phys_to_spectral(C_yzf,C_yzs,0)
   call phys_to_spectral(C_zzf,C_zzs,0)
  
  
end subroutine

  
subroutine assemble_eigen

integer :: i,j,k

double precision :: CTOT(3,3),Eigenvectors(3,3),Eigenvalues(3)

!GET PSI^n+1 in physical
call spectral_to_phys(C_xxs,C_xxf,0)
call spectral_to_phys(C_xys,C_xyf,0)
call spectral_to_phys(C_xzs,C_xzf,0)
call spectral_to_phys(C_yys,C_yyf,0)
call spectral_to_phys(C_yzs,C_yzf,0)
call spectral_to_phys(C_zzs,C_zzf,0)


!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz 
     do i=1,nx 
        CTOT = reshape([C_xxf(i,k,j), C_xyf(i,k,j), C_xzf(i,k,j), &
                        C_xyf(i,k,j), C_yyf(i,k,j), C_yzf(i,k,j), &
                        C_xzf(i,k,j), C_yzf(i,k,j), C_zzf(i,k,j)], [3,3])
    
        call Diagonalization(Eigenvalues,Eigenvectors,CTOT,3) 
        !call Diagonalization(CTOT,Eigenvectors,Eigenvalues) 

        diag11(i,k,j)=Eigenvalues(1)
        diag22(i,k,j)=Eigenvalues(2)
        diag33(i,k,j)=Eigenvalues(3)
        erre11(i,k,j)=Eigenvectors(1,1)
        erre12(i,k,j)=Eigenvectors(1,2)
        erre13(i,k,j)=Eigenvectors(1,3)
        erre21(i,k,j)=Eigenvectors(2,1)
        erre22(i,k,j)=Eigenvectors(2,2)
        erre23(i,k,j)=Eigenvectors(2,3)
        erre31(i,k,j)=Eigenvectors(3,1)
        erre32(i,k,j)=Eigenvectors(3,2)
        erre33(i,k,j)=Eigenvectors(3,3)
     enddo
   enddo
enddo  

!$acc kernels
diag11=exp(diag11)
diag22=exp(diag22)
diag33=exp(diag33)
!$acc end kernels

end subroutine

end module assemble



