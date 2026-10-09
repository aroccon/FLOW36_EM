subroutine calculate_pol(hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz)
  use commondata
  use par_size
  use polymer
  use sim_par
  use wavenumber
  use dctz_bwd_module
  use velocity
  
  double precision, dimension(spx,spy) :: beta2c
  double precision, allocatable :: r1(:,:,:),r2(:,:,:)
  double precision,intent(inout), dimension(spx,nz,spy,2) :: hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz  
  integer :: i,j,c,indx,indy

  indx=cstart(1)
  indy=cstart(3)
  
  ! r term of the boundary conditions
  !$acc parallel loop collapse(2)
  do j=1,spy
    do i=1,spx
      beta2c(i,j)=1.0d0/gammaconf+k2(i+indx,j+indy)
    enddo
  enddo  
  
  !$acc kernels
   hconfxx=-hconfxx/gammaconf
   hconfxy=-hconfxy/gammaconf
   hconfxz=-hconfxz/gammaconf
   hconfyy=-hconfyy/gammaconf 
   hconfyz=-hconfyz/gammaconf
   hconfzz=-hconfzz/gammaconf
  !$acc end kernels



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Global diffusion is implemented for each mode of the boundary condition
  ! keep in mind that diffusion is applied for the real and the immaginary part 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! one solve per tensor component for all the (i,j) columns and both real/imaginary parts
  ! (helmholtz_cols, one GPU kernel) instead of a host loop over the columns calling
  ! helmholtz_rred: same algorithm, operations and boundary values, bitwise identical on CPU.
  ! r1: boundary value at z(1)=+1 (bc_c..(i,nz,j,c)), r2: at z(2)=-1 (sIdiag for the diagonal
  ! components, 0 otherwise)
  allocate(r1(spx,spy,2),r2(spx,spy,2))
  !!!!!!!!!!!!!!! Cxx
  !$acc parallel loop collapse(3)
  do c=1,2
   do j=1,spy
    do i=1,spx
    r1(i,j,c)=bc_cxx(i,nz,j,c)
    r2(i,j,c)=sIdiag(i,1,j,c)
    enddo
   enddo
  enddo
  call helmholtz_cols(hconfxx,beta2c,[1.0d0,1.0d0],[0.0d0,0.0d0],r1,r2,zp)
  !!!!!!!!!!!!!!! Cxy
  !$acc parallel loop collapse(3)
  do c=1,2
   do j=1,spy
    do i=1,spx
    r1(i,j,c)=bc_cxy(i,nz,j,c)
    r2(i,j,c)=0.0d0
    enddo
   enddo
  enddo
  call helmholtz_cols(hconfxy,beta2c,[1.0d0,1.0d0],[0.0d0,0.0d0],r1,r2,zp)
  !!!!!!!!!!!!!!! Cxz
   ! imaginary part: same boundary value as the real part (bc_cxz(i,nz,j,1)), as in the
   ! previous per-column version where s_bc(1)=bc_cxz(i,nz,j,2) was commented out
  !$acc parallel loop collapse(3)
  do c=1,2
   do j=1,spy
    do i=1,spx
    r1(i,j,c)=bc_cxz(i,nz,j,1)
    r2(i,j,c)=0.0d0
    enddo
   enddo
  enddo
  call helmholtz_cols(hconfxz,beta2c,[1.0d0,1.0d0],[0.0d0,0.0d0],r1,r2,zp)
  !!!!!!!!!!!!!!! Cyy
  !$acc parallel loop collapse(3)
  do c=1,2
   do j=1,spy
    do i=1,spx
    r1(i,j,c)=bc_cyy(i,nz,j,c)
    r2(i,j,c)=sIdiag(i,1,j,c)
    enddo
   enddo
  enddo
  call helmholtz_cols(hconfyy,beta2c,[1.0d0,1.0d0],[0.0d0,0.0d0],r1,r2,zp)
  !!!!!!!!!!!!!!! Cyz
  !$acc parallel loop collapse(3)
  do c=1,2
   do j=1,spy
    do i=1,spx
    r1(i,j,c)=bc_cyz(i,nz,j,c)
    r2(i,j,c)=0.0d0
    enddo
   enddo
  enddo
  call helmholtz_cols(hconfyz,beta2c,[1.0d0,1.0d0],[0.0d0,0.0d0],r1,r2,zp)
  !!!!!!!!!!!!!!! Czz
  !$acc parallel loop collapse(3)
  do c=1,2
   do j=1,spy
    do i=1,spx
    r1(i,j,c)=bc_czz(i,nz,j,c)
    r2(i,j,c)=sIdiag(i,1,j,c)
    enddo
   enddo
  enddo
  call helmholtz_cols(hconfzz,beta2c,[1.0d0,1.0d0],[0.0d0,0.0d0],r1,r2,zp)
  deallocate(r1,r2)
   
  return
end subroutine calculate_pol 