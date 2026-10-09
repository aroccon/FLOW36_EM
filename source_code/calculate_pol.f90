subroutine calculate_pol(hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz)
  use commondata
  use par_size
  use polymer
  use sim_par
  use wavenumber
  use dctz_bwd_module
  use velocity
  
  double precision, dimension(spx,spy) :: beta2c
  double precision, dimension(2) :: s_bc
  double precision, dimension(nz) :: s_force
  double precision,intent(inout), dimension(spx,nz,spy,2) :: hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz  
  integer :: i,j,indx,indy

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

   do j=1,spy
    do i=1,spx 
      !!!!!!!!!!!!!!! Cxx
      s_force=hconfxx(i,:,j,1)
      s_bc(1)=bc_cxx(i,nz,j,1)
                                   !s_bc(2)=bc_cxx(i,1,j,1)
      s_bc(2)=sIdiag(i,1,j,1)
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfxx(i,:,j,1)=s_force
      s_force=hconfxx(i,:,j,2)
      s_bc(1)=bc_cxx(i,nz,j,2)
                                    !s_bc(2)=bc_cxx(i,1,j,2)
      s_bc(2)=sIdiag(i,1,j,2)
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfxx(i,:,j,2)=s_force
      !!!!!!!!!!!!!!! Cxy
      s_force=hconfxy(i,:,j,1)
      s_bc(1)=bc_cxy(i,nz,j,1)
                                    !s_bc(2)=bc_cxy(i,1,j,1)
      s_bc(2)=0.d0
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfxy(i,:,j,1)=s_force
      s_force=hconfxy(i,:,j,2)
      s_bc(1)=bc_cxy(i,nz,j,2)
                                     !s_bc(2)=bc_cxy(i,1,j,2)
      s_bc(2)=0.d0
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfxy(i,:,j,2)=s_force
      !!!!!!!!!!!!!!! Cxz
      s_force=hconfxz(i,:,j,1)
      s_bc(1)=bc_cxz(i,nz,j,1)
                                     !s_bc(2)=bc_cxz(i,1,j,1)
      s_bc(2)=0.d0
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfxz(i,:,j,1)=s_force
      s_force=hconfxz(i,:,j,2)
                                     !s_bc(1)=bc_cxz(i,nz,j,2)
      s_bc(2)=0.d0
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfxz(i,:,j,2)=s_force
      !!!!!!!!!!!!!!! Cyy
      s_force=hconfyy(i,:,j,1)
      s_bc(1)=bc_cyy(i,nz,j,1)
                                     !s_bc(2)=bc_cyy(i,1,j,1)
      s_bc(2)=sIdiag(i,1,j,1)
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfyy(i,:,j,1)=s_force
      s_force=hconfyy(i,:,j,2)
      s_bc(1)=bc_cyy(i,nz,j,2)
                                     !s_bc(2)=bc_cyy(i,1,j,2)
      s_bc(2)=sIdiag(i,1,j,2)
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfyy(i,:,j,2)=s_force
      !!!!!!!!!!!!!!! Cyz
      s_force=hconfyz(i,:,j,1)
      s_bc(1)=bc_cyz(i,nz,j,1)
                                     !s_bc(2)=bc_cyz(i,1,j,1)
      s_bc(2)=0.d0
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfyz(i,:,j,1)=s_force
      s_force=hconfyz(i,:,j,2)
      s_bc(1)=bc_cyz(i,nz,j,2)
                                     !s_bc(2)=bc_cyz(i,1,j,2)
      s_bc(2)=0.d0
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfyz(i,:,j,2)=s_force
      !!!!!!!!!!!!!!! Czz
      s_force=hconfzz(i,:,j,1)
      s_bc(1)=bc_czz(i,nz,j,1)
                                     !s_bc(2)=bc_czz(i,1,j,1)
      s_bc(2)=sIdiag(i,1,j,1)
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfzz(i,:,j,1)=s_force
      s_force=hconfzz(i,:,j,2)
      s_bc(1)=bc_czz(i,nz,j,2)
                                     !s_bc(2)=bc_czz(i,1,j,2)
      s_bc(2)=sIdiag(i,1,j,2)
      call helmholtz_rred(s_force,beta2c(i,j),[1.0d0,1.0d0],[0.0d0,0.0d0],s_bc,zp)
      hconfzz(i,:,j,2)=s_force                           
   enddo
  enddo
   
  return
end subroutine calculate_pol 