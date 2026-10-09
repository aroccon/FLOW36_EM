subroutine hist_term_pol(hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz)

    use commondata
    use par_size
    use polymer
    use wavenumber
    use sim_par
    
    double precision, dimension(spx,nz,spy,2) :: hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz
    double precision, dimension(spx,nz,spy,2) :: dC_xxs,ddC_xxs,dC_xys,ddC_xys,dC_xzs,ddC_xzs
    double precision, dimension(spx,nz,spy,2) :: dC_yys,ddC_yys,dC_yzs,ddC_yzs,dC_zzs,ddC_zzs

    
    integer :: indx,indy
    integer :: i,j,k
    
    
    indx=cstart(1)
    indy=cstart(3)
    
  !$acc kernels
   bc_cxx=hconfxx+C_xxs
   bc_cxy=hconfxy+C_xys
   bc_cxz=hconfxz+C_xzs
   bc_cyy=hconfyy+C_yys
   bc_cyz=hconfyz+C_yzs
   bc_czz=hconfzz+C_zzs
  !$acc end kernels

   call spectral_to_phys_bc(bc_cxx,bc_cxx,0)
   call spectral_to_phys_bc(bc_cxy,bc_cxy,0)
   call spectral_to_phys_bc(bc_cxz,bc_cxz,0)
   call spectral_to_phys_bc(bc_cyy,bc_cyy,0)
   call spectral_to_phys_bc(bc_cyz,bc_cyz,0)
   call spectral_to_phys_bc(bc_czz,bc_czz,0)

    ! hconfxx
call dz(C_xxs,dC_xxs)
call dz(dC_xxs,ddC_xxs)

!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    hconfxx(i,:,j,1)=hconfxx(i,:,j,1)+gammaconf*ddC_xxs(i,:,j,1)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_xxs(i,:,j,1) 
    hconfxx(i,:,j,2)=hconfxx(i,:,j,2)+gammaconf*ddC_xxs(i,:,j,2)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_xxs(i,:,j,2)
  enddo 
enddo 

    ! hconfxy
call dz(C_xys,dC_xys)
call dz(dC_xys,ddC_xys)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
   hconfxy(i,:,j,1)=hconfxy(i,:,j,1)+gammaconf*ddC_xys(i,:,j,1)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_xys(i,:,j,1) 
   hconfxy(i,:,j,2)=hconfxy(i,:,j,2)+gammaconf*ddC_xys(i,:,j,2)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_xys(i,:,j,2)
   enddo 
enddo 

    ! hconfxz
call dz(C_xzs,dC_xzs)
call dz(dC_xzs,ddC_xzs)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
   hconfxz(i,:,j,1)=hconfxz(i,:,j,1)+gammaconf*ddC_xzs(i,:,j,1)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_xzs(i,:,j,1) 
   hconfxz(i,:,j,2)=hconfxz(i,:,j,2)+gammaconf*ddC_xzs(i,:,j,2)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_xzs(i,:,j,2)
   enddo 
enddo 

    ! hconfyy
call dz(C_yys,dC_yys)
call dz(dC_yys,ddC_yys)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
   hconfyy(i,:,j,1)=hconfyy(i,:,j,1)+gammaconf*ddC_yys(i,:,j,1)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_yys(i,:,j,1) 
   hconfyy(i,:,j,2)=hconfyy(i,:,j,2)+gammaconf*ddC_yys(i,:,j,2)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_yys(i,:,j,2)
   enddo 
enddo 

    ! hconfyz
call dz(C_yzs,dC_yzs)
call dz(dC_yzs,ddC_yzs)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
   hconfyz(i,:,j,1)=hconfyz(i,:,j,1)+gammaconf*ddC_yzs(i,:,j,1)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_yzs(i,:,j,1) 
   hconfyz(i,:,j,2)=hconfyz(i,:,j,2)+gammaconf*ddC_yzs(i,:,j,2)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_yzs(i,:,j,2)
   enddo 
enddo 

    ! hconfzz
call dz(C_zzs,dC_zzs)
call dz(dC_zzs,ddC_zzs)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
   hconfzz(i,:,j,1)=hconfzz(i,:,j,1)+gammaconf*ddC_zzs(i,:,j,1)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_zzs(i,:,j,1) 
   hconfzz(i,:,j,2)=hconfzz(i,:,j,2)+gammaconf*ddC_zzs(i,:,j,2)+(1.0d0-k2(i+indx,j+indy)*gammaconf)*C_zzs(i,:,j,2)
   enddo 
enddo 
   
return
end