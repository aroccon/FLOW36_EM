subroutine pol_stress_tensor(pols1,pols2,pols3)

  use commondata
  use par_size
  use phase_field
  use wavenumber
  use velocity
  use sim_par
  use polymer
 
  double precision, dimension(spx,nz,spy,2) :: pols1,pols2,pols3
  double precision, allocatable, dimension(:,:,:,:) :: a4,a5,a6,a7
  double precision, allocatable, dimension(:,:,:) :: a4f,a5f,a6f,a7f
  double precision, allocatable, dimension(:,:,:) :: pC_xxf,pC_xyf,pC_xzf,pC_yyf,pC_yzf,pC_zzf
  double precision, allocatable, dimension(:,:,:,:) :: pC_xxs,pC_xys,pC_xzs,pC_yys,pC_yzs,pC_zzs
  double precision :: phif,alpha,alpha_smooth

  integer :: i,j,k
 
 
#define phiflag phicompflag
#define match_dens matched_density
#define polflag polcompflag 
#define fenepflag fenepcompflag
#define carflag carnatureflag
 
#if (phiflag == 1 && polflag == 1)
 call spectral_to_phys(phic,phi,1)
#endif
 
 allocate(a4f(nx,fpz,fpy))
 allocate(a5f(nx,fpz,fpy))
 allocate(a6f(nx,fpz,fpy))
 allocate(a7f(nx,fpz,fpy))
 allocate(a4(spx,nz,spy,2))
 allocate(a5(spx,nz,spy,2))
 allocate(a6(spx,nz,spy,2))
 allocate(a7(spx,nz,spy,2))
 
 
 call spectral_to_phys(C_xxs,C_xxf,0)
 call spectral_to_phys(C_xys,C_xyf,0)
 call spectral_to_phys(C_xzs,C_xzf,0)
 call spectral_to_phys(C_yys,C_yyf,0)
 call spectral_to_phys(C_yzs,C_yzf,0)
 call spectral_to_phys(C_zzs,C_zzf,0)

#if (fenepflag == 1 || fenepflag == 2)
#if fenepflag == 1  
 zitaVel=(L**2.0d0-3.0d0)/(L**2.0d0-(C_xxf+C_yyf+C_zzf))
 allocate(pC_xxf(nx,fpz,fpy))
 allocate(pC_xyf(nx,fpz,fpy))
 allocate(pC_xzf(nx,fpz,fpy))
 allocate(pC_yyf(nx,fpz,fpy))
 allocate(pC_yzf(nx,fpz,fpy))
 allocate(pC_zzf(nx,fpz,fpy))
 !$acc parallel loop collapse(3)
 do j=1,fpy
   do k=1,fpz
     do i=1,nx
       pC_xxf(i,k,j)=C_xxf(i,k,j)*zitaVel(i,k,j)
       pC_xyf(i,k,j)=C_xyf(i,k,j)*zitaVel(i,k,j)
       pC_xzf(i,k,j)=C_xzf(i,k,j)*zitaVel(i,k,j)
       pC_yyf(i,k,j)=C_yyf(i,k,j)*zitaVel(i,k,j)
       pC_yzf(i,k,j)=C_yzf(i,k,j)*zitaVel(i,k,j)
       pC_zzf(i,k,j)=C_zzf(i,k,j)*zitaVel(i,k,j)
     enddo
   enddo
 enddo
 allocate(pC_xxs(spx,nz,spy,2))
 allocate(pC_xys(spx,nz,spy,2))
 allocate(pC_xzs(spx,nz,spy,2))
 allocate(pC_yys(spx,nz,spy,2))
 allocate(pC_yzs(spx,nz,spy,2))
 allocate(pC_zzs(spx,nz,spy,2))
 call phys_to_spectral(pC_xxf,pC_xxs,0)
 call phys_to_spectral(pC_xyf,pC_xys,0)
 call phys_to_spectral(pC_xzf,pC_xzs,0)
 call phys_to_spectral(pC_yyf,pC_yys,0)
 call phys_to_spectral(pC_yzf,pC_yzs,0)
 call phys_to_spectral(pC_zzf,pC_zzs,0)
 deallocate(pC_xxf)
 deallocate(pC_xyf)
 deallocate(pC_xzf)
 deallocate(pC_yyf)
 deallocate(pC_yzf)
 deallocate(pC_zzf)
#elif fenepflag == 2
  zitaVel=(L**2.0d0)/(L**2.0d0-(C_xxf+C_yyf+C_zzf))
  allocate(pC_xxf(nx,fpz,fpy))
  allocate(pC_xyf(nx,fpz,fpy))
  allocate(pC_xzf(nx,fpz,fpy))
  allocate(pC_yyf(nx,fpz,fpy))
  allocate(pC_yzf(nx,fpz,fpy))
  allocate(pC_zzf(nx,fpz,fpy))
  !$acc parallel loop collapse(3)
  do j=1,fpy
    do k=1,fpz
      do i=1,nx
        pC_xxf(i,k,j)=(C_xxf(i,k,j)-1.d0)*zitaVel(i,k,j)
        pC_xyf(i,k,j)=C_xyf(i,k,j)*zitaVel(i,k,j)
        pC_xzf(i,k,j)=C_xzf(i,k,j)*zitaVel(i,k,j)
        pC_yyf(i,k,j)=(C_yyf(i,k,j)-1.d0)*zitaVel(i,k,j)
        pC_yzf(i,k,j)=C_yzf(i,k,j)*zitaVel(i,k,j)
        pC_zzf(i,k,j)=(C_zzf(i,k,j)-1.d0)*zitaVel(i,k,j)
      enddo
    enddo
  enddo
  allocate(pC_xxs(spx,nz,spy,2))
  allocate(pC_xys(spx,nz,spy,2))
  allocate(pC_xzs(spx,nz,spy,2))
  allocate(pC_yys(spx,nz,spy,2))
  allocate(pC_yzs(spx,nz,spy,2))
  allocate(pC_zzs(spx,nz,spy,2))
  call phys_to_spectral(pC_xxf,pC_xxs,0)
  call phys_to_spectral(pC_xyf,pC_xys,0)
  call phys_to_spectral(pC_xzf,pC_xzs,0)
  call phys_to_spectral(pC_yyf,pC_yys,0)
  call phys_to_spectral(pC_yzf,pC_yzs,0)
  call phys_to_spectral(pC_zzf,pC_zzs,0)
  deallocate(pC_xxf)
  deallocate(pC_xyf)
  deallocate(pC_xzf)
  deallocate(pC_yyf)
  deallocate(pC_yzf)
  deallocate(pC_zzf)
#endif
 !!!!!!!!!!!!!!!!!!!! FENE-P
 ! assemble x component of the polymer stress
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a4(i,k,j,1)=-kx(i+cstart(1))*pC_xxs(i,k,j,2)-ky(j+cstart(3))*pC_xys(i,k,j,2)
       a4(i,k,j,2)=+kx(i+cstart(1))*pC_xxs(i,k,j,1)+ky(j+cstart(3))*pC_xys(i,k,j,1)
     enddo
   enddo
 enddo
 
 call dz(pC_xzs,a5)
 
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a4(i,k,j,1)=a4(i,k,j,1)+a5(i,k,j,1)
       a4(i,k,j,2)=a4(i,k,j,2)+a5(i,k,j,2)
     enddo
   enddo
 enddo
 
 
!assemble y component of the polymer stress
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a5(i,k,j,1)=-kx(i+cstart(1))*pC_xys(i,k,j,2)-ky(j+cstart(3))*pC_yys(i,k,j,2)
       a5(i,k,j,2)=+kx(i+cstart(1))*pC_xys(i,k,j,1)+ky(j+cstart(3))*pC_yys(i,k,j,1)
     enddo
   enddo
 enddo
 
 call dz(pC_yzs,a6)
 
 
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a5(i,k,j,1)=a5(i,k,j,1)+a6(i,k,j,1)
       a5(i,k,j,2)=a5(i,k,j,2)+a6(i,k,j,2)
     enddo
   enddo
 enddo
 
 
!assemble z component of the polymer stress
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a6(i,k,j,1)=-kx(i+cstart(1))*pC_xzs(i,k,j,2)-ky(j+cstart(3))*pC_yzs(i,k,j,2)
       a6(i,k,j,2)=+kx(i+cstart(1))*pC_xzs(i,k,j,1)+ky(j+cstart(3))*pC_yzs(i,k,j,1)
     enddo
   enddo
 enddo
 
 call dz(pC_zzs,a7)

 deallocate(pC_xxs)
 deallocate(pC_xys)
 deallocate(pC_xzs)
 deallocate(pC_yys)
 deallocate(pC_yzs)
 deallocate(pC_zzs)
 
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a6(i,k,j,1)=a6(i,k,j,1)+a7(i,k,j,1)
       a6(i,k,j,2)=a6(i,k,j,2)+a7(i,k,j,2)
     enddo
   enddo
 enddo
 
#else
if(rank.eq.0)write(*,*) "Pol stress Saramito"
#if fenepflag == 3
  allocate(pC_xxf(nx,fpz,fpy))
  allocate(pC_xyf(nx,fpz,fpy))
  allocate(pC_xzf(nx,fpz,fpy))
  allocate(pC_yyf(nx,fpz,fpy))
  allocate(pC_yzf(nx,fpz,fpy))
  allocate(pC_zzf(nx,fpz,fpy))
do j=1,fpy
  do k=1,fpz
    do i=1,nx
    phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
    pC_xxf(i,k,j)=(C_xxf(i,k,j)-1.d0)*(1.d0-phif)/2.d0 
    pC_yyf(i,k,j)=(C_yyf(i,k,j)-1.d0)*(1.d0-phif)/2.d0 
    pC_zzf(i,k,j)=(C_zzf(i,k,j)-1.d0)*(1.d0-phif)/2.d0 
    pC_xyf(i,k,j)=C_xyf(i,k,j)*(1.d0-phif)/2.d0
    pC_xzf(i,k,j)=C_xzf(i,k,j)*(1.d0-phif)/2.d0
    pC_yzf(i,k,j)=C_yzf(i,k,j)*(1.d0-phif)/2.d0
  enddo
 enddo
enddo

  allocate(pC_xxs(spx,nz,spy,2))
  allocate(pC_xys(spx,nz,spy,2))
  allocate(pC_xzs(spx,nz,spy,2))
  allocate(pC_yys(spx,nz,spy,2))
  allocate(pC_yzs(spx,nz,spy,2))
  allocate(pC_zzs(spx,nz,spy,2))
  call phys_to_spectral(pC_xxf,pC_xxs,0)
  call phys_to_spectral(pC_xyf,pC_xys,0)
  call phys_to_spectral(pC_xzf,pC_xzs,0)
  call phys_to_spectral(pC_yyf,pC_yys,0)
  call phys_to_spectral(pC_yzf,pC_yzs,0)
  call phys_to_spectral(pC_zzf,pC_zzs,0)
  deallocate(pC_xxf)
  deallocate(pC_xyf)
  deallocate(pC_xzf)
  deallocate(pC_yyf)
  deallocate(pC_yzf)
  deallocate(pC_zzf)
 
#endif

 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a4(i,k,j,1)=-kx(i+cstart(1))*pC_xxs(i,k,j,2)-ky(j+cstart(3))*pC_xys(i,k,j,2)
       a4(i,k,j,2)=+kx(i+cstart(1))*pC_xxs(i,k,j,1)+ky(j+cstart(3))*pC_xys(i,k,j,1)
     enddo
   enddo
 enddo
 
 call dz(pC_xzs,a5)
 
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a4(i,k,j,1)=a4(i,k,j,1)+a5(i,k,j,1)
       a4(i,k,j,2)=a4(i,k,j,2)+a5(i,k,j,2)
     enddo
   enddo
 enddo
 
 
!assemble y component of the polymer stress
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a5(i,k,j,1)=-kx(i+cstart(1))*pC_xys(i,k,j,2)-ky(j+cstart(3))*pC_yys(i,k,j,2)
       a5(i,k,j,2)=+kx(i+cstart(1))*pC_xys(i,k,j,1)+ky(j+cstart(3))*pC_yys(i,k,j,1)
     enddo
   enddo
 enddo
 
 call dz(pC_yzs,a6)
 
 
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a5(i,k,j,1)=a5(i,k,j,1)+a6(i,k,j,1)
       a5(i,k,j,2)=a5(i,k,j,2)+a6(i,k,j,2)
     enddo
   enddo
 enddo
 
 
!assemble z component of the polymer stress
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a6(i,k,j,1)=-kx(i+cstart(1))*pC_xzs(i,k,j,2)-ky(j+cstart(3))*pC_yzs(i,k,j,2)
       a6(i,k,j,2)=+kx(i+cstart(1))*pC_xzs(i,k,j,1)+ky(j+cstart(3))*pC_yzs(i,k,j,1)
     enddo
   enddo
 enddo
 
 call dz(pC_zzs,a7)

 deallocate(pC_xxs)
 deallocate(pC_xys)
 deallocate(pC_xzs)
 deallocate(pC_yys)
 deallocate(pC_yzs)
 deallocate(pC_zzs)
 
 !$acc parallel loop collapse(3)
 do j=1,spy
   do k=1,nz
     do i=1,spx
       a6(i,k,j,1)=a6(i,k,j,1)+a7(i,k,j,1)
       a6(i,k,j,2)=a6(i,k,j,2)+a7(i,k,j,2)
     enddo
   enddo
 enddo

#endif

! add to pols1/pols2/pols3 terms 
#if match_dens == 2
!$acc kernels
! rescale NS if rhor > 1 for improved stability
pols1=((1.0d0-xo)*a4)/(re*Wi*rhor)
pols2=((1.0d0-xo)*a5)/(re*Wi*rhor)
pols3=((1.0d0-xo)*a6)/(re*Wi*rhor)
!$acc end kernels
#else
!$acc kernels
pols1=((1.0d0-xo)*a4)/(re*Wi)
pols2=((1.0d0-xo)*a5)/(re*Wi)
pols3=((1.0d0-xo)*a6)/(re*Wi)
!$acc end kernels

#if (phiflag == 1 && carflag == 1)
pols1=pols1*visr
pols2=pols2*visr
pols3=pols3*visr
#endif

#endif
    
deallocate(a4f)
deallocate(a5f)
deallocate(a6f)
deallocate(a7f)
deallocate(a4)
deallocate(a5)
deallocate(a6)
deallocate(a7)

return
end
    