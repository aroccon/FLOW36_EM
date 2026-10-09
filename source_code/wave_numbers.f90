subroutine wave_numbers

  use commondata
  use wavenumber
  use sim_par
  use phase_field
  use dual_grid
  use polymer
  
  integer :: i,j
  
#define match_dens matched_density
#define match_visc matched_viscosity
#define polflag polcompflag 
#define fenepflag fenepcompflag
#define phiflag phicompflag
#define carflag carnatureflag
#define non_new_flag non_newtonian
#define car_carreau_flag carriercarreauflag
  
  kx(1)=0.0d0
  do i=2,nx/2+1
    kx(i)=dble(i-1)*2.0d0*pi/xl
  enddo
  
  ky(1)=0.0d0
  do i=2,ny/2+1
    ky(ny-i+2)=-dble(i-1)*2.0d0*pi/yl
    ky(i)=dble(i-1)*2.0d0*pi/yl
  enddo
  
  do j=1,ny
    do i=1,nx/2+1
      k2(i,j)=kx(i)*kx(i)+ky(j)*ky(j)
    enddo
  enddo
  
  allocate(kxpsi(npsix/2+1))
  allocate(kypsi(npsiy))
  allocate(k2psi(npsix/2+1,npsiy))
  
  kxpsi(1)=0.0d0
  do i=2,npsix/2+1
    kxpsi(i)=dble(i-1)*2.0d0*pi/xl
  enddo
  
  kypsi(1)=0.0d0
  do i=2,npsiy/2+1
    kypsi(npsiy-i+2)=-dble(i-1)*2.0d0*pi/yl
    kypsi(i)=dble(i-1)*2.0d0*pi/yl
  enddo
  
  do j=1,npsiy
    do i=1,npsix/2+1
      k2psi(i,j)=kxpsi(i)*kxpsi(i)+kypsi(j)*kypsi(j)
    enddo
  enddo
  
#if (polflag == 0 && non_new_flag == 0)
#if match_dens == 2
gamma=dt/(2.0d0*re*rhor)
#else
gamma=dt/(2.0d0*re)
#endif
!!!!!!
#elif (phiflag == 0 && polflag == 1)
gammaconf=dt/(2.0d0*re*Sc)
#if fenepflag == 1
CP=-dt/(2.0d0*Wi) 
#elif fenepflag == 2
CP=(dt/(2.d0*Wi))-(3.d0*dt)/(2.d0*Wi*L**2.d0)
#endif
#if match_dens == 2
gamma=(dt*xo)/(2.0d0*re*rhor)
#else
gamma=(dt*xo)/(2.0d0*re)
#endif
!!!!!!
#elif (phiflag == 1 && polflag == 1)

#if fenepflag == 1
CP=-dt/(2.0d0*Wi) 
#elif fenepflag == 2
CP=(dt/(2.d0*Wi))-(3.d0*dt)/(2.d0*Wi*L**2.d0)
#endif
#if carflag == 0 
gammaconf=dt/(2.0d0*re*Sc)
#if match_dens == 2
if (visr .le. xo) then
gamma=(dt*xo)/(2.0d0*re*rhor)
else if (xo < visr) then
gamma=(dt*visr)/(2.0d0*re*rhor)
endif 
#else
if (visr .le. xo) then
gamma=(dt*xo)/(2.0d0*re)
else if (xo < visr) then
gamma=(dt*visr)/(2.0d0*re)
endif
#endif
!!!!!!
#elif carflag == 1
gammaconf=(dt*visr)/(2.0d0*re*Sc)
#if match_dens == 2
gamma=dt/(2.0d0*re*rhor)
#else
gamma=dt/(2.0d0*re)
#endif
#endif
#endif

#if polflag == 1
#if (match_visc == 2 && carflag == 0)
gamma=gamma
#elif (match_visc == 2 && carflag == 1)
gamma=gamma*(visr*xo)
#endif
#else
gamma=gamma*visr
#endif

#if non_newtonian == 1
#if car_carreau_flag == 0 
gamma=(dt*visr)/(2.0d0*re)
#elif car_carreau_flag == 1
gamma=dt/(2.0d0*re)
#endif
#endif

return
end