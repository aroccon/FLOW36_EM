subroutine solver(ntime)

  use commondata
  use par_size
  use sim_par
  use velocity
  use sterms
  use grid
  use phase_field
  use surfactant
  use temperature
  use particle
  use dual_grid
  use polymer
  use derivatives
  use assemble

double precision, allocatable, dimension(:,:,:,:) :: s1,s2,s3,sphi,hphi,spsi,hpsi,stheta,htheta
double precision, allocatable, dimension(:,:,:,:) :: pols1,pols2,pols3,polh1,polh2,polh3
double precision, allocatable, dimension(:,:,:,:) :: sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz 
double precision, allocatable, dimension(:,:,:,:) :: hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz
double precision, dimension(spx,nz,spy,2) :: h1,h2,h3,h,omega

integer :: ntime

#define openaccflag openacccompflag
#define phiflag phicompflag
#define phicorflag phicorcompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define polflag polcompflag 
#define artpolflag polymerdifflag
#define solvpolflag polysolvflag 
#define partflag particlecompflag
#define twowayc twowaycflag
#define boussinnesq boussinnesqcompflag
#define match_dens matched_density
#define expx expansionx
#define expy expansiony
#define expz expansionz
#define cpiflag cpicompflag

! flow part (not particle)
if(rank.lt.flow_comm_lim)then
  
! Cahn-Hilliard equation solution
#if phiflag == 1
allocate(sphi(spx,nz,spy,2))
allocate(hphi(spx,nz,spy,2))
! calculate non-linear terms of Cahn-Hilliard equation
call sterm_ch(sphi)
!print *, "sphi", sphi

if (ntime.eq.nstart+1)then
call euler_phi(sphi,hphi)
else
call adams_bashforth_phi(sphi,hphi)
endif

!$acc kernels
sphi_o=sphi
!$acc end kernels

deallocate(sphi)

!history term
!$acc kernels
hphi=hphi+phic
! save phase field at current value
!phi^(n)
 phicp=phic
!$acc end kernels

! Solving Cahn-Hilliard equation, 4th order (phicor 0 to 6)
#if phicorflag != 7 
!phase field helmoltz equations
call calculate_phi(hphi)
!now phic=phi^(n+1)
#endif 
!! Solving Allen-Cahn equation, 2nd order (phicor 7 or 8)
#if (phicorflag == 7 || phicorflag == 8)
call calculate_phi_ac(hphi) 
#endif 

deallocate(hphi)
             
call chop_modes(phic)
#endif  

! transform variables back to physical space and perform dealiasing
call spectral_to_phys(uc,u,1)
call spectral_to_phys(vc,v,1)
call spectral_to_phys(wc,w,1)

!conformation tensor equation
#if polflag == 1  
allocate(sconfxx(spx,nz,spy,2))
allocate(sconfxy(spx,nz,spy,2))
allocate(sconfxz(spx,nz,spy,2))
allocate(sconfyy(spx,nz,spy,2))
allocate(sconfyz(spx,nz,spy,2))
allocate(sconfzz(spx,nz,spy,2))


call velocity_derivatives


#if solvpolflag == 0
!compute the non linear terms of the 6 conformation equation
call conf_non_linear(sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz)
call destroy_derivatives
#endif


#if solvpolflag == 1 

! get matrix omega and B matrix and logC at t^n
call assemble_matrices

call conf_non_linear_pol(sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz)

call destroy_derivatives
#endif

allocate(hconfxx(spx,nz,spy,2))
allocate(hconfxy(spx,nz,spy,2))
allocate(hconfxz(spx,nz,spy,2))
allocate(hconfyy(spx,nz,spy,2))
allocate(hconfyz(spx,nz,spy,2))
allocate(hconfzz(spx,nz,spy,2))

! time integration of the non-linear terms
! for first time step the code uses an explicit Euler algorithm,
! while from the second time step on an Adams-Bashforth algorithm
if (ntime.eq.nstart+1)then
 call euler_pol(sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz,hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz)
else
  call adams_bashforth_pol(sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz,hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz)
endif

  !$acc kernels
  sconfxx_o=sconfxx + 0.0d0
  sconfxy_o=sconfxy + 0.0d0
  sconfxz_o=sconfxz + 0.0d0
  sconfyy_o=sconfyy + 0.0d0
  sconfyz_o=sconfyz + 0.0d0
  sconfzz_o=sconfzz + 0.0d0
  !$acc end kernels
  
deallocate(sconfxx)
deallocate(sconfxy)
deallocate(sconfxz)
deallocate(sconfyy)
deallocate(sconfyz)
deallocate(sconfzz)


#if artpolflag == 0 
!if(rank.eq.0)write(*,*) 'adams bashforth'

!$acc kernels
 hconfxx=hconfxx+C_xxs
 hconfxy=hconfxy+C_xys
 hconfxz=hconfxz+C_xzs
 hconfyy=hconfyy+C_yys
 hconfyz=hconfyz+C_yzs
 hconfzz=hconfzz+C_zzs
!$acc end kernels
 
#elif artpolflag == 1 

allocate(bc_cxx(spx,nz,spy,2))
allocate(bc_cxy(spx,nz,spy,2))
allocate(bc_cxz(spx,nz,spy,2))
allocate(bc_cyy(spx,nz,spy,2))
allocate(bc_cyz(spx,nz,spy,2))
allocate(bc_czz(spx,nz,spy,2))
  
 ! add linear contribution to history term
 call hist_term_pol(hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz)
 
! helmoltz equations for C
 call calculate_pol(hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz)
 
deallocate(bc_cxx)
deallocate(bc_cxy)
deallocate(bc_cxz)
deallocate(bc_cyy)
deallocate(bc_cyz)
deallocate(bc_czz)

#endif 

  !now C_xxs,C_xys,C_xzs,C_yys,C_yzs,C_zzs are the unknowns at the n+1 time step
  !$acc kernels
  C_xxs=hconfxx
  C_xys=hconfxy
  C_xzs=hconfxz
  C_yys=hconfyy
  C_yzs=hconfyz
  C_zzs=hconfzz
  !$acc end kernels 

#if(solvpolflag == 1)
  !hconf is psi^(n+1)=log(C)^(n+1)
  !get Eigenvalues and Eigenvectors
  call assemble_eigen
#endif

  deallocate(hconfxx)
  deallocate(hconfxy)
  deallocate(hconfxz)
  deallocate(hconfyy)
  deallocate(hconfyz)
  deallocate(hconfzz)
  
  call chop_modes(C_xxs)
  call chop_modes(C_xys)
  call chop_modes(C_xzs)
  call chop_modes(C_yys)
  call chop_modes(C_yzs)
  call chop_modes(C_zzs)
  
#endif  
 
!! Navier stokes
  allocate(s1(spx,nz,spy,2))
  allocate(s2(spx,nz,spy,2))
  allocate(s3(spx,nz,spy,2))
  
  ! why not adding some polymer stress
#if polflag == 1
  allocate(pols1(spx,nz,spy,2))
  allocate(pols2(spx,nz,spy,2))
  allocate(pols3(spx,nz,spy,2))
#endif
  
  !$acc kernels
  s1=0.0d0
  s2=0.0d0
  s3=0.0d0
#if polflag == 1
  pols1=0.0d0
  pols2=0.0d0
  pols3=0.0d0 
#endif
  !$acc end kernels
  
  ! calculate convective terms of N-S equation and store them in s1,s2,s3
  call convective_ns(s1,s2,s3)

! add mean pressure gradient to S term
#if cpiflag == 0
#if match_dens == 2
  ! rescale NS equation if rhor > 1 for improved stability
  !$acc kernels
    s1=s1-sgradpx/rhor
    s2=s2-sgradpy/rhor
   !$acc end kernels
#else
   !$acc kernels
      s1=s1-sgradpx
      s2=s2-sgradpy
   !$acc end kernels
#endif
#elif cpiflag == 1
#if match_dens == 2
  ! rescale NS equation if rhor > 1 for improved stability
  !$acc kernels
    s1=s1-sgradpx*dabs(gradpx)/rhor
    s2=s2-sgradpy*dabs(gradpy)/rhor
  !$acc end kernels
#else
  !$acc kernels
    s1=s1-sgradpx*dabs(gradpx)
    s2=s2-sgradpy*dabs(gradpy)
  !$acc end kernels
#endif
#endif

! add non-linear part of phase field to N-S non-linear terms
#if phiflag == 1
  call phi_non_linear(s1,s2,s3)
#endif

#if polflag == 1
  ! calculate polymer stress terms and store them in pols1,pols2,pols3
  call pol_stress_tensor(pols1,pols2,pols3)
#endif 

! add temperature contribution to N-S non-linear terms (Boussinesq)
#if tempflag == 1
#if boussinnesq == 1
  ! write(*,*) grav*Ra/(16.0d0*Pr*Re**2)
  !$acc kernels
    s1=s1-grav(1)*Ra/(16.0d0*Pr*Re**2)*thetac
    s2=s2-grav(3)*Ra/(16.0d0*Pr*Re**2)*thetac
    s3=s3-grav(2)*Ra/(16.0d0*Pr*Re**2)*thetac
  !$acc end kernels
#endif
#endif

#if polflag == 1
allocate(polh1(spx,nz,spy,2))
allocate(polh2(spx,nz,spy,2))
allocate(polh3(spx,nz,spy,2))
#endif

! time integration of the non-linear terms
! for first time step the code uses an explicit Euler algorithm,
! while from the second time step on an Adams-Bashforth algorithm
if (ntime.eq.nstart+1)then
call euler(s1,s2,s3,h1,h2,h3)
#if polflag == 1
call euler(pols1,pols2,pols3,polh1,polh2,polh3)
#endif
else
call adams_bashforth(s1,s2,s3,h1,h2,h3)
#if polflag == 1
call crank_nicholson_pol(pols1,pols2,pols3,polh1,polh2,polh3)
#endif
endif

!openACC: adding the useless +0.0d0, Memory movement on GPUs and there is performance gain
!Ignore the +0.0d0 if you are reading the code, it is only a computationl trick
!$acc kernels
s1_o=s1 + 0.0d0
s2_o=s2 + 0.0d0
s3_o=s3 + 0.0d0
#if polflag == 1
pols1_n=pols1 + 0.0d0
pols2_n=pols2 + 0.0d0
pols3_n=pols3 + 0.0d0
h1 = h1 + polh1
h2 = h2 + polh2
h3 = h3 + polh3
#endif
!$acc end kernels

deallocate(s1)
deallocate(s2)
deallocate(s3)

#if polflag == 1
deallocate(pols1)
deallocate(pols2)
deallocate(pols3)
deallocate(polh1)
deallocate(polh2)
deallocate(polh3)
#endif

! add linear contribution to history term
call hist_term(h1,h2,h3,h)

! solve Helmholtz equation for w
call calculate_w(h)

! solve Helmholtz equation for omega_z and store in h1
call calculate_omega(h1,h2,omega)

! calculate u,v from continuity and vorticity definition
call calculate_uv(omega,h1,h2)

call chop_modes(uc)
call chop_modes(vc)
call chop_modes(wc)
! only if surfactant calculated on finer grid than phase field
#if expx != 1 || expy != 1 || expz != 1
  ! for phase field calculation
  call spectral_to_phys(uc,u,1)
  call spectral_to_phys(vc,v,1)
  call spectral_to_phys(wc,w,1)
  ! for surfactant calculation
  call coarse2fine(uc,uc_fg)
  call spectral_to_phys_fg(uc_fg,u_fg,1)
  call coarse2fine(vc,vc_fg)
  call spectral_to_phys_fg(vc_fg,v_fg,1)
  call coarse2fine(wc,wc_fg)
  call spectral_to_phys_fg(wc_fg,w_fg,1)
#else
  call spectral_to_phys(uc,u,1)
  call spectral_to_phys(vc,v,1)
  call spectral_to_phys(wc,w,1)
  !This is only a backup copy not used with ex* = 1, gpu code skip this part (strong slow down)
#if openaccflag == 0
    uc_fg=uc
    vc_fg=vc
    wc_fg=wc
    u_fg=u
    v_fg=v
    w_fg=w
#endif
#endif

#if phiflag == 1
! Cahn-Hilliard equation for surfactant
#if psiflag == 1
 allocate(spsi(spxpsi,npsiz,spypsi,2))
 allocate(hpsi(spxpsi,npsiz,spypsi,2))
  ! calculate non-linear terms of surfactant Cahn-Hilliard equation
 call sterm_surf(spsi)

 if (ntime.eq.nstart+1)then
 call euler_psi(spsi,hpsi)
 else
 call adams_bashforth_psi(spsi,hpsi)
 endif

 spsi_o=spsi

deallocate(spsi)
hpsi=hpsi+psic_fg
call calculate_psi(hpsi)
deallocate(hpsi)
call chop_modes_fg(psic_fg)
#endif
! surfactant part executed only iff the phase field is activated
#endif


! Temperature transport equation
#if tempflag == 1
allocate(stheta(spx,nz,spy,2))
allocate(htheta(spx,nz,spy,2))

! calculate non-linear terms of temperature equation
call sterm_temp(stheta)

if (ntime.eq.nstart+1)then
call euler_theta(stheta,htheta)
else
call adams_bashforth_theta(stheta,htheta)
endif

!$acc kernels
stheta_o=stheta
!$acc end kernels

deallocate(stheta)

! assemble history term
call hist_term_temp(htheta)

call calculate_theta(htheta)

deallocate(htheta)

call chop_modes(thetac)
#endif
endif

!particle part
#if partflag == 1
if(rank.ge.leader)then
call lagrangian_tracker
endif
   
call get_velocity

#if twowayc == 1

call get_2WCforces
#endif
#endif


return
end
