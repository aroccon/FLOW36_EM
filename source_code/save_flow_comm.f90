subroutine save_flow_comm(i)

use commondata
use sim_par
use velocity
use phase_field
use stats
use surfactant
use temperature
use polymer

#define machine machineflag
#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define polflag polcompflag 
#define fdump physical_dump_frequency
#define sdump spectral_dump_frequency

integer :: i
character(len=12) :: namevar


#if fdump > 0
  ! save fields at the end of the timestep according to ndump value
  if(mod(i, ndump) == 0) then
    if(rank == 0) write(*,*) 'saving solution in physical space' 
    call spectral_to_phys(uc, u, 0)
    namevar = 'u'
    call write_output(u, i, namevar)
    call spectral_to_phys(vc, v, 0)
    namevar = 'v'
    call write_output(v, i, namevar)
    call spectral_to_phys(wc, w, 0)
    namevar = 'w'
    call write_output(w, i, namevar)
#    if phiflag == 1
      call spectral_to_phys(phic, phi, 0)
      namevar = 'phi'
      call write_output(phi, i, namevar)
#      if psiflag == 1
        call spectral_to_phys_fg(psic_fg, psi_fg, 0)
        namevar = 'psi'
        call write_output_fg(psi_fg, i, namevar)
#      endif
#    endif
#    if tempflag == 1
      call spectral_to_phys(thetac, theta, 0)
      namevar = 'T'
      call write_output(theta, i, namevar)
#    endif
#    if polflag == 1
      call spectral_to_phys(C_xxs, C_xxf, 0)
      call spectral_to_phys(C_xys, C_xyf, 0)
      call spectral_to_phys(C_xzs, C_xzf, 0)
      call spectral_to_phys(C_yys, C_yyf, 0)
      call spectral_to_phys(C_yzs, C_yzf, 0)
      call spectral_to_phys(C_zzs, C_zzf, 0)
      namevar = 'Cxxf'
      call write_output(C_xxf, i, namevar)   
      namevar = 'Cxyf'
      call write_output(C_xyf, i, namevar)    
      namevar = 'Cxzf'
      call write_output(C_xzf, i, namevar)    
      namevar = 'Cyyf'
      call write_output(C_yyf, i, namevar)  
      namevar = 'Cyzf'
      call write_output(C_yzf, i, namevar)    
      namevar = 'Czzf'
      call write_output(C_zzf, i, namevar)
#    endif
  endif
#endif

#if sdump > 0
#  if machine != 2 && machine != 5
    ! save fields at the end of the timestep according to sdump value
    if(mod(i, sdump) == 0) then
      if(rank == 0) write(*,*) 'saving solution in spectral space'
      namevar = 'uc'
      call write_output_spectral(uc, i, namevar)
      namevar = 'vc'
      call write_output_spectral(vc, i, namevar)
      namevar = 'wc'
      call write_output_spectral(wc, i, namevar)
#      if phiflag == 1
        namevar = 'phic'
        call write_output_spectral(phic, i, namevar)
#        if psiflag == 1
          namevar = 'psic'
          call write_output_spectral_fg(psic_fg, i, namevar)
#        endif
#      endif
#      if tempflag == 1
        namevar = 'Tc'
        call write_output_spectral(thetac, i, namevar)
#      endif
#      if polflag == 1
       namevar = 'Cxxs'
       call write_output_spectral(C_xxs, i, namevar)
       namevar = 'Cxys'
       call write_output_spectral(C_xys, i, namevar)
       namevar = 'Cxzs'
       call write_output_spectral(C_xzs, i, namevar)
       namevar = 'Cyys'
       call write_output_spectral(C_yys, i, namevar) 
       namevar = 'Cyzs'
       call write_output_spectral(C_yzs, i, namevar)
       namevar = 'Czzs'
       call write_output_spectral(C_zzs, i, namevar)
#      endif
    endif
#  endif
#endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine save_flow_comm_final(i)

use commondata
use sim_par
use velocity
use phase_field
use stats
use surfactant
use temperature
use polymer

#define machine machineflag
#define phiflag phicompflag
#define psiflag psicompflag
#define tempflag tempcompflag
#define polflag polcompflag 
#define fdump physical_dump_frequency
#define sdump spectral_dump_frequency

integer :: i
character(len=5) :: namevar


#if fdump > 0
  ! Save fields at the end of timestep according to ndump value
  if (mod(i, ndump) .ne. 0) then
    if (rank == 0) write(*,*) 'Writing final fields in physical space' 
    call spectral_to_phys(uc, u, 0)
    namevar = 'u'
    call write_output(u, i, namevar) 
    call spectral_to_phys(vc, v, 0)
    namevar = 'v'
    call write_output(v, i, namevar) 
    call spectral_to_phys(wc, w, 0)
    namevar = 'w'
    call write_output(w, i, namevar)
#    if phiflag == 1
      call spectral_to_phys(phic, phi, 0)
      namevar = 'phi'
      call write_output(phi, i, namevar)
#      if psiflag == 1
        call spectral_to_phys_fg(psic_fg, psi_fg, 0)
        namevar = 'psi'
        call write_output_fg(psi_fg, i, namevar)
#      endif
#    endif
#    if tempflag == 1
      call spectral_to_phys(thetac, theta, 0)
      namevar = 'T'
      call write_output(theta, i, namevar)
#    endif
#    if polflag == 1
     call spectral_to_phys(C_xxs, C_xxf, 0)
     call spectral_to_phys(C_xys, C_xyf, 0)
     call spectral_to_phys(C_xzs, C_xzf, 0)
     call spectral_to_phys(C_yys, C_yyf, 0)
     call spectral_to_phys(C_yzs, C_yzf, 0)
     call spectral_to_phys(C_zzs, C_zzf, 0)
     namevar = 'Cxxf'
     call write_output(C_xxf, i, namevar)   
     namevar = 'Cxyf'
     call write_output(C_xyf, i, namevar)    
     namevar = 'Cxzf'
     call write_output(C_xzf, i, namevar)    
     namevar = 'Cyyf'
     call write_output(C_yyf, i, namevar)  
     namevar = 'Cyzf'
     call write_output(C_yzf, i, namevar)    
     namevar = 'Czzf'
     call write_output(C_zzf, i, namevar)
#    endif
  endif
#endif

#if sdump > 0
#  if machine != 2 && machine != 5
    ! Save fields at the end of timestep according to ndump value
    if (mod(i, sdump) .ne. 0) then
      if (rank == 0) write(*,*) 'Writing final fields in spectral space'   
      namevar = 'uc'
      call write_output_spectral(uc, i, namevar)  
      namevar = 'vc'
      call write_output_spectral(vc, i, namevar)   
      namevar = 'wc'
      call write_output_spectral(wc, i, namevar)
#      if phiflag == 1
        namevar = 'phic'
        call write_output_spectral(phic, i, namevar)       
#        if psiflag == 1
          namevar = 'psic'
          call write_output_spectral_fg(psic_fg, i, namevar)
#        endif
#      endif
#      if tempflag == 1
        namevar = 'Tc'
        call write_output_spectral(thetac, i, namevar)
#      endif
#      if polflag == 1
       namevar = 'Cxxs'
       call write_output_spectral(C_xxs, i, namevar)
       namevar = 'Cxys'
       call write_output_spectral(C_xys, i, namevar)
       namevar = 'Cxzs'
       call write_output_spectral(C_xzs, i, namevar)
       namevar = 'Cyys'
       call write_output_spectral(C_yys, i, namevar) 
       namevar = 'Cyzs'
       call write_output_spectral(C_yzs, i, namevar)
       namevar = 'Czzs'
       call write_output_spectral(C_zzs, i, namevar)
#      endif
    endif
#  endif
#endif

return
end
