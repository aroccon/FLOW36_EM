subroutine initialize_pol

  use commondata
  use par_size
  use grid
  use polymer
  use sterms
  use sim_par
  use velocity_old
  use polymer
  use derivatives
  use assemble
  use phase_field
!  use Eigen3x3

#define fenepflag fenepcompflag
#define solvpolflag polysolvflag 
#define carflag carnatureflag
#define phiflag phicompflag

  
  logical :: checkC_xxf,checkC_xyf,checkC_xzf,checkC_yyf,checkC_yzf,checkC_zzf
  logical :: checkC_xxs,checkC_xys,checkC_xzs,checkC_yys,checkC_yzs,checkC_zzs

  double precision, allocatable, dimension(:,:,:) :: Kappa,jota,omicron
  double precision :: CTOT(3,3),Eigenvectors(3,3),Eigenvalues(3)
  double precision :: phif



  character(len=8) :: time

  allocate(C_xxf(nx,fpz,fpy))
  allocate(C_xyf(nx,fpz,fpy))
  allocate(C_xzf(nx,fpz,fpy))
  allocate(C_yyf(nx,fpz,fpy))
  allocate(C_yzf(nx,fpz,fpy))
  allocate(C_zzf(nx,fpz,fpy))
  allocate(C_xxs(spx,nz,spy,2))
  allocate(C_xys(spx,nz,spy,2))
  allocate(C_xzs(spx,nz,spy,2))
  allocate(C_yys(spx,nz,spy,2))
  allocate(C_yzs(spx,nz,spy,2))
  allocate(C_zzs(spx,nz,spy,2))    
  allocate(sconfxx_o(spx,nz,spy,2))
  allocate(sconfxy_o(spx,nz,spy,2))
  allocate(sconfxz_o(spx,nz,spy,2))
  allocate(sconfyy_o(spx,nz,spy,2))
  allocate(sconfyz_o(spx,nz,spy,2))
  allocate(sconfzz_o(spx,nz,spy,2))    
  allocate(pols1_n(spx,nz,spy,2))    
  allocate(pols2_n(spx,nz,spy,2))    
  allocate(pols3_n(spx,nz,spy,2))    
#if (fenepflag == 1 || fenepflag == 2)
  allocate(zitaVel(nx,fpz,fpy))
  allocate(rkkn(nx,fpz,fpy))
#endif
#if fenepflag == 3
  allocate(modtauD(nx,fpz,fpy))   
#endif


#if solvpolflag == 1
  allocate(diag11(nx,fpz,fpy))
  allocate(diag22(nx,fpz,fpy))
  allocate(diag33(nx,fpz,fpy))
  allocate(erre11(nx,fpz,fpy))
  allocate(erre12(nx,fpz,fpy))
  allocate(erre13(nx,fpz,fpy))
  allocate(erre21(nx,fpz,fpy))
  allocate(erre22(nx,fpz,fpy))
  allocate(erre23(nx,fpz,fpy))
  allocate(erre31(nx,fpz,fpy))
  allocate(erre32(nx,fpz,fpy))
  allocate(erre33(nx,fpz,fpy))
  diag11=1.d0
  diag22=1.d0 
  diag33=1.d0
  erre11=1.d0
  erre12=0.d0 
  erre13=0.d0
  erre21=0.d0
  erre22=1.d0 
  erre23=0.d0
  erre31=0.d0
  erre32=0.d0
  erre33=1.d0
#endif

  if(in_cond_pol.eq.0)then
      if(rank.eq.0)write(*,*) 'Initializing polymer coiled state C=I'
      C_xxf=1.0d0
      C_xyf=0.0d0
      C_xzf=0.0d0
      C_yyf=1.0d0 
      C_yzf=0.0d0
      C_zzf=1.0d0
      call phys_to_spectral(C_xxf,C_xxs,0)
      call phys_to_spectral(C_xyf,C_xys,0)
      call phys_to_spectral(C_xzf,C_xzs,0)
      call phys_to_spectral(C_yyf,C_yys,0)
      call phys_to_spectral(C_yzf,C_yzs,0)
      call phys_to_spectral(C_zzf,C_zzs,0)

#if fenepflag == 3
if(rank.eq.0)write(*,*) ' init tau Deviatorico'
    modtauD=0.0d0
    shear_stress=C_yzf(1,1,1)/Wi
    normal=(C_yyf(1,1,1)-C_zzf(1,1,1))/Wi
#endif

  elseif(in_cond_pol.eq.1)then
          if(rank.eq.0)write(*,*) 'Initializing conformation tensor fields from data file (parallel read)'
          write(time,'(I8.8)') nt_restart
          if(restart.eq.1)then
            inquire(file=trim(folder)//'/Cxxf_'//time//'.dat',exist=checkC_xxf)
            inquire(file=trim(folder)//'/Cxyf_'//time//'.dat',exist=checkC_xyf)
            inquire(file=trim(folder)//'/Cxzf_'//time//'.dat',exist=checkC_xzf)
            inquire(file=trim(folder)//'/Cyyf_'//time//'.dat',exist=checkC_yyf)
            inquire(file=trim(folder)//'/Cyzf_'//time//'.dat',exist=checkC_yzf)
            inquire(file=trim(folder)//'/Czzf_'//time//'.dat',exist=checkC_zzf)  
            inquire(file=trim(folder)//'/Cxxs_'//time//'.dat',exist=checkC_xxs)
            inquire(file=trim(folder)//'/Cxys_'//time//'.dat',exist=checkC_xys)
            inquire(file=trim(folder)//'/Cxzs_'//time//'.dat',exist=checkC_xzs)
            inquire(file=trim(folder)//'/Cyys_'//time//'.dat',exist=checkC_yys)
            inquire(file=trim(folder)//'/Cyzs_'//time//'.dat',exist=checkC_yzs)
            inquire(file=trim(folder)//'/Czzs_'//time//'.dat',exist=checkC_zzs)          
          else
            checkC_xxf=.true.
            checkC_xyf=.true. 
            checkC_xzf=.true. 
            checkC_yyf=.true. 
            checkC_yzf=.true. 
            checkC_zzf=.true.
          endif
          if((checkC_xxf.eqv..true.).and.(checkC_xyf.eqv..true.).and.(checkC_xzf.eqv..true.).and. &
           & (checkC_yyf.eqv..true.).and.(checkC_yzf.eqv..true.).and.(checkC_zzf.eqv..true.))then
              call read_fields(C_xxf,nt_restart,'Cxxf ',restart)
              call read_fields(C_xyf,nt_restart,'Cxyf ',restart)
              call read_fields(C_xzf,nt_restart,'Cxzf ',restart)
              call read_fields(C_yyf,nt_restart,'Cyyf ',restart)
              call read_fields(C_yzf,nt_restart,'Cyzf ',restart)
              call read_fields(C_zzf,nt_restart,'Czzf ',restart)
            if(rank.eq.0) write(*,'(1x,a,a,a)') 'Conformation tensor MASK '
           !$acc parallel loop collapse(3)
           do j=1,fpy
             do k=1,fpz
               do i=1,nx
                 phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
                 C_xxf(i,k,j)=1.0d0*0.5d0*(1.0d0+phif)+C_xxf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_xyf(i,k,j)=0.0d0*0.5d0*(1.0d0+phif)+C_xyf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_xzf(i,k,j)=0.0d0*0.5d0*(1.0d0+phif)+C_xzf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_yyf(i,k,j)=1.0d0*0.5d0*(1.0d0+phif)+C_yyf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_yzf(i,k,j)=0.0d0*0.5d0*(1.0d0+phif)+C_yzf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_zzf(i,k,j)=1.0d0*0.5d0*(1.0d0+phif)+C_zzf(i,k,j)*0.5d0*(1.0d0-phif)
               enddo
             enddo
           enddo


              ! transform physical variable to spectral space
            call phys_to_spectral(C_xxf,C_xxs,0)
            call phys_to_spectral(C_xyf,C_xys,0)
            call phys_to_spectral(C_xzf,C_xzs,0)
            call phys_to_spectral(C_yyf,C_yys,0)
            call phys_to_spectral(C_yzf,C_yzs,0)
            call phys_to_spectral(C_zzf,C_zzs,0)
#if solvpolflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz 
     do i=1,nx 
        CTOT = reshape([C_xxf(i,k,j), C_xyf(i,k,j), C_xzf(i,k,j), &
                        C_xyf(i,k,j), C_yyf(i,k,j), C_yzf(i,k,j), &
                        C_xzf(i,k,j), C_yzf(i,k,j), C_zzf(i,k,j)], [3,3])
    
        call Diagonalization(CTOT,Eigenvectors,Eigenvalues) 

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
#endif
          elseif((checkC_xxs.eqv..true.).and.(checkC_xys.eqv..true.).and.(checkC_xzs.eqv..true.).and. &
               & (checkC_yys.eqv..true.).and.(checkC_yzs.eqv..true.).and.(checkC_zzs.eqv..true.))then
              call read_fields_s(C_xxs,nt_restart,'Cxxs ',restart)
              call read_fields_s(C_xys,nt_restart,'Cxys ',restart)
              call read_fields_s(C_xzs,nt_restart,'Cxzs ',restart)
              call read_fields_s(C_yys,nt_restart,'Cyys ',restart)
              call read_fields_s(C_yzs,nt_restart,'Cyzs ',restart)
              call read_fields_s(C_zzs,nt_restart,'Czzs ',restart)
              ! transform to physical space
              call spectral_to_phys(C_xxs,C_xxf,0)
              call spectral_to_phys(C_xys,C_xyf,0)
              call spectral_to_phys(C_xzs,C_xzf,0)
              call spectral_to_phys(C_yys,C_yyf,0)
              call spectral_to_phys(C_yzs,C_yzf,0)
              call spectral_to_phys(C_zzs,C_zzf,0)

              call assemble_eigen
          else
            if(rank.eq.0) write(*,'(1x,a,a,a)') 'Missing conformation tensor field input file ',time,' , stopping simulation'
            call exit(0)
          endif
  elseif(in_cond_pol.eq.2)then
          if(rank.eq.0) write(*,*) 'Initializing conformation tensor fields from data file (serial read)'
          call read_fields_serial(C_xxf,nt_restart,'Cxxf ',restart)
          call read_fields_serial(C_xyf,nt_restart,'Cxyf ',restart)
          call read_fields_serial(C_xzf,nt_restart,'Cxzf ',restart)
          call read_fields_serial(C_yyf,nt_restart,'Cyyf ',restart)
          call read_fields_serial(C_yzf,nt_restart,'Cyzf ',restart)
          call read_fields_serial(C_zzf,nt_restart,'Czzf ',restart)
          ! transform physical variable to spectral space
          call phys_to_spectral(C_xxf,C_xxs,0)
          call phys_to_spectral(C_xyf,C_xys,0)
          call phys_to_spectral(C_xzf,C_xzs,0)
          call phys_to_spectral(C_yyf,C_yys,0)
          call phys_to_spectral(C_yzf,C_yzs,0)
          call phys_to_spectral(C_zzf,C_zzs,0)
#if solvpolflag == 1
          !$acc parallel loop collapse(3)
          do j=1,fpy
            do k=1,fpz 
               do i=1,nx 
                  CTOT = reshape([C_xxf(i,k,j), C_xyf(i,k,j), C_xzf(i,k,j), &
                                  C_xyf(i,k,j), C_yyf(i,k,j), C_yzf(i,k,j), &
                                  C_xzf(i,k,j), C_yzf(i,k,j), C_zzf(i,k,j)], [3,3])
              
                  call Diagonalization(CTOT,Eigenvectors,Eigenvalues) 
          
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
#endif
   elseif(in_cond_pol.eq.3)then
            call velocity_derivatives
       if ((fenepflag.eq.1).or.(fenepflag.eq.2)) then
           allocate(Kappa(nx,fpz,fpy))
           allocate(jota(nx,fpz,fpy))
           allocate(omicron(nx,fpz,fpy))

           rkkn=0.0d0
           Kappa = (Wi*2.0d0**0.5d0)*dudz/L
           jota = asinh((3.0d0*(3.0d0**0.5d0)/2.0d0)*Kappa)
           omicron = (3.0d0**0.5d0)*Kappa/(2.0d0*sinh(jota/3.0d0))
 
           C_xxf = (1.0d0+2.0d0*(Wi*dudz/omicron)**2.0d0)/omicron
           C_xyf = 0.0d0
           C_xzf = Wi*dudz/omicron**2.0d0
           C_yyf = 1.0d0/omicron
           C_yzf = 0.0d0
           C_zzf = 1.0d0/omicron
 
           deallocate(Kappa)
           deallocate(jota)
           deallocate(omicron)
#if (phiflag == 1 && carflag == 0)
           !$acc parallel loop collapse(3)
           do j=1,fpy
             do k=1,fpz
               do i=1,nx
                 phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
                 C_xxf(i,k,j)=1.0d0*0.5d0*(1.0d0+phif)+C_xxf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_xyf(i,k,j)=0.0d0*0.5d0*(1.0d0+phif)+C_xyf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_xzf(i,k,j)=0.0d0*0.5d0*(1.0d0+phif)+C_xzf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_yyf(i,k,j)=1.0d0*0.5d0*(1.0d0+phif)+C_yyf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_yzf(i,k,j)=0.0d0*0.5d0*(1.0d0+phif)+C_yzf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_zzf(i,k,j)=1.0d0*0.5d0*(1.0d0+phif)+C_zzf(i,k,j)*0.5d0*(1.0d0-phif)
               enddo
             enddo
           enddo
#elif (phiflag == 1 && carflag == 1)
           !$acc parallel loop collapse(3)
           do j=1,fpy
             do k=1,fpz
               do i=1,nx
                 phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
                 a6f(i,k,j)=a6f(i,k,j)*0.5d0*(1.0d0+phif)
                 C_xxf(i,k,j)=1.0d0*0.5d0*(1.0d0-phif)+C_xxf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_xyf(i,k,j)=0.0d0*0.5d0*(1.0d0-phif)+C_xyf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_xzf(i,k,j)=0.0d0*0.5d0*(1.0d0-phif)+C_xzf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_yyf(i,k,j)=1.0d0*0.5d0*(1.0d0-phif)+C_yyf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_yzf(i,k,j)=0.0d0*0.5d0*(1.0d0-phif)+C_yzf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_zzf(i,k,j)=1.0d0*0.5d0*(1.0d0-phif)+C_zzf(i,k,j)*0.5d0*(1.0d0+phif)
               enddo
             enddo
           enddo
#endif
       else
           C_xxf = (1.0d0+2.0d0*(Wi*dudz)**2.0d0)
           C_xyf = 0.0d0
           C_xzf = Wi*dudz
           C_yyf = 1.0d0
           C_yzf = 0.0d0
           C_zzf = 1.0d0
#if (phiflag == 1 && carflag == 0)
           !$acc parallel loop collapse(3)
           do j=1,fpy
             do k=1,fpz
               do i=1,nx
                 phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
                 C_xxf(i,k,j)=1.0d0*0.5d0*(1.0d0+phif)+C_xxf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_xyf(i,k,j)=0.0d0*0.5d0*(1.0d0+phif)+C_xyf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_xzf(i,k,j)=0.0d0*0.5d0*(1.0d0+phif)+C_xzf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_yyf(i,k,j)=1.0d0*0.5d0*(1.0d0+phif)+C_yyf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_yzf(i,k,j)=0.0d0*0.5d0*(1.0d0+phif)+C_yzf(i,k,j)*0.5d0*(1.0d0-phif)
                 C_zzf(i,k,j)=1.0d0*0.5d0*(1.0d0+phif)+C_zzf(i,k,j)*0.5d0*(1.0d0-phif)
               enddo
             enddo
           enddo
#elif (phiflag == 1 && carflag == 1)
           !$acc parallel loop collapse(3)
           do j=1,fpy
             do k=1,fpz
               do i=1,nx
                 phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
                 a6f(i,k,j)=a6f(i,k,j)*0.5d0*(1.0d0+phif)
                 C_xxf(i,k,j)=1.0d0*0.5d0*(1.0d0-phif)+C_xxf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_xyf(i,k,j)=0.0d0*0.5d0*(1.0d0-phif)+C_xyf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_xzf(i,k,j)=0.0d0*0.5d0*(1.0d0-phif)+C_xzf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_yyf(i,k,j)=1.0d0*0.5d0*(1.0d0-phif)+C_yyf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_yzf(i,k,j)=0.0d0*0.5d0*(1.0d0-phif)+C_yzf(i,k,j)*0.5d0*(1.0d0+phif)
                 C_zzf(i,k,j)=1.0d0*0.5d0*(1.0d0-phif)+C_zzf(i,k,j)*0.5d0*(1.0d0+phif)
               enddo
             enddo
           enddo
#endif
       endif
 
       call destroy_derivatives

#if solvpolflag == 1
       !$acc parallel loop collapse(3)
       do j=1,fpy
         do k=1,fpz 
            do i=1,nx 
               CTOT = reshape([C_xxf(i,k,j), C_xyf(i,k,j), C_xzf(i,k,j), &
                               C_xyf(i,k,j), C_yyf(i,k,j), C_yzf(i,k,j), &
                               C_xzf(i,k,j), C_yzf(i,k,j), C_zzf(i,k,j)], [3,3])
           
               call Diagonalization(CTOT,Eigenvectors,Eigenvalues) 
       
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
#endif
 
       call phys_to_spectral(C_xxf,C_xxs,0)
       call phys_to_spectral(C_xyf,C_xys,0)
       call phys_to_spectral(C_xzf,C_xzs,0)
       call phys_to_spectral(C_yyf,C_yys,0)
       call phys_to_spectral(C_yzf,C_yzs,0)
       call phys_to_spectral(C_zzf,C_zzs,0)
 
  endif

  return
end
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine destroy_pol
  
  use polymer
  use velocity_old
  use sterms
  use assemble
      
  deallocate(C_xxf)
  deallocate(C_xyf)
  deallocate(C_xzf)
  deallocate(C_yyf)
  deallocate(C_yzf)
  deallocate(C_zzf)
  deallocate(C_xxs)
  deallocate(C_xys)
  deallocate(C_xzs)
  deallocate(C_yys)
  deallocate(C_yzs)
  deallocate(C_zzs)  
  deallocate(sconfxx_o)  
  deallocate(sconfxy_o)  
  deallocate(sconfxz_o)  
  deallocate(sconfyy_o)  
  deallocate(sconfyz_o)  
  deallocate(sconfzz_o) 
  deallocate(pols1_n)    
  deallocate(pols2_n)    
  deallocate(pols3_n)   
#if (fenepflag == 1 || fenepflag == 2)
  deallocate(zitaVel)
  deallocate(rkkn)
#endif
#if fenepflag == 3
  deallocate(modtauD)
#endif

#if solvpolflag == 1
  deallocate(diag11)
  deallocate(diag22)
  deallocate(diag33)
  deallocate(erre11)
  deallocate(erre12)
  deallocate(erre13)
  deallocate(erre21)
  deallocate(erre22)
  deallocate(erre23)
  deallocate(erre31)
  deallocate(erre32)
  deallocate(erre33)
#endif
  
  return
  end
