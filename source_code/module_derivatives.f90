module derivatives

  use commondata
  use velocity
  use par_size
  use wavenumber

  implicit none 
  
 double precision, allocatable, dimension(:,:,:) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
 double precision, allocatable, dimension(:,:,:,:) :: sdudx,sdudy,sdudz,sdvdx,sdvdy,sdvdz,sdwdx,sdwdy,sdwdz
 private 
 public dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,velocity_derivatives,destroy_derivatives
      
 contains

subroutine velocity_derivatives
 integer :: indx,indy
 integer :: i,j,k
 
 ! x derivatives
 allocate(sdudx(spx,nz,spy,2))
 allocate(sdvdx(spx,nz,spy,2))
 allocate(sdwdx(spx,nz,spy,2))

 indx=cstart(1)
 !$acc parallel loop collapse(2)
do j=1,spy 
  do i=1,spx
    sdudx(i,:,j,1)=-uc(i,:,j,2)*kx(indx+i)
    sdudx(i,:,j,2)=uc(i,:,j,1)*kx(indx+i)
    sdvdx(i,:,j,1)=-vc(i,:,j,2)*kx(indx+i)
    sdvdx(i,:,j,2)=vc(i,:,j,1)*kx(indx+i)
    sdwdx(i,:,j,1)=-wc(i,:,j,2)*kx(indx+i)
    sdwdx(i,:,j,2)=wc(i,:,j,1)*kx(indx+i)
  enddo
enddo

allocate(dudx(nx,fpz,fpy))
allocate(dvdx(nx,fpz,fpy))
allocate(dwdx(nx,fpz,fpy))

call spectral_to_phys(sdudx,dudx,0)
call spectral_to_phys(sdvdx,dvdx,0)
call spectral_to_phys(sdwdx,dwdx,0)
 
deallocate(sdudx)
deallocate(sdvdx)
deallocate(sdwdx)

! y derivatives
allocate(sdudy(spx,nz,spy,2))
allocate(sdvdy(spx,nz,spy,2))
allocate(sdwdy(spx,nz,spy,2))

indy=cstart(3)
!$acc parallel loop collapse(2)
do j=1,spy
  do i=1,spx
    sdudy(i,:,j,1)=-uc(i,:,j,2)*ky(indy+j)
    sdudy(i,:,j,2)=uc(i,:,j,1)*ky(indy+j)
    sdvdy(i,:,j,1)=-vc(i,:,j,2)*ky(indy+j)
    sdvdy(i,:,j,2)=vc(i,:,j,1)*ky(indy+j)
    sdwdy(i,:,j,1)=-wc(i,:,j,2)*ky(indy+j)
    sdwdy(i,:,j,2)=wc(i,:,j,1)*ky(indy+j)
  enddo
enddo

allocate(dudy(nx,fpz,fpy))
allocate(dvdy(nx,fpz,fpy))
allocate(dwdy(nx,fpz,fpy))

call spectral_to_phys(sdudy,dudy,0)
call spectral_to_phys(sdvdy,dvdy,0)
call spectral_to_phys(sdwdy,dwdy,0)

deallocate(sdudy)
deallocate(sdvdy)
deallocate(sdwdy)

! z derivatives
allocate(sdudz(spx,nz,spy,2))
allocate(sdvdz(spx,nz,spy,2))
allocate(sdwdz(spx,nz,spy,2))

 call dz(uc,sdudz)
 call dz(vc,sdvdz)
 call dz(wc,sdwdz)
 
allocate(dudz(nx,fpz,fpy))
allocate(dvdz(nx,fpz,fpy))
allocate(dwdz(nx,fpz,fpy))

 call spectral_to_phys(sdudz,dudz,1)
 call spectral_to_phys(sdvdz,dvdz,1)
 call spectral_to_phys(sdwdz,dwdz,1)
 
deallocate(sdudz)
deallocate(sdvdz)
deallocate(sdwdz)


 ! if (rank.eq.0) then
 !   print *,'dudx', dudx(1,1,1)
 !   print *,'dwdx', dvdx(1,1,1)
 !   print *,'dwdx', dwdx(1,1,1)
 !   print *,'dudy', dudy(1,1,1)
 !   print *,'dwdy', dvdy(1,1,1)
 !   print *,'dwdy', dwdy(1,1,1)
 !   print *,'dudz', dudz(1,1,1)
 !   print *,'dwdz', dvdz(1,1,1)
 !   print *,'dwdz', dwdz(1,1,1)
 ! endif

return
      
end subroutine velocity_derivatives
    
subroutine destroy_derivatives

        
  deallocate(dudx)
  deallocate(dudy)
  deallocate(dudz)
  deallocate(dvdx)
  deallocate(dvdy)
  deallocate(dvdz)
  deallocate(dwdx)
  deallocate(dwdy)
  deallocate(dwdz)

return
end subroutine destroy_derivatives
end module derivatives