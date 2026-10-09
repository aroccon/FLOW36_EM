subroutine save_interface(ntime)

use commondata
use grid
use par_size
use mpi
use velocity
use phase_field
use sim_par
use polymer
use stats


implicit none

integer, intent(in) :: ntime
integer :: i,j,k,k_local
character(len=256) :: filename

double precision, allocatable :: xi(:,:), xi_loc(:,:)
integer, allocatable :: k_interface(:,:)

double precision :: stencil

double precision, allocatable :: vel_int(:,:,:), vel_int_loc(:,:,:)
double precision, allocatable :: C_int(:,:,:), C_int_loc(:,:,:)

!--------------------------------------------------
! Convert phase field to physical space
!--------------------------------------------------
call spectral_to_phys(phic,phi,0)

! Allocate 2D array for integral
allocate(xi_loc(nx,ny))
xi_loc=0.0d0

! Compute z-integral of phi for local domain
do j=1,fpy
  do i=1,nx
    do k=2,fpz-1
     xi_loc(i,fstart(3)+j) = xi_loc(i,fstart(3)+j) + &
       &0.5d0*(0.5d0*(phi(i,k,j)+1.0d0))* &
       &(z(fstart(2)+k+1)-z(fstart(2)+k-1))
    enddo
  enddo
enddo

! handling boundaries (bottom)
k=1
if(fstart(2).eq.0)then
  do j=1,fpy
    do i=1,nx
        xi_loc(i,fstart(3)+j) = xi_loc(i,fstart(3)+j)+ &
         &0.5d0*(0.5d0*(phi(i,k,j)+1.0d0))* &
         &(z(fstart(2)+k+1)-z(fstart(2)+k))
    enddo
  enddo
else
  do j=1,fpy
    do i=1,nx
        xi_loc(i,fstart(3)+j) = xi_loc(i,fstart(3)+j)+ &
         &0.5d0*(0.5d0*(phi(i,k,j)+1.0d0))* &
         &(z(fstart(2)+k+1)-z(fstart(2)+k-1))
    enddo
  enddo
endif
  
! handling boundaries (top)
k=fpz
if(fstart(2)+fpz.eq.nz)then 
  do j=1,fpy
    do i=1,nx
     xi_loc(i,fstart(3)+j)=xi_loc(i,fstart(3)+j) + &
       0.5d0*(0.5d0*(phi(i,k,j)+1.0d0))*&
       (z(fstart(2)+k)-z(fstart(2)+k-1))
    enddo
  enddo
else
  do j=1,fpy
    do i=1,nx
     xi_loc(i,fstart(3)+j)=xi_loc(i,fstart(3)+j) + &
       0.5d0*(0.5d0*(phi(i,k,j)+1.0d0))*&
       (z(fstart(2)+k+1)-z(fstart(2)+k-1))
    enddo
  enddo
endif


! Reduce
allocate(xi(nx,ny))
call mpi_allreduce(xi_loc,xi,nx*ny,mpi_double_precision,mpi_sum,flow_comm,ierr)
deallocate(xi_loc)

xi = xi +1.d0

! Save to file
if(rank.eq.0)then
  write(filename,'(a,"/xi_",i8.8,".dat")') trim(folder),ntime
  open(80,file=filename,form='unformatted',status='new')
  write(80) xi
  close(80)
endif

!--------------------------------------------------
! Locate interface index k
!--------------------------------------------------
allocate(k_interface(nx,ny))
k_interface = 0

do j = 1, ny
  do i = 1, nx
    do k = 2, nz
      if ((z(k-1)-xi(i,j))*(z(k)-xi(i,j)) <= 0.0d0) then
        k_interface(i,j) = k
        exit
      endif
    enddo
  enddo
enddo

!--------------------------------------------------
! Allocate interpolated quantities
!--------------------------------------------------
allocate(vel_int_loc(nx,ny,3))
allocate(vel_int(nx,ny,3))

allocate(C_int_loc(nx,ny,6))
allocate(C_int(nx,ny,6))

vel_int_loc = 0.0d0
C_int_loc   = 0.0d0

!--------------------------------------------------
! Interpolate at interface
!--------------------------------------------------
do j = 1, fpy
  do i = 1, nx

    k_local = k_interface(i,fstart(3)+j) - fstart(2)

    ! -------------------------------
    ! Contribution from k_local
    ! -------------------------------
    if (k_local > 0 .and. k_local <= fpz) then

      stencil = (xi(i,fstart(3)+j) - z(k_interface(i,fstart(3)+j)-1)) / &
                (z(k_interface(i,fstart(3)+j)) - z(k_interface(i,fstart(3)+j)-1))

      vel_int_loc(i,fstart(3)+j,1) = vel_int_loc(i,fstart(3)+j,1) + &
                                     (1.d0 - stencil) * u(i,k_local,j)
      vel_int_loc(i,fstart(3)+j,2) = vel_int_loc(i,fstart(3)+j,2) + &
                                     (1.d0 - stencil) * v(i,k_local,j)
      vel_int_loc(i,fstart(3)+j,3) = vel_int_loc(i,fstart(3)+j,3) + &
                                     (1.d0 - stencil) * w(i,k_local,j)

      C_int_loc(i,fstart(3)+j,1) = C_int_loc(i,fstart(3)+j,1) + &
                                   (1.d0 - stencil) * C_xxf(i,k_local,j)
      C_int_loc(i,fstart(3)+j,2) = C_int_loc(i,fstart(3)+j,2) + &
                                   (1.d0 - stencil) * C_yyf(i,k_local,j)
      C_int_loc(i,fstart(3)+j,3) = C_int_loc(i,fstart(3)+j,3) + &
                                   (1.d0 - stencil) * C_zzf(i,k_local,j)
      C_int_loc(i,fstart(3)+j,4) = C_int_loc(i,fstart(3)+j,4) + &
                                   (1.d0 - stencil) * C_xyf(i,k_local,j)
      C_int_loc(i,fstart(3)+j,5) = C_int_loc(i,fstart(3)+j,5) + &
                                   (1.d0 - stencil) * C_xzf(i,k_local,j)
      C_int_loc(i,fstart(3)+j,6) = C_int_loc(i,fstart(3)+j,6) + &
                                   (1.d0 - stencil) * C_yzf(i,k_local,j)
    endif

    ! -------------------------------
    ! Contribution from k_local-1
    ! -------------------------------
    if (k_local-1 > 0 .and. k_local-1 <= fpz) then

      stencil = 1.d0 - (xi(i,fstart(3)+j) - z(k_interface(i,fstart(3)+j)-1)) / &
                (z(k_interface(i,fstart(3)+j)) - z(k_interface(i,fstart(3)+j)-1))

      vel_int_loc(i,fstart(3)+j,1) = vel_int_loc(i,fstart(3)+j,1) + &
                                     (1.d0 - stencil) * u(i,k_local-1,j)
      vel_int_loc(i,fstart(3)+j,2) = vel_int_loc(i,fstart(3)+j,2) + &
                                     (1.d0 - stencil) * v(i,k_local-1,j)
      vel_int_loc(i,fstart(3)+j,3) = vel_int_loc(i,fstart(3)+j,3) + &
                                     (1.d0 - stencil) * w(i,k_local-1,j)

      C_int_loc(i,fstart(3)+j,1) = C_int_loc(i,fstart(3)+j,1) + &
                                   (1.d0 - stencil) * C_xxf(i,k_local-1,j)
      C_int_loc(i,fstart(3)+j,2) = C_int_loc(i,fstart(3)+j,2) + &
                                   (1.d0 - stencil) * C_yyf(i,k_local-1,j)
      C_int_loc(i,fstart(3)+j,3) = C_int_loc(i,fstart(3)+j,3) + &
                                   (1.d0 - stencil) * C_zzf(i,k_local-1,j)
      C_int_loc(i,fstart(3)+j,4) = C_int_loc(i,fstart(3)+j,4) + &
                                   (1.d0 - stencil) * C_xyf(i,k_local-1,j)
      C_int_loc(i,fstart(3)+j,5) = C_int_loc(i,fstart(3)+j,5) + &
                                   (1.d0 - stencil) * C_xzf(i,k_local-1,j)
      C_int_loc(i,fstart(3)+j,6) = C_int_loc(i,fstart(3)+j,6) + &
                                   (1.d0 - stencil) * C_yzf(i,k_local-1,j)
    endif

  enddo
enddo

!--------------------------------------------------
! MPI reduction
!--------------------------------------------------
call mpi_allreduce(vel_int_loc,vel_int,nx*ny*3,mpi_double_precision,mpi_sum,flow_comm,ierr)
call mpi_allreduce(C_int_loc,C_int,nx*ny*6,mpi_double_precision,mpi_sum,flow_comm,ierr)

!--------------------------------------------------
! Save interface fields
!--------------------------------------------------
if(rank == 0) then
  write(filename,'(a,"/interface_fields_",i8.8,".dat")') trim(folder), ntime
  open(20,file=filename,form='unformatted',status='replace')
  write(20) vel_int
  write(20) C_int
  close(20)
endif

!--------------------------------------------------
! Cleanup
!--------------------------------------------------
deallocate(xi, k_interface)
deallocate(vel_int, vel_int_loc)
deallocate(C_int, C_int_loc)

return
end subroutine save_interface