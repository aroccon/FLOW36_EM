subroutine yz2xz(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,npzz,ny,2),wa(nx/2+1,npzz,npyy,2)

call yz2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,nx,ny)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine yz2xz_fg(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use dual_grid

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,npzz,npsiy,2),wa(npsix/2+1,npzz,npyy,2)

call yz2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,npsix,npsiy)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! y-z to x-z pencils (direction 1, nycpu ranks with the same z coordinate), gnx x gny grid
! wa: x complete (gnx/2+1 modes), y split -> uc: y complete, x split
! All the blocks are packed with one kernel, exchanged with one mpi_alltoall on cart_comm_dir(1)
! and unpacked with one kernel. Blocks have the padded size ngxx*ngzz*ngyy*2 as in the previous
! pairwise exchange; block j goes to / comes from the rank with y coordinate j.
subroutine yz2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,gnx,gny)

use mpi
use commondata

integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz,gnx,gny
integer :: i,ky,kz,c,j,rx,ry,numel,i0,cnt
double precision :: uc(nsxx,npzz,gny,2),wa(gnx/2+1,npzz,npyy,2)
double precision, allocatable :: bufs(:,:,:,:,:),bufr(:,:,:,:,:)

rx=mod(gnx/2+1,nycpu)
ry=mod(gny,nycpu)

allocate(bufs(ngxx,ngzz,ngyy,2,0:nycpu-1))
allocate(bufr(ngxx,ngzz,ngyy,2,0:nycpu-1))
numel=ngxx*ngyy*ngzz*2

! pack: block j holds the x slab of rank j (offset i0, cnt points), all local y and z
!$acc parallel loop collapse(5) private(i0,cnt)
do j=0,nycpu-1
 do c=1,2
  do ky=1,npyy
   do kz=1,npzz
    do i=1,ngxx
     if(j.lt.rx .or. rx.eq.0)then
      i0=j*ngxx
      cnt=ngxx
     else
      i0=(j-rx)*(ngxx-1)+rx*ngxx
      cnt=ngxx-1
     endif
     if(i.le.cnt) bufs(i,kz,ky,c,j)=wa(i0+i,kz,ky,c)
    enddo
   enddo
  enddo
 enddo
enddo

!CUDA-aware MPI GPU-GPU communicaton (by default hpc-sdk is CUDA-aware)
!$acc host_data use_device(bufs,bufr)
call mpi_alltoall(bufs,numel,mpi_double_precision,bufr,numel,mpi_double_precision,cart_comm_dir(1),ierr)
!$acc end host_data

! unpack: block j holds the local x slab with the y points of rank j (offset i0, cnt points)
!$acc parallel loop collapse(5) private(i0,cnt)
do j=0,nycpu-1
 do c=1,2
  do ky=1,ngyy
   do kz=1,npzz
    do i=1,nsxx
     if(j.lt.ry .or. ry.eq.0)then
      i0=j*ngyy
      cnt=ngyy
     else
      i0=(j-ry)*(ngyy-1)+ry*ngyy
      cnt=ngyy-1
     endif
     if(ky.le.cnt) uc(i,kz,i0+ky,c)=bufr(i,kz,ky,c,j)
    enddo
   enddo
  enddo
 enddo
enddo

deallocate(bufs)
deallocate(bufr)

return
end
