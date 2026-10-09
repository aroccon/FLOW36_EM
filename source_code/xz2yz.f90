subroutine xz2yz(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nx/2+1,npzz,npyy,2),wa(nsxx,npzz,ny,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nycpu,8))
call xz2yz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,nx,ny,a2a_send,a2a_recv)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine xz2yz_fg(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use dual_grid
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(npsix/2+1,npzz,npyy,2),wa(nsxx,npzz,npsiy,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nycpu,8))
call xz2yz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,npsix,npsiy,a2a_send,a2a_recv)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! x-z to y-z pencils (direction 1, nycpu ranks with the same z coordinate), gnx x gny grid
! wa: y complete, x split -> uc: x complete (gnx/2+1 modes), y split
! Buffers bufs/bufr come from module a2a_buffers (allocated once).
! All the blocks are packed with one kernel, exchanged with one mpi_alltoall on cart_comm_dir(1)
! and unpacked with one kernel. Blocks have the padded size ngxx*ngzz*ngyy*2 as in the previous
! pairwise exchange; block j goes to / comes from the rank with y coordinate j.
subroutine xz2yz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,gnx,gny,bufs,bufr)

use mpi
use commondata

integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz,gnx,gny
integer :: i,ky,kz,c,j,rx,ry,numel,i0,cnt
double precision :: uc(gnx/2+1,npzz,npyy,2),wa(nsxx,npzz,gny,2)
double precision :: bufs(ngxx,ngzz,ngyy,2,0:nycpu-1),bufr(ngxx,ngzz,ngyy,2,0:nycpu-1)

ry=mod(gny,nycpu)
rx=mod(gnx/2+1,nycpu)

numel=ngxx*ngyy*ngzz*2

! pack: block j holds the y slab of rank j (offset i0, cnt points), all local x and z
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
     if(ky.le.cnt) bufs(i,kz,ky,c,j)=wa(i,kz,i0+ky,c)
    enddo
   enddo
  enddo
 enddo
enddo

!CUDA-aware MPI GPU-GPU communicaton (by default hpc-sdk is CUDA-aware)
!$acc host_data use_device(bufs,bufr)
call mpi_alltoall(bufs,numel,mpi_double_precision,bufr,numel,mpi_double_precision,cart_comm_dir(1),ierr)
!$acc end host_data

! unpack: block j holds the x slab of rank j (offset i0, cnt points) with the local y and z
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
     if(i.le.cnt) uc(i0+i,kz,ky,c)=bufr(i,kz,ky,c,j)
    enddo
   enddo
  enddo
 enddo
enddo

return
end
