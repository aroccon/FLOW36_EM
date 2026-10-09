subroutine yz2xz(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,npzz,ny,2),wa(nx/2+1,npzz,npyy,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nycpu,8))
call yz2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,nx,ny,a2a_send,a2a_recv)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine yz2xz_fg(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use dual_grid
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,npzz,npsiy,2),wa(npsix/2+1,npzz,npyy,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nycpu,8))
call yz2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,npsix,npsiy,a2a_send,a2a_recv)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! y-z to x-z pencils (direction 1, nycpu ranks with the same z coordinate), gnx x gny grid
! wa: x complete (gnx/2+1 modes), y split -> uc: y complete, x split
! Buffers bufs/bufr come from module a2a_buffers (allocated once).
! All the blocks are packed with one kernel, exchanged with non-blocking messages on cart_comm_dir(1)
! and unpacked with one kernel. Blocks have the padded size ngxx*ngzz*ngyy*2 as in the previous
! pairwise exchange; block j goes to / comes from the rank with y coordinate j.
subroutine yz2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,gnx,gny,bufs,bufr)

use mpi
use commondata

integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz,gnx,gny
integer :: me,nreq,req(2*nycpu)
integer :: i,ky,kz,c,j,rx,ry,numel,i0,cnt
double precision :: uc(nsxx,npzz,gny,2),wa(gnx/2+1,npzz,npyy,2)
double precision :: bufs(ngxx,ngzz,ngyy,2,0:nycpu-1),bufr(ngxx,ngzz,ngyy,2,0:nycpu-1)

rx=mod(gnx/2+1,nycpu)
ry=mod(gny,nycpu)

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

! exchange with all the other ranks of cart_comm_dir(1) at once: all the receives and sends
! are posted together (no sequential rounds); the own block (j=me) is not sent, the unpack
! kernel takes it directly from bufs. Point-to-point instead of mpi_alltoall: the collective
! would copy the own block (and possibly stage blocks) on the host with managed memory.
call mpi_comm_rank(cart_comm_dir(1),me,ierr)
nreq=0
!CUDA-aware MPI GPU-GPU communicaton (by default hpc-sdk is CUDA-aware)
!$acc host_data use_device(bufs,bufr)
do j=0,nycpu-1
 if(j.ne.me)then
  nreq=nreq+1
  call mpi_irecv(bufr(1,1,1,1,j),numel,mpi_double_precision,j,0,cart_comm_dir(1),req(nreq),ierr)
 endif
enddo
do j=0,nycpu-1
 if(j.ne.me)then
  nreq=nreq+1
  call mpi_isend(bufs(1,1,1,1,j),numel,mpi_double_precision,j,0,cart_comm_dir(1),req(nreq),ierr)
 endif
enddo
call mpi_waitall(nreq,req,mpi_statuses_ignore,ierr)
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
     if(ky.le.cnt)then
      if(j.eq.me)then
       uc(i,kz,i0+ky,c)=bufs(i,kz,ky,c,j)
      else
       uc(i,kz,i0+ky,c)=bufr(i,kz,ky,c,j)
      endif
     endif
    enddo
   enddo
  enddo
 enddo
enddo

return
end
