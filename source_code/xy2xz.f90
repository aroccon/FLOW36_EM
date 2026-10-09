subroutine xy2xz(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,npzz,ny,2),wa(nsxx,nz,npyy,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nzcpu,8))
call xy2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,ny,nz,a2a_send,a2a_recv)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine xy2xz_fg(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use dual_grid
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,npzz,npsiy,2),wa(nsxx,npsiz,npyy,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nzcpu,8))
call xy2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,npsiy,npsiz,a2a_send,a2a_recv)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! x-y to x-z pencils (direction 0, nzcpu ranks with the same y coordinate), gny x gnz grid
! wa: z complete, y split -> uc: y complete, z split
! Buffers bufs/bufr come from module a2a_buffers (allocated once).
! All the blocks are packed with one kernel, exchanged with non-blocking messages on cart_comm_dir(0)
! and unpacked with one kernel. Blocks have the padded size ngxx*ngzz*ngyy*2 as in the previous
! pairwise exchange; block j goes to / comes from the rank with z coordinate j.
subroutine xy2xz_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,gny,gnz,bufs,bufr)

use mpi
use commondata

integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz,gny,gnz
integer :: me,nreq,req(2*nzcpu)
integer :: i,ky,kz,c,j,ry,rz,numel,i0,cnt
double precision :: uc(nsxx,npzz,gny,2),wa(nsxx,gnz,npyy,2)
double precision :: bufs(ngxx,ngzz,ngyy,2,0:nzcpu-1),bufr(ngxx,ngzz,ngyy,2,0:nzcpu-1)

ry=mod(gny,nzcpu)
rz=mod(gnz,nzcpu)

numel=ngxx*ngyy*ngzz*2

! pack: block j holds the z slab of rank j (offset i0, cnt points), all local y
!$acc parallel loop collapse(5) private(i0,cnt)
do j=0,nzcpu-1
 do c=1,2
  do ky=1,npyy
   do kz=1,ngzz
    do i=1,nsxx
     if(j.lt.rz .or. rz.eq.0)then
      i0=j*ngzz
      cnt=ngzz
     else
      i0=(j-rz)*(ngzz-1)+rz*ngzz
      cnt=ngzz-1
     endif
     if(kz.le.cnt) bufs(i,kz,ky,c,j)=wa(i,i0+kz,ky,c)
    enddo
   enddo
  enddo
 enddo
enddo

! exchange with all the other ranks of cart_comm_dir(0) at once: all the receives and sends
! are posted together (no sequential rounds); the own block (j=me) is not sent, the unpack
! kernel takes it directly from bufs. Point-to-point instead of mpi_alltoall: the collective
! would copy the own block (and possibly stage blocks) on the host with managed memory.
call mpi_comm_rank(cart_comm_dir(0),me,ierr)
nreq=0
!CUDA-aware MPI GPU-GPU communicaton (by default hpc-sdk is CUDA-aware)
!$acc host_data use_device(bufs,bufr)
do j=0,nzcpu-1
 if(j.ne.me)then
  nreq=nreq+1
  call mpi_irecv(bufr(1,1,1,1,j),numel,mpi_double_precision,j,0,cart_comm_dir(0),req(nreq),ierr)
 endif
enddo
do j=0,nzcpu-1
 if(j.ne.me)then
  nreq=nreq+1
  call mpi_isend(bufs(1,1,1,1,j),numel,mpi_double_precision,j,0,cart_comm_dir(0),req(nreq),ierr)
 endif
enddo
call mpi_waitall(nreq,req,mpi_statuses_ignore,ierr)
!$acc end host_data

! unpack: block j holds the local z slab with the y points of rank j (offset i0, cnt points)
!$acc parallel loop collapse(5) private(i0,cnt)
do j=0,nzcpu-1
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
