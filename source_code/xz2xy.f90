subroutine xz2xy(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,nz,npyy,2),wa(nsxx,npzz,ny,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nzcpu,8))
call xz2xy_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,ny,nz,a2a_send,a2a_recv)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine xz2xy_fg(wa,uc,dims,ngxx,nsxx,ngyy,npyy,ngzz,npzz)

use commondata
use dual_grid
use a2a_buffers

integer :: dims(2)
integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz
double precision :: uc(nsxx,npsiz,npyy,2),wa(nsxx,npzz,npsiy,2)

call a2a_reserve(int(ngxx,8)*int(ngzz,8)*int(ngyy,8)*2_8*int(nzcpu,8))
call xz2xy_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,npsiy,npsiz,a2a_send,a2a_recv)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! x-z to x-y pencils (direction 0, nzcpu ranks with the same y coordinate), gny x gnz grid
! wa: y complete, z split -> uc: z complete, y split
! Buffers bufs/bufr come from module a2a_buffers (allocated once).
! All the blocks are packed with one kernel, exchanged with one mpi_alltoall on cart_comm_dir(0)
! and unpacked with one kernel. Blocks have the padded size ngxx*ngzz*ngyy*2 as in the previous
! pairwise exchange; block j goes to / comes from the rank with z coordinate j.
subroutine xz2xy_a2a(wa,uc,ngxx,nsxx,ngyy,npyy,ngzz,npzz,gny,gnz,bufs,bufr)

use mpi
use commondata

integer :: nsxx,npyy,npzz,ngxx,ngyy,ngzz,gny,gnz
integer :: i,ky,kz,c,j,ry,rz,numel,i0,cnt
double precision :: uc(nsxx,gnz,npyy,2),wa(nsxx,npzz,gny,2)
double precision :: bufs(ngxx,ngzz,ngyy,2,0:nzcpu-1),bufr(ngxx,ngzz,ngyy,2,0:nzcpu-1)

ry=mod(gny,nzcpu)
rz=mod(gnz,nzcpu)

numel=ngxx*ngyy*ngzz*2

! pack: block j holds the y slab of rank j (offset i0, cnt points), all local z
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
     if(ky.le.cnt) bufs(i,kz,ky,c,j)=wa(i,kz,i0+ky,c)
    enddo
   enddo
  enddo
 enddo
enddo

!CUDA-aware MPI GPU-GPU communicaton (by default hpc-sdk is CUDA-aware)
!$acc host_data use_device(bufs,bufr)
call mpi_alltoall(bufs,numel,mpi_double_precision,bufr,numel,mpi_double_precision,cart_comm_dir(0),ierr)
!$acc end host_data

! unpack: block j holds the local y points with the z slab of rank j (offset i0, cnt points)
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
     if(kz.le.cnt) uc(i,i0+kz,ky,c)=bufr(i,kz,ky,c,j)
    enddo
   enddo
  enddo
 enddo
enddo

return
end
