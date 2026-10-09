subroutine conf_non_linear_pol(sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz)

    use commondata
    use par_size
    use velocity
    use wavenumber
    use polymer
    use mpi
    use phase_field
    use assemble
    use derivatives
    
#define gieflag giecompflag
#define fenepflag fenepcompflag
#define phiflag phicompflag
#define carflag carnatureflag
    
    double precision :: sconfxx(spx,nz,spy,2)
    double precision, allocatable, dimension(:,:,:,:) :: ucxxs,vcxxs,wcxxs,a4s
    double precision, allocatable, dimension(:,:,:) :: ucxx,vcxx,wcxx,a4f
    
    double precision :: sconfxy(spx,nz,spy,2)
    double precision, allocatable, dimension(:,:,:,:) :: ucxys,vcxys,wcxys
    double precision, allocatable, dimension(:,:,:) :: ucxy,vcxy,wcxy
    
    double precision :: sconfxz(spx,nz,spy,2)
    double precision, allocatable, dimension(:,:,:,:) :: ucxzs,vcxzs,wcxzs
    double precision, allocatable, dimension(:,:,:) :: ucxz,vcxz,wcxz
    
    double precision :: sconfyy(spx,nz,spy,2)
    double precision, allocatable, dimension(:,:,:,:) :: ucyys,vcyys,wcyys
    double precision, allocatable, dimension(:,:,:) :: ucyy,vcyy,wcyy
    
    double precision :: sconfyz(spx,nz,spy,2)
    double precision, allocatable, dimension(:,:,:,:) :: ucyzs,vcyzs,wcyzs
    double precision, allocatable, dimension(:,:,:) :: ucyz,vcyz,wcyz
    
    double precision :: sconfzz(spx,nz,spy,2)
    double precision, allocatable, dimension(:,:,:,:) :: uczzs,vczzs,wczzs
    double precision, allocatable, dimension(:,:,:) :: uczz,vczz,wczz

    double precision :: phif,zita,mask,epsnum
    
    
    integer :: indx,indy
    integer :: i,j,k
    
    epsnum=0.0000000001d0 

    indx=cstart(1)
    indy=cstart(3)
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! NON LINEAR XX 

allocate(ucxx(nx,fpz,fpy))
allocate(vcxx(nx,fpz,fpy))
allocate(wcxx(nx,fpz,fpy))

allocate(ucxxs(spx,nz,spy,2))
allocate(vcxxs(spx,nz,spy,2))
allocate(wcxxs(spx,nz,spy,2))

! form products of the advective terms, ...
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ucxx(i,k,j)=u(i,k,j)*C_xxf(i,k,j)
      vcxx(i,k,j)=v(i,k,j)*C_xxf(i,k,j)
      wcxx(i,k,j)=w(i,k,j)*C_xxf(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(ucxx,ucxxs,1)
call phys_to_spectral(vcxx,vcxxs,1)
call phys_to_spectral(wcxx,wcxxs,1)

deallocate(ucxx)
deallocate(vcxx)
deallocate(wcxx)

call dz(wcxxs,sconfxx)

!$acc parallel loop collapse(2) 
do j=1,spy
  do i=1,spx
    sconfxx(i,:,j,1)=-(sconfxx(i,:,j,1)-kx(i+indx)*ucxxs(i,:,j,2)-ky(j+indy)*vcxxs(i,:,j,2))
    sconfxx(i,:,j,2)=-(sconfxx(i,:,j,2)+kx(i+indx)*ucxxs(i,:,j,1)+ky(j+indy)*vcxxs(i,:,j,1))
  enddo
enddo

deallocate(ucxxs)
deallocate(vcxxs)
deallocate(wcxxs)
!now we have the advective terms
!polymer isotropic terms 
allocate(a4f(nx,fpz,fpy))
allocate(a4s(spx,nz,spy,2))


!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      mask=merge(1.0d0, 0.0d0, &
      (dabs(diag22(i,k,j)-diag11(i,k,j))<=epsnum).and.&
      (dabs(diag33(i,k,j)-diag11(i,k,j))<=epsnum).and.&
      (dabs(diag33(i,k,j)-diag22(i,k,j))<=epsnum))
      !!!!!!!!!!!!!!!!
      a4f(i,k,j)=2.0d0*(C_xxf(i,k,j)*omega1(i,k,j)+C_xzf(i,k,j)*omega2(i,k,j)+ &
              &         (1.d0-mask)*(emme11(i,k,j)*erre11(i,k,j)**2.d0+emme22(i,k,j)*erre12(i,k,j)**2.d0 + &
              &                      emme33(i,k,j)*erre13(i,k,j)**2.d0)+mask*(dudx(i,k,j)))
    enddo
  enddo
enddo


!additional non linear term for shear thinning
#if gieflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)-mob*((diag11(i,k,j)+1.d0/diag11(i,k,j))*erre11(i,k,j)**2.d0+   &
        &                        (diag22(i,k,j)+1.d0/diag22(i,k,j))*erre12(i,k,j)**2.d0+   &
        &                        (diag33(i,k,j)+1.d0/diag33(i,k,j))*erre13(i,k,j)**2.d0-2.d0)/Wi  
    enddo
  enddo
enddo
#endif

#if fenepflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
      zitaVel(i,k,j)=(L**2.0d0-3.0d0)/(L**2.0d0-(diag11(i,k,j)+diag22(i,k,j)+diag33(i,k,j)))
      !!!!!!!!!!!!!!!!!!!!!!!!
      a4f(i,k,j)=a4f(i,k,j)+((erre11(i,k,j)**2.d0)/diag11(i,k,j)+(erre12(i,k,j)**2.d0)/diag22(i,k,j)+ &
            &                (erre13(i,k,j)**2.d0)/diag33(i,k,j)-zitaVel(i,k,j))/Wi
   enddo
 enddo
#elif fenepflag != 1 
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
      do i=1,nx
        a4f(i,k,j)=a4f(i,k,j)+((erre11(i,k,j)**2.d0)/diag11(i,k,j)+ &
        &                      (erre12(i,k,j)**2.d0)/diag22(i,k,j)+ &
        &                      (erre13(i,k,j)**2.d0)/diag33(i,k,j)-1.d0)/Wi
    enddo
  enddo
enddo
#endif


#if (phiflag == 1 && carflag == 0)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0-phif)
    enddo
  enddo
enddo
#elif (phiflag == 1 && carflag == 1)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0+phif)
    enddo
  enddo
enddo
#endif

call phys_to_spectral(a4f,a4s,1)

deallocate(a4f)

!$acc  kernels
sconfxx=sconfxx+a4s
!$acc end kernels

deallocate(a4s) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! NON LINEAR XY 
allocate(ucxy(nx,fpz,fpy))
allocate(vcxy(nx,fpz,fpy))
allocate(wcxy(nx,fpz,fpy))

allocate(ucxys(spx,nz,spy,2))
allocate(vcxys(spx,nz,spy,2))
allocate(wcxys(spx,nz,spy,2))

! form products of the advective terms, ...
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ucxy(i,k,j)=u(i,k,j)*C_xyf(i,k,j)
      vcxy(i,k,j)=v(i,k,j)*C_xyf(i,k,j)
      wcxy(i,k,j)=w(i,k,j)*C_xyf(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(ucxy,ucxys,1)
call phys_to_spectral(vcxy,vcxys,1)
call phys_to_spectral(wcxy,wcxys,1)

deallocate(ucxy)
deallocate(vcxy)
deallocate(wcxy)

call dz(wcxys,sconfxy)

!$acc parallel loop collapse(2) 
do j=1,spy
  do i=1,spx
    sconfxy(i,:,j,1)=-(sconfxy(i,:,j,1)-kx(i+indx)*ucxys(i,:,j,2)-ky(j+indy)*vcxys(i,:,j,2))
    sconfxy(i,:,j,2)=-(sconfxy(i,:,j,2)+kx(i+indx)*ucxys(i,:,j,1)+ky(j+indy)*vcxys(i,:,j,1))
  enddo
enddo

deallocate(ucxys)
deallocate(vcxys)
deallocate(wcxys)

!now we have the advective terms
!polymer isotropic terms 

allocate(a4f(nx,fpz,fpy))
allocate(a4s(spx,nz,spy,2))

!$acc parallel loop collapse(3)
do j=1,fpy
   do k=1,fpz  
      do i=1,nx
        mask=merge(1.0d0, 0.0d0, &
                  (dabs(diag22(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                  (dabs(diag33(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                  (dabs(diag33(i,k,j)-diag22(i,k,j))<=epsnum))
            !!!!!!!!!!!!!!!!        
        a4f(i,k,j)=omega1(i,k,j)*(C_yyf(i,k,j)-C_xxf(i,k,j))+omega2(i,k,j)*C_yzf(i,k,j)+omega3(i,k,j)*C_xzf(i,k,j) + &
          &         (1.d0-mask)*2.0d0*(emme11(i,k,j)*erre11(i,k,j)*erre21(i,k,j)+emme22(i,k,j)*erre12(i,k,j)*erre22(i,k,j) + &
          &                            emme33(i,k,j)*erre13(i,k,j)*erre23(i,k,j))+mask*(dudy(i,k,j)+dvdx(i,k,j))
      enddo
   enddo
enddo

 
!additional non linear term for shear thinning
#if gieflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)-mob*((diag11(i,k,j)+1.d0/diag11(i,k,j))*erre11(i,k,j)*erre21(i,k,j)+   &
        &                        (diag22(i,k,j)+1.d0/diag22(i,k,j))*erre12(i,k,j)*erre22(i,k,j)+   &
        &                        (diag33(i,k,j)+1.d0/diag33(i,k,j))*erre13(i,k,j)*erre23(i,k,j))/Wi  
    enddo
  enddo
enddo
#endif

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
      do i=1,nx
        a4f(i,k,j)=a4f(i,k,j)+((erre11(i,k,j)*erre21(i,k,j))/diag11(i,k,j)+(erre12(i,k,j)*erre22(i,k,j))/diag22(i,k,j)+ &
        &                      (erre13(i,k,j)*erre23(i,k,j))/diag33(i,k,j))/Wi
    enddo
  enddo
enddo


#if (phiflag == 1 && carflag == 0)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
     phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
     a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0-phif)
    enddo
  enddo
enddo
#elif (phiflag == 1 && carflag == 1)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0+phif)
    enddo
  enddo
enddo
#endif

call phys_to_spectral(a4f,a4s,1)

deallocate(a4f)

!$acc  kernels
sconfxy=sconfxy+a4s
!$acc end kernels
!
deallocate(a4s) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! NON LINEAR XZ 

!maybe here we should add a check of some type to avoid hadamard instabilities

allocate(ucxz(nx,fpz,fpy))
allocate(vcxz(nx,fpz,fpy))
allocate(wcxz(nx,fpz,fpy))

allocate(ucxzs(spx,nz,spy,2))
allocate(vcxzs(spx,nz,spy,2))
allocate(wcxzs(spx,nz,spy,2))

! form products of the advective terms, ...
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ucxz(i,k,j)=u(i,k,j)*C_xzf(i,k,j)
      vcxz(i,k,j)=v(i,k,j)*C_xzf(i,k,j)
      wcxz(i,k,j)=w(i,k,j)*C_xzf(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(ucxz,ucxzs,1)
call phys_to_spectral(vcxz,vcxzs,1)
call phys_to_spectral(wcxz,wcxzs,1)

deallocate(ucxz)
deallocate(vcxz)
deallocate(wcxz)

call dz(wcxzs,sconfxz)

!$acc parallel loop collapse(2) 
do j=1,spy
  do i=1,spx
    sconfxz(i,:,j,1)=-(sconfxz(i,:,j,1)-kx(i+indx)*ucxzs(i,:,j,2)-ky(j+indy)*vcxzs(i,:,j,2))
    sconfxz(i,:,j,2)=-(sconfxz(i,:,j,2)+kx(i+indx)*ucxzs(i,:,j,1)+ky(j+indy)*vcxzs(i,:,j,1))
  enddo
enddo

deallocate(ucxzs)
deallocate(vcxzs)
deallocate(wcxzs)

!now we have the advective terms

!polymer isotropic terms 

allocate(a4f(nx,fpz,fpy))
allocate(a4s(spx,nz,spy,2))


!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      mask=merge(1.0d0, 0.0d0, &
          (dabs(diag22(i,k,j)-diag11(i,k,j))<=epsnum).and.&
          (dabs(diag33(i,k,j)-diag11(i,k,j))<=epsnum).and.&
          (dabs(diag33(i,k,j)-diag22(i,k,j))<=epsnum))
            !!!!!!!!!!!!!!!! 
       a4f(i,k,j)=omega1(i,k,j)*C_yzf(i,k,j)+omega2(i,k,j)*(C_zzf(i,k,j)-C_xxf(i,k,j))-omega3(i,k,j)*C_xyf(i,k,j)+ &
         &        2.0d0*(1.d0-mask)*(emme11(i,k,j)*erre11(i,k,j)*erre31(i,k,j)+emme22(i,k,j)*erre12(i,k,j)*erre32(i,k,j) + &
         &               emme33(i,k,j)*erre13(i,k,j)*erre33(i,k,j))+mask*(dudz(i,k,j)+dwdx(i,k,j))
     enddo
  enddo
enddo

!if(rank.eq.0)then
!print *, 'a4f-xz: ', a4f(1,1,1)
!endif
!additional non linear term for shear thinning
#if gieflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)-mob*((diag11(i,k,j)+1.d0/diag11(i,k,j))*erre11(i,k,j)*erre31(i,k,j)+   &
        &                        (diag22(i,k,j)+1.d0/diag22(i,k,j))*erre12(i,k,j)*erre32(i,k,j)+   &
        &                        (diag33(i,k,j)+1.d0/diag33(i,k,j))*erre13(i,k,j)*erre33(i,k,j))/Wi  
    enddo
  enddo
enddo
#endif

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
      do i=1,nx
        a4f(i,k,j)=a4f(i,k,j)+((erre11(i,k,j)*erre31(i,k,j))/diag11(i,k,j)+&
        &                      (erre12(i,k,j)*erre32(i,k,j))/diag22(i,k,j)+ &
        &                      (erre13(i,k,j)*erre33(i,k,j))/diag33(i,k,j))/Wi
    enddo
  enddo
enddo


#if (phiflag == 1 && carflag == 0)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0-phif)
    enddo
  enddo
enddo
#elif (phiflag == 1 && carflag == 1)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0+phif)
    enddo
  enddo
enddo
#endif

call phys_to_spectral(a4f,a4s,1)

deallocate(a4f)

!$acc  kernels
sconfxz=sconfxz+a4s
!$acc end kernels

deallocate(a4s)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! NON LINEAR YY

allocate(ucyy(nx,fpz,fpy))
allocate(vcyy(nx,fpz,fpy))
allocate(wcyy(nx,fpz,fpy))

allocate(ucyys(spx,nz,spy,2))
allocate(vcyys(spx,nz,spy,2))
allocate(wcyys(spx,nz,spy,2))

! form products of the advective terms, ...
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ucyy(i,k,j)=u(i,k,j)*C_yyf(i,k,j)
      vcyy(i,k,j)=v(i,k,j)*C_yyf(i,k,j)
      wcyy(i,k,j)=w(i,k,j)*C_yyf(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(ucyy,ucyys,1)
call phys_to_spectral(vcyy,vcyys,1)
call phys_to_spectral(wcyy,wcyys,1)

deallocate(ucyy)
deallocate(vcyy)
deallocate(wcyy)

call dz(wcyys,sconfyy)

!$acc parallel loop collapse(2) 
do j=1,spy
  do i=1,spx
    sconfyy(i,:,j,1)=-(sconfyy(i,:,j,1)-kx(i+indx)*ucyys(i,:,j,2)-ky(j+indy)*vcyys(i,:,j,2))
    sconfyy(i,:,j,2)=-(sconfyy(i,:,j,2)+kx(i+indx)*ucyys(i,:,j,1)+ky(j+indy)*vcyys(i,:,j,1))
  enddo
enddo

deallocate(ucyys)
deallocate(vcyys)
deallocate(wcyys)


!now we have the advective terms

!polymer isotropic terms 

allocate(a4f(nx,fpz,fpy))
allocate(a4s(spx,nz,spy,2))

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      mask=merge(1.0d0, 0.0d0, &
                (dabs(diag22(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                (dabs(diag33(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                (dabs(diag33(i,k,j)-diag22(i,k,j))<=epsnum))
      !!!!!!!!!!!!!!!!
      a4f(i,k,j)=2.0d0*(-C_xyf(i,k,j)*omega1(i,k,j)+C_yzf(i,k,j)*omega3(i,k,j)+ &
              &         (1.d0-mask)*(emme11(i,k,j)*erre21(i,k,j)**2.d0+emme22(i,k,j)*erre22(i,k,j)**2.d0 + &
              &              emme33(i,k,j)*erre23(i,k,j)**2.d0)+mask*dvdy(i,k,j))
    enddo
  enddo
enddo


!additional non linear term for shear thinning
#if gieflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)-mob*((diag11(i,k,j)+1.d0/diag11(i,k,j))*erre21(i,k,j)**2.d0+   &
        &                        (diag22(i,k,j)+1.d0/diag22(i,k,j))*erre22(i,k,j)**2.d0+   &
        &                        (diag33(i,k,j)+1.d0/diag33(i,k,j))*erre23(i,k,j)**2.d0-2.d0)/Wi  
    enddo
  enddo
enddo
#endif

#if fenepflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)+((erre21(i,k,j)**2.d0)/diag11(i,k,j)+(erre22(i,k,j)**2.d0)/diag22(i,k,j)+ &
            &                (erre23(i,k,j)**2.d0)/diag33(i,k,j)-zitaVel(i,k,j))/Wi
   enddo
 enddo
enddo
#elif fenepflag != 1 
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)+((erre21(i,k,j)**2.d0)/diag11(i,k,j)+(erre22(i,k,j)**2.d0)/diag22(i,k,j)+ &
           &                 (erre23(i,k,j)**2.d0)/diag33(i,k,j)-1.d0)/Wi
   enddo
  enddo
enddo
#endif



#if (phiflag == 1 && carflag == 0)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0-phif)
    enddo
  enddo
enddo
#elif (phiflag == 1 && carflag == 1)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0+phif)
    enddo
  enddo
enddo
#endif

call phys_to_spectral(a4f,a4s,1)

deallocate(a4f)

!$acc  kernels
sconfyy=sconfyy+a4s
!$acc end kernels

deallocate(a4s)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! NON LINEAR YZ

allocate(ucyz(nx,fpz,fpy))
allocate(vcyz(nx,fpz,fpy))
allocate(wcyz(nx,fpz,fpy))

allocate(ucyzs(spx,nz,spy,2))
allocate(vcyzs(spx,nz,spy,2))
allocate(wcyzs(spx,nz,spy,2))

! form products of the advective terms, ...
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      ucyz(i,k,j)=u(i,k,j)*C_yzf(i,k,j)
      vcyz(i,k,j)=v(i,k,j)*C_yzf(i,k,j)
      wcyz(i,k,j)=w(i,k,j)*C_yzf(i,k,j)
    enddo
  enddo
enddo

call phys_to_spectral(ucyz,ucyzs,1)
call phys_to_spectral(vcyz,vcyzs,1)
call phys_to_spectral(wcyz,wcyzs,1)

deallocate(ucyz)
deallocate(vcyz)
deallocate(wcyz)

call dz(wcyzs,sconfyz)

!$acc parallel loop collapse(2) 
do j=1,spy
  do i=1,spx
    sconfyz(i,:,j,1)=-(sconfyz(i,:,j,1)-kx(i+indx)*ucyzs(i,:,j,2)-ky(j+indy)*vcyzs(i,:,j,2))
    sconfyz(i,:,j,2)=-(sconfyz(i,:,j,2)+kx(i+indx)*ucyzs(i,:,j,1)+ky(j+indy)*vcyzs(i,:,j,1))
  enddo
enddo

deallocate(ucyzs)
deallocate(vcyzs)
deallocate(wcyzs)

!now we have the advective terms

!polymer isotropic terms 

allocate(a4f(nx,fpz,fpy))
allocate(a4s(spx,nz,spy,2))

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz  
     do i=1,nx
      mask=merge(1.0d0, 0.0d0, &
                (dabs(diag22(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                (dabs(diag33(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                (dabs(diag33(i,k,j)-diag22(i,k,j))<=epsnum))
           !!!!!!!!!!!!!!!!
           a4f(i,k,j)=-omega1(i,k,j)*C_xzf(i,k,j)-omega2(i,k,j)*C_xyf(i,k,j)+omega3(i,k,j)*(C_zzf(i,k,j)-C_yyf(i,k,j))+ &
             &         2.0d0*(1.d0-mask)*(emme11(i,k,j)*erre21(i,k,j)*erre31(i,k,j)+emme22(i,k,j)*erre22(i,k,j)*erre32(i,k,j) + &
             &                emme33(i,k,j)*erre23(i,k,j)*erre33(i,k,j))+mask*(dvdz(i,k,j)+dwdy(i,k,j))
     enddo
  enddo
enddo


deallocate(omega1)

!additional non linear term for shear thinning
#if gieflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)-mob*((diag11(i,k,j)+1.d0/diag11(i,k,j))*erre21(i,k,j)*erre31(i,k,j)+   &
        &                        (diag22(i,k,j)+1.d0/diag22(i,k,j))*erre22(i,k,j)*erre32(i,k,j)+   &
        &                        (diag33(i,k,j)+1.d0/diag33(i,k,j))*erre23(i,k,j)*erre33(i,k,j))/Wi  
    enddo
  enddo
enddo
#endif

!!!!!!!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
      do i=1,nx
        a4f(i,k,j)=a4f(i,k,j)+((erre21(i,k,j)*erre31(i,k,j))/diag11(i,k,j)+(erre22(i,k,j)*erre32(i,k,j))/diag22(i,k,j)+ &
        &                      (erre23(i,k,j)*erre33(i,k,j))/diag33(i,k,j))/Wi
    enddo
  enddo
enddo


#if (phiflag == 1 && carflag == 0)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0-phif)
    enddo
  enddo
enddo
#elif (phiflag == 1 && carflag == 1)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0+phif)
    enddo
  enddo
enddo
#endif

call phys_to_spectral(a4f,a4s,1)

deallocate(a4f)

!$acc  kernels
sconfyz=sconfyz+a4s
!$acc end kernels

deallocate(a4s)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
! NON LINEAR ZZ


allocate(uczz(nx,fpz,fpy))
allocate(vczz(nx,fpz,fpy))
allocate(wczz(nx,fpz,fpy))

allocate(uczzs(spx,nz,spy,2))
allocate(vczzs(spx,nz,spy,2))
allocate(wczzs(spx,nz,spy,2))

! form products of the advective terms, ...
!$acc parallel loop collapse(3)
do j=1,fpy
    do k=1,fpz
      do i=1,nx
        uczz(i,k,j)=u(i,k,j)*C_zzf(i,k,j)
        vczz(i,k,j)=v(i,k,j)*C_zzf(i,k,j)
        wczz(i,k,j)=w(i,k,j)*C_zzf(i,k,j)
      enddo
    enddo
enddo

call phys_to_spectral(uczz,uczzs,1)
call phys_to_spectral(vczz,vczzs,1)
call phys_to_spectral(wczz,wczzs,1)

deallocate(uczz)
deallocate(vczz)
deallocate(wczz)

call dz(wczzs,sconfzz)

!$acc parallel loop collapse(2) 
do j=1,spy
    do i=1,spx
      sconfzz(i,:,j,1)=-(sconfzz(i,:,j,1)-kx(i+indx)*uczzs(i,:,j,2)-ky(j+indy)*vczzs(i,:,j,2))
      sconfzz(i,:,j,2)=-(sconfzz(i,:,j,2)+kx(i+indx)*uczzs(i,:,j,1)+ky(j+indy)*vczzs(i,:,j,1))
    enddo
enddo

deallocate(uczzs)
deallocate(vczzs)
deallocate(wczzs)

!now we have the advective terms
!polymer isotropic terms 

allocate(a4f(nx,fpz,fpy))
allocate(a4s(spx,nz,spy,2))

!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      mask=merge(1.0d0, 0.0d0, &
                (dabs(diag22(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                (dabs(diag33(i,k,j)-diag11(i,k,j))<=epsnum).and.&
                (dabs(diag33(i,k,j)-diag22(i,k,j))<=epsnum))
             !!!!!!!!!!!!!!!!
      a4f(i,k,j)=2.0d0*(-C_xzf(i,k,j)*omega2(i,k,j)-C_yzf(i,k,j)*omega3(i,k,j)+ &
              &         (1.d0-mask)*(emme11(i,k,j)*erre31(i,k,j)**2.d0+emme22(i,k,j)*erre32(i,k,j)**2.d0 + &
              &         emme33(i,k,j)*erre33(i,k,j)**2.d0)+mask*dwdz(i,k,j))
    enddo
  enddo
enddo



deallocate(omega2)
deallocate(omega3)
deallocate(emme11)
deallocate(emme22)
deallocate(emme33)

!additional non linear term for shear thinning
#if gieflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)-mob*((diag11(i,k,j)+1.d0/diag11(i,k,j))*erre31(i,k,j)**2.d0  +   &
        &                        (diag22(i,k,j)+1.d0/diag22(i,k,j))*erre32(i,k,j)**2.d0+   &
        &                        (diag33(i,k,j)+1.d0/diag33(i,k,j))*erre33(i,k,j)**2.d0-2.d0)/Wi  
    enddo
  enddo
enddo
#endif
!
#if fenepflag == 1
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)+((erre31(i,k,j)**2.d0)/diag11(i,k,j)+(erre32(i,k,j)**2.d0)/diag22(i,k,j)+ &
            &                (erre33(i,k,j)**2.d0)/diag33(i,k,j)-zitaVel(i,k,j))/Wi
   enddo
 enddo
enddo
#elif fenepflag != 1 
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      a4f(i,k,j)=a4f(i,k,j)+((erre31(i,k,j)**2.d0)/diag11(i,k,j)+&
           &                 (erre32(i,k,j)**2.d0)/diag22(i,k,j)+ &
           &                 (erre33(i,k,j)**2.d0)/diag33(i,k,j)-1.d0)/Wi
   enddo
  enddo
enddo
#endif



#if (phiflag == 1 && carflag == 0)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0-phif)
    enddo
  enddo
enddo
#elif (phiflag == 1 && carflag == 1)
!$acc parallel loop collapse(3)
do j=1,fpy
  do k=1,fpz
    do i=1,nx
      phif=min(1.0d0,max(-1.0d0,phi(i,k,j)))
      a4f(i,k,j)=a4f(i,k,j)*0.5d0*(1.0d0+phif)
    enddo
  enddo
enddo
#endif


call phys_to_spectral(a4f,a4s,1)

deallocate(a4f)

!$acc  kernels
sconfzz=sconfzz+a4s
!$acc end kernels

deallocate(a4s)


!!!!!!!!!!!!!!!!!!!!!  
 return
end subroutine