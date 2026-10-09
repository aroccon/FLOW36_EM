subroutine euler(s1,s2,s3,h1,h2,h3)

use commondata
use sim_par
use par_size

double precision :: s1(spx,nz,spy,2), s2(spx,nz,spy,2), s3(spx,nz,spy,2)
double precision :: h1(spx,nz,spy,2), h2(spx,nz,spy,2), h3(spx,nz,spy,2)

!$acc kernels
h1=dt*s1
h2=dt*s2
h3=dt*s3
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adams_bashforth(s1,s2,s3,h1,h2,h3)

use commondata
use sim_par
use par_size
use sterms

double precision :: s1(spx,nz,spy,2), s2(spx,nz,spy,2), s3(spx,nz,spy,2)
double precision :: h1(spx,nz,spy,2), h2(spx,nz,spy,2), h3(spx,nz,spy,2)

!$acc kernels
h1=0.5d0*dt*(3.0d0*s1-s1_o)
h2=0.5d0*dt*(3.0d0*s2-s2_o)
h3=0.5d0*dt*(3.0d0*s3-s3_o)
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine euler_phi(sphi,hphi)

use commondata
use sim_par
use par_size

double precision :: sphi(spx,nz,spy,2), hphi(spx,nz,spy,2)

!$acc kernels
hphi=dt*sphi
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adams_bashforth_phi(sphi,hphi)

use commondata
use sim_par
use par_size
use sterms

double precision :: sphi(spx,nz,spy,2), hphi(spx,nz,spy,2)

!$acc kernels
hphi=0.5d0*dt*(3.0d0*sphi-sphi_o)
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine euler_psi(spsi,hpsi)

use commondata
use sim_par
use par_size
use dual_grid

double precision, dimension(spxpsi,npsiz,spypsi,2) :: spsi, hpsi

hpsi=dt*spsi

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adams_bashforth_psi(spsi,hpsi)

use commondata
use sim_par
use par_size
use sterms
use dual_grid

double precision, dimension(spxpsi,npsiz,spypsi,2) :: spsi, hpsi

hpsi=0.5d0*dt*(3.0d0*spsi-spsi_o)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine euler_theta(stheta,htheta)

use commondata
use sim_par
use par_size

double precision :: stheta(spx,nz,spy,2), htheta(spx,nz,spy,2)

!$acc kernels
htheta=dt*stheta
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine adams_bashforth_theta(stheta,htheta)

use commondata
use sim_par
use par_size
use sterms

double precision :: stheta(spx,nz,spy,2), htheta(spx,nz,spy,2)

!$acc kernels
htheta=0.5d0*dt*(3.0d0*stheta-stheta_o)
!$acc end kernels

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine crank_nicholson_pol(pols1,pols2,pols3,polh1,polh2,polh3)

    use commondata
    use sim_par
    use par_size
    use sterms
    
    double precision :: pols1(spx,nz,spy,2), pols2(spx,nz,spy,2), pols3(spx,nz,spy,2)
    double precision :: polh1(spx,nz,spy,2), polh2(spx,nz,spy,2), polh3(spx,nz,spy,2)
    
    !$acc kernels
    polh1=0.5d0*dt*(pols1+pols1_n)
    polh2=0.5d0*dt*(pols2+pols2_n)
    polh3=0.5d0*dt*(pols3+pols3_n)
    !$acc end kernels
    
    return
 end

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine euler_pol(sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz,hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz)
!subroutine euler_pol(sconfxx,sconfxz,sconfzz,hconfxx,hconfxz,hconfzz)

    use commondata
    use sim_par
    use par_size
    
    double precision,dimension(spx,nz,spy,2)   :: sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz
    double precision,dimension(spx,nz,spy,2)   :: hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz
    
    !$acc kernels
    hconfxx=dt*sconfxx
    hconfxy=dt*sconfxy
    hconfxz=dt*sconfxz
    hconfyy=dt*sconfyy
    hconfyz=dt*sconfyz
    hconfzz=dt*sconfzz
    !$acc end kernels
    
 return
end subroutine euler_pol
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine adams_bashforth_pol(sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz,hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz)
!subroutine adams_bashforth_pol(sconfxx,sconfxz,sconfzz,hconfxx,hconfxz,hconfzz)
    
    use commondata
    use sim_par
    use par_size
    use sterms
    
   double precision,dimension(spx,nz,spy,2)   :: sconfxx,sconfxy,sconfxz,sconfyy,sconfyz,sconfzz
   double precision,dimension(spx,nz,spy,2)   :: hconfxx,hconfxy,hconfxz,hconfyy,hconfyz,hconfzz
    
   !$acc kernels
    hconfxx=0.5d0*dt*(3.0d0*sconfxx-sconfxx_o)
    hconfxy=0.5d0*dt*(3.0d0*sconfxy-sconfxy_o)
    hconfxz=0.5d0*dt*(3.0d0*sconfxz-sconfxz_o)
    hconfyy=0.5d0*dt*(3.0d0*sconfyy-sconfyy_o)
    hconfyz=0.5d0*dt*(3.0d0*sconfyz-sconfyz_o)
    hconfzz=0.5d0*dt*(3.0d0*sconfzz-sconfzz_o)
   !$acc end kernels
    
 return
end subroutine adams_bashforth_pol