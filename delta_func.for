!=======================================================================
!			Pablo Ouro Barba
!			Cardiff 2013-2014
!=======================================================================
!######################################################################
      Real function phi_r1smth(r)
!######################################################################
      implicit none 
      REAL, intent(in) :: r
      REAL :: PI,abr
       PI = 4.D0*DATAN(1.D0)
       abr=SQRT(r*r)
      IF (abr.ge.1.5) THEN
        phi_r1smth = 0.0
      ELSE IF ((abr.lt.1.5).and.(abr.ge.0.5)) THEN
        phi_r1smth = 9./8.-3.*abr/2+abr**2/2
      ELSE IF ((abr.lt.0.5).and.(abr.ge.0.0)) THEN
        phi_r1smth = 3./4.-abr**2
      END IF
      RETURN
      End Function
!######################################################################
      Real function phi_r2smth(r)
!######################################################################
      implicit none 
      REAL, intent(in) :: r
      REAL :: PI
       PI = 4.D0*DATAN(1.D0)
      IF (r.le.-2.5) THEN
        phi_r2smth = 0.0
      ELSE IF ((r.ge.-2.5).and.(r.le.-1.5)) THEN
        phi_r2smth= -1./8./PI*(-5.*PI-2.*PI*r+4.*sin(PI/4.*(-2.*r-1.)))
      ELSE IF ((r.ge.-1.5).and.(r.le.0.0)) THEN
        phi_r2smth = 1./4./PI*(PI+2.*sin(PI/4.*(-2.*r+1.))
     &                           -2.*sin(PI/4.*(-2.*r-1.)))
      ELSE IF ((r.ge.0.0).and.(r.le.1.5)) THEN
        phi_r2smth = 1./4./PI*(PI+2.*sin(PI/4.*(2.*r+1.))
     &                           -2.*sin(PI/4.*(2.*r-1.)))
      ELSE IF ((r.ge.1.5).and.(r.le.2.5)) THEN
        phi_r2smth= -1./8./PI*(-5.*PI+2.*PI*r+4.*sin(PI/4.*(2.*r-1.)))
      ELSE IF (r.ge.2.5) THEN
        phi_r2smth = 0.0
      END IF

      RETURN
      End Function
!######################################################################
      Real function phi_r3(r)
!######################################################################
      implicit none 
      REAL, intent(in) :: r
      IF (r.le.-1.5) THEN
        phi_r3 = 0.0
      ELSE IF ((r.ge.-1.5).and.(r.le.-0.5)) THEN
        phi_r3 = 1.0/6.0*(5.0+3.0*r-sqrt(-3.0*(1.0+r)**2+1.0))
      ELSE IF ((r.ge.-0.5).and.(r.le.0.0)) THEN
        phi_r3 = 1.0/3.0*(1.0+sqrt(-3.0*r**2+1.0))
      ELSE IF ((r.ge.0.0).and.(r.le.0.5)) THEN
        phi_r3 = 1.0/3.0*(1.0+sqrt(-3.0*r**2+1.0))
      ELSE IF ((r.ge.0.5).and.(r.le.1.5)) THEN
        phi_r3 = 1.0/6.0*(5.0-3.0*r-sqrt(-3.0*(1.0-r)**2+1.0))
      ELSE IF (r.ge.1.5) THEN
        phi_r3 = 0.0
      END IF

      RETURN
      End Function
!######################################################################
      Real function phi_r3smth2(r)
!######################################################################
      implicit none 
      REAL, intent(in) :: r
      REAL :: PI
       PI = 4.D0*DATAN(1.D0)
      IF (r.le.-2.0) THEN
        phi_r3smth2 = 0.0
      ELSE IF ((r.ge.-2.0).and.(r.le.-1.0)) THEN
        phi_r3smth2 = 55.0/48.0 - sqrt(3.0)*pi/108.0 + 13.0*r/12.0
     &+ r**2/4.0 + (-2.0*r-3.0)/48.0*sqrt(-12.0*r**2-36.0*r-23.0)
     &+ sqrt(3.0)/36.0*ASIN(sqrt(3.0)/2.0*(-2.0*r-3.0))
      ELSE IF ((r.ge.-1.0).and.(r.le.0.0)) THEN
        phi_r3smth2 = 17.0/48.0 + sqrt(3.0)*pi/108.0 - r/4.0
     &- r**2/4.0 + (2.0*r+1.0)/16.0*sqrt(-12.0*r**2-12.0*r+1.0)
     &- sqrt(3.0)/12.0*ASIN(sqrt(3.0)/2.0*(-2.0*r-1.0))
      ELSE IF ((r.ge.0.0).and.(r.le.1.0)) THEN
        phi_r3smth2 = 17.0/48.0 + sqrt(3.0)*pi/108.0 + r/4.0
     &- r**2/4.0 + (-2.0*r+1.0)/16.0*sqrt(-12.0*r**2+12.0*r+1.0)
     &- sqrt(3.0)/12.0*ASIN(sqrt(3.0)/2.0*(2.0*r-1.0))
      ELSE IF ((r.ge.1.0).and.(r.le.2.0)) THEN
        phi_r3smth2 = 55.0/48.0 - sqrt(3.0)*pi/108.0 - 13.0*r/12.0
     &+ r**2/4.0 + (2.0*r-3.0)/48.0*sqrt(-12.0*r**2+36.0*r-23.0)
     &+ sqrt(3.0)/36.0*ASIN(sqrt(3.0)/2.0*(2.0*r-3.0))
      ELSE IF (r.ge.2.0) THEN
        phi_r3smth2 = 0.0
      END IF

      RETURN
      End Function
!######################################################################
      Real function phi_r3smth(r)
!######################################################################
      implicit none 
      Double precision, intent(in) :: r
!      REAL :: PI
      REAL :: scal1,scal2,scal3,scal4,scal5,scal6

      scal1 = 1.095450017660160615261  ! 55.0/48.0 - sqrt(3.0)*pi/108.0
      scal2 = 1.083333333333333333333  ! 13.0/12.0
      scal3 = 0.4045499823398393847387 ! 17.0/48.0 + sqrt(3.0)*pi/108.0
      scal4 = 0.0481125224324688137091 ! sqrt(3.0)/36.0
      scal5 = 0.1443375672974064411273 ! sqrt(3.0)/12.0
      scal6 = 0.8660254037844386467637 ! sqrt(3.0)/2.0

!       PI = 4.D0*DATAN(1.D0)

      IF (r.le.-2.0) THEN
        phi_r3smth = 0.0
      ELSE IF ((r.ge.-2.0).and.(r.le.-1.0)) THEN
        phi_r3smth = scal1 + scal2*r + scal4*ASIN(scal6*(-2.0*r-3.0))
     &+ 0.25*r**2 + (-2.0*r-3.0)/48.0*sqrt(-12.0*r**2-36.0*r-23.0)
      ELSE IF ((r.ge.-1.0).and.(r.le.0.0)) THEN
        phi_r3smth = scal3 - r/4.0 - scal5*ASIN(scal6*(-2.0*r-1.0))
     &- 0.25*r**2 + (2.0*r+1.0)/16.0*sqrt(-12.0*r**2-12.0*r+1.0)
      ELSE IF ((r.ge.0.0).and.(r.le.1.0)) THEN

        phi_r3smth = scal3 + r/4.0 - scal5*ASIN(scal6*(2.0*r-1.0))
     &- 0.25*r**2 + (-2.0*r+1.0)/16.0*sqrt(-12.0*r**2+12.0*r+1.0)
      ELSE IF ((r.ge.1.0).and.(r.le.2.0)) THEN
        phi_r3smth = scal1 - scal2*r + scal4*ASIN(scal6*(2.0*r-3.0))
     &+ 0.25*r**2 + (2.0*r-3.0)/48.0*sqrt(-12.0*r**2+36.0*r-23.0)
      ELSE IF (r.ge.2.0) THEN
        phi_r3smth = 0.0
      END IF

      RETURN
      End Function
!######################################################################
      Real function phi_r4(r)
!######################################################################
      implicit none 
      REAL, intent(in) :: r
      IF (r.le.-2.0) THEN
        phi_r4 = 0.0
      ELSE IF ((r.ge.-2.0).and.(r.le.-1.0)) THEN
        phi_r4 = 1.0/8.0*(5.0+2.0*r-sqrt(-7.0-12.0*r-4.0*r**2))
      ELSE IF ((r.ge.-1.0).and.(r.le.0.0)) THEN
        phi_r4 = 1.0/8.0*(3.0+2.0*r+sqrt(1.0-4.0*r-4.0*r**2))
      ELSE IF ((r.ge.0.0).and.(r.le.1.0)) THEN
        phi_r4 = 1.0/8.0*(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*r**2))
      ELSE IF ((r.ge.1.0).and.(r.le.2.0)) THEN
        phi_r4 = 1.0/8.0*(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*r**2))
      ELSE IF (r.ge.2.0) THEN
        phi_r4 = 0.0
      END IF

      RETURN
      End Function
!######################################################################
      Real function phi_r4smth(r)
!######################################################################
      implicit none
      REAL, intent(in) :: r
      REAL :: PI,ar
       PI = 4.D0*DATAN(1.D0)
	ar=abs(r)
      IF ((ar.ge.0.0).and.(ar.le.0.5)) THEN
	phi_r4smth = 3.0/8.0+ PI/32.0 - r*r/4.0
      ELSE IF ((ar.ge.0.5).and.(ar.le.1.5)) THEN
	phi_r4smth = 1.0/4.0 + (1.0-ar)/8.0*SQRT(-2.0+8.0*ar -4.0*r*r)
     & -1.0/8.0*ASIN(sqrt(2.0)*(ar-1.0))
      ELSE IF ((ar.ge.1.5).and.(ar.le.2.5)) THEN
	phi_r4smth = 17.0/16.0 - PI/64.0 - 3.0*ar/4.0+ r*r/8.0
     & + (ar-2.0)/16.0 * SQRT(-14.0+ 16.0*ar - 4.0*r*r)
     & + 1.0/16.0*ASIN(sqrt(2.0)*(ar-2.0))
      ELSE IF (ar.ge.2.5) THEN
        phi_r4smth = 0.0
      END IF
      RETURN
      End Function
!######################################################################
      double precision function dh(dx,dy,dz,xij,yij,zij,Xl,Yl,Zl,order)
!######################################################################
! ..... The salient properties of the kernels DIRAC DELTA,dh are the following:
! * dh is a continuously diﬀerentiable function and therefore yields 
!   a smoother transfer than e.g. linear interpolation.
! * Interpolation using the kernels dh is second-order accurate 
!   for smooth ﬁelds (Uhlmann(2005) cf. Section 5.1.1).
! * The support of the regularized delta function is small, which makes
!   the evaluation of the sums in Eq. (9) relatively cheap. In particular,
!   we use the expression for dh deﬁned by Roma et al., involving only
!   three grid points in each coordinate direction.
! * (ORDER=3) A. Roma, C. Peskin, M. Berger, An adaptive version of the
!   immersed boundary method, J. Comput. Phys. 153 (1999)
! * (ORDER=4) C. Peskin, The immersed boundary method, 
!   Acta Numerica 11 (2002) 1–39
!
      implicit none 

      REAL,    intent(in) :: dx,dy,dz,xij,yij,zij,Xl,Yl,Zl
      INTEGER, intent(in) :: order
      REAL*8 ::phi_r2smth,phi_r3smth,phi_r4,phi_r3,phi_r1smth,phi_r4smth

      SELECT CASE (order)

         CASE (1)
      dh =  phi_r1smth((xij-Xl)/dx) 
     & * phi_r1smth((yij-Yl)/dy) * phi_r1smth((zij-Zl)/dz)
         CASE (2)
      dh =   phi_r2smth((xij-Xl)/dx) 
     & * phi_r2smth((yij-Yl)/dy) * phi_r2smth((zij-Zl)/dz)
         CASE (3)
       dh =   phi_r3smth((xij-Xl)/dx) 
     &   * phi_r3smth((yij-Yl)/dy) * phi_r3smth((zij-Zl)/dz)
         CASE (4)
       dh =   phi_r4smth((xij-Xl)/dx) 
     &   * phi_r4smth((yij-Yl)/dy) * phi_r4smth((zij-Zl)/dz)
         CASE (5)
       dh =   phi_r3((xij-Xl)/dx) 
     &   * phi_r3((yij-Yl)/dy) * phi_r3((zij-Zl)/dz)
         CASE (6)
      dh =   phi_r4((xij-Xl)/dx)
     &   * phi_r4((yij-Yl)/dy) * phi_r4((zij-Zl)/dz)

      Case default

         print*, '===ERROR==='
         print*, ' order of delta function is not selected '
          STOP

      End Select

      RETURN
      End Function
