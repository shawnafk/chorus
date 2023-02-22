subroutine test01
!
!******************************************************************************
!
!! TEST01 tests QAG.
!
!  QAG is an adaptive automatic integrator using a gauss-kronrod rule.
!
!  integrate cos(100*sin(x)) from 0 to pi.
!  exact answer is pi*j0(100), or roughly 0.06278740.
!
!  KEY chooses the order of the integration rule, from 1 to 6.
!
  implicit none
!
  real(kind=8), parameter :: a = -0.1D0
  real(kind=8) :: abserr
  real(kind=8) :: b
  real(kind=8), parameter :: epsabs = 1.0D-16
  real(kind=8), parameter :: epsrel = 1.0D-10
  real(kind=8), external :: f02
  integer :: ier
  integer, parameter :: npts2 = 4
  real(kind=8) :: points(npts2)
  integer, parameter :: key = 6
  integer :: neval
  real(kind=8) :: pi=3.1415926D0
  real(kind=8) :: result
  real(kind=8), parameter :: true = 0.06278740D+00
  integer :: limit = 1000
  integer :: leniw = 2004
  integer :: lenw = 4004
  integer :: last 
  integer :: iwork(1000)
  real(kind=8) :: work(2004)
!
  b = 1.1D0
  points(1) = 0.0D0
  points(2) = 1.0D0
  call dqag(f02, -0.1D0, 1.1D0, epsabs, epsrel, key, result, &
       abserr, neval, ier, limit, lenw, last, iwork, work)
  !call dqagp(f02, a, b, npts2, points, epsabs, epsrel, &
  !     result, abserr, neval, ier, &
  !     leniw, lenw, last, iwork, work)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test QAG'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is COS(100*SIN(X))'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end subroutine test01
 


subroutine test07
!
!******************************************************************************
!
!! TEST07 tests QAWO.
!
!
!  QAWO integrates functions of the form f(x)*sin(omega*x)
!  or f(x)*cos(omega*x)
!
!  integrate log(x)*sin(10*pi*x)dx from 0 to 1
!  exact answer is -(gamma+log(10*pi)-ci(10*pi))/(10*pi)
!  =-0.1281316
!
!  gamma is euler's constant.
!  ci is the cosine integral ci(x)=integral(x to infinity) cos(v)/v dv.
!
!  specify integr=1 for cosine integral, 2 for sine integral
!
  implicit none
!
  real(kind=8), parameter :: a = 0.0D+00
  real(kind=8) :: abserr
  real(kind=8), parameter :: b = 1.0D+00
  real(kind=8), parameter :: ci = -0.000842D+00
  real(kind=8), parameter :: epsabs = 1.0D-16
  real(kind=8), parameter :: epsrel = 1.0D-10
  real(kind=8), external :: f06
  real(kind=8), parameter :: gamma = 0.5772156649D+00
  integer :: ier
  integer, parameter :: integr = 2
  integer :: neval
  real(kind=8) :: omega
  real(kind=8) :: pi = 3.1415926D0
  real(kind=8) :: result
  real(kind=8) :: true
  integer :: leniw = 2000
  integer :: maxp1 = 1000
  integer :: lenw = 30000
  integer :: last
  integer :: iwork(2000)
  real(kind=8) :: work(30000)
 
!
!  set argument of the sine or cosine function
!
  omega = 10.0D+00 * pi

  call dqawo(f06, a, b, omega, integr, epsabs, epsrel, &
       result, abserr, neval, ier, leniw, maxp1, lenw, &
       last, iwork, work)
!
!  I can't find an evaluation of ci(10*pi)
!  the following is an estimate
!
  true = - ( gamma + log ( 10.0D+00 * pi) - ci ) / ( 10.0D+00 * pi )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Test QAWO'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)*sin(10*pi*x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end subroutine test07


function f02 ( x )
!
!******************************************************************************
!
!! F02 is the integrand function COS(100*SIN(X)).
!
  implicit none
!
  real(kind=8) :: f02
  real(kind=8) :: x
  if (x>=0.0D0 .and. x<=1.0D0) then
     f02 = 1.0D0/sqrt(x*(1.0D0-x))
  else
     f02 = 0.0D0
  end if
  return
end function f02



function f06 ( x )
!
!******************************************************************************
!
!! F06 is the integrand function LOG(X).
!
  implicit none
!
  real(kind=8) :: f06
  real(kind=8) :: x

  if ( x <= 0.0D+00 ) then
    f06 = 0.0D+00
  else
    f06 = log ( x )
  end if

  return
end function f06



program test
  implicit none
  call test01()
  call test07()
#ifndef INTEL
  return
#endif
end program test
