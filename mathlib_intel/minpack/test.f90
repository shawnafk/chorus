  subroutine test_hybrd()
    implicit none
    integer :: J,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR,NWRITE
    real(kind=8) :: XTOL,EPSFCN,FACTOR,FNORM
    real(kind=8) :: X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9), &
         WA1(9),WA2(9),WA3(9),WA4(9)
    real(kind=8) :: ENORM,DPMPAR
    external :: FCN1
    NWRITE = 6
    N = 9
    DO J = 1, 9
       X(J) = -1.D0
    end DO
    LDFJAC = 9
    LR = 45
    XTOL = DSQRT(DPMPAR(1))
    MAXFEV = 2000
    ML = 1
    MU = 1
    EPSFCN = 0.D0
    MODE = 2
    DO  J = 1, 9
       DIAG(J) = 1.D0
    end DO
    FACTOR = 1.D2
    NPRINT = 0
    CALL HYBRD(FCN1,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,DIAG, &
         MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC, &
         R,LR,QTF,WA1,WA2,WA3,WA4)
    FNORM = ENORM(N,FVEC)
    WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N)
1000 FORMAT (5X,'FINAL L2 NORM OF THE RESIDUALS',D15.7 // &
         5X,'NUMBER OF FUNCTION EVALUATIONS',I10 // &
         5X,'EXIT PARAMETER',16X,I10 // &
         5X,'FINAL APPROXIMATE SOLUTION' // (5X,3D15.7))
    return
  end subroutine test_hybrd

  SUBROUTINE FCN1(N,X,FVEC,IFLAG)
    implicit none
    integer :: N,IFLAG
    real(kind=8) :: X(N),FVEC(N)
    integer :: K
    real(kind=8) :: ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
    DATA ZERO,ONE,TWO,THREE /0.D0,1.D0,2.D0,3.D0/
    IF (IFLAG .NE. 0) GO TO 5
    RETURN
5   CONTINUE
    DO K = 1, N
       TEMP = (THREE - TWO*X(K))*X(K)
       TEMP1 = ZERO
       IF (K .NE. 1) TEMP1 = X(K-1)
       TEMP2 = ZERO
       IF (K .NE. N) TEMP2 = X(K+1)
       FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
    end do
    RETURN
  end SUBROUTINE FCN1



  subroutine test_hybrd1()
    implicit none
    integer :: J,N,INFO,LWA,NWRITE
    real(kind=8) :: TOL,FNORM
    real(kind=8) :: X(1),FVEC(1),WA(180)
    real(kind=8) :: ENORM,DPMPAR
    EXTERNAL :: FCN2
    NWRITE = 6
    N = 1
    DO J = 1, 1
       X(J) = 1.D0
    end DO
    LWA = 180
    TOL = DSQRT(DPMPAR(1))
    CALL HYBRD1(FCN2,N,X,FVEC,TOL,INFO,WA,LWA)
    FNORM = ENORM(N,FVEC)
    WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
1000 FORMAT (5X,'FINAL L2 NORM OF THE RESIDUALS',D15.7 // &
         5X,'EXIT PARAMETER',16X,I10 // &
         5X,'FINAL APPROXIMATE SOLUTION' // (5X,3D15.7))
    return
  END subroutine test_hybrd1


  SUBROUTINE FCN2(N,X,FVEC,IFLAG)
    implicit none
    integer :: N,IFLAG
    real(kind=8) :: X(N),FVEC(N)
    FVEC(1) = X(1)**2+1.0D0
    RETURN
  END SUBROUTINE FCN2


  subroutine test_hybrj()
    implicit none
    integer :: J,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV,LR,NWRITE
    real(kind=8) :: XTOL,FACTOR,FNORM
    real(kind=8) :: X(9),FVEC(9),FJAC(9,9),DIAG(9),R(45),QTF(9), &
         WA1(9),WA2(9),WA3(9),WA4(9)
    real(kind=8) :: ENORM,DPMPAR
    EXTERNAL :: FCN3
    NWRITE = 6
    N = 9
    DO J = 1, 9
       X(J) = -1.D0
    end DO
    LDFJAC = 9
    LR = 45
    XTOL = DSQRT(DPMPAR(1))
    MAXFEV = 1000
    MODE = 2
    DO J = 1, 9
       DIAG(J) = 1.D0
    end DO
    FACTOR = 1.D2
    NPRINT = 0
    CALL HYBRJ(FCN3,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,DIAG, &
         MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,R,LR,QTF, &
         WA1,WA2,WA3,WA4)
    FNORM = ENORM(N,FVEC)
    WRITE (NWRITE,1000) FNORM,NFEV,NJEV,INFO,(X(J),J=1,N)
1000 FORMAT (5X,'FINAL L2 NORM OF THE RESIDUALS',D15.7 // &
         5X,'NUMBER OF FUNCTION EVALUATIONS',I10 // &
         5X,'NUMBER OF JACOBIAN EVALUATIONS',I10 // &
         5X,'EXIT PARAMETER',16X,I10 // &
         5X,'FINAL APPROXIMATE SOLUTION' // (5X,3D15.7))
    return
  END subroutine test_hybrj

  SUBROUTINE FCN3(N,X,FVEC,FJAC,LDFJAC,IFLAG)
    implicit none
    integer :: N,LDFJAC,IFLAG
    real(kind=8) :: X(N),FVEC(N),FJAC(LDFJAC,N)
    integer :: J,K
    real(kind=8) :: ONE,TEMP,TEMP1,TEMP2,FOUR,THREE,TWO,ZERO
    DATA ZERO,ONE,TWO,THREE,FOUR /0.D0,1.D0,2.D0,3.D0,4.D0/
    IF (IFLAG .NE. 0) GO TO 5
    RETURN
5   CONTINUE
    IF (IFLAG .EQ. 2) GO TO 20
    DO K = 1, N
       TEMP = (THREE - TWO*X(K))*X(K)
       TEMP1 = ZERO
       IF (K .NE. 1) TEMP1 = X(K-1)
       TEMP2 = ZERO
       IF (K .NE. N) TEMP2 = X(K+1)
       FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
    end do
    GO TO 50
20  CONTINUE
    DO K = 1, N
       DO J = 1, N
          FJAC(K,J) = ZERO
       end do
       FJAC(K,K) = THREE - FOUR*X(K)
       IF (K .NE. 1) FJAC(K,K-1) = -ONE
       IF (K .NE. N) FJAC(K,K+1) = -TWO
    end do
50  CONTINUE
    RETURN
  END SUBROUTINE FCN3


  subroutine test_hybrj1()
    implicit none
    integer :: J,N,LDFJAC,INFO,LWA,NWRITE
    real(kind=8) :: TOL,FNORM
    real(kind=8) :: X(9),FVEC(9),FJAC(9,9),WA(99)
    real(kind=8) :: ENORM,DPMPAR
    EXTERNAL :: FCN4
    NWRITE = 6
    N = 9
    DO J = 1, 9
       X(J) = -1.D0
    end do
    LDFJAC = 9
    LWA = 99
    TOL = DSQRT(DPMPAR(1))
    CALL HYBRJ1(FCN4,N,X,FVEC,FJAC,LDFJAC,TOL,INFO,WA,LWA)
    FNORM = ENORM(N,FVEC)
    WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
1000 FORMAT (5X,'FINAL L2 NORM OF THE RESIDUALS',D15.7 // &
         5X,'EXIT PARAMETER',16X,I10 // &
         5X,'FINAL APPROXIMATE SOLUTION' // (5X,3D15.7))
    return
  END subroutine test_hybrj1

  SUBROUTINE FCN4(N,X,FVEC,FJAC,LDFJAC,IFLAG)
    implicit none
    integer :: N,LDFJAC,IFLAG
    real(kind=8) :: X(N),FVEC(N),FJAC(LDFJAC,N)
    INTEGER :: J,K
    real(kind=8) :: ONE,TEMP,TEMP1,TEMP2,FOUR,THREE,TWO,ZERO
    DATA ZERO,ONE,TWO,THREE,FOUR /0.D0,1.D0,2.D0,3.D0,4.D0/
    IF (IFLAG .EQ. 2) GO TO 20
    DO K = 1, N
       TEMP = (THREE - TWO*X(K))*X(K)
       TEMP1 = ZERO
       IF (K .NE. 1) TEMP1 = X(K-1)
       TEMP2 = ZERO
       IF (K .NE. N) TEMP2 = X(K+1)
       FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
    end do
    GO TO 50
20  CONTINUE
    DO K = 1, N
       DO  J = 1, N
          FJAC(K,J) = ZERO
       end do
       FJAC(K,K) = THREE - FOUR*X(K)
       IF (K .NE. 1) FJAC(K,K-1) = -ONE
       IF (K .NE. N) FJAC(K,K+1) = -TWO
    end do
50  CONTINUE
    RETURN
  END SUBROUTINE FCN4



program test
  implicit none
  call test_hybrd()
  call test_hybrd1()
  call test_hybrj()
  call test_hybrj1()
#ifndef INTEL
  return
#endif
end program test



