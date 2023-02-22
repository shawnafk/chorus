module chorus_LIB
  use chorus_IO
  use chorus_TYPE
  implicit none

  interface inverse
     module procedure inverse_2D
     module procedure inverse_all
  end interface inverse

  interface trapez
     module procedure ctrapez
     module procedure rtrapez
  end interface trapez

contains
  subroutine  lin_solver(MAT, RHS, SOL, rank)
    ! The routine uses partial pivoting to solve a linear system using Gauss elimination
    ! Input parameters :
    !      integer, rank : the order of the matrix
    !      complex, MAT(rank,rank) :  the coefficient matrix
    !      complex, RHS(rank) : the right hand side
    ! Output, parameter :
    !      complex, SOL(rank) : the solution
    implicit none
    integer, intent(in) :: rank
    complex(fp), intent(in) :: MAT(rank,rank), RHS(rank)
    complex(fp), intent(out) :: SOL(rank)
    complex(fp), allocatable :: pivot(:,:), row(:)
    complex(fp) :: tmp
    integer :: i, j, k, p
    allocate(pivot(rank,rank), row(rank))
    pivot = MAT
    SOL = RHS
    do k = 1, rank
       !  Find the maximum element in column I.
       p = k
       do i = k+1, rank
          if (abs(pivot(p,k)) < abs(pivot(i,k))) then
             p = i
          end if
       end do
       if (pivot(p,k) == 0.0_fp) then
          call write_log('zero pivot on step', k)
          stop 
       end if
       !  Switch rows k and k.
       if (k /= p) then
          row = pivot(k,:)
          pivot(k,:) = pivot(p,:)
          pivot(p,:) = row
          tmp = SOL(k)
          SOL(k) = SOL(p)
          SOL(p) = tmp
       end if
       !  Scale the pivot row.
       pivot(k,k+1:rank) = pivot(k,k+1:rank)/pivot(k,k)
       SOL(k) = SOL(k)/pivot(k,k)
       pivot(k,k) = 1.0_fp
       !  Use the pivot row to eliminate lower entries in that column.
       do i = k+1, rank
          if (pivot(i,k) /= 0.0_fp) then
             tmp = -pivot(i,k)
             pivot(i,k) = 0.0_fp
             pivot(i,k+1:rank) = pivot(i,k+1:rank) + tmp*pivot(k,k+1:rank)
             SOL(i) = SOL(i) + tmp*SOL(k)
          end if
       end do
    end do
    !  Back solve.
    do j = rank, 2, -1
       SOL(1:j-1) = SOL(1:j-1) - pivot(1:j-1,j)*SOL(j)
    end do
    deallocate(pivot, row)
    return
  end subroutine lin_solver

  
  subroutine eigenvalue_2D(MAT)
    implicit none
    complex(fp), intent(inout) :: MAT(2,2)
    complex(fp) :: a, b, c, d, det, tree
    a = MAT(1,1)
    b = MAT(1,2)
    c = MAT(2,1)
    d = MAT(2,2)
    det = a*d - b*c
    tree = a + d
    MAT(1,1) = 0.5_fp*tree + sqrt(0.25_fp*tree**2-det)
    MAT(1,2) = (0.0_fp, 0.0_fp)
    MAT(2,1) = (0.0_fp, 0.0_fp)
    MAT(2,2) = 0.5_fp*tree - sqrt(0.25_fp*tree**2-det)
    return
  end subroutine eigenvalue_2D
    

  subroutine inverse_2D(MAT)
    implicit none
    complex(fp), intent(inout) :: MAT(2,2)
    complex(fp) :: a, b, c, d, det
    a = MAT(1,1)
    b = MAT(1,2)
    c = MAT(2,1)
    d = MAT(2,2)
    det = a*d - b*c
    MAT(1,1) = d/det
    MAT(1,2) = -b/det
    MAT(2,1) = -c/det
    MAT(2,2) = a/det
    return
  end subroutine inverse_2D


  subroutine Dinverse(DMAT, MAT, type)
    ! type 0:
    ! calculate MAT*dMAT^{-1}/dx
    ! type 1:
    ! calculate dMAT^{-1}/dx
    implicit none
    complex(fp), intent(inout) :: DMAT(2,2)
    complex(fp), intent(in) :: MAT(2,2)
    integer, intent(in) :: type
    complex(fp) :: a, b, c, d, da, db, dc, dd, det
    a = MAT(1,1)
    b = MAT(1,2)
    c = MAT(2,1)
    d = MAT(2,2)
    da = DMAT(1,1)
    db = DMAT(1,2)
    dc = DMAT(2,1)
    dd = DMAT(2,2)
    det = a*d - b*c
    if (type == 0) then
       DMAT(1,1) = (c*db-d*da)/det
       DMAT(1,2) = (b*da-a*db)/det
       DMAT(2,1) = (c*dd-d*dc)/det
       DMAT(2,2) = (b*dc-a*dd)/det
    elseif (type == 1) then
       DMAT(1,1) = (-d*d*da+d*c*db+d*b*dc-b*c*dd)/(det*det)
       DMAT(1,2) = (-a*d*db-b*b*dc+b*d*da+b*a*dd)/(det*det)
       DMAT(2,1) = (-c*c*db-a*d*dc+c*d*da+c*a*dd)/(det*det)
       DMAT(2,2) = (-b*c*da+b*a*dc+a*c*db-a*a*dd)/(det*det)
    end if
    return 
  end subroutine Dinverse
  
  
  subroutine inverse_all(MAT, rank)
    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! Alex G. December 2009
    !-----------------------------------------------------------
    ! input ...
    ! MAT(rank,rank) - array of coefficients for matrix MAT
    ! rank      - dimension
    ! output ...
    ! MAT(rank,rank) - inverse matrix of MAT
    !===========================================================
    implicit none 
    integer, intent(in) :: rank
    complex(fp), intent(inout) :: MAT(rank,rank)
    complex(fp), allocatable :: pivot(:,:), L(:,:), U(:,:)
    complex(fp), allocatable :: b(:), d(:), sol(:)
    complex(fp) :: coeff
    integer :: i, j, k
    allocate(pivot(rank,rank), L(rank,rank), U(rank,rank))
    allocate(b(rank), d(rank), sol(rank))
    ! step 0: initialization for matrices L and U and b
    pivot = MAT
    L = (0.0_fp, 0.0_fp)
    U = (0.0_fp, 0.0_fp)
    b = (0.0_fp, 0.0_fp)
    ! step 1: forward elimination
    do k = 1, rank-1
       do i = k+1, rank
          coeff = pivot(i,k) / pivot(k,k)
          L(i,k) = coeff
          do j = k+1, rank
             pivot(i,j) = pivot(i,j) - coeff*pivot(k,j)
          end do
       end do
    end do
    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient the diagonal elements are 1.0
    do i = 1, rank
       L(i,i) = (1.0_fp, 0.0_fp)
    end do
    ! U matrix is the upper triangular part of MAT
    do j = 1, rank
       do i = 1, j
          U(i,j) = pivot(i,j)
       end do
    end do
    ! Step 3: compute columns of the inverse matrix 
    do k = 1, rank
       b(k) = (1.0_fp, 0.0_fp)
       d(1) = b(1)
       ! Step 3a: Solve Ld=b using the forward substitution
       do i = 2, rank
          d(i) = b(i)
          do j = 1, i-1
             d(i) = d(i) - L(i,j)*d(j)
          end do
       end do
       ! Step 3b: Solve Uxsol=d using the back substitution
       sol(rank) = d(rank) / U(rank,rank)
       do i = rank-1, 1, -1
          sol(i) = d(i)
          do j = rank, i+1, -1
             sol(i) = sol(i) - U(i,j)*sol(j)
          end do
          sol(i) = sol(i) / U(i,i)
       end do
       ! Step 3c: fill the solutions sol(rank) into column k of inverse matrix
       do i = 1, rank
          MAT(i,k) = sol(i)
       end do
       b(k) = (0.0_fp, 0.0_fp)
    end do
    deallocate(pivot, L, U, b, d, sol)
    return
  end subroutine inverse_all


  function DET(MAT, rank)
    !Function to find the determinant of a complex square matrix
    !Description: The subroutine is based on two key points:
    !1. A determinant is unaltered when row operations are performed:
    !   Hence, using this principle,
    !   row operations (column operations would work as well) are used
    !   to convert the matrix into upper traingular form
    !2. The determinant of a triangular matrix is obtained by 
    !   finding the product of the diagonal elements
    implicit none
    integer, intent(in) :: rank
    complex(fp), intent(in) :: MAT(rank,rank)
    complex(fp) :: DET
    complex(fp) :: m, temp
    complex(fp), allocatable :: matrix(:,:)
    integer :: i, j, k, l
    logical :: DetExists = .TRUE.
    allocate(matrix(rank,rank))
    l = 1
    matrix = MAT
    !Convert to upper triangular form
    do k = 1, rank-1
       if (matrix(k,k) == 0) then
          DetExists = .FALSE.
          do i = k+1, rank
             if (matrix(i,k) /= 0) then
                do j = 1, rank
                   temp = matrix(i,j)
                   matrix(i,j)= matrix(k,j)
                   matrix(k,j) = temp
                end do
                DetExists = .TRUE.
                l = -l
                exit
             end if
          end do
          if (DetExists .EQV. .FALSE.) then
             DET = 0
             return
          end if
       end if
       do j = k+1, rank
          m = matrix(j,k)/matrix(k,k)
          do i = k+1, rank
             matrix(j,i) = matrix(j,i) - m*matrix(k,i)
          end do
       end do
    end do  
    !Calculate determinant by finding product of diagonal elements
    DET = l
    do i = 1, rank
       DET = DET * matrix(i,i)
    end do
    deallocate(matrix)
    return
  end function DET


  subroutine interp(yinp, xinp, ninp, y, x, narray, deriv)
    implicit none
    !*************************************************************************************
    ! Piecewise liner interpolation.
    ! Given input arrays x (independent variable) and y (dependent variable),
    ! both of dimension narray, this routine finds, by linear interpolation, 
    ! the value of yinp(x=xinp). 
    ! Array x must be in ascending order.
    ! The flag ierr is returned as -1 if xinp is below the low end of x (an error), 
    ! +1 if xinp is above the high end of x (also an error), or 0 if there was no error.
    !***************************************************************************************
    integer, intent(in) :: ninp
    complex(fp), intent(out) :: yinp(ninp)
    real(fp), intent(in) :: xinp(ninp)
    integer, intent(in) :: narray
    complex(fp), intent(in) :: y(narray)
    real(fp), intent(in) :: x(narray)
    integer, intent(in) :: deriv
    integer :: i, j
    if (deriv == 0) then
       do  i = 1, ninp
          if (xinp(i) < x(1)) then
             yinp(i) = y(1)
          else if (xinp(i) > x(narray)) then
             yinp(i) = y(narray)
          else
             do j = 2, narray
                if (x(j) > xinp(i)) exit
             end do
             yinp(i) = (y(j)-y(j-1))/(x(j)-x(j-1))*(xinp(i)-x(j-1)) + y(j-1)
          end if
       end do
    elseif (deriv == 1) then
       do  i = 1, ninp
          if (xinp(i) < x(1)) then
             yinp(i) = (y(2)-y(1))/(x(2)-x(1))
          else if (xinp(i) > x(narray)) then
             yinp(i) = (y(narray)-y(narray-1))/(x(narray)-x(narray-1))
          else
             do j = 2, narray
                if (x(j) > xinp(i)) exit
             end do
             yinp(i) = (y(j)-y(j-1))/(x(j)-x(j-1))
          end if
       end do
    end if
    return
  end subroutine interp


  pure subroutine ispline(yintp, xintp, nintp, y, x, narray, deriv)
    ! output : yintp = interpolated value(s) at point xintp
    implicit none
    integer, intent(in) :: nintp, narray
    real(fp), intent(in) :: xintp(nintp),  x(narray)
    complex(fp), intent(in) :: y(narray)
    complex(fp), intent(out) :: yintp(nintp)
    integer, intent(in) :: deriv
    complex(fp), allocatable :: b(:), c(:), d(:)
    real(fp) :: incx
    integer :: idx, i, j, k
    allocate(b(narray), c(narray), d(narray))
    call spline(x, y, b, c, d, narray)
    if (deriv == 0) then
       do idx = 1, nintp
          ! if u is ouside the x interval take a boundary value (left or right)
          if(xintp(idx) < x(1)) then
             yintp(idx) = y(1)
          elseif (xintp(idx) > x(narray)) then
             yintp(idx) = y(narray)
          else
             ! binary search for for i, such that x(i) <= u <= x(i+1)
             i = 1
             j = narray+1
             do while (j > i+1)
                k = (i+j)/2
                if(xintp(idx) < x(k)) then
                   j = k
                else
                   i = k
                end if
             end do
             ! evaluate spline interpolation
             incx = xintp(idx) - x(i)
             yintp(idx) = y(i) + incx*(b(i) + incx*(c(i) + incx*d(i)))
          end if
       end do
    elseif (deriv == 1) then
       do idx = 1, nintp
          ! if u is ouside the x interval take a boundary value (left or right)
          if(xintp(idx) < x(1)) then
             yintp(idx) = b(1)
          elseif (xintp(idx) > x(narray)) then
             yintp(idx) = b(narray)
          else
             ! binary search for for i, such that x(i) <= u <= x(i+1)
             i = 1
             j = narray+1
             do while (j > i+1)
                k = (i+j)/2
                if(xintp(idx) < x(k)) then
                   j = k
                else
                   i = k
                end if
             end do
             ! evaluate spline interpolation
             incx = xintp(idx) - x(i)
             yintp(idx) = b(i)+2.0_fp*incx*c(i)+3.0_fp*incx*incx*d(i)
          end if
       end do
    end if
    deallocate(b, c, d)
    return
  end subroutine ispline
    
  
  pure subroutine spline (x, y, b, c, d, narray)
    !======================================================================
    !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,narray
    !  for cubic spline interpolation
    !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
    !  for  x(i) <= x <= x(i+1)
    !----------------------------------------------------------------------
    !  input..
    !  x = the arrays of data point (in strictly increasing order)
    !  y = the arrays of data ordinates
    !  narray = size of the arrays xi and yi (narray>=2)
    !  output..
    !  b, c, d  = arrays of spline coefficients
    !---------------------------------------------------------------------- 
    ! comments ...
    !  using p to denote differentiation
    !  y(i) = s(x(i))
    !  b(i) = sp(x(i))
    !  c(i) = spp(x(i))/2
    !  d(i) = sppp(x(i))/6 (derivative from the right)
    !======================================================================
    implicit none
    integer, intent(in) :: narray
    real(fp), intent(in) :: x(narray)
    complex(fp), intent(in) :: y(narray)
    complex(fp), intent(out) :: b(narray), c(narray), d(narray)
    integer :: i, j, gap
    real(fp) :: h
    gap = narray-1
    ! check input
    if ( narray < 2 ) return
    if ( narray < 3 ) then
       b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
       c(1) = (0.0_fp, 0.0_fp)
       d(1) = (0.0_fp, 0.0_fp)
       b(2) = b(1)
       c(2) = (0.0_fp, 0.0_fp)
       d(2) = (0.0_fp, 0.0_fp)
       return
    end if
    ! step 1: preparation
    ! b = diagonal, d = offdiagonal, c = right hand side
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, gap
       d(i) = x(i+1) - x(i)
       b(i) = 2.0_fp*(d(i-1) + d(i))
       c(i+1) = (y(i+1) - y(i))/d(i)
       c(i) = c(i+1) - c(i)
    end do
    ! step 2: end conditions 
    ! third derivatives at x(1) and x(narray) obtained from divided differences
    b(1) = -d(1)
    b(narray) = -d(narray-1)
    c(1) = (0.0_fp, 0.0_fp)
    c(narray) = (0.0_fp, 0.0_fp)
    if(narray .ne. 3) then
       c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
       c(narray) = c(narray-1)/(x(narray)-x(narray-2)) - c(narray-2)/(x(narray-1)-x(narray-3))
       c(1) = c(1)*d(1)**2/(x(4)-x(1))
       c(narray) = -c(narray)*d(narray-1)**2/(x(narray)-x(narray-3))
    end if
    ! step 3: forward elimination 
    do i = 2, narray
       h = d(i-1)/b(i-1)
       b(i) = b(i) - h*d(i-1)
       c(i) = c(i) - h*c(i-1)
    end do
    ! step 4: back substitution
    c(narray) = c(narray)/b(narray)
    do j = 1, gap
       i = narray-j
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do
    ! step 5: compute spline coefficients
    b(narray) = (y(narray) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0_fp*c(narray))
    do i = 1, gap
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0_fp*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3.0_fp*c(i)
    end do
    c(narray) = 3.0_fp*c(narray)
    d(narray) = d(narray-1)
    return
  end subroutine spline


  pure function line_intg(f0, f1, f2, g0, g1, g2, incx)
  !! integral \int_{xi}^{xi+incx} f(x)*g(x) dx is calculated, where
  !! f0 = f(xi), f1 = f(xi+0.5*incx),               f2 = f(xi+incx),
  !! g0 = g(xi), g1 = \int_{xi}^{xi+incx} g(x) dx,  g2 = g(xi+incx).
    implicit none
    real(fp) :: line_intg
    real(fp), intent(in) :: f0, f1, f2, g0, g1, g2, incx
    line_intg = incx/30.0_fp * ( (4.0_fp*f0-3.0_fp*f1-f2)*g0   &
                           - 3.0_fp*(f0-12.0_fp*f1+f2)*g1/incx &
                           + (-f0-3.0_fp*f1+4.0_fp*f2)*g2 )                        
    return
  end function line_intg


  pure function area_intg(f, g)
  !! integral \int_{qi}^{qi+dq} \int_{pj}^{pj+dp} f(q,p)*g(q,p) dqdp
  !! f(1,1) = f(qi,pj), f(1,2) = f(qi,pj+dp), f(2,1) = f(qi+dq,pj), f(2,2) = f(qi+dq,pj+dp)
  !! g(1,1)=fep(qi,pj,4), g(1,2)=fep(qi,pj+dp,4), g(2,1)=fep(qi+dq,pj,4), g(2,2)=fep(qi+dq,pj+dp,4)  
    implicit none
    real(fp) :: area_intg
    real(fp), intent(in) :: f(2,2), g(2,2,4)
    real(fp) :: f0, f1, f2, g0, g1, g2
    f0 = f(1,1)+f(1,2)
    f2 = f(2,1)+f(2,2)
    f1 = 0.5_fp*(f0+f2)
    g0 = g(1,1,1)
    g1 = g(1,1,2)
    g2 = g(2,1,1)
    area_intg = line_intg(f0, f1, f2, g0, g1, g2, dq)*dp/2.0_fp
    f0 = f(1,1)+2.0_fp*f(1,2)
    f2 = f(2,1)+2.0_fp*f(2,2)
    f1 = 0.5_fp*(f0+f2)
    g0 = 2.0_fp*g(1,1,1)-3.0_fp*g(1,1,3)/dp+g(1,2,1)
    g1 = 2.0_fp*g(1,1,2)-3.0_fp*g(1,1,4)/dp+g(1,2,2)
    g2 = 2.0_fp*g(2,1,1)-3.0_fp*g(2,1,3)/dp+g(2,2,1)
    area_intg = area_intg - line_intg(f0, f1, f2, g0, g1, g2, dq)*dp/3.0_fp
    f0 = f(1,1)+3.0_fp*f(1,2)
    f2 = f(2,1)+3.0_fp*f(2,2)
    f1 = 0.5_fp*(f0+f2)
    g0 = g(1,1,1)-2.0_fp*g(1,1,3)/dp+g(1,2,1)
    g1 = g(1,1,2)-2.0_fp*g(1,1,4)/dp+g(1,2,2)
    g2 = g(2,1,1)-2.0_fp*g(2,1,3)/dp+g(2,2,1)
    area_intg = area_intg + line_intg(f0, f1, f2, g0, g1, g2, dq)*dp/4.0_fp
    return
  end function area_intg
    

  pure function upwind_deriv(f0, f1, f2, g0, g1, flow, incx)
  !! upwind interpolation is used to calculate the derivative df/dx, where
  !! f0 = f(xi-incx), f1 = f(xi), f2 = f(xi+incx),
  !! g0 = \int_{xi-incx}^{xi} f(x) dx, g1 = \int_{xi}^{xi+incx} f(x) dx.
    implicit none
    real(fp) :: upwind_deriv
    real(fp), intent(in) :: f0, f1, f2, g0, g1, flow, incx
    real(fp) :: up_flow, down_flow
    real(fp) :: up_deriv, down_deriv
    up_flow = max(flow, 0.0_fp)
    down_flow = min(flow, 0.0_fp)   
    select case (upwind)
    case (2) 
       up_deriv   = 2.0_fp*(f0+2.0_fp*f1-3.0_fp*g0/incx)/incx
       down_deriv = 2.0_fp*(-f2-2.0_fp*f1+3.0_fp*g1/incx)/incx
    case (3)
       up_deriv = (5.0_fp*f0+8.0_fp*f1-f2)/(6.0_fp*incx)                     &
                + (-3.0_fp*g0+g1)/(incx**2)
       down_deriv = (f0-8.0_fp*f1-5.0_fp*f2)/(6.0_fp*incx)                   &
                  + (-g0+3.0_fp*g1)/(incx**2)
    end select
    if (flow .ne. 0.0_fp) then 
       upwind_deriv = (up_flow*up_deriv + down_flow*down_deriv)/flow
    else
       upwind_deriv = 0.5_fp*(up_deriv + down_deriv)
    end if
    return
  end function upwind_deriv  

  
  subroutine central_diff(d1f, d2f, f0, f1, f2, g0, g1, incx)
  !! first d1f and second d2f central difference of distribution is calculated, where,
  !! output:
  !!          d1f = df/dx,    d2f = d^2f/dx^2
  !! input:
  !!          f0 = f(xi-incx), f1 = f(xi), f2 = f(xi+incx)
  !!          g0  = \int_{xi-incx}^{xi} f(x) dx, g1 = \int_{xi}^{xi+incx} f(x) dx.
    implicit none
    real(fp), intent(out) :: d1f, d2f
    real(fp), intent(in) :: f0, f1, f2, g0, g1, incx
    d1f = -(f2-f0)/(2.0_fp*incx)+2.0_fp*(g1-g0)/(incx*incx)
    d2f = -3.0_fp*(f0+8.0_fp*f1+f2)/(2.0_fp*incx*incx) &
                 + 15.0_fp*(g0+g1)/(2.0_fp*incx*incx*incx)
    return
  end subroutine central_diff


  pure function ctrapez(x, y, narray)
    implicit none
    complex(fp) :: ctrapez
    integer, intent(in) :: narray
    real(fp), intent(in) :: x(narray)
    complex(fp), intent(in) :: y(narray)
    integer :: i
    ctrapez = (0.0_fp, 0.0_fp)
    do i = 1, narray-1
       ctrapez = ctrapez + 0.5_fp*(y(i+1)+y(i))*(x(i+1)-x(i))
    end do
    return
  end function ctrapez


  pure function rtrapez(x, y, narray)
    implicit none
    real(fp) :: rtrapez
    integer, intent(in) :: narray
    real(fp), intent(in) :: x(narray)
    real(fp), intent(in) :: y(narray)
    integer :: i
    rtrapez = 0.0_fp
    do i = 1, narray-1
       rtrapez = rtrapez + 0.5_fp*(y(i+1)+y(i))*(x(i+1)-x(i))
    end do
    return
  end function rtrapez


  function norm_rand(mean, std)
    implicit none
    real(fp), intent(in) :: mean, std
    real(fp) :: norm_rand
    real(fp) :: r1, r2
    call random_number(r1)
    call random_number(r2)
    norm_rand = sqrt(-2.0_fp*dlog(r1))*cos(2.0_fp*pi*r2)
    norm_rand = mean + std*norm_rand
    return
  end function norm_rand

  
  subroutine set_win(win, xlft0, xlft1, xrht1, xrht0, x, narray)
    implicit none
    integer, intent(in) :: narray
    real(fp), intent(out) :: win(narray)
    real(fp), intent(in) :: xlft0, xlft1, xrht1, xrht0
    real(fp), intent(in) :: x(narray)
    real(fp) :: lft(narray), rht(narray)
    real(fp) :: a=35.0_fp, b=84.0_fp, c=70.0_fp, d=20.0_fp
    lft = (x-xlft1) / (xlft1-xlft0)
    rht = (x-xrht1) / (xrht0-xrht1)
    where (x < xlft0 .or. x > xrht0)
       win = 0.0_fp
    elsewhere (x >= xlft0 .and. x <= xlft1) 
       win = 1.0_fp-a*lft**4-b*lft**5-c*lft**6-d*lft**7       
    elsewhere (x <= xrht0 .and. x >= xrht1) 
       win = 1.0_fp-a*rht**4+b*rht**5-c*rht**6+d*rht**7
    elsewhere
       win = 1.0_fp
    end where
    return
  end subroutine set_win

end module chorus_LIB
