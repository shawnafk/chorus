!! Source conforms to 
!! Fortran 95 (ISO/IEC 1539-1:1997) "Base language"  + 
!! TR-15580: Floating-point exception handling


!! If your f95 compiler does not support "TR-15580: Floating-point exception
!! handling", i.e. warns about (the lack of) the ieee_exceptions module,
!! please uncomment the following module
!!
MODULE ieee_exceptions
  IMPLICIT NONE

  INTEGER, PARAMETER :: IEEE_OVERFLOW       = 1
  INTEGER, PARAMETER :: IEEE_DIVIDE_BY_ZERO = 2
  INTEGER, PARAMETER :: IEEE_INVALID        = 3
  INTEGER, PARAMETER :: IEEE_UNDERFLOW      = 4
  INTEGER, PARAMETER :: IEEE_INEXACT        = 5

CONTAINS
  ELEMENTAL SUBROUTINE IEEE_SET_FLAG(flag, flag_value)
    IMPLICIT NONE
    INTEGER, INTENT (in) :: flag
    LOGICAL, INTENT (in) :: flag_value

    REAL :: x, y

    IF (flag_value .EQV. .FALSE.) RETURN

    SELECT CASE (flag)
    CASE (IEEE_INVALID)
       x = -1.0
       y = SQRT(x)
    CASE (IEEE_DIVIDE_BY_ZERO)
       x = 0.0
       y = 1.0/x
    END SELECT
  END SUBROUTINE IEEE_SET_FLAG

END MODULE ieee_exceptions


MODULE ad_types
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 1         ! number of independent variables

  !----------------------------------------------------------------------------
  !                You SHOULD NOT NEED to change anything below.
  !----------------------------------------------------------------------------

  ! ***************** WARNING: dpk, spk MUST BE NOT EQUAL. ******************

  INTEGER, PARAMETER :: dpk = KIND(1.d0)  
  ! kind for real for dependent and independent variables

  INTEGER, PARAMETER :: spk = KIND(1.0)   
  ! other kind for real variables in mixed-mode arithmetic

  INTEGER, PARAMETER :: ik  = KIND(1)     
  ! kind for integer variables in mixed-mode arithmetic

  ! the above to support expressions like: f = 2 * x + 3.0 * y + 4.d0

  !----------------------------------------------------------------------------
  !                    You MUST NOT change anything below.
  !----------------------------------------------------------------------------

  INTEGER, PARAMETER :: nhes = (n * (n + 1)) / 2 ! dimension of the Hessian pack

  LOGICAL :: order_is_1or2, order_is_2

  TYPE func 
     SEQUENCE  ! for use with common and equivalence blocks
     REAL (dpk) :: value    = 0.0_dpk
     REAL (dpk) :: x(n)     = 0.0_dpk
     REAL (dpk) :: xx(nhes) = 0.0_dpk
  END TYPE func

END MODULE ad_types

MODULE ad_utilities
  USE ad_types

  IMPLICIT NONE

  INTERFACE independent
     MODULE PROCEDURE indep_scalar, indep_vector
  END INTERFACE independent

CONTAINS 
  SUBROUTINE derivative(order)
    IMPLICIT NONE
    INTEGER, INTENT (in) :: order

    order_is_2    =  order == 2
    order_is_1or2 = (order == 1) .OR. (order_is_2)
  END SUBROUTINE derivative

  SUBROUTINE indep_scalar(i, x, val)
    IMPLICIT NONE
    INTEGER,     INTENT (in)  :: i
    TYPE (func), INTENT (out) :: x
    REAL (dpk),  INTENT (in)  :: val

    IF ((i < 1) .OR. (i > n)) STOP "error in auto_deriv: indep_scalar"

    x%value = val

    IF (order_is_1or2) THEN 
       x%x    = 0.0_dpk
       x%x(i) = 1.0_dpk
    ENDIF

    IF (order_is_2) x%xx = 0.0_dpk
  END SUBROUTINE indep_scalar

  SUBROUTINE indep_vector(x, val)
    IMPLICIT NONE
    TYPE (func), DIMENSION(:), INTENT (out) :: x
    REAL (dpk),  DIMENSION(:), INTENT (in)  :: val

    INTEGER :: i

    IF (SIZE(x) /= n) STOP "error in auto_deriv: indep_vector"

    DO i=1,n
       CALL indep_scalar(i, x(i), val(i))
    ENDDO
  END SUBROUTINE indep_vector

  PURE SUBROUTINE extract(x, val, Dx, DDx)
    IMPLICIT NONE
    TYPE (func), INTENT (in) :: x
    REAL (dpk), INTENT (out) :: val
    REAL (dpk), INTENT (out), OPTIONAL :: Dx(n), DDx(nhes)

    val = x%value
    IF (PRESENT(Dx) .AND. order_is_1or2) Dx  = x%x
    IF (PRESENT(DDx) .AND. order_is_2  ) DDx = x%xx
  END SUBROUTINE extract
END MODULE ad_utilities

MODULE ad_auxiliary
  USE ad_types
  IMPLICIT NONE

  INTEGER, PRIVATE :: i,j
  INTEGER, PARAMETER :: row(nhes) = (/ ( (i, i=j,n), j=1,n) /)
  INTEGER, PARAMETER :: col(nhes) = (/ ( (j, i=j,n), j=1,n) /)

CONTAINS 
  PURE FUNCTION tensor(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in)     :: a, b
    REAL (dpk), DIMENSION (nhes) :: res

    res = a%x(col) * b%x(row)
  END FUNCTION tensor

  ! is_small is used for checking denominators in fractions.
  ELEMENTAL FUNCTION is_small(x) RESULT (res)
    IMPLICIT NONE

    REAL (dpk), INTENT (in) :: x
    LOGICAL :: res

    res = ABS(x) < TINY(x)
  END FUNCTION is_small
END MODULE ad_auxiliary

MODULE ad_assign
  USE ad_types

  IMPLICIT NONE

  INTERFACE ASSIGNMENT (=)
     MODULE PROCEDURE assig_FF, assig_FR, assig_FS, assig_FI
  END INTERFACE ASSIGNMENT (=)

  PRIVATE

  PUBLIC :: ASSIGNMENT (=)

CONTAINS
  ELEMENTAL SUBROUTINE assig_FF(res, a) ! res = a
    IMPLICIT NONE

    TYPE (func), INTENT (out) :: res
    TYPE (func), INTENT (in)  :: a

    res%value = a%value
    IF (order_is_1or2) res%x  = a%x
    IF (order_is_2)    res%xx = a%xx
  END SUBROUTINE assig_FF

  ELEMENTAL SUBROUTINE assig_FR(res, lambda)  ! res = lambda
    IMPLICIT NONE

    TYPE (func), INTENT (inout) :: res
    REAL (dpk),  INTENT (in)    :: lambda

    res%value = lambda  ! what is the correct value for the derivatives?
  END SUBROUTINE assig_FR

  ELEMENTAL SUBROUTINE assig_FS(res, lambda) ! res = lambda
    IMPLICIT NONE

    TYPE (func), INTENT (inout) :: res
    REAL (spk),  INTENT (in)    :: lambda

    res = REAL(lambda, dpk)
  END SUBROUTINE assig_FS

  ELEMENTAL SUBROUTINE assig_FI(res, lambda) ! res = lambda
    IMPLICIT NONE

    TYPE (func),  INTENT (inout) :: res
    INTEGER (ik), INTENT (in)    :: lambda

    res = REAL(lambda, dpk)
  END SUBROUTINE assig_FI

END MODULE ad_assign

MODULE ad_relational
  USE ad_types

  IMPLICIT NONE

  INTERFACE OPERATOR (<)
     MODULE PROCEDURE less_FF
     MODULE PROCEDURE less_FR, less_FS, less_FI, less_RF, less_SF, less_IF
  END INTERFACE

  INTERFACE OPERATOR (<=)
     MODULE PROCEDURE less_equal_FF
     MODULE PROCEDURE less_equal_FR, less_equal_RF
     MODULE PROCEDURE less_equal_FI, less_equal_IF
     MODULE PROCEDURE less_equal_FS, less_equal_SF
  END INTERFACE

  INTERFACE OPERATOR (>)
     MODULE PROCEDURE greater_FF 
     MODULE PROCEDURE greater_FR, greater_FS, greater_FI
     MODULE PROCEDURE greater_RF, greater_SF, greater_IF
  END INTERFACE

  INTERFACE OPERATOR (>=)
     MODULE PROCEDURE greater_equal_FR, greater_equal_RF
     MODULE PROCEDURE greater_equal_FI, greater_equal_IF
     MODULE PROCEDURE greater_equal_FF, greater_equal_SF, greater_equal_FS
  END INTERFACE

  INTERFACE OPERATOR (==)
     MODULE PROCEDURE equal_FR, equal_FI, equal_RF, equal_IF
     MODULE PROCEDURE equal_FF, equal_SF, equal_FS
  END INTERFACE

  INTERFACE OPERATOR (/=)
     MODULE PROCEDURE not_equal_FR, not_equal_FI, not_equal_RF
     MODULE PROCEDURE not_equal_IF, not_equal_FF, not_equal_SF, not_equal_FS
  END INTERFACE

  PRIVATE 
  PUBLIC :: OPERATOR (<), OPERATOR (<=), OPERATOR (>), OPERATOR (>=)
  PUBLIC :: OPERATOR (/=), OPERATOR (==)

CONTAINS
  ELEMENTAL FUNCTION less_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value < lambda
  END FUNCTION less_FR

  ELEMENTAL FUNCTION less_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value < REAL(lambda, dpk)
  END FUNCTION less_FS

  ELEMENTAL FUNCTION less_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value < REAL(lambda, dpk)
  END FUNCTION less_FI

  ELEMENTAL FUNCTION less_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) < a%value
  END FUNCTION less_RF

  ELEMENTAL FUNCTION less_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) < a%value
  END FUNCTION less_SF

  ELEMENTAL FUNCTION less_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) < a%value
  END FUNCTION less_IF

  ELEMENTAL FUNCTION less_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    LOGICAL :: res

    res = a%value < b%value
  END FUNCTION less_FF

  ELEMENTAL FUNCTION less_equal_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value <= lambda
  END FUNCTION less_equal_FR

  ELEMENTAL FUNCTION less_equal_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value <= REAL(lambda,dpk)
  END FUNCTION less_equal_FS

  ELEMENTAL FUNCTION less_equal_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value <= REAL(lambda,dpk)
  END FUNCTION less_equal_FI

  ELEMENTAL FUNCTION less_equal_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = lambda <= a%value
  END FUNCTION less_equal_RF

  ELEMENTAL FUNCTION less_equal_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) <= a%value
  END FUNCTION less_equal_SF

  ELEMENTAL FUNCTION less_equal_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) <= a%value
  END FUNCTION less_equal_IF

  ELEMENTAL FUNCTION less_equal_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    LOGICAL :: res

    res = a%value <= b%value
  END FUNCTION less_equal_FF

  ELEMENTAL FUNCTION greater_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value > lambda
  END FUNCTION greater_FR

  ELEMENTAL FUNCTION greater_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value > REAL(lambda,dpk)
  END FUNCTION greater_FS

  ELEMENTAL FUNCTION greater_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value > REAL(lambda,dpk)
  END FUNCTION greater_FI

  ELEMENTAL FUNCTION greater_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = lambda > a%value
  END FUNCTION greater_RF

  ELEMENTAL FUNCTION greater_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) > a%value
  END FUNCTION greater_SF

  ELEMENTAL FUNCTION greater_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) > a%value
  END FUNCTION greater_IF

  ELEMENTAL FUNCTION greater_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    LOGICAL :: res

    res = a%value > b%value
  END FUNCTION greater_FF

  ELEMENTAL FUNCTION greater_equal_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value >= lambda
  END FUNCTION greater_equal_FR

  ELEMENTAL FUNCTION greater_equal_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value >= REAL(lambda,dpk)
  END FUNCTION greater_equal_FS

  ELEMENTAL FUNCTION greater_equal_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value >= REAL(lambda,dpk)
  END FUNCTION greater_equal_FI

  ELEMENTAL FUNCTION greater_equal_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = lambda >= a%value
  END FUNCTION greater_equal_RF

  ELEMENTAL FUNCTION greater_equal_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) >= a%value
  END FUNCTION greater_equal_SF

  ELEMENTAL FUNCTION greater_equal_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) >= a%value
  END FUNCTION greater_equal_IF

  ELEMENTAL FUNCTION greater_equal_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    LOGICAL :: res

    res = a%value >= b%value
  END FUNCTION greater_equal_FF

  ELEMENTAL FUNCTION equal_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value == lambda
  END FUNCTION equal_FR

  ELEMENTAL FUNCTION equal_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value == REAL(lambda,dpk)
  END FUNCTION equal_FS

  ELEMENTAL FUNCTION equal_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value == REAL(lambda,dpk)
  END FUNCTION equal_FI

  ELEMENTAL FUNCTION equal_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = lambda == a%value
  END FUNCTION equal_RF

  ELEMENTAL FUNCTION equal_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) == a%value
  END FUNCTION equal_SF

  ELEMENTAL FUNCTION equal_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) == a%value
  END FUNCTION equal_IF

  ELEMENTAL FUNCTION equal_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    LOGICAL :: res

    res = a%value == b%value
  END FUNCTION equal_FF

  ELEMENTAL FUNCTION not_equal_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value /= lambda
  END FUNCTION not_equal_FR

  ELEMENTAL FUNCTION not_equal_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value /= REAL(lambda,dpk)
  END FUNCTION not_equal_FS

  ELEMENTAL FUNCTION not_equal_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = a%value /= REAL(lambda,dpk)
  END FUNCTION not_equal_FI

  ELEMENTAL FUNCTION not_equal_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = lambda /= a%value
  END FUNCTION not_equal_RF

  ELEMENTAL FUNCTION not_equal_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) /= a%value
  END FUNCTION not_equal_SF

  ELEMENTAL FUNCTION not_equal_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    LOGICAL :: res

    res = REAL(lambda,dpk) /= a%value
  END FUNCTION not_equal_IF

  ELEMENTAL FUNCTION not_equal_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    LOGICAL :: res

    res = a%value /= b%value
  END FUNCTION not_equal_FF

END MODULE ad_relational


MODULE ad_operator_plus
  USE ad_types

  IMPLICIT NONE

  INTERFACE OPERATOR (+)
     MODULE PROCEDURE unary_plus
     MODULE PROCEDURE add_FF, add_FR, add_FS, add_FI, add_RF, add_SF, add_IF
  END INTERFACE

  PRIVATE 
  PUBLIC :: OPERATOR (+)

CONTAINS
  ELEMENTAL FUNCTION unary_plus(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = a
  END FUNCTION unary_plus

  ELEMENTAL FUNCTION add_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = a
    res%value = lambda + res%value
  END FUNCTION add_RF

  ELEMENTAL FUNCTION add_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = REAL(lambda, dpk) + a
  END FUNCTION add_SF

  ELEMENTAL FUNCTION add_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func),  INTENT (in) :: a
    TYPE (func) :: res

    res = REAL(lambda, dpk) + a
  END FUNCTION add_IF

  ELEMENTAL FUNCTION add_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = lambda + a
  END FUNCTION add_FR

  ELEMENTAL FUNCTION add_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = lambda + a 
  END FUNCTION add_FS

  ELEMENTAL FUNCTION add_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func),  INTENT (in) :: a
    TYPE (func) :: res

    res = lambda + a 
  END FUNCTION add_FI

  ELEMENTAL FUNCTION add_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    TYPE (func) :: res

    res%value = a%value + b%value
    IF (order_is_1or2) res%x = a%x + b%x
    IF (order_is_2)    res%xx = a%xx + b%xx
  END FUNCTION add_FF

END MODULE ad_operator_plus

MODULE ad_operator_minus
  USE ad_types
  USE ad_operator_plus

  IMPLICIT NONE

  INTERFACE OPERATOR (-)
     MODULE PROCEDURE negate
     MODULE PROCEDURE sub_FF, sub_FR, sub_FS, sub_FI, sub_RF, sub_SF, sub_IF
  END INTERFACE

  PRIVATE 
  PUBLIC :: OPERATOR (-)

CONTAINS
  ELEMENTAL FUNCTION negate(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res%value = -a%value
    IF (order_is_1or2) res%x  = -a%x
    IF (order_is_2)    res%xx = -a%xx
  END FUNCTION negate

  ELEMENTAL FUNCTION sub_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = lambda + (-a)
  END FUNCTION sub_RF


  ELEMENTAL FUNCTION sub_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = lambda + (-a)
  END FUNCTION sub_SF


  ELEMENTAL FUNCTION sub_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func),  INTENT (in) :: a
    TYPE (func) :: res

    res = lambda + (-a)
  END FUNCTION sub_IF

  ELEMENTAL FUNCTION sub_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func) :: res

    res = a + (-lambda)
  END FUNCTION sub_FR

  ELEMENTAL FUNCTION sub_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    TYPE (func) :: res

    res = a + (-lambda)
  END FUNCTION sub_FS

  ELEMENTAL FUNCTION sub_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func) :: res

    res = a + (-lambda)
  END FUNCTION sub_FI

  ELEMENTAL FUNCTION sub_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    TYPE (func) :: res

    res = a + (-b)
  END FUNCTION sub_FF

END MODULE ad_operator_minus

MODULE ad_operator_star
  USE ad_types
  USE ad_auxiliary

  IMPLICIT NONE

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE mul_FF, mul_FR, mul_FS, mul_FI, mul_RF, mul_SF, mul_IF
  END INTERFACE

  PRIVATE 
  PUBLIC :: OPERATOR (*)

CONTAINS
  ELEMENTAL FUNCTION mul_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res%value = lambda * a%value
    IF (order_is_1or2) res%x  = lambda * a%x
    IF (order_is_2)    res%xx = lambda * a%xx
  END FUNCTION mul_RF

  ELEMENTAL FUNCTION mul_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = REAL(lambda, dpk) * a
  END FUNCTION mul_SF

  ELEMENTAL FUNCTION mul_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func),  INTENT (in) :: a
    TYPE (func) :: res

    res = REAL(lambda, dpk) * a
  END FUNCTION mul_IF

  ELEMENTAL FUNCTION mul_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func) :: res

    res = lambda * a
  END FUNCTION mul_FR

  ELEMENTAL FUNCTION mul_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    TYPE (func) :: res

    res = lambda * a
  END FUNCTION mul_FS

  ELEMENTAL FUNCTION mul_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func),  INTENT (in) :: a
    TYPE (func) :: res

    res = lambda * a
  END FUNCTION mul_FI

  ELEMENTAL FUNCTION mul_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    TYPE (func) :: res

    res%value = a%value * b%value

    IF (order_is_1or2) res%x = a%value * b%x + b%value * a%x

    IF (order_is_2) THEN
       res%xx = a%value * b%xx + tensor(a,b) + tensor(b,a) + b%value * a%xx
    ENDIF
  END FUNCTION mul_FF

END MODULE ad_operator_star

MODULE ad_operator_slash
  USE ad_types
  USE ad_operator_star
  USE ad_auxiliary

  IMPLICIT NONE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE div_FF, div_FR, div_FS, div_FI, div_RF, div_SF, div_IF 
  END INTERFACE

  PRIVATE 
  PUBLIC :: OPERATOR (/)

CONTAINS
  ELEMENTAL FUNCTION div_FR(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func) :: res

    res = (1.0_dpk / lambda) * a
  END FUNCTION div_FR

  ELEMENTAL FUNCTION div_FS(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    TYPE (func) :: res

    res = (1.0_dpk / REAL(lambda,dpk)) * a  
  END FUNCTION div_FS

  ELEMENTAL FUNCTION div_FI(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func) :: res

    res = (1.0_dpk / REAL(lambda,dpk)) * a  
  END FUNCTION div_FI

  ELEMENTAL FUNCTION div_RF(lambda, a) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res%value = lambda / a%value
    IF (order_is_1or2) res%x = (-res%value / a%value) * a%x

    IF (order_is_2) THEN
       res%xx = -( res%value * a%xx + 2.0_dpk * tensor(res,a) ) / a%value
    ENDIF
  END FUNCTION div_RF

  ELEMENTAL FUNCTION div_SF(lambda, a) RESULT (res)
    IMPLICIT NONE

    REAL (spk), INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = REAL(lambda, dpk) / a
  END FUNCTION div_SF

  ELEMENTAL FUNCTION div_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = REAL(lambda, dpk) / a
  END FUNCTION div_IF

  ELEMENTAL FUNCTION div_FF(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    TYPE (func) :: res

    res = a * (1.0_dpk / b)
  END FUNCTION div_FF

END MODULE ad_operator_slash

MODULE ad_fortran_library
  USE ad_types
  USE ad_relational
  USE ad_assign
  USE ad_operator_plus
  USE ad_operator_minus
  USE ad_operator_star
  USE ad_operator_slash
  USE ad_auxiliary

  USE ieee_exceptions

  IMPLICIT NONE

  INTERFACE abs
     MODULE PROCEDURE abs_
  END INTERFACE

  INTERFACE acos
     MODULE PROCEDURE acos_
  END INTERFACE

  INTERFACE aint
     MODULE PROCEDURE aint_
  END INTERFACE

  INTERFACE anint
     MODULE PROCEDURE anint_
  END INTERFACE

  INTERFACE asin
     MODULE PROCEDURE asin_
  END INTERFACE

  INTERFACE atan
     MODULE PROCEDURE atan_
  END INTERFACE

  INTERFACE atan2
     MODULE PROCEDURE atan2_FF_, atan2_RF_, atan2_FR_, atan2_SF_, atan2_FS_
  END INTERFACE

  INTERFACE ceiling
     MODULE PROCEDURE ceiling_
  END INTERFACE

  INTERFACE cos
     MODULE PROCEDURE cos_
  END INTERFACE

  INTERFACE cosh
     MODULE PROCEDURE cosh_
  END INTERFACE

  INTERFACE digits
     MODULE PROCEDURE digits_
  END INTERFACE

  INTERFACE dim
     MODULE PROCEDURE dim_
  END INTERFACE

  INTERFACE dot_product
     MODULE PROCEDURE dot_product_FF_
     MODULE PROCEDURE dot_product_RF_, dot_product_SF_, dot_product_IF_
     MODULE PROCEDURE dot_product_FR_, dot_product_FS_, dot_product_FI_
  END INTERFACE

  INTERFACE epsilon
     MODULE PROCEDURE epsilon_
  END INTERFACE

  INTERFACE exp
     MODULE PROCEDURE exp_
  END INTERFACE

  INTERFACE exponent
     MODULE PROCEDURE exponent_
  END INTERFACE

  INTERFACE floor
     MODULE PROCEDURE floor_
  END INTERFACE

  INTERFACE fraction
     MODULE PROCEDURE fraction_
  END INTERFACE

  INTERFACE huge
     MODULE PROCEDURE huge_
  END INTERFACE

  INTERFACE int
     MODULE PROCEDURE int_
  END INTERFACE

  INTERFACE kind
     MODULE PROCEDURE kind_
  END INTERFACE

  INTERFACE log
     MODULE PROCEDURE log_
  END INTERFACE

  INTERFACE log10
     MODULE PROCEDURE log10_
  END INTERFACE

  INTERFACE matmul
     MODULE PROCEDURE matmul_FF_12_
     MODULE PROCEDURE matmul_RF_12_, matmul_SF_12_, matmul_IF_12_ 
     MODULE PROCEDURE matmul_FR_12_, matmul_FS_12_, matmul_FI_12_

     MODULE PROCEDURE matmul_FF_21_
     MODULE PROCEDURE matmul_RF_21_, matmul_SF_21_, matmul_IF_21_
     MODULE PROCEDURE matmul_FR_21_, matmul_FS_21_, matmul_FI_21_

     MODULE PROCEDURE matmul_FF_22_
     MODULE PROCEDURE matmul_RF_22_, matmul_SF_22_, matmul_IF_22_
     MODULE PROCEDURE matmul_FR_22_, matmul_FS_22_, matmul_FI_22_
  END INTERFACE

  INTERFACE max
     MODULE PROCEDURE max2_FF_, max2_RF_, max2_FR_, max2_SF_, max2_FS_, max3_
  END INTERFACE

  INTERFACE maxexponent
     MODULE PROCEDURE maxexponent_
  END INTERFACE

  INTERFACE maxloc
     MODULE PROCEDURE maxloc_1, maxloc__dim_1, maxloc__mask_1
     MODULE PROCEDURE maxloc__dim_mask_1
  END INTERFACE

  INTERFACE maxval
     MODULE PROCEDURE maxval_1
  END INTERFACE

  INTERFACE min
     MODULE PROCEDURE min2_FF_, min2_RF_, min2_FR_, min2_SF_, min2_FS_, min3_
  END INTERFACE

  INTERFACE minexponent
     MODULE PROCEDURE minexponent_
  END INTERFACE

  INTERFACE minloc
     MODULE PROCEDURE minloc_1, minloc__dim_1, minloc__mask_1
     MODULE PROCEDURE minloc__dim_mask_1
  END INTERFACE

  INTERFACE minval
     MODULE PROCEDURE minval_1
  END INTERFACE

  INTERFACE mod
     MODULE PROCEDURE mod_FF_, mod_RF_, mod_FR_
  END INTERFACE

  INTERFACE modulo
     MODULE PROCEDURE modulo_FR_, modulo_RF_, modulo_FF_
  END INTERFACE

  INTERFACE nearest
     MODULE PROCEDURE nearest_FF_, nearest_FR_, nearest_RF_
     MODULE PROCEDURE nearest_FS_, nearest_SF_
  END INTERFACE

  INTERFACE nint
     MODULE PROCEDURE nint_
  END INTERFACE

  INTERFACE PRECISION
     MODULE PROCEDURE precision_
  END INTERFACE

  INTERFACE product
     MODULE PROCEDURE product_1
  END INTERFACE

  INTERFACE radix
     MODULE PROCEDURE radix_
  END INTERFACE

  INTERFACE range
     MODULE PROCEDURE range_
  END INTERFACE

  INTERFACE rrspacing
     MODULE PROCEDURE rrspacing_
  END INTERFACE

  INTERFACE scale
     MODULE PROCEDURE scale_
  END INTERFACE

  INTERFACE set_exponent
     MODULE PROCEDURE set_exponent_
  END INTERFACE

  INTERFACE sign
     MODULE PROCEDURE sign_RF_, sign_FR_, sign_FF_
  END INTERFACE

  INTERFACE sin
     MODULE PROCEDURE sin_
  END INTERFACE

  INTERFACE sinh
     MODULE PROCEDURE sinh_
  END INTERFACE

  INTERFACE spacing
     MODULE PROCEDURE spacing_
  END INTERFACE

  INTERFACE sqrt
     MODULE PROCEDURE sqrt_
  END INTERFACE

  INTERFACE sum
     MODULE PROCEDURE sum_1
  END INTERFACE

  INTERFACE tan
     MODULE PROCEDURE tan_
  END INTERFACE

  INTERFACE tanh
     MODULE PROCEDURE tanh_
  END INTERFACE

  INTERFACE tiny
     MODULE PROCEDURE tiny_
  END INTERFACE

  PRIVATE 
  PUBLIC :: abs, acos, aint, anint, asin, atan, atan2, ceiling, cos, cosh
  PUBLIC :: digits, dim, dot_product, epsilon, exp, exponent, floor, fraction
  PUBLIC :: huge, int, log, log10, matmul, max, maxexponent, maxloc, maxval
  PUBLIC :: min, minexponent, minloc, minval, mod, modulo, nearest, nint
  PUBLIC :: precision, product, radix, range, rrspacing, scale, set_exponent
  PUBLIC :: sign, sin, sinh, spacing, sqrt, sum, tan, tanh, tiny, kind

CONTAINS
  ELEMENTAL FUNCTION abs_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = a
    IF (a%value < 0.0_dpk) res = -res

    IF (a%value == 0.0_dpk) CALL IEEE_SET_FLAG(IEEE_INVALID, .TRUE.)
    ! discontinuous function.
  END FUNCTION abs_

  ELEMENTAL FUNCTION acos_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    REAL (dpk)  :: sin_

    res%value = ACOS(a%value)

    IF (order_is_1or2) THEN 
       sin_ = SIN(res%value)
       IF (is_small(sin_)) CALL IEEE_SET_FLAG(IEEE_DIVIDE_BY_ZERO, .TRUE.)

       res%x = -a%x / sin_
    END IF

    IF (order_is_2) res%xx = -(a%xx + a%value * tensor(res, res)) / sin_
  END FUNCTION acos_

  ELEMENTAL FUNCTION aint_(a) RESULT (res) 
    ! optional argument kind cannot be implemented
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk) :: res

    res = AINT(a%value)
  END FUNCTION aint_

  ELEMENTAL FUNCTION anint_(a) RESULT (res) 
    ! optional argument kind cannot be implemented
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk) :: res

    res = ANINT(a%value)
  END FUNCTION anint_

  ELEMENTAL FUNCTION asin_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    REAL (dpk) :: cos_

    res%value = ASIN(a%value)

    IF (order_is_1or2) THEN
       cos_ = COS(res%value)
       IF (is_small(cos_)) CALL IEEE_SET_FLAG(IEEE_DIVIDE_BY_ZERO, .TRUE.)

       res%x = a%x / cos_
    END IF

    IF (order_is_2) res%xx = (a%xx + a%value * tensor(res, res)) / cos_
  END FUNCTION asin_

  ELEMENTAL FUNCTION atan_(a) RESULT (res)
    IMPLICIT NONE
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    REAL (dpk) :: cos2

    res%value = ATAN(a%value)
    IF (order_is_1or2) THEN
       cos2  = 1.0_dpk / (1.0_dpk + a%value**2) ! COS(res%value)**2  
       res%x = cos2 * a%x
    ENDIF

    IF (order_is_2) THEN
!    res%xx = -SIN(2.0_dpk * res%value) * tensor(a,res) + cos2 * a%xx
!    res%xx = -2.0_dpk * SIN(res%value) * COS(res%value) * tensor(a,res) + cos2 * a%xx
!    res%xx = cos2 * (-2.0_dpk * TAN(res%value) * tensor(a,res) + a%xx)
!    res%xx = -2.0_dpk * a%value * cos2 * tensor(a,res) + cos2 * a%xx
       res%xx = -2.0_dpk * a%value * tensor(res,res) + cos2 * a%xx
    ENDIF
  END FUNCTION atan_

  ELEMENTAL FUNCTION atan2_FF_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    TYPE (func) :: res

    IF (is_small(b%value)) CALL IEEE_SET_FLAG(IEEE_DIVIDE_BY_ZERO, .TRUE.)

    res = ATAN(a / b)
    res%value = ATAN2(a%value, b%value)
  END FUNCTION atan2_FF_

  ELEMENTAL FUNCTION atan2_RF_(a_, b) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: a_
    TYPE (func), INTENT (in) :: b
    TYPE (func) :: res

    TYPE (func) :: a

    a = func(a_, 0.0_dpk, 0.0_dpk)
    res = ATAN2(a,b)
  END FUNCTION atan2_RF_

  ELEMENTAL FUNCTION atan2_FR_(a, b_) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: b_
    TYPE (func) :: res

    TYPE (func) :: b 

    b = func(b_, 0.0_dpk, 0.0_dpk)
    res = ATAN2(a,b)
  END FUNCTION atan2_FR_

  ELEMENTAL FUNCTION atan2_SF_(a, b) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  INTENT (in) :: a
    TYPE (func), INTENT (in) :: b
    TYPE (func) :: res

    res = ATAN2(REAL(a,dpk), b)
  END FUNCTION atan2_SF_

  ELEMENTAL FUNCTION atan2_FS_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: b
    TYPE (func) :: res

    res = ATAN2(a, REAL(b,dpk))
  END FUNCTION atan2_FS_

  ELEMENTAL FUNCTION ceiling_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = CEILING(a%value)
  END FUNCTION ceiling_

  ELEMENTAL FUNCTION cos_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    REAL (dpk) :: minus_sin

    res%value = COS(a%value)
    IF (order_is_1or2) THEN
       minus_sin = -SIN(a%value)
       res%x = minus_sin * a%x
    ENDIF

    IF (order_is_2) res%xx = minus_sin * a%xx - res%value * tensor(a,a)
  END FUNCTION cos_

  ELEMENTAL FUNCTION cosh_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    REAL (dpk) :: sinha 

    res%value = COSH(a%value)

    IF (order_is_1or2) THEN 
       sinha = SINH(a%value)
       res%x  = sinha * a%x
    ENDIF
    IF (order_is_2)  res%xx = res%value * tensor(a,a) + sinha * a%xx
  END FUNCTION cosh_

  ELEMENTAL FUNCTION digits_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = DIGITS(a%value)
  END FUNCTION digits_

  ELEMENTAL FUNCTION dim_(x, y) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: x, y
    TYPE (func) :: res

    res = func(0.0_dpk, 0.0_dpk, 0.0_dpk)
    IF (x > y) res = x - y
  END FUNCTION dim_

  PURE FUNCTION dot_product_FF_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a, b
    TYPE (func) :: res

    INTEGER :: i

    res = a(1) * b(1)
    DO i = 2, SIZE(a)
       res = res + a(i) * b(i)
    ENDDO
  END FUNCTION dot_product_FF_

  PURE FUNCTION dot_product_RF_(lambda, b) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  DIMENSION (:), INTENT (in) :: lambda
    TYPE (func), DIMENSION (:), INTENT (in) :: b
    TYPE (func) :: res

    INTEGER :: i

    res%value = DOT_PRODUCT(lambda, b%value)
    IF (order_is_1or2) FORALL (i=1:n)    res%x(i) =DOT_PRODUCT(lambda, b%x(i))
    IF (order_is_2)    FORALL (i=1:nhes) res%xx(i)=DOT_PRODUCT(lambda, b%xx(i))

!!%
!!%    res = 0.0_dpk
!!%    DO i = 1, SIZE(lambda)
!!%       res = res + lambda(i) * b(i)
!!%    ENDDO
  END FUNCTION dot_product_RF_

  PURE FUNCTION dot_product_SF_(a, b) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  DIMENSION (:), INTENT (in) :: a
    TYPE (func), DIMENSION (:), INTENT (in) :: b
    TYPE (func) :: res

    res = DOT_PRODUCT(REAL(a,dpk), b)
  END FUNCTION dot_product_SF_

  PURE FUNCTION dot_product_IF_(a, b) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), DIMENSION (:), INTENT (in) :: a
    TYPE (func),  DIMENSION (:), INTENT (in) :: b
    TYPE (func) :: res

    res = DOT_PRODUCT(REAL(a,dpk), b)
  END FUNCTION dot_product_IF_

  PURE FUNCTION dot_product_FR_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    REAL (dpk),  DIMENSION (:), INTENT (in) :: b
    TYPE (func) :: res

    res = DOT_PRODUCT(b, a)
  END FUNCTION dot_product_FR_

  PURE FUNCTION dot_product_FS_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    REAL (spk),  DIMENSION (:), INTENT (in) :: b
    TYPE (func) :: res

    res = DOT_PRODUCT(b, a)
  END FUNCTION dot_product_FS_

  PURE FUNCTION dot_product_FI_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  DIMENSION (:), INTENT (in) :: a
    INTEGER (ik), DIMENSION (:), INTENT (in) :: b
    TYPE (func) :: res

    res = DOT_PRODUCT(b, a)
  END FUNCTION dot_product_FI_

  ELEMENTAL FUNCTION epsilon_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk) :: res

    res = EPSILON(a%value)
  END FUNCTION epsilon_

  ELEMENTAL FUNCTION exp_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res%value = EXP(a%value)
    IF (order_is_1or2) res%x = res%value * a%x
    IF (order_is_2)   res%xx = res%value * a%xx + tensor(a,res)
  END FUNCTION exp_

  ELEMENTAL FUNCTION exponent_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = EXPONENT(a%value)
  END FUNCTION exponent_

  ELEMENTAL FUNCTION floor_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = FLOOR(a%value)
  END FUNCTION floor_

  ELEMENTAL FUNCTION fraction_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk) :: res

    res = FRACTION(a%value)
  END FUNCTION fraction_

  ELEMENTAL FUNCTION huge_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk) :: res

    res = HUGE(a%value)
  END FUNCTION huge_

  ELEMENTAL FUNCTION int_(a) RESULT (res) 
    ! optional argument kind cannot be implemented
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = INT(a%value)
  END FUNCTION int_

  ELEMENTAL FUNCTION kind_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = KIND(a%value)     ! by default = dpk
  END FUNCTION kind_

  ELEMENTAL FUNCTION log_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res%value = LOG(a%value)
    IF (order_is_1or2) res%x  = a%x / a%value
    IF (order_is_2)    res%xx = a%xx / a%value - tensor(res, res)
  END FUNCTION log_

  ELEMENTAL FUNCTION log10_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = LOG(a) / LOG(10.0_dpk)
  END FUNCTION log10_

  PURE FUNCTION matmul_FF_12_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:),    INTENT (in) :: a
    TYPE (func), DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(b, 2)) :: res

    INTEGER :: i

    res%value = MATMUL(a%value, b%value)

    IF (order_is_1or2) THEN
       FORALL (i = 1:n) &
            res%x(i) = MATMUL(a%x(i), b%value) + MATMUL(a%value, b%x(i)) 
    END IF

    IF (order_is_2) THEN
       FORALL (i=1:nhes) res%xx(i) = &
            MATMUL(a%xx(i), b%value) + &
            MATMUL(a%x(col(i)), b%x(row(i))) + & ! tensor
            MATMUL(a%value, b%xx(i))
    END IF
  END FUNCTION matmul_FF_12_

  PURE FUNCTION matmul_RF_12_(lambda, b) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  DIMENSION (:),    INTENT (in) :: lambda
    TYPE (func), DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(b, 2)) :: res

    INTEGER :: i

    res%value = MATMUL(lambda, b%value)

    IF (order_is_1or2) FORALL(i=1:n)     res%x(i)  = MATMUL(lambda, b%x(i)) 
    IF (order_is_2)    FORALL (i=1:nhes) res%xx(i) = MATMUL(lambda, b%xx(i))
  END FUNCTION matmul_RF_12_

  PURE FUNCTION matmul_SF_12_(lambda, b) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  DIMENSION (:),    INTENT (in) :: lambda
    TYPE (func), DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(b, 2)) :: res

    res = MATMUL(REAL(lambda,dpk), b)
  END FUNCTION matmul_SF_12_

  PURE FUNCTION matmul_IF_12_(lambda, b) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), DIMENSION (:),    INTENT (in) :: lambda
    TYPE (func),  DIMENSION (:, :), INTENT (in) :: b
    TYPE (func),  DIMENSION (SIZE(b, 2)) :: res

    res = MATMUL(REAL(lambda,dpk), b)
  END FUNCTION matmul_IF_12_

  PURE FUNCTION matmul_FR_12_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:),    INTENT (in) :: a
    REAL (dpk),  DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(b, 2)) :: res

    INTEGER :: i

    res%value = MATMUL(a%value, b)

    IF (order_is_1or2) FORALL(i=1:n)    res%x(i)  = MATMUL(a%x(i), b)
    IF (order_is_2)    FORALL(i=1:nhes) res%xx(i) = MATMUL(a%xx(i), b)
  END FUNCTION matmul_FR_12_

  PURE FUNCTION matmul_FS_12_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:),    INTENT (in) :: a
    REAL (spk),  DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(b, 2)) :: res

    res = MATMUL(a, REAL(b,dpk))
  END FUNCTION matmul_FS_12_

  PURE FUNCTION matmul_FI_12_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:),     INTENT (in) :: a
    INTEGER (ik), DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(b, 2)) :: res

    res = MATMUL(a, REAL(b,dpk))
  END FUNCTION matmul_FI_12_

  PURE FUNCTION matmul_FF_21_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:, :), INTENT (in) :: a
    TYPE (func), DIMENSION (:), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1)) :: res

    INTEGER :: i

    res%value = MATMUL(a%value, b%value)

    IF (order_is_1or2) THEN
       FORALL (i = 1:n) &
            res%x(i) = MATMUL(a%x(i), b%value) + MATMUL(a%value, b%x(i)) 
    END IF

    IF (order_is_2) THEN
       FORALL (i=1:nhes) res%xx(i) = &
            MATMUL(a%xx(i), b%value) + &
            MATMUL(a%x(col(i)), b%x(row(i))) + & ! tensor
            MATMUL(a%value, b%xx(i))
    END IF
  END FUNCTION matmul_FF_21_

  PURE FUNCTION matmul_RF_21_(a, b) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  DIMENSION (:, :), INTENT (in) :: a
    TYPE (func), DIMENSION (:),    INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1)) :: res

    INTEGER :: i

    res%value = MATMUL(a, b%value)

    IF (order_is_1or2) FORALL(i=1:n)    res%x(i)  = MATMUL(a, b%x(i))  
    IF (order_is_2)    FORALL(i=1:nhes) res%xx(i) = MATMUL(a, b%xx(i))
  END FUNCTION matmul_RF_21_

  PURE FUNCTION matmul_SF_21_(a, b) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  DIMENSION (:, :), INTENT (in) :: a
    TYPE (func), DIMENSION (:),    INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1)) :: res

    res = MATMUL(REAL(a,dpk), b)
  END FUNCTION matmul_SF_21_

  PURE FUNCTION matmul_IF_21_(a, b) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), DIMENSION (:, :), INTENT (in) :: a
    TYPE (func),  DIMENSION (:),    INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1)) :: res

    res = MATMUL(REAL(a,dpk), b)
  END FUNCTION matmul_IF_21_

  PURE FUNCTION matmul_FR_21_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:, :), INTENT (in) :: a
    REAL (dpk),  DIMENSION (:),    INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1)) :: res

    INTEGER :: i

    res%value = MATMUL(a%value, b)

    IF (order_is_1or2) FORALL (i=1:n)    res%x(i)  = MATMUL(a%x(i), b) 
    IF (order_is_2)    FORALL (i=1:nhes) res%xx(i) = MATMUL(a%xx(i), b) 
  END FUNCTION matmul_FR_21_

  PURE FUNCTION matmul_FS_21_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:, :), INTENT (in) :: a
    REAL (spk),  DIMENSION (:),    INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1)) :: res

    res = MATMUL(a, REAL(b,dpk))
  END FUNCTION matmul_FS_21_

  PURE FUNCTION matmul_FI_21_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  DIMENSION (:, :), INTENT (in) :: a
    INTEGER (ik), DIMENSION (:),    INTENT (in) :: b
    TYPE (func),  DIMENSION (SIZE(a, 1)) :: res

    res = MATMUL(a, REAL(b,dpk))
  END FUNCTION matmul_FI_21_

  PURE FUNCTION matmul_FF_22_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:, :), INTENT (in) :: a
    TYPE (func), DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1), SIZE(b, 2)) :: res

    INTEGER :: i

    res%value = MATMUL(a%value, b%value)

    IF (order_is_1or2) THEN
       FORALL (i=1:n) &
            res%x(i) = MATMUL(a%x(i), b%value) + MATMUL(a%value, b%x(i)) 
    END IF

    IF (order_is_2) THEN
       FORALL (i=1:nhes) res%xx(i) = &
            MATMUL(a%xx(i), b%value) + &
            MATMUL(a%x(col(i)), b%x(row(i))) + & ! tensor
            MATMUL(a%value, b%xx(i))
    END IF
  END FUNCTION matmul_FF_22_

  PURE FUNCTION matmul_RF_22_(a, b) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  DIMENSION (:, :), INTENT (in) :: a
    TYPE (func), DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1), SIZE(b, 2)) :: res

    INTEGER :: i

    res%value = MATMUL(a, b%value)

    IF (order_is_1or2) FORALL (i=1:n)    res%x(i)  = MATMUL(a, b%x(i)) 
    IF (order_is_2)    FORALL (i=1:nhes) res%xx(i) = MATMUL(a, b%xx(i))

  END FUNCTION matmul_RF_22_

  PURE FUNCTION matmul_SF_22_(a, b) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  DIMENSION (:, :), INTENT (in) :: a
    TYPE (func), DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1), SIZE(b, 2)) :: res

    res = MATMUL(REAL(a, dpk), b)
  END FUNCTION matmul_SF_22_

  PURE FUNCTION matmul_IF_22_(a, b) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), DIMENSION (:, :), INTENT (in) :: a
    TYPE (func),  DIMENSION (:, :), INTENT (in) :: b
    TYPE (func),  DIMENSION (SIZE(a, 1), SIZE(b, 2)) :: res

    res = MATMUL(REAL(a, dpk), b)
  END FUNCTION matmul_IF_22_

  PURE FUNCTION matmul_FR_22_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:, :), INTENT (in) :: a
    REAL (dpk),  DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1), SIZE(b, 2)) :: res

    INTEGER :: i

    res%value = MATMUL(a%value, b)

    IF (order_is_1or2) FORALL (i=1:n)    res%x(i)  = MATMUL(a%x(i),  b)
    IF (order_is_2)    FORALL (i=1:nhes) res%xx(i) = MATMUL(a%xx(i), b)
  END FUNCTION matmul_FR_22_

  PURE FUNCTION matmul_FS_22_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:, :), INTENT (in) :: a
    REAL (spk),  DIMENSION (:, :), INTENT (in) :: b
    TYPE (func), DIMENSION (SIZE(a, 1), SIZE(b, 2)) :: res

    res = MATMUL(a, REAL(b,dpk))

  END FUNCTION matmul_FS_22_

  PURE FUNCTION matmul_FI_22_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func),  DIMENSION (:, :), INTENT (in) :: a
    INTEGER (ik), DIMENSION (:, :), INTENT (in) :: b
    TYPE (func),  DIMENSION (SIZE(a, 1), SIZE(b, 2)) :: res

    res = MATMUL(a, REAL(b,dpk))

  END FUNCTION matmul_FI_22_

  ELEMENTAL FUNCTION max2_FF_(a1, a2) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a1
    TYPE (func), INTENT (in) :: a2
    TYPE (func) :: res

    IF (a1 >= a2) THEN 
       res = a1
    ELSE
       res = a2
    END IF
  END FUNCTION max2_FF_

  ELEMENTAL FUNCTION max2_FR_(a1, a2) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a1
    REAL (dpk),  INTENT (in) :: a2
    TYPE (func) :: res

    res = MAX(a1, func(a2, 0.0_dpk, 0.0_dpk))
  END FUNCTION max2_FR_

  ELEMENTAL FUNCTION max2_RF_(a1, a2) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: a1
    TYPE (func), INTENT (in) :: a2
    TYPE (func) :: res

    res = MAX(a2, a1)
  END FUNCTION max2_RF_

  ELEMENTAL FUNCTION max2_FS_(a1, a2) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a1
    REAL (spk),  INTENT (in) :: a2
    TYPE (func) :: res

    res = MAX(a1, REAL(a2,dpk))
  END FUNCTION max2_FS_

  ELEMENTAL FUNCTION max2_SF_(a1, a2) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  INTENT (in) :: a1
    TYPE (func), INTENT (in) :: a2
    TYPE (func) :: res

    res = MAX(a2, a1)
  END FUNCTION max2_SF_

  ELEMENTAL FUNCTION max3_(a1, a2, a3) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a1, a2, a3
    TYPE (func) :: res

    res = MAX(a1, MAX(a2,a3))
  END FUNCTION max3_

  ELEMENTAL FUNCTION maxexponent_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = MAXEXPONENT(a%value)
  END FUNCTION maxexponent_

  PURE FUNCTION maxloc_1(a) RESULT (res) 
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    INTEGER, DIMENSION (SIZE(SHAPE(a))) :: res

    res = MAXLOC(a%value)
  END FUNCTION maxloc_1

  PURE FUNCTION maxloc__dim_1(a, dim) RESULT (res) ! F95 interface
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    INTEGER, INTENT (in) :: dim   ! not optional
    INTEGER, DIMENSION (SIZE(SHAPE(a))) :: res

    res = MAXLOC(a%value, dim)
  END FUNCTION maxloc__dim_1

  PURE FUNCTION maxloc__mask_1(a, mask) RESULT (res) ! F90 interface
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    LOGICAL, DIMENSION (:), INTENT (in) :: mask   ! not optional
    INTEGER, DIMENSION (SIZE(SHAPE(a))) :: res

    res = MAXLOC(a%value, mask)
  END FUNCTION maxloc__mask_1

  PURE FUNCTION maxloc__dim_mask_1(a, dim, mask) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    INTEGER, INTENT (in) :: dim
    LOGICAL, DIMENSION (:), INTENT (in) :: mask
    INTEGER, DIMENSION (SIZE(SHAPE(a))) :: res

    res = MAXLOC(a%value, dim, mask)  
  END FUNCTION maxloc__dim_mask_1

  PURE FUNCTION maxval_1(a, dim, mask) RESULT (res) 
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a 
    INTEGER, INTENT (in), OPTIONAL :: dim
    LOGICAL, DIMENSION (:), INTENT (in), OPTIONAL :: mask ! same rank with a
    TYPE (func) :: res   ! res is scalar

    INTEGER :: ind(SIZE(SHAPE(a)))
    INTEGER :: mydim
    LOGICAL, DIMENSION (SIZE(a)) :: mymask

    mydim = 1
    IF (PRESENT(dim)) mydim = dim

    mymask = .TRUE.
    IF (PRESENT(mask)) mymask = mask

    ind = MAXLOC(a, mydim, mymask)

    res = a(ind(1))
  END FUNCTION maxval_1

  ELEMENTAL FUNCTION min2_FF_(a1, a2) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a1, a2
    TYPE (func) :: res

    IF (a1 <= a2) res = a1
    IF (a2 <= a1) res = a2
  END FUNCTION min2_FF_

  ELEMENTAL FUNCTION min2_FR_(a1, a2) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a1
    REAL (dpk),  INTENT (in) :: a2
    TYPE (func) :: res

    res = MIN(a1, func(a2, 0.0_dpk, 0.0_dpk))
  END FUNCTION min2_FR_

  ELEMENTAL FUNCTION min2_RF_(a1, a2) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: a1
    TYPE (func), INTENT (in) :: a2
    TYPE (func) :: res

    res = MIN(a2, a1)
  END FUNCTION min2_RF_

  ELEMENTAL FUNCTION min2_FS_(a1, a2) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a1
    REAL (spk),  INTENT (in) :: a2
    TYPE (func) :: res

    res = MIN(a1, REAL(a2,dpk))
  END FUNCTION min2_FS_

  ELEMENTAL FUNCTION min2_SF_(a1, a2) RESULT (res)
    IMPLICIT NONE

    REAL (spk),  INTENT (in) :: a1
    TYPE (func), INTENT (in) :: a2
    TYPE (func) :: res

    res = MIN(a2, a1)
  END FUNCTION min2_SF_

  ELEMENTAL FUNCTION min3_(a1, a2, a3) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a1, a2, a3
    TYPE (func) :: res

    res = MIN(a1, MIN(a2,a3))
  END FUNCTION min3_

  ELEMENTAL FUNCTION minexponent_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = MINEXPONENT(a%value)
  END FUNCTION minexponent_

  PURE FUNCTION minloc_1(a) RESULT (res) 
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    INTEGER, DIMENSION (SIZE(SHAPE(a))) :: res

    res = MINLOC(a%value)
  END FUNCTION minloc_1

  PURE FUNCTION minloc__dim_1(a, dim) RESULT (res) ! F95 interface
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    INTEGER, INTENT (in) :: dim   ! not optional
    INTEGER, DIMENSION (SIZE(SHAPE(a))) :: res

    res = MINLOC(a%value, dim)
  END FUNCTION minloc__dim_1

  PURE FUNCTION minloc__mask_1(a, mask) RESULT (res) ! F90 interface
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    LOGICAL, DIMENSION (:), INTENT (in) :: mask   ! not optional
    INTEGER, DIMENSION (SIZE(SHAPE(a))) :: res

    res = MINLOC(a%value, mask)
  END FUNCTION minloc__mask_1

  PURE FUNCTION minloc__dim_mask_1(a, dim, mask) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a
    INTEGER, INTENT (in) :: dim
    LOGICAL, DIMENSION (:), INTENT (in) :: mask
    INTEGER, DIMENSION (SIZE(SHAPE(a))) :: res

    res = MINLOC(a%value, dim, mask)  
  END FUNCTION minloc__dim_mask_1

  PURE FUNCTION minval_1(a, dim, mask) RESULT (res) 
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a 
    INTEGER, INTENT (in), OPTIONAL :: dim
    LOGICAL, DIMENSION (:), INTENT (in), OPTIONAL :: mask ! same rank with a
    TYPE (func) :: res   ! res is scalar

    INTEGER :: ind(SIZE(SHAPE(a)))
    INTEGER :: mydim
    LOGICAL, DIMENSION (SIZE(a)) :: mymask

    mydim = 1
    IF (PRESENT(dim)) mydim = dim

    mymask = .TRUE.
    IF (PRESENT(mask)) mymask = mask

    ind = MINLOC(a, mydim, mymask)

    res = a(ind(1))
  END FUNCTION minval_1

  ELEMENTAL FUNCTION mod_FF_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b

    TYPE (func) :: res
    INTEGER :: div

    res%value = MOD(a%value, b%value)

    IF (order_is_1or2) THEN
       div = (a%value - res%value) / b%value
       res%x = a%x - div * b%x
    END IF
    IF (order_is_2) res%xx = a%xx - div * b%xx
  END FUNCTION mod_FF_

  ELEMENTAL FUNCTION mod_FR_(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda

    TYPE (func) :: res

    res = a
    res%value = MOD(a%value, lambda)
  END FUNCTION mod_FR_

  ELEMENTAL FUNCTION mod_RF_(lambda, b) RESULT (res)  
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: b

    REAL (dpk) :: res

    res = MOD(lambda, b%value)
  END FUNCTION mod_RF_
  ! no mod_sf, mod_fs; a, lambda must be of the same kind

  ELEMENTAL FUNCTION modulo_FF_(a, b) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    TYPE (func) :: res

    INTEGER :: div

    res%value = MODULO(a%value, b%value)

    IF (order_is_1or2) THEN
       div = (a%value - res%value) / b%value
       res%x = a%x - div * b%x
    END IF
    IF (order_is_2) res%xx = a%xx - div * b%xx
  END FUNCTION modulo_FF_

  ELEMENTAL FUNCTION modulo_FR_(a, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda

    TYPE (func) :: res

    res = a
    res%value = MODULO(a%value, lambda)
  END FUNCTION modulo_FR_

  ! no modulo_sf, modulo_fs; a, lambda must be of the same kind 
  ELEMENTAL FUNCTION modulo_RF_(lambd, a) RESULT (res)
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: lambd
    TYPE (func), INTENT (in) :: a

    REAL (dpk) :: res

    res = MODULO(lambd, a%value)
  END FUNCTION modulo_RF_

  ELEMENTAL FUNCTION nearest_FF_(a, b) RESULT (res) 
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    TYPE (func) :: res

    ! the derivatives should be kept. Otherwise, don't use nearest().
    res = a 

    res = NEAREST(a%value, b%value)
  END FUNCTION nearest_FF_

  ELEMENTAL FUNCTION nearest_RF_(a, b) RESULT (res) 
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: a
    TYPE (func), INTENT (in) :: b
    REAL (dpk) :: res

    res = NEAREST(a, b%value)
  END FUNCTION nearest_RF_

  ELEMENTAL FUNCTION nearest_FR_(a, b) RESULT (res) 
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: b
    TYPE (func) :: res

    res = NEAREST(a, func(b, 0.0_dpk, 0.0_dpk))
  END FUNCTION nearest_FR_

  ELEMENTAL FUNCTION nearest_SF_(a, b) RESULT (res) 
    IMPLICIT NONE

    REAL (spk),  INTENT (in) :: a
    TYPE (func), INTENT (in) :: b
    REAL (spk) :: res

    res = NEAREST(a, b%value)
  END FUNCTION nearest_SF_

  ELEMENTAL FUNCTION nearest_FS_(a, b) RESULT (res) 
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: b
    TYPE (func) :: res

    res = NEAREST(a, REAL(b, dpk))
  END FUNCTION nearest_FS_

  ELEMENTAL FUNCTION nint_(a) RESULT (res) 
    ! optional argument kind cannot be implemented
    IMPLICIT NONE
    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = NINT(a%value)
  END FUNCTION nint_

  ELEMENTAL FUNCTION precision_(a) RESULT (res)
    IMPLICIT NONE
    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = PRECISION(a%value)
  END FUNCTION precision_

  PURE FUNCTION product_1(a, dim, mask) RESULT (res)
    IMPLICIT NONE

    TYPE (func), DIMENSION (:), INTENT (in) :: a  
    INTEGER, INTENT (in), OPTIONAL :: dim  
    LOGICAL, DIMENSION (:), INTENT (in), OPTIONAL :: mask
    TYPE (func) :: res      ! res is array with rank = rank(a) - 1

    INTEGER :: i, j, k
    REAL (dpk) :: old_a, old_a_2
    REAL (dpk), DIMENSION(SIZE(a,1)) :: a_  ! copy of a
    INTEGER :: mydim
    LOGICAL, DIMENSION (SIZE(a)) :: mymask

    mydim = 1
    IF (PRESENT(dim)) mydim = dim

    mymask = .TRUE.
    IF (PRESENT(mask)) mymask = mask

    a_ = a%value

    res%value = PRODUCT(a_, mydim, mymask)

    IF (order_is_1or2) THEN 
       DO j = 1, n
          res%x(j) = 0.0_dpk
          DO i = 1, SIZE(a)
             old_a = a_(i)
             a_(i) = a(i)%x(j)
             res%x(j) = res%x(j) + PRODUCT(a_, mydim, mymask)
             a_(i) = old_a
          ENDDO
       ENDDO
    ENDIF

    IF (order_is_2) THEN 
       ! First add to the sum the terms : 
       ! a_1%xx() * a_2 * a_3 ... + a_1 * a_2%xx() * a_3 ... + ...
       DO j = 1, nhes
          res%xx(j) = 0.0_dpk
          DO i = 1, SIZE(a)
             old_a = a_(i)
             a_(i) = a(i)%xx(j)
             res%xx(j) = res%xx(j) + PRODUCT(a_, mydim, mymask)
             a_(i) = old_a
          ENDDO
       ENDDO

       ! Now compute the tensor products
       DO j = 1, nhes
          DO k = 1, SIZE(a)
             old_a = a_(k)
             a_(k) = a(k)%x(col(j))
             
             DO i = 1, SIZE(a)
                IF (i == k) CYCLE
                old_a_2 = a_(i)
                a_(i) = a(i)%x(row(j))
                res%xx(j) = res%xx(j) + PRODUCT(a_, mydim, mymask)
                a_(i) = old_a_2
             ENDDO
             
             a_(k) = old_a
          ENDDO
       ENDDO
    ENDIF
  END FUNCTION product_1

  ELEMENTAL FUNCTION radix_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = RADIX(a%value)
  END FUNCTION radix_

  ELEMENTAL FUNCTION range_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER :: res

    res = RANGE(a%value)
  END FUNCTION range_

  ELEMENTAL FUNCTION rrspacing_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk) :: res

    res = RRSPACING(a%value)
  END FUNCTION rrspacing_

  ELEMENTAL FUNCTION scale_(a, i) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER,     INTENT (in) :: i
    TYPE (func) :: res

    res%value                 = SCALE(a%value, i)
    IF (order_is_1or2) res%x  = SCALE(a%x, i)
    IF (order_is_2)    res%xx = SCALE(a%xx, i)
  END FUNCTION scale_

  ELEMENTAL FUNCTION set_exponent_(a, i) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    INTEGER,     INTENT (in) :: i
    REAL (dpk) :: res

    res = SET_EXPONENT(a%value, i)
  END FUNCTION set_exponent_

  ELEMENTAL FUNCTION sign_FF_(a, b) RESULT (res)
    IMPLICIT NONE
    TYPE (func), INTENT (in) :: a
    TYPE (func), INTENT (in) :: b
    TYPE (func) :: res

    res = ABS(a)
    IF (b < 0.0_dpk) res = -res
  END FUNCTION sign_FF_

  ELEMENTAL FUNCTION sign_RF_(lambda, a) RESULT (res)
    IMPLICIT NONE

    REAL (dpk) , INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    REAL (dpk) :: res

    res = SIGN(lambda, a%value)
  END FUNCTION sign_RF_

  ELEMENTAL FUNCTION sign_FR_(b, lambda) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: b
    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func) :: res

    res =  ABS(b)
    IF (lambda < 0.0_dpk) res = -res
  END FUNCTION sign_FR_
  ! no sign_sf_, sign_fs_; a, lambda must be of the same kind 

  ELEMENTAL FUNCTION sin_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    REAL (dpk) :: cosa 

    res%value = SIN(a%value)
    IF (order_is_1or2) THEN 
       cosa  = COS(a%value) 
       res%x = cosa * a%x
    ENDIF

    IF (order_is_2) res%xx = -res%value * tensor(a,a) + cosa * a%xx
  END FUNCTION sin_

  ELEMENTAL FUNCTION sinh_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    REAL (dpk) :: cosha 

    res%value = SINH(a%value)
    IF (order_is_1or2) THEN 
       cosha = COSH(a%value)
       res%x = cosha * a%x
    ENDIF

    IF (order_is_2) res%xx = res%value * tensor(a,a) + cosha * a%xx
  END FUNCTION sinh_

  ELEMENTAL FUNCTION spacing_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk) :: res

    res = SPACING(a%value)
  END FUNCTION spacing_

  ELEMENTAL FUNCTION sqrt_(a) RESULT (res)
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    IF (a%value == 0.0_dpk) CALL IEEE_SET_FLAG(IEEE_INVALID, .TRUE.)
    ! discontinuous function.

    res%value = SQRT(a%value)
    IF (order_is_1or2) res%x = 0.5_dpk / res%value * a%x
    IF  (order_is_2)  res%xx = (0.5_dpk * a%xx - tensor(res, res)) / res%value
  END FUNCTION sqrt_

  PURE FUNCTION sum_1(a, dim, mask) RESULT (res) 
    IMPLICIT NONE
    TYPE (func), DIMENSION (:), INTENT (in)           :: a
    INTEGER,                    INTENT (in), OPTIONAL :: dim
    LOGICAL,     DIMENSION (:), INTENT (in), OPTIONAL :: mask 
    ! mask has the shape of a
    TYPE (func) :: res  ! res is array with rank = rank(a) - 1

    INTEGER :: i, mydim
    LOGICAL, DIMENSION (SIZE(a)) :: mymask

    mydim = 1
    IF (PRESENT(dim)) mydim = dim

    mymask = .TRUE.
    IF (PRESENT(mask)) mymask = mask

    res%value = SUM(a%value, mydim, mymask)

    IF (order_is_1or2) FORALL(i=1:n)  res%x(i) = SUM(a%x(i),  mydim, mymask)
    IF (order_is_2) FORALL(i=1:nhes) res%xx(i) = SUM(a%xx(i), mydim, mymask)
  END FUNCTION sum_1

  ELEMENTAL FUNCTION tan_(a) RESULT (res)
    IMPLICIT NONE
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    REAL (dpk)  :: sec2

    res%value = TAN(a%value)
    IF (order_is_1or2) THEN
       sec2 = 1.0_dpk + res%value**2        ! 1.0_dpk / COS(a%value)**2
       res%x = sec2 * a%x
    ENDIF

    IF (order_is_2) res%xx = sec2 * a%xx + 2.0_dpk * res%value * tensor(a,res)
  END FUNCTION tan_

  ELEMENTAL FUNCTION tanh_(a) RESULT (res)
    IMPLICIT NONE
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    REAL (dpk)  :: sech2

    res%value = TANH(a%value)
    IF (order_is_1or2) THEN
       sech2 = 1.0_dpk - res%value**2         ! 1.0_dpk / COSH(a%value)**2
       res%x = sech2 * a%x
    ENDIF

    IF (order_is_2) res%xx = sech2 * a%xx - 2.0_dpk * res%value * tensor(a,res)
  END FUNCTION tanh_

  ELEMENTAL FUNCTION tiny_(a) RESULT (res)
    IMPLICIT NONE
    TYPE (func), INTENT (in) :: a
    REAL (dpk) :: res

    res = TINY(a%value)
  END FUNCTION tiny_
END MODULE ad_fortran_library

MODULE ad_operator_power
  USE ad_types
  USE ad_auxiliary
  USE ad_assign

  USE ieee_exceptions

  IMPLICIT NONE

  INTERFACE OPERATOR (**)
     MODULE PROCEDURE power_FF, power_FR, power_FS, power_FI
     MODULE PROCEDURE power_RF, power_SF, power_IF
  END INTERFACE

  PRIVATE 
  PUBLIC :: OPERATOR (**)

CONTAINS
  ELEMENTAL FUNCTION power_FR(a, lambda) RESULT (res)   ! res = a**lambda
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func) :: res

    IF (lambda == 0.0_dpk) THEN
       IF (a%value == 0.0_dpk) CALL IEEE_SET_FLAG(IEEE_INVALID, .TRUE.)

       res%value = 1.0_dpk
       IF (order_is_1or2) res%x  = 0.0_dpk
       IF (order_is_2)    res%xx = 0.0_dpk

       RETURN
    END IF

    IF (lambda == 1.0_dpk) THEN 
       res = a
       RETURN
    ENDIF

    IF (a%value == 0.0_dpk) THEN
       res%value = 0.0_dpk
       IF (order_is_1or2) res%x  = 0.0_dpk
       IF (order_is_2)    res%xx = 0.0_dpk

       RETURN
    ENDIF

    res%value = a%value**lambda
    
    IF (order_is_1or2) res%x  = lambda * res%value / a%value * a%x
    
    IF (order_is_2)    res%xx = ((lambda - 1.0_dpk) * tensor(res,a) &
         + lambda * res%value * a%xx) / a%value

  END FUNCTION power_FR

  ELEMENTAL FUNCTION power_FS(a, lambda) RESULT (res) ! res = a**lambda
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a
    REAL (spk),  INTENT (in) :: lambda
    TYPE (func) :: res

    res = a**REAL(lambda, dpk)
  END FUNCTION power_FS

  ELEMENTAL FUNCTION power_FI(a, lambda) RESULT (res) ! res = a**lambda
    IMPLICIT NONE

    TYPE (func),  INTENT (in) :: a
    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func) :: res

    res = a**REAL(lambda, dpk)
  END FUNCTION power_FI

  ELEMENTAL FUNCTION power_RF(lambda, a) RESULT (res)    ! res = lambda**a
    IMPLICIT NONE

    REAL (dpk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    IF (lambda <= 0.0_dpk) CALL IEEE_SET_FLAG(IEEE_INVALID, .TRUE.)

    res%value = lambda**a%value

    IF (order_is_1or2) res%x  = LOG(lambda) * res%value * a%x
    IF (order_is_2)    res%xx = LOG(lambda) * (tensor(res,a) + res%value * a%xx)
  END FUNCTION power_RF

  ELEMENTAL FUNCTION power_SF(lambda, a) RESULT (res)     ! res = lambda**a
    IMPLICIT NONE

    REAL (spk),  INTENT (in) :: lambda
    TYPE (func), INTENT (in) :: a
    TYPE (func) :: res

    res = REAL(lambda, dpk)**a
  END FUNCTION power_SF

  ELEMENTAL FUNCTION power_IF(lambda, a) RESULT (res)
    IMPLICIT NONE

    INTEGER (ik), INTENT (in) :: lambda
    TYPE (func),  INTENT (in) :: a
    TYPE (func) :: res

    res = REAL(lambda, dpk)**a
  END FUNCTION power_IF

  ELEMENTAL FUNCTION power_FF(a, b) RESULT (res)     ! res = a**b
    USE ad_operator_star
    USE ad_fortran_library, ONLY : exp, log
    IMPLICIT NONE

    TYPE (func), INTENT (in) :: a, b
    TYPE (func) :: res

    res = EXP(LOG(a) * b)
  END FUNCTION power_FF

END MODULE ad_operator_power

MODULE deriv_class
  USE ad_types
  USE ad_utilities
  USE ad_fortran_library
  USE ad_assign
  USE ad_operator_plus
  USE ad_operator_minus
  USE ad_operator_star
  USE ad_operator_slash
  USE ad_operator_power
  USE ad_relational
  IMPLICIT NONE

  PRIVATE :: n, nhes, order_is_2, order_is_1or2, dpk, spk, ik
END MODULE deriv_class
