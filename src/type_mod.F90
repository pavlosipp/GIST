
MODULE type_mod


!Constants

  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197d0

  INTEGER, PARAMETER :: ni = 4

!------------ NOT USED ACTUALLY !!! ----------------------
  TYPE POINT_3D
    REAL (DP) :: x
    REAL (DP) :: y
    REAL (DP) :: z
  END TYPE

END MODULE type_mod
