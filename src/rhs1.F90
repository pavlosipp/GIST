
SUBROUTINE rhs1(ppp,y,yp)
  USE specifications_mod
  USE field_mod
  IMPLICIT NONE

  !Input arguments:
  !    
  !phi     
  !y(1)=r
  !y(2)=z

  REAL(DP), dimension(:), INTENT(OUT) :: yp
  REAL(DP), dimension(:), INTENT(IN)  :: y
  REAL(DP), INTENT(IN) :: ppp
 
  REAL(DP) :: Bpm1,df1_dy1,df2_dy1,df1_dy2,df2_dy2
  REAL(DP) :: df3_dy1,df3_dy2,df1_dy3,df2_dy3,df3_dy3
  REAL(DP) :: rrr,zzz

!Coords 
  
  rrr = y(1)
  zzz = y(2)

  CALL field(rrr,ppp,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,&
            &dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ,device)

  Bpm1 = 1.d0/Bp

  yp(1) = Br*rrr*Bpm1
  yp(2) = Bz*rrr*Bpm1

!For D3D/NSTX/JET we invert the toroidal angle

!  IF (device == 3.OR. device==5.OR.device==7)  yp=-yp


END SUBROUTINE rhs1




