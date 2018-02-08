
SUBROUTINE rhs(ppp,y,yp)
  USE specifications_mod
  USE field_mod
  IMPLICIT NONE

  !Input arguments:
  !    
  !phi     
  !y(1)=r
  !y(2)=z
  !y(3)=c(1,1)
  !y(4)=c(1,2)
  !y(5)=c(1,3)
  !y(6)=c(2,1)
  !y(7)=c(2,2)
  !y(8)=c(3,3)

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

!Auxilliaries

  df1_dy1 = dBrdR*rrr * Bpm1
  df2_dy1 = dBzdR*rrr * Bpm1
  df1_dy2 = rrr*dBrdZ * Bpm1
  df2_dy2 = rrr*dBzdZ * Bpm1
  df3_dy1 = Bpm1*dBpdR - 1.d0/rrr
  df3_dy2 = Bpm1*dBpdZ
  df1_dy3 = rrr*Bpm1*dBrdp
  df2_dy3 = rrr*Bpm1*dBzdp
  df3_dy3 = dBpdp *Bpm1
 
!RHS

  yp(3) = -(  y(3)*df1_dy1 + y(4)*df2_dy1 + y(5)*df3_dy1 )
  yp(4) = -(  y(3)*df1_dy2 + y(4)*df2_dy2 + y(5)*df3_dy2 )
  yp(5) = -(  y(3)*df1_dy3 + y(4)*df2_dy3 + y(5)*df3_dy3 )
  yp(6) = -(  y(6)*df1_dy1 + y(7)*df2_dy1 + y(8)*df3_dy1 )
  yp(7) = -(  y(6)*df1_dy2 + y(7)*df2_dy2 + y(8)*df3_dy2 )
  yp(8) = -(  y(6)*df1_dy3 + y(7)*df2_dy3 + y(8)*df3_dy3 )

!For D3D/NSTX/JET we invert the toroidal angle

!  IF (device == 3.OR. device==5.OR.device==7)  yp=-yp


END SUBROUTINE rhs




