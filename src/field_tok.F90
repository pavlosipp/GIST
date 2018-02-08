
SUBROUTINE field_tok(rrr,ppp,zzz,brad,bphi,bzet,&
     &dbrdr,dbrdz,dbpdr,dbpdz,dbzdr,dbzdz)
  USE lagrange_mod, ONLY: indef, plag2d
  USE specifications_mod, ONLY : rmin,rmax,zmin,zmax,rmag
  USE type_mod
  IMPLICIT NONE

  REAL(DP),INTENT(IN) :: rrr,ppp,zzz
  REAL(DP),INTENT(OUT) :: brad,bphi,bzet,dbrdr,dbrdz,dbpdr,dbpdz,dbzdr,dbzdz 

  REAL(DP) :: dr,dz,overdr,overdz,r,z
  INTEGER, PARAMETER :: nr=200,nz=200,values=40000,mp=4
  INTEGER :: l,ii,jj,kk,i,j,m,k,ir,iz,count=0,first_time=0,ierr,icp
  REAL(DP), DIMENSION(nr,nz) :: br,bp,bz
  REAL(DP), DIMENSION(mp,mp) :: Brsten,Bpsten,Bzsten 
  REAL(DP), DIMENSION(mp) :: rsten,zsten 
  REAL(DP) :: polylag,poly1r,poly1z
  INTEGER, DIMENSION(mp) :: indr,indz

  SAVE 

  first_reading: if (first_time == 0) then
     first_time=first_time+1


     OPEN(UNIT=11, FILE='../mesh/mesh_tok', ACTION="READ", STATUS='OLD', ERR=1000)
     READ(11,*) rmag
     DO kk=1,16
        READ(11,*) 
     ENDDO

     DO ir=1,nr
        DO iz=1,nz
           count=count+1
           READ(11,*) r,z,Br(ir,iz),Bp(ir,iz),Bz(ir,iz)
           IF (count == 1) THEN
              rmin=r
              zmin=z
           ELSEIF (count == values) THEN
              rmax=r
              zmax=z
           ENDIF
        ENDDO
        READ(11,*)
     ENDDO

     CLOSE(11)

     dr=(rmax-rmin)/REAL(nr-1)
     dz=(zmax-zmin)/REAL(nz-1)
     
  ENDIF FIRST_READING

!---------------------------- Interpolation --------------------------------!

  dr=(rmax-rmin)/REAL(nr-1)
  dz=(zmax-zmin)/REAL(nz-1)
  overdr=1./dr
  overdz=1./dz

!4-point stencil for each (r,p,z) 

  CALL indef(rrr,rmin,overdr,nr,indr)
  CALL indef(zzz,zmin,overdz,nz,indz)

!Coordinates of stencil points

  DO i=1,mp
     rsten(i)=rmin+(indr(i)-1)*dr
     zsten(i)=zmin+(indz(i)-1)*dz
  ENDDO

  DO l=1,mp
     DO m=1,mp
        Brsten(m,l)=Br(indr(m),indz(l))
        Bpsten(m,l)=Bp(indr(m),indz(l)) 
        Bzsten(m,l)=Bz(indr(m),indz(l))
     ENDDO
  ENDDO

  CALL plag2d(rrr,zzz,Brsten,overdr,overdz,&
       rsten,zsten,polylag,poly1r,poly1z)
  Brad=polylag
  dBrdr=poly1r
  dBrdz=poly1z

  CALL plag2d(rrr,zzz,Bpsten,overdr,overdz,&
       rsten,zsten,polylag,poly1r,poly1z)
  Bphi=polylag
  dBpdr=poly1r
  dBpdz=poly1z

  CALL plag2d(rrr,zzz,Bzsten,overdr,overdz,&
       rsten,zsten,polylag,poly1r,poly1z)
  Bzet=polylag
  dBzdr=poly1r
  dBzdz=poly1z

  RETURN
1000 STOP 'Fatal error(field_tok.f90): missing S-A file'


END SUBROUTINE field_tok







