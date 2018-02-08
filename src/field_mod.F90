
MODULE field_mod
USE type_mod
USE vmec_mod
USE specifications_mod, ONLY: vmec
IMPLICIT NONE


!=============================================================================!
CONTAINS               !                MODULE SUBPROGRAMS                    !
!=============================================================================!

  SUBROUTINE field(r,p,z,br,bp,bz,dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,&
                  &dbpdz,dbzdr,dbzdp,dbzdz,device)

    INTERFACE

       SUBROUTINE field_w7x(r,p,z,brad,bphi,bzet,&
            &dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,dbpdz,dbzdr,dbzdp,dbzdz)
         USE type_mod
         IMPLICIT NONE
         REAL(DP),INTENT(IN) :: r,p,z
         REAL(DP),INTENT(OUT) :: brad,bphi,bzet,dbrdr,dbrdp,dbrdz
         REAL(DP),INTENT(OUT) :: dbpdr,dbpdp,dbpdz,dbzdr,dbzdp,dbzdz 
       END SUBROUTINE field_w7x

       SUBROUTINE field_ncsx(r,p,z,brad,bphi,bzet,&
            &dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,dbpdz,dbzdr,dbzdp,dbzdz)
         USE type_mod
         IMPLICIT NONE
         REAL(DP),INTENT(IN) :: r,p,z
         REAL(DP),INTENT(OUT) :: brad,bphi,bzet,dbrdr,dbrdp,dbrdz
         REAL(DP),INTENT(OUT) :: dbpdr,dbpdp,dbpdz,dbzdr,dbzdp,dbzdz 
       END SUBROUTINE field_ncsx

    END INTERFACE

    REAL(DP), INTENT(IN):: r,p,z
    REAL(DP), INTENT(OUT) :: br,bp,bz,dbrdr,dbrdz,dbpdr,dbpdz,dbzdr,dbzdz
    REAL(DP), INTENT(OUT) :: dbrdp,dbpdp,dbzdp
    INTEGER, INTENT(IN) :: device

         CALL field_vmec(r,p,z,br,bp,bz,dbrdr,dbrdp,dbrdz,dbpdr,dbpdp,&
                        &dbpdz,dbzdr,dbzdp,dbzdz)

  END SUBROUTINE field

!------------------------------------------------------------------------!

END MODULE field_mod



