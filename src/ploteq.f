C-----------------------------------------------------------------------
C  22.06.88        LAST MODIFICATION 23.11.05     UHS: PLOTEQ DATA D
C-----------------------------------------------------------------------
C
      subroutine PLOTEQ(LMAP)
C 
      use specifications_mod, ONLY:ns_vmec,mpnt_vmec
      use tprb_param
      use triggs
      use vmcfou
      use boofou
C.. Implicits .. 
      implicit none
C
C     WRITES EQUILIBRIUM COEFFICIENTS TO UNIT 17
C     IF LMAP =1: ORIGINAL DATA ARE USED,
C     IF LMAP =2: BC AFTER THE MAPPING ARE USED
C
!      include 'parbal.inc'
C 
C.. Formal Arguments .. 
C.. In/Out Status: Read, Not Written ..
      integer LMAP
C 
C.. Local Scalars .. 
      integer I,L,NI1, MLMNV1
!      include 'tprcom.bal'
C 
C ... Executable Statements ...
C 
C

       NI1=ns_vmec-1 
       MLMNV1=mpnt_vmec

      if (LMAP .eq. 1) then
C
 1000   format (
     &  /' TERPSICHORE EQUILIBRIUM DATA'//,' LMAP = ',i2,' NI1 = ',i3,
     &  ' LMN =',i3,//)
C
!        write (17,1000) LMAP, NI1, LMNV !
                                       !
C
C
        do I = 0,NI1
C
          do L = 1,LMNV
 1001       format (1x,4i5,1p2e25.14)
C
!            write (17,1001) I, L, MX(L), NX(L), FRV(L,I), FZV(L,I)
C
          end do
C
        end do
      else
C
!        write (17,1000) LMAP, NI1, LMNB
C
        do I = 0,NI1
C
          do L = 1,LMNB
C
!            write (17,1001) I, L, MB(L), NB(L), FRI(L,I), FZI(L,I)
C
          end do
C
        end do
      end if
      end
