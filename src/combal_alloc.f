C-----------------------------------------------------------------------
C  29.11.05        LAST MODIFICATION 29.11.05   UHS: tprall_bal DATA   D
C-----------------------------------------------------------------------
C
        subroutine combal_alloc
C
C     FORMERLY  TPRCOM DATA D
C     COMMON BLOCKS FOR TERPSICHORE BALLOONING --LAST MODIFIED 20.11.05
c--------0---------0---------0---------0---------0---------0---------0-c
                                                                        COM0020
        use vbal_param
        use cmodes
        use cfoura
        use ccoef
!
C..   ALLOCATE BLOCKS
!
        allocate  (       mb(ntmp)    ,nb(ntmp) )                       COM0220

                                                                        COM0240
        allocate ( a0mn(ntmp)  ,a1mn(ntmp)  ,a2mn(ntmp)  ,rgmn(ntmp)    COM0250
     &                         ,d0mn(ntmp)  ,d1mn(ntmp)  ,d2mn(ntmp),
     &                          gssmn(ntmp), gstmn(ntmp) ,bsmn(ntmp))

        allocate ( s(mxx)  ,c0(mxx),pax1(mxx),c1(mxx),c2(mxx),c3(mxx)       COM0350
     &  ,fs(mxx),c111(mxx),fq(mxx),y(mxx), gss(mxx), gst(mxx),bs(mxx)
     &  ,cbel(mxx))

        end subroutine combal_alloc
