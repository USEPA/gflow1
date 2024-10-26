C     Last change:  HMH  18 Dec 2006    3:00 pm
C     This file contains the following routines and functions
C
C     SUBROUTINE TWNEAR  check for streamline trace near well and slows down trace or ends trace
C
c
c  ------------------------------------------------------------------------------------------------------
c
      SUBROUTINE TWNEAR (CZ,CZNEW,RZ0,RZNEW,lredo)
c
c  ------------------------------------------------------------------------------------------------------
c
C
C     Routine check for streamline trace near Theis well and slows down trace or ends trace
C
      IMPLICIT NONE
      INTEGER(4) I
      LOGICAL L3DEND,L3DREV,LSRCE,lredo
      REAL(8) RZ0,RZNEW,RDIS,RZONE,RSD0,RSTEP,RADW,
     &     RDISI(3),RTEST(3),RSCALP,RFSP3D,RDISEND
      COMPLEX(8) CZ,CZNEW,CDIS,CDISEND
      INCLUDE 'twcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
C      write (ilume,1000) ntw
C 1000 format (' twnear0: ntw ',i3)
      IF (NTW.EQ.0.or.lredo) RETURN
      CALL GETSTEP (RSD0,RSTEP,L3DEND,L3DREV)
      if (l3dend) RETURN ! streamline ended elsewhere
      DO I=1,NTW
      IF (RTWCT.GT.RTWST(I)) THEN ! well is active
      LSRCE=RTWQ(I).LT.0.0
C ----------------------only end of streamline check when moving toward well.
      IF (L3DREV.AND.LSRCE.OR..NOT.L3DREV.AND..NOT.LSRCE) THEN
        CDIS=CZ-CTWZ(I)
        RDIS=ABS(CDIS)
        RADW=RTWR0(I)
        RZONE=RADW+RSD0
C        write (ilume,1001) RDIS,RSTEP
C 1001   format (' twnear1: RDIS,RSTEP',2(E14.7))
        IF (RDIS.LT.RZONE) THEN  ! starting point CZ near well, check end point
          CDISEND=CZNEW-CTWZ(I)
          RDISEND=ABS(CDISEND)
            IF (RDISEND.LT.RADW) THEN  ! CZNEW inside well radius: stop trace
C            write (ilume,1002) rsd0,rstep,rdis,rzone
C 1002       format (' twnear2: rsd0,rstep,rdis,rzone',4(e14.7))
            L3DEND=.TRUE.
            iFlag=2
            iElementType=5
            IF (LEN(ATWLAB(I)).GT.0)
     &      aElementLabel=ATWLAB(I)
            CALL SETSTEP (RSTEP,L3DEND)
            RETURN
            ENDIF
C ---------------------------------- check if CZNEW overshot the well
            RDISI(1)=REAL(CDIS)
            RDISI(2)=AIMAG(CDIS)
            RDISI(3)=0.0
            RTEST(1)=REAL(CDISEND)
            RTEST(2)=AIMAG(CDISEND)
            RTEST(3)=0.0
            RSCALP=RFSP3D(RDISI,RTEST)
            IF (RSCALP.LT.0.0) THEN     ! CZ and CZNEW on different sides of well: stop trace
C            write (ilume,1003) rscalp
C 1003       format (' twnear3: rscalp',e14.7)
              L3DEND=.TRUE.
              iFlag=2
              iElementType=5
              IF (LEN(ATWLAB(I)).GT.0)
     &        aElementLabel=ATWLAB(I)
              CALL SETSTEP (RSTEP,L3DEND)
              RETURN
            END IF
C            write (ilume,1005) RDIS,RDISEND,RSTEP
C 1005       format (' twnear5: RDIS,RDISEND,RSTEP',3(e14.7))
        ENDIF
      ENDIF
      ENDIF
      END DO
      RETURN
      END SUBROUTINE

