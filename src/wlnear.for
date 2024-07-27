C     Last change:  HMH  18 Dec 2006    3:05 pm
C     This file contains the following routines and functions
C
C     SUBROUTINE WLNEAR  check for streamline trace near well and slows down trace or ends trace
C
C
C -----------------------------------------------------------------------
C
      SUBROUTINE WLNEAR (CZ,CZNEW,RZ0,RZNEW,lredo)
C
C -----------------------------------------------------------------------
C
C
C     Routine check for streamline trace near well and slows down trace or ends trace
C     Trace must start outside the RZONE around the well.
C
      IMPLICIT NONE
      INTEGER(4) I
      LOGICAL L3DEND,L3DREV,LSRCE,lredo
      REAL(8) RZ0,RZNEW,RDIS,RZONE,RSD0,RSTEP,RADW,
     &     RDISI(3),RTEST(3),RSCALP,RFSP3D,RDISEND
      COMPLEX(8) CZ,CZNEW,CDIS,CDISEND
      INCLUDE 'wlcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      IF (NWL.EQ.0.or.lredo) RETURN
      CALL GETSTEP (RSD0,RSTEP,L3DEND,L3DREV)
      if (l3dend) RETURN ! streamline ended elsewhere
      DO I=1,NWL
C -----------------------no end of streamline check when moving away from well.
      LSRCE=RWLQ(I).LT.0.0
      IF (L3DREV.AND.LSRCE.OR..NOT.L3DREV.AND..NOT.LSRCE) THEN
        CDIS=CZ-CWLZ(I)
        RDIS=ABS(CDIS)
        RADW=RWLR(I)
        RZONE=RADW+RSD0
        IF (RDIS.LT.RZONE) THEN  ! starting point CZ near well, check point
           CDISEND=CZNEW-CWLZ(I)
           RDISEND=ABS(CDISEND)
            IF (RDISEND.LT.RADW) THEN !  CZNEW inside well radius: stop trace
            L3DEND=.TRUE.
            iFlag=2
            iElementType=1
            IF (LEN(AWLAB(I)).GT.0)
     &      aElementLabel=AWLAB(I)
            CALL SETSTEP (RSTEP,L3DEND)
            RETURN
            ENDIF
C ----------------------- check if CZNEW overshot the well
            RDISI(1)=REAL(CDIS)
            RDISI(2)=AIMAG(CDIS)
            RDISI(3)=0.0
            RTEST(1)=REAL(CDISEND)
            RTEST(2)=AIMAG(CDISEND)
            RTEST(3)=0.0
            RSCALP=RFSP3D(RDISI,RTEST)
            IF (RSCALP.LT.0.0) THEN ! CZ and CZNEW on different sides of well: stop trace
              L3DEND=.TRUE.
              iFlag=2
              iElementType=1
              IF (LEN(AWLAB(I)).GT.0)
     &        aElementLabel=AWLAB(I)
              CALL SETSTEP (RSTEP,L3DEND)
              RETURN
            END IF
        ENDIF
      ENDIF
      END DO
      RETURN
      END SUBROUTINE
