C     Last change:  HMH  18 Dec 2006    2:59 pm
c     This file contoins the following routines and functions
c
c     SUBROUTINE PDNEAR      aborts streamline tracing when at sink disc
c
c
c
c
c --------------------------------------------------------------------
c
      SUBROUTINE PDNEAR (CZ,CZNEW,R3DZ,RZNEW,LREDO)
c
c --------------------------------------------------------------------
c
C
C     Routine ends streamline at sink disc.
C
C     INPUT:
C
C     CZ         starting position of current step
C     CZNEW      end position of current step
C     R3DZ       starting elevation of current step
C     RZNEW      end elevation of current step
C
C     OUTPUT:
C
C     CZNEW      (corrected?) end position
C     RZNEW      (corrected?) end elevation
C     RSTEP      (corrected?) stepsize
C     LREDO      TRUE if new point should be calculated with smaller stepsize
C
      IMPLICIT NONE
      INTEGER(4) I,ID,IDCLOSE,iredo
      LOGICAL LB,LT,LSOURCE,LSINK,L3DEND,L3DREV,LREDO
      REAL(8) R3DZ,RZNEW,RHGHT,RFHGHT,RBOT,RTOP,RTAU,RAD,
     &        RS,RSTOP,RSBOT,RDS0,RSTEP,RFBASE
      COMPLEX(8) CZ,CZNEW,CZT
      INCLUDE 'pdcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      DATA IDCLOSE,iredo /0,0/
C     
      IF (NPD.EQ.0.or.lredo) RETURN
      CALL GETSTEP (RDS0,RSTEP,L3DEND,L3DREV)
      if (l3dend)  return ! streamline ended elsewhere
      RHGHT=RFHGHT(CZ)
      RBOT=RFBASE(CZ)
      RTOP=RBOT+RHGHT
      IF (IDCLOSE.NE.0) THEN ! already close to a disc  -----------------(1)
        CZT=CZ-CPDZ(IDCLOSE)
        RTAU=ABS(CZT)
        RAD=RPDR(IDCLOSE)
        IF (RTAU.LT.RAD+RDS0) THEN ! close to same disc --------------(3)
          GOTO 15
        ENDIF
      ENDIF
      IDCLOSE=0
      ID=1
  5   IF (ID.GT.NPD) RETURN
      DO 10 I=ID,NPD
      CZT=CZ-CPDZ(I)
      RTAU=ABS(CZT)
      RAD=RPDR(I)
      IF (RTAU.LT.RAD+RDS0) THEN ! close to another disc -----------(2)
        IDCLOSE=I
        ID=I+1
        GOTO 15
      ENDIF
  10  CONTINUE 
      IDCLOSE=0       
      RETURN      
  15  IF (ABS(R3DZ-RZNEW).GT.0.5*RHGHT) THEN  ! large vertical step, slow down
        RSTEP=0.3*RHGHT
        if (.not.l3dend) LREDO=.TRUE.
        GOTO 20
      ENDIF
      RS=RPDS(IDCLOSE)
      RSBOT=RPDSB(IDCLOSE)
      RSTOP=RS-RSBOT
      LSOURCE=RSTOP.LE.0.0.AND.RSBOT.LE.0.0
      LSINK=RSTOP.GE.0.0.AND.RSBOT.GE.0.0
      IF (L3DREV.AND.LSINK) GOTO 5 ! ------------------------------------(4)
      IF (.NOT.L3DREV.AND.LSOURCE) GOTO 5
      LB=(.NOT.L3DREV.AND.RSBOT.GT.0.0).OR.(L3DREV.AND.RSBOT.LT.0.0)
      IF (LB.AND.RZNEW.LT.RBOT+0.2*RHGHT) THEN ! --------------------------(5)
        IF (RZNEW.LT.RBOT) THEN ! -----------------------------------------(7)
          RSTEP=R3DZ-RBOT
          RZNEW=RBOT
          CZNEW=CZ
          L3DEND=.TRUE.
          iFlag=2
          iElementType=3
          IF (LEN(APDLAB(IDCLOSE)).GT.0)
     &    aElementLabel=APDLAB(IDCLOSE)
        ELSE
           RSTEP=0.1*RHGHT
           if (.not.l3dend) LREDO=.TRUE.
        ENDIF
        GOTO 20
      ELSE
        LT=(.NOT.L3DREV.AND.RSTOP.GT.0.0).OR.(L3DREV.AND.RSTOP.LT.0.0)
        IF (LT.AND.RZNEW.GT.RTOP-0.2*RHGHT) THEN ! -------------------------(6)
          IF (RZNEW.GT.RTOP) THEN ! ----------------------------------------(8)
            RSTEP=RTOP-R3DZ
            RZNEW=RTOP
            CZNEW=CZ
            L3DEND=.TRUE.
            iFlag=2
            iElementType=3
            IF (LEN(APDLAB(IDCLOSE)).GT.0)
     &      aElementLabel=APDLAB(IDCLOSE)
          ELSE
           RSTEP=0.1*RHGHT
           if (.not.l3dend) LREDO=.TRUE.
          ENDIF
          GOTO 20
        ENDIF
        GOTO 5
      ENDIF
  20  CALL SETSTEP (RSTEP,L3DEND)      
      RETURN
      END
