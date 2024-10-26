C     Last change:  HMH  18 Dec 2006    3:00 pm
C     This file contains the following routines and functions
C
C     REAL FUNCTION RFTWPT    contribution to the discharge potential due to all transient wells
C     SUBROUTINE TWQI         returns the discharge vector components due to all transient wells
C
C -------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFTWPT(CZ)
C
C -------------------------------------------------------------------------
C
C
c     Routine returns the contribution to the discharge potential due to all Theis wells.
c
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RTOLD,RPTOLD,RT,R,RCOEF,RU2,RE1
      COMPLEX(8) CZOLD,CZ
      INCLUDE 'lusys.inc'
      INCLUDE 'twcom.inc'      
      DATA CZOLD,RTOLD,RPTOLD /(1.E30,1.E30),1.E30,0.0/
      RFTWPT=RPTOLD
      IF (CZ.EQ.CZOLD.AND.RTWCT.EQ.RTOLD) RETURN
      CZOLD=CZ
      RTOLD=RTWCT
      RFTWPT=0.0
      DO 10 I=1,NTW
      IF (RTWCT.LE.RTWST(I)) GOTO 10
      IF (RTWTRANS(I).LE.0.0D0) GOTO 10
      RT=RTWCT-RTWST(I)
      R=ABS(CZ-CTWZ(I))
      R=MAX(R,RTWR0(I))
      RCOEF=0.25*RTWSTO(I)/(RTWTRANS(I)*RT)
      RU2=R*R*RCOEF
      RFTWPT=RFTWPT-RTWQ(I)*RE1(RU2)/RPI4
   10 CONTINUE
      RPTOLD=RFTWPT
      RETURN
C
      END
C
C -------------------------------------------------------------------------
C
      SUBROUTINE TWQI (CZ,RQI)
C
C -------------------------------------------------------------------------
C
C
C     Routine returns the discharge vector components due to all
C     transient wells
C      
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RQI,R,RT,RCOEF,RU2
      COMPLEX(8) CZ,CZZ,CQ,CQADD
      INCLUDE 'lusys.inc'
      INCLUDE 'twcom.inc'      
      INCLUDE 'tracom.inc'
      DIMENSION RQI(3)
      IF (NTW.EQ.0) RETURN
      CQ=(0.0,0.0)
      DO 10 I=1,NTW
      IF (RTWCT.LE.RTWST(I)) GOTO 10
      IF (RTWTRANS(I).LE.0.0D0) GOTO 10
      RT=RTWCT-RTWST(I)
      CZZ=CZ-CTWZ(I)
      R=ABS(CZZ)
      R=MAX(R,RTWR0(I))
      RCOEF=0.25*RTWSTO(I)/(RTWTRANS(I)*RT)
      RU2=R*R*RCOEF
      CQADD=-CZZ/R*2.0*RTWQ(I)/(RPI4*R)*EXP(-RU2)
      CQ=CQ+CQADD
  10  CONTINUE
      RQI(1)=RQI(1)+REAL(CQ)
      RQI(2)=RQI(2)+AIMAG(CQ)
      RETURN
      END




