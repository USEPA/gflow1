C     Last change:  HMH  18 Dec 2006    3:01 pm
C     This file contains the following routines and functions
C
C     SUBROUTINE TWTIME      input routine for current time
C     SUBROUTINE TIME        called by TWTIME
C     REAL FUNCTION RTWTIM   reports the current time setting
C     LOGICAL FUNCTION LTWL  TRUE if transient wells are present and time is larger than zero
C
C
C ------------------------------------------------------------------------------
C
      SUBROUTINE TWTIME (CZ)
C
C ------------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      REAL(8) RXC,RYC,RDUM,RVAR
      COMPLEX(8) CZ
      INCLUDE 'dwmn.inc'
      INCLUDE 'twcom.inc'      
      INCLUDE 'tracom.inc'
      INCLUDE 'match.inc'
      INCLUDE 'lusys.inc'
      IF(NTW.EQ.0) RETURN
      RXC=REAL(CZ)
      RYC=AIMAG(CZ)
      IF (LSINGL.AND.LOPEN) THEN
      if (lucon) WRITE (ILUME,1000) RTWCT
C      CALL CURGET (RXC,RYC)         ! NOT AVAILABLE IN BATCH MODE
      ELSE
      if (lucon) WRITE (ILUME,1050) RTWCT
      CALL INLINE
      ENDIF
      RDUM=RVAR(1)
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      ELSE
      RTWCT=RDUM
      WRITE (ILUME,2000) RTWCT
      ENDIF
      RETURN
 1000 FORMAT ('+Time =',G11.4,' Give new time or press <Enter> ')
 1050 FORMAT (' Time =',G11.4,' Give new time or press <Enter> ')
 2000 FORMAT ('+Time =',G11.4,'                                      ')
      END
C
C ------------------------------------------------------------------------------
C
      SUBROUTINE TIME (RT)
C
C ------------------------------------------------------------------------------
C
C
C     Set time for transient solutions, call in main program by TIME command.
C     New routine 12/3/99
C
      IMPLICIT none
      REAL(8) RT
      INCLUDE 'twcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      RTWCT=RT
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RTWTIM ()
C
C ------------------------------------------------------------------------------
C
C
      IMPLICIT NONE
      INCLUDE 'twcom.inc'      
      RTWTIM=RTWCT
      RETURN
      END
C
C ------------------------------------------------------------------------------
C
      LOGICAL FUNCTION LTWL()
C
C ------------------------------------------------------------------------------
C
C
C     Function is TRUE when transient wells are introduced and time is larger than 0
C
      IMPLICIT none
      INCLUDE 'twcom.inc'      
      LTWL = (NTW.GT.0.AND.RTWCT.GT.0.0)
      RETURN
      END
C
C -----------------------------------------------------------------------------
C
      SUBROUTINE SETLOCALTRANSMISSIVITY
C
C
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RH,RK,RFPERM,RFBASE,RTOP,RFTOP
      COMPLEX(8) CZ
      INCLUDE 'lusys.inc'
      INCLUDE 'twcom.inc'
      IF (NTW.GT.0) THEN
      DO 10 I=1,NTW
      CZ=CTWZ(I)
      RK=RFPERM(CZ)
      RTOP=RFTOP(CZ)
      RH=MIN(RTWH0(I),RTOP)
      RTWTRANS(I)=(RH-RFBASE(CZ))*RK
      IF (RTWTRANS(I).LE.0.0D0) THEN
        WRITE (ILUER,1000) ATWLAB(I)
        RTWTRANS(I)=0.0D0
      END IF
   10 CONTINUE
      ENDIF
      RETURN
 1000 FORMAT (' ***ERROR in SETLOCALTRANSMISSIVITY: Zero or negative'
     &        ' transmissivity at transient well: ',A16,/,
     &        ' The well will be ignored.')
      END SUBROUTINE



