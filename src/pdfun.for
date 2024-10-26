C     Last change:  HMH  18 Dec 2006    2:58 pm
c     This file contains the following routines and functions
c
c     COMPLEX FUNCTION CFPDOM     returns complex potential when outside disc
c     COMPLEX FUNCTION CFPDCO     coefficient function for complex potential when outside disc
c     REAL FUNCTION RFPDPT        real potential for all discs
c     REAL FUNCTION RFPDCO        coefficient function for discharge potential when inside disc
c     SUBROUTINE PDQI             discharge vector for all discs
c
c
c
C
C ----------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFPDOM (CZ)
C
C ----------------------------------------------------------------------
C
C
C     Function returns the complex potential outside discsinks.
C
      IMPLICIT NONE
      INTEGER(4) I
      COMPLEX(8) CZ,CFPDCO
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      CFPDOM=(0.0,0.0)
      IF (NPD.EQ.0) RETURN
      DO 10 I=1,NPD
      CFPDOM=CFPDOM+RPDS(I)*CFPDCO(CZ,I)
  10  CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFPDCO (CZ,IPD)
C
C ----------------------------------------------------------------------
C
C
C     Function returns the complex potential due to discsink IPD and strength 1.
C     CFPDCO is set to (0.0,0.0) when CZ is INSIDE the discsink.
C
      IMPLICIT NONE
      INTEGER(4) IPD
      REAL(8) RAD,RL,R,rgvrefdist
      COMPLEX(8) CZ,CR
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      CFPDCO=(0.0,0.0)
      RAD=RPDR(IPD)
      CR=CZ-CPDZ(IPD)
      rl=rgvrefdist(cpdz(ipd))
      R=ABS(CR)
      IF (R.GE.RAD) CFPDCO=0.5*RAD*RAD*LOG(CR/rl)
      RETURN
      END
C
C ----------------------------------------------------------------------
C      
      
      REAL(8) FUNCTION RFPDPT (CZ)
C
C ----------------------------------------------------------------------
C
C
C     Function return the real potential inside discsinks.
C      
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RFPDCO
      COMPLEX(8) CZ
      INCLUDE 'pdcom.inc'
      INCLUDE 'tracom.inc'
      RFPDPT=0.0
      IF (NPD.EQ.0) RETURN
      DO 10 I=1,NPD
      RFPDPT=RFPDPT+RPDS(I)*RFPDCO(CZ,I)
  10  CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------------------
C
      REAL(8) FUNCTION RFPDCO (CZ,IPD)
C
C ----------------------------------------------------------------------
C
C
C     Function returns the real potential inside discsink IPD.
C     RFPDCO is set to zero if CZ is NOT INSIDE the discsink.
C      
      IMPLICIT NONE
      INTEGER(4) IPD
      REAL(8) R,RAD,rgvrefdist,RL
      COMPLEX(8) CZ,CR
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      RFPDCO=0.0
      CR=CZ-CPDZ(IPD)
      R=ABS(CR)
      RAD=RPDR(IPD)
      rl=rgvrefdist(cpdz(ipd))
      RAD=RPDR(IPD)
      IF (R.LT.RAD) RFPDCO=0.25*(R*R-RAD*RAD)+0.5*rad*rad*LOG(rad/rl)
      RETURN
      END
C
C ----------------------------------------------------------------------
C     
      SUBROUTINE PDQI (CZ,RQI)
C
C ----------------------------------------------------------------------
C
C
C     Discharge vector RQI at CZ, R3DZ for all (2D) discsinks.
C     RQI(3) is set to the vertical spec. discharge component.
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RQI,RQTAU,RQ3,RHGHT,RFHGHT,R3DZ,
     &        RTAU,RX,RY,RADP,RS,R3DZMINBASE,RFBASE
      COMPLEX(8) CZ,CZT
      INCLUDE 'pdcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      DIMENSION RQI(3)
      IF (NPD.EQ.0) RETURN
      RQTAU=0.0
      RQ3=0.0
      RHGHT=RFHGHT(CZ)
      IF (RHGHT.LT.0.001) THEN
       WRITE (ILUER,1000)
       RETURN
      ENDIF
      CALL ZVAL (R3DZ)
      DO 10 I=1,NPD
      CZT=CZ-CPDZ(I)
      RTAU=ABS(CZT)
      RX=REAL(CZT)
      RY=AIMAG(CZT)
      RADP=RPDR(I)
      RS=RPDS(I)
C ------------------------calculate discharge.          
      IF (RTAU.GT.RADP) THEN
C ---------------------------CZ outside discsink.
      RQTAU=-RS*RADP*RADP/(2.0*RTAU)
      RQ3=0.0
      ELSE
C ---------------------------CZ inside discsink.    
      RQTAU=-RS*RTAU/2.0
      R3DZMINBASE=MAX(R3DZ-RFBASE(CZ),0.0)
      RQ3=RS*R3DZMINBASE/RHGHT-RPDSB(I)
      ENDIF
      RTAU=MAX(RTAU,1.E-5)
      RQI(1)=RQI(1)+RX*RQTAU/RTAU
      RQI(2)=RQI(2)+RY*RQTAU/RTAU
      RQI(3)=RQI(3)+RQ3
  10  CONTINUE
      RETURN      
 1000 FORMAT (' ***ERROR in DISCSINK module: zero aquifer heigth.')      
      END
