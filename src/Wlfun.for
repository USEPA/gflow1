C     Last change:  HMH  18 Dec 2006    3:04 pm
C     This file contains the following routines and functions
C
C     SUBROUTINE WELQI         returns the discharge vector due to all (2D) wells
C     COMPLEX FUNCTION CFWLQC  returns the "W-function" due to well i
C     COMPLEX FUNCTION CFWEOM  returns the complex potential due to all (2D) wells   
C     COMPLEX FUNCTION CFWLOMC returns the complex potential for well i with strength 1.  
C
C -------------------------------------------------------------------------
C
      SUBROUTINE WELQI(CZ,RQI)
C
C -------------------------------------------------------------------------
C
C
C     Discharge vector RQI at CZ for all 2D wells.
C     RQI(3) is set to zero.
C
      IMPLICIT NONE
      INTEGER I
      REAL(8) RQI
      COMPLEX(8) CQW,CFWLQC,CZ
      INCLUDE 'WLCOM.INC'
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
      DIMENSION RQI(3)
      IF (NWL.EQ.0) RETURN
      CQW=(0.0,0.0)
      DO I=1,NWL
      CQW=CQW+RWLQ(I)*CFWLQC(CZ,I)
      END DO
      RQI(1)=RQI(1)+REAL(CQW)
      RQI(2)=RQI(2)-AIMAG(CQW)
      RETURN
      END
C
C -------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFWLQC (CZ,IW)
C
C -------------------------------------------------------------------------
C
C
C     Function returns the W-function at CZ due to well IW with strength 1.
C     Note: When exactly at the center of the well, the well is ignored (no
C     contribution to the discharge).
      IMPLICIT NONE
      INTEGER(4) IW
      REAL(8) RDIS
      COMPLEX(8) CZ,CDIS
      INCLUDE 'WLCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      CFWLQC=(0.0,0.0)
      CDIS=CZ-CWLZ(IW)
      RDIS=ABS(CDIS)
      IF (RDIS.lt.1.0E-10) RETURN       ! at the well center, ignore the well.
      IF (RDIS.LT.RWLR(IW)) THEN
      CDIS=CDIS/RDIS*RWLR(IW)
      ENDIF
      CFWLQC=-1.0/CDIS/RPI2
      RETURN
      END
C
C -------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFWEOM (CZ)
C
C -------------------------------------------------------------------------
C
C     Function returns the complex potential due to all (2D) wells
      IMPLICIT NONE
      INTEGER(4) IW
      COMPLEX(8) CZ,CFWLOMC
      INCLUDE 'WLCOM.INC'
      INCLUDE 'TRACOM.INC'
      CFWEOM=(0.0,0.0)
      IF (NWL.EQ.0) RETURN      
      DO 10 IW=1,NWL
      CFWEOM=CFWEOM+RWLQ(IW)*CFWLOMC(CZ,IW)
  10  CONTINUE
      RETURN
      END
C
C -------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFWLOMC (CZ,IW)
C
C -------------------------------------------------------------------------
C
C
C     Function returns the complex potential at CZ for well IW with strength 1.
C
      IMPLICIT NONE
      INTEGER(4) IW
      REAL(8) RDIS,RADW,RGVREFDIST
      COMPLEX(8) CZ,CDIS,CZR
      INCLUDE 'WLCOM.INC'
      INCLUDE 'TRACOM.INC'
      CFWLOMC=(0.0,0.0)
      RADW=RWLR(IW)      
      CDIS=CZ-CWLZ(IW)
      RDIS=ABS(CDIS)
      IF (RDIS.LT.1.0E-10) THEN
      CDIS=CMPLX(RADW,0.0)
      ELSE
      IF (RDIS.LT.RADW) CDIS=CDIS/RDIS*RADW
      ENDIF
      CZR=CDIS/RGVREFDIST(CWLZ(IW))    ! to force contribution to reference point equal to zero
      CFWLOMC=LOG(CZR)/RPI2
      RETURN
      END

