C     Last change:  HMH  18 Dec 2006    2:22 pm
c     This file contains the following routines and functions
c
c     REAL FUNCTION RFNFDI    returns flow across line due to all 3D discs 
c     REAL FUNCTION RFNFDICO  coefficient functions for RFNFDI
c     REAL  FUNCTION RFDPDIO  coefficient function when line is outside the disc
c     REAL  FUNCTION RFDPDII  coefficient function when line is inside the disc 
c
c
c
c
c
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFNFDI(CZ1,CZ2)
c
c --------------------------------------------------------------------------
c
C
C     Function returns the total flow, integrated over the aquifer height, 
c     across a line between CZ1 and CZ2 due to all three-dimensional sink discs.
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RFNFDICO
      COMPLEX(8) CZ1,CZ2
      INCLUDE 'dicom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      RFNFDI=0.0
      IF (NDIS.EQ.0) RETURN
      DO I=1,NDIS
      RFNFDI=RFNFDI+RDIS(I)*RFNFDICO(I,CZ1,CZ2)
      END DO
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFNFDICO(I,CZ1,CZ2)
c
c --------------------------------------------------------------------------
c
C
C     Function returns the total flow, integrated over the aquifer height, 
c     across the line between CZ1 and CZ2 due to sink disc I with sink density 1.
C
      IMPLICIT NONE
      INTEGER(4) I
      LOGICAL LINSECTCIRCLE
      REAL(8) RFDPDII,RFDPDIO,RAD,RR1,RR2,RA,RB,RC,RD,RFSCALAR
      COMPLEX(8) CZ0,CZ1,CZ2,CZ3,CZ4,CR1,CR2,CV12,CV21
      INCLUDE 'dicom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      RFNFDICO = 0.0
      RAD=RDIR(I)
      CZ0=CMPLX(RDIZ(1,I),RDIZ(2,I))
      CR1=CZ1-CZ0
      RR1=ABS(CR1)
      CR2=CZ2-CZ0
      RR2=ABS(CR2)
      IF ( RR1.LE.RAD .AND. RR2.LE.RAD ) THEN ! CZ1 and CZ2 inside
          RFNFDICO=RFDPDII(I,CZ1,CZ2)
          RETURN
      ELSEIF ( RR1.GT.RAD .AND. RR2.LT.RAD) THEN  ! CZ1 outside, CZ2 inside
          CV12=CZ2-CZ1 
          RA=ABS(RFSCALAR(CR1,CV12)/ABS(CV12))
          RB=ABS(CR1)
          RC=SQRT(RB*RB-RA*RA)
          RD=SQRT(RAD*RAD-RC*RC)
          CZ3=(RA-RD)*CV12/ABS(CV12)+CZ1 ! intersection of CZ1 - CZ2 with boundary of disc
          RFNFDICO=RFDPDIO(I,CZ1,CZ3)+RFDPDII(I,CZ3,CZ2)
          RETURN
      ELSEIF ( RR2.GT.RAD .AND. RR1.LT.RAD) THEN  ! CZ1 inside, CZ2 outside
          CV21=CZ1-CZ2
          RA=ABS(RFSCALAR(CR2,CV21)/ABS(CV21))
          RB=ABS(CR2)
          RC=SQRT(RB*RB-RA*RA)
          RD=SQRT(RAD*RAD-RC*RC)
          CZ3=(RA-RD)*CV21/ABS(CV21)+CZ2 ! intersection of CZ1 - CZ2 with boundary of disc
          RFNFDICO=RFDPDII(I,CZ1,CZ3)+RFDPDIO(I,CZ3,CZ2)
          RETURN
      ELSEIF ( RR1.GT.RAD .AND. RR2.GT.RAD ) THEN ! CZ1 and CZ2 outside
          IF (LINSECTCIRCLE(CZ0,RAD,CZ1,CZ2,CZ3,CZ4)) THEN ! line segment intersect disc boundary
            RFNFDICO=
     &      RFDPDIO(I,CZ1,CZ3)+RFDPDII(I,CZ3,CZ4)+RFDPDIO(I,CZ4,CZ2)
          ELSE  ! line segment entirely outside disc
            RFNFDICO=RFDPDIO(I,CZ1,CZ2)
          END IF
      ENDIF
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      REAL(8)  FUNCTION RFDPDIO(I,CZ1,CZ2)
c
c --------------------------------------------------------------------------
c
C
C     Function returns total flow, integrated over the aquifer height, 
c     across the line between CZ1 and CZ2 due to the 3D sink disc I with sink density 1. 
c     The line is entirely OUTSIDE the disc.
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RAD,RFBRANCH
      COMPLEX(8) CZ0,CZ1,CZ2,CR1,CR2
      INCLUDE 'dicom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
c     Compute the delta-psi for the sink disc being on the outside.
C     NOTE: The coefficient function CDICOMC is not used to avoid a point just inside
C     the disc (numerical inaccuracies) to miss a PSI allocation.
      RAD=RDIR(I)
      CZ0=CMPLX(RDIZ(1,I),RDIZ(2,I))
      CR1=CZ1-CZ0
      CR2=CZ2-CZ0
      RFDPDIO=0.5*RAD*RAD*AIMAG(LOG(CR1)-LOG(CR2))
C     Add flow in branch cut of PSI, in case CZ1 - CZ2 crosses branch cut.
      RFDPDIO=RFDPDIO+3.1415926*RAD*RAD*RFBRANCH(CZ0,CZ1,CZ2)
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDPDII(I,CZ1,CZ2)
c
c --------------------------------------------------------------------------
c
C
C     Function returns total flow, integrated over the aquifer height, 
c     across the line between CZ1 and CZ2 due to the 3D sink disc I with sink density 1. 
c     The line is entirely INSIDE the disc.
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RFAREATRIANGLE
      COMPLEX(8) CZ0,CZ1,CZ2
      INCLUDE 'dicom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      CZ0=CMPLX(RDIZ(1,I),RDIZ(2,I))
      RFDPDII=RFAREATRIANGLE (CZ1,CZ2,CZ0)
      RETURN
      END

