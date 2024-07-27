C     Last change:  HMH  18 Dec 2006    2:20 pm
c     This file contains the following routines and functions
c
c     SUBROUTINE DICZC     generates control points
c     SUBROUTINE DIMAT     generates matrix coefficients
c     SUBROUTINE DIKNO     generates known vector
c     SUBROUTINE DISUB     substitutes solution vector
c
c
c ------------------------------------------------------------------------
c
      SUBROUTINE DICZC(CZI,N,RFAC,CALPH,ITYPE)
c
c ------------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER(4) N,ITYPE,I,ICPT
      REAL(8) RFAC,RTOL
      COMPLEX(8) CZI,CALPH,CZ0
      INCLUDE 'dicom.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION CZI(1),RFAC(4,1),CALPH(1),ITYPE(1)
      IF (NDIS.EQ.0) RETURN
      RTOL=1.0E31
      DO 3 I=1,NDIS
      RTOL=MIN(RTOL,RDIR(I))
   3  CONTINUE
      RTOL=0.1*RTOL      ! 10% of smallest sinkdisc radius
      DO 10 I=1,NDIS
      IF (LDIH(I)) THEN
      N=N+1
      CZ0=CMPLX(RDIZ(1,I),RDIZ(2,I))
      CZI(N)=CZ0+RDICPT(I)
      R3DZA(N)=RDIZ(3,I)
      ITYPE(N)=1
      RFAC(1,N)=1.0D0
      RFAC(4,N)=1.0D0
      CALPH(N)=(0.0,0.0)
C --------------------------------  Check for nearby control points
C                                   RTOL is 10% of smallest sinkdisc radius
      DO 5 ICPT=1,N-1
      IF (ABS(CZI(N)-CZI(ICPT)).LT.RTOL) 
     & WRITE (ILUER,1000) I,ADILAB(I),ICPT,CZI(ICPT)
   5  CONTINUE
      ENDIF
  10  CONTINUE
 1000 FORMAT (' ***WARNING: 3Dsinkdisc ',I3,' with label ',A16,/
     &' may be too close to control point # ',I3,' =',2G11.4)  
      RETURN
      END
c
c ------------------------------------------------------------------------
c
      SUBROUTINE DIMAT (RA,CZI,N,J,RFAC,CALPH,ITYPE)
c
c ------------------------------------------------------------------------
c
C
C     NOTE: discharges normal and parallel to a line are not implemented!
C
      IMPLICIT NONE
      INTEGER(4) N,J,ITYPE,ID,I,IEQS,IEQ
      LOGICAL LNEG
      REAL(8) RA,RFAC,RFDIPC,RFNFDICO
      COMPLEX(8) CZI,CALPH,CZ,CZA,CDICOMC
      INCLUDE 'dicom.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION RA(N,1),CZI(1),CALPH(1),ITYPE(1),RFAC(4,1)
      IF (NDIS.EQ.0) RETURN
      DO 20 ID=1,NDIS
      IF (LDIH(ID)) THEN
      J=J+1
!      WRITE (ILUME,1000) N,J
      DO 10 I=1,N
      CZ=CZI(I)
      R3DZ=R3DZA(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between Cz & CZA, positive to the left when at CZ
C       +3 discharge component normal to a line  (NOT IMPLEMENTED!!!!!!)
C       +4 discharge component parallel to a line (NOT IMPLEMENTED!!!!!!)
C       +5 continuity equation: provide total discharge
C        6 request for zero matrix coefficient
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      GOTO (1,2,3,4,5,6),IEQ
  1   RA(I,J)=RA(I,J)+REAL(CDICOMC(CZ,ID))+RFDIPC(CZ,ID) ! provide potential to CZ
      IF (LNEG) RA(I,J)=RA(I,J)-REAL(CDICOMC(CZA,ID))-RFDIPC(CZA,ID) ! subtract potential at CZA
      GOTO 9
  2   IF (LNEG) THEN
      RA(I,J)=RA(I,J)+RFNFDICO(ID,CZ,CZA) ! provide flow across CZ & CZA
      ELSE
      RA(I,J)=RA(I,J)+AIMAG(CDICOMC(CZ,ID)) ! provide PSI at CZ
      END IF
      GOTO 9
C ----------------NOTE: discharges normal and parallel to CZA not implemented.
  3   CONTINUE
      GOTO 9
  4   CONTINUE
      GOTO 9
C ---------------------------------------------------------------------------
  5   RA(I,J)=RA(I,J)+0.5*RPI2*RDIR(ID)*RDIR(ID)  ! provide total discharge
      GOTO 9
  6   RA(I,J)=0.0
      GOTO 9
  9   RA(I,J)=RA(I,J)*RFAC(1,I)
  10  CONTINUE
      ENDIF
  20  CONTINUE
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)      
      END
c
c ------------------------------------------------------------------------
c
      SUBROUTINE DIKNO (RB,J,CZI)
c
c ------------------------------------------------------------------------
c
C-------------------------------     Note: for ITYPE=1 only!!      
      IMPLICIT NONE
      INTEGER(4) J,ID
      REAL(8) RB,RFPOT,RFPOTH
      COMPLEX(8) CZI,CZ
      INCLUDE 'dicom.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION RB(*),CZI(*)
      IF (NDIS.EQ.0) RETURN
      DO 10 ID=1,NDIS
      IF (LDIH(ID)) THEN
      J=J+1
      CZ=CZI(J)
      R3DZ=R3DZA(J)
      RB(J)=RB(J)+RFPOTH(RDIH(ID),CZ)-RFPOT(CZ)
      ENDIF      
  10  CONTINUE
      RETURN
      END
c
c ------------------------------------------------------------------------
c
      SUBROUTINE DISUB (RB,J)
c
c ------------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER(4) J,ID
      REAL(8) RB
      INCLUDE 'dicom.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'tracom.inc'
      DIMENSION RB(1)
      IF (NDIS.EQ.0) RETURN
      DO 10 ID=1,NDIS
      IF (LDIH(ID)) THEN
      J=J+1
      RDIS(ID)=RDIS(ID)+RB(J)
      ENDIF
  10  CONTINUE
      RETURN
      END

