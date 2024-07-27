C     Last change:  HMH  18 Dec 2006    2:19 pm
c     This file contains the following routines and functions
c
c     SUBROUTINE DIERROR      returns max. error at control points
c
c
c
c
c
c ------------------------------------------------------------------------
c
      SUBROUTINE DIERROR (RERMAX)
c
c ------------------------------------------------------------------------
c
C
C     Routine calculates and reports the maximum error in the
C     boundary conditions specified at sinkdiscs.
C
      IMPLICIT NONE
      INTEGER(4) I,IHED
      LOGICAL LNOTHING
      REAL(8) RERMAX,RERRH,RS,RFBASE,RF,RFHEAD,RFFF,RDUM
      COMPLEX(8) CZ,CZ0
      INCLUDE 'DICOM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      IF (NDIS.EQ.0) RETURN
      RERRH=0.0
      IHED=0
      LNOTHING=.TRUE.
      DO 10 I=1,NDIS
      IF (LDIH(I)) THEN
      LNOTHING=.FALSE.
      CZ0=CMPLX(RDIZ(1,I),RDIZ(2,I))
      CZ=CZ0+RDICPT(I)
      R3DZ=RDIZ(3,I)
      RS=RDIH(I)-rfbase(cz)! 9/2/99 use heads with respect to aquifer base (independent from terrain elevation)
      RF=RFHEAD(CZ)-rfbase(cz)! 9/2/99 use heads with respect to aquifer base (independent from terrain elevation)
      RFFF=MAX(ABS(RF+RS),1.0E-10)     ! avoid division by zero error
      RDUM=ABS(RF-RS)/(RFFF/2.0)*100
      RDIERR(I)=RDUM     ! save for use in extract routine
      IF (RDUM.GT.RERRH) RERRH=RDUM
      IHED=IHED+1
      ENDIF
  10  CONTINUE
      IF (LNOTHING) RETURN
      IF (IHED.GT.0) then   ! 9/2/99 modified
      IF (.not.lucon) WRITE (ILUME,1000) IHED,RERRH
      write (*,1000) ihed,rerrh
      endif
      RERMAX=MAX(RERMAX,RERRH)
      RETURN
 1000 FORMAT (' ',I3,'      three-dimensional sinkdiscs:  max. error=',
     &        E11.4,' %')
      END

 