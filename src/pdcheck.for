C     Last change:  HMH  18 Dec 2006    2:57 pm
c     This file contains the following routines and functions
c
c     SUBROUTINE PDERROR     determines max. error at control points
c
c
c
c ---------------------------------------------------------------------------
c
      SUBROUTINE PDERROR (RERMAX)
c
c ---------------------------------------------------------------------------
c
C
C     Routine calculates and reports the maximum error in the
C     boundary conditions specified at discsinks.
C
      IMPLICIT NONE
      INTEGER(4) IHED,IRES,I,K
      REAL(8) RERMAX,RERRH,RERRS,RBASE,RFBASE,RS,RF,
     &        RFHEAD,RC,RFFF,RDUM,RSIG,RSIGC,RSSS
      COMPLEX(8) CZ
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      IF (NPDHD.EQ.0) RETURN
      RERRH=0.0
      RERRS=0.0
      IHED=0
      IRES=0
      DO 10 I=1,NPDHD
      K=IPDHD(I)
      IF (LPDGIV(K)) GOTO 10
      CZ=CPDZ(K)
      rbase=rfbase(cz)
      RS=RPDH(K)-rbase
      RF=RFHEAD(CZ)-rbase
      RC=RPDRES(K)
      IF (RC.EQ.0.0) THEN
C ---------Head specified discsinks without resistance
      RFFF=MAX(ABS(RF+RS),1.0E-10)     ! avoid division by zero error
      RDUM=ABS(RF-RS)/(RFFF/2.0)*100
      RERRH=MAX(RDUM,RERRH)
      IHED=IHED+1
      ELSE
C ---------Head specified discsinks with resistance      
      RSIG=RPDS(K)
      RSIGC=(RF-RS)/RC
      RSSS=MAX(ABS(RSIG+RSIGC),1.0E-10)      ! avoid division by zero error
      RDUM=ABS((RSIG-RSIGC)/(RSSS/2.0))*100
      RERRS=MAX(RDUM,RERRS)
      IRES=IRES+1
      ENDIF
  10  CONTINUE
      IF (IHED.GT.0) then
      if (.not.lucon) WRITE (ILUME,1000) IHED,RERRH
      write (*,1000) ihed,rerrh
      endif
      IF (IRES.GT.0) then
      if (.not.lucon) WRITE (ILUME,2000) IRES,RERRS
      write (*,2000) ires,rerrs
      endif
      RERMAX=MAX(RERMAX,RERRH)
      RERMAX=MAX(RERMAX,RERRS)
      RETURN
 1000 FORMAT (' ',I3,'    discsinks without resistance:   max. error=',
     &        G11.4,' %')
 2000 FORMAT (' ',I3,'    discsinks with resistance:      max. error=',
     &        G11.4,' %')     
      END