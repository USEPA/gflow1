C     Last change:  HMH  18 Dec 2006    2:58 pm
c     This file contains the following routines and functions
c
c     SUBROUTINE PDCZC       generates control points
c     SUBROUTINE PDMAT       generates matrix coefficients
c     SUBROUTINE PDKNO       generates known vector
c     SUBROUTINE PDSUB       substitutes solution
c
c
c ---------------------------------------------------------------------------
c
      SUBROUTINE PDCZC(CZI,N,RFAC,CALPH,ITYPE)
c
c ---------------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER(4) N,ITYPE,I,IAD,ICPT
      REAL(8) RFAC,RTOL,RDIS,RC,RS,RFPERM
      COMPLEX(8) CZI,CALPH
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION CZI(1),RFAC(4,1),CALPH(1),ITYPE(1)
      IF (NPDHD.EQ.0) RETURN
      RTOL=1.0E31
      DO 3 I=1,NPDHD
      IAD=IPDHD(I)      
      RTOL=MIN(RTOL,RPDR(IAD))
   3  CONTINUE
      RTOL=0.1*RTOL      ! 10% of smallest discsink radius
      DO 10 I=1,NPDHD
      IAD=IPDHD(I)
      LPDGIV(IAD)=.FALSE.
      RDIS=RPDEP(IAD)
      RC=RPDRES(IAD)
      IF (RPDS(IAD).NE.0.0) THEN
C -------------------------------If depth or resistance are not zero,
C                                check for unsaturated zone underneath discsink.
C                                Except if at first iteration (RSL=0).
        IF (RDIS*RC.NE.0.0) THEN
          RS=-RDIS/RC
        ELSE
          RS=-RFPERM(CPDZ(IAD))
        ENDIF
        IF (RPDS(IAD).LE.RS) THEN
C -------------------------------Limit infiltration when an unsaturated
C                                zone develops underneath the discsink.
C                                Remove discsink from system of equations.
          RPDS(IAD)=RS
          LPDGIV(IAD)=.TRUE.
        ENDIF
      ENDIF
      IF (LPDGIV(IAD)) GOTO 10
      N=N+1
      CZI(N)=CPDZ(IAD)
      ITYPE(N)=1
      RFAC(1,N)=1.0D0
      RFAC(4,N)=1.0D0
      CALPH(N)=(0.0,0.0)
C --------------------------------  Check for nearby control points
C                                   RTOL is 10% of smallest discsink radius
      DO 5 ICPT=1,N-1
      IF (ABS(CZI(N)-CZI(ICPT)).LT.RTOL) 
     & WRITE (ILUER,1000) IAD,APDLAB(IAD),ICPT,CZI(ICPT)
   5  CONTINUE
  10  CONTINUE
 1000 FORMAT (' ***WARNING: discsink ',I3,' with label ',A16,/
     &' may be too close to control point # ',I3,' =',2G11.4)  
      RETURN
      END
c
c ---------------------------------------------------------------------------
c
      SUBROUTINE PDMAT (RA,CZI,N,J,RFAC,CALPH,ITYPE)
c
c ---------------------------------------------------------------------------
c
C
C     NOTE: discharges normal and parallel to a line are not implemented!
C
      IMPLICIT NONE
      INTEGER(4) N,J,ITYPE,II,IAD,I,IEQS,IEQ
      LOGICAL LNEG
      REAL(8) RA,RFAC,RFPDCO,RK,RFPERM,RTHICK,
     &        RFBASE,RFTOP,RC,RSIG,RH,RFHEAD,RFNFPDCO
      COMPLEX(8) CZI,CALPH,CZ,CZA,CFPDCO
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION RA(N,1),CZI(1),CALPH(1),ITYPE(1),RFAC(4,1)
      IF (NPDHD.EQ.0) RETURN
      DO 20 II=1,NPDHD
      IAD=IPDHD(II)
      IF (LPDGIV(IAD)) GOTO 20
      J=J+1
      DO 10 I=1,N
      CZ=CZI(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ and CZA, positive from left to right when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request for zero matrix coefficient
C
      LNEG=IEQS.LT.0
      IEQ=ABS(IEQS)
      GOTO (1,2,3,4,5,6),IEQ
  1   RA(I,J)=RA(I,J)+REAL(CFPDCO(CZ,IAD))+RFPDCO(CZ,IAD) ! provide potential at CZ
      IF (.NOT.LNEG.AND.CZ.EQ.CPDZ(IAD)) THEN ! provide correction term when at own collocation point
C
C  At sink disc there is resistance between surfacewater and aquifer:
C  phi(cz)-phi(specified)/c=qz (qz=specific discharge into surface water)
C      
            RK=RFPERM(CZ)
            RTHICK=RFTOP(CZ)-RFBASE(CZ)
            RC=RPDRES(IAD)
            RSIG=RPDS(IAD)
            IF (RSIG.EQ.0.0) THEN
            RH=RPDH(IAD)-RFBASE(CZ)
            ELSE
            RH=(RFHEAD(CZ)+RPDH(IAD))*0.5-RFBASE(CZ)
            ENDIF
            RH=MIN(RH,RTHICK)
            RA(I,J)=RA(I,J)-RK*RC*RH
      ENDIF
      IF (LNEG) RA(I,J)=RA(I,J)-REAL(CFPDCO(CZA,IAD))-RFPDCO(CZA,IAD)
      GOTO 9
  2   IF (LNEG) THEN
      RA(I,J)=RA(I,J)+RFNFPDCO(IAD,CZ,CZA)  ! provide flow across CZ & CZA (sink density 1)
      ELSE
      RA(I,J)=RA(I,J)+AIMAG(CFPDCO(CZ,IAD)) ! provide PSI at CZ
      END IF
      GOTO 9
C ----------------NOTE: discharges normal and paralleL to unit vector CZA not implemented.
  3   CONTINUE
      GOTO 9
  4   CONTINUE
      GOTO 9
C ---------------------------------------------------------------------------
  5   RA(I,J)=RA(I,J)+0.5*RPI2*RPDR(IAD)*RPDR(IAD)  ! provide total discharge for continuity equation
      GOTO 9
  6   RA(I,J)=0.0
      GOTO 9
  9   RA(I,J)=RA(I,J)*RFAC(1,I)
  10  CONTINUE
  20  CONTINUE
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)      
      END
c
c ---------------------------------------------------------------------------
c
      SUBROUTINE PDKNO (RB,J,CZI)
c
c ---------------------------------------------------------------------------
c
C-------------------------------     Note: for ITYPE=1 only!!      
      IMPLICIT NONE
      INTEGER(4) J,I,IAD
      REAL(8) RB,RFPOT,RFPOTH,RK,RFPERM,RTHICK,
     &        RFBASE,RFTOP,RC,RSIG,RH,RFHEAD
      COMPLEX(8) CZI,CZ
      INCLUDE 'pdcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION RB(*),CZI(*)
      IF (NPDHD.EQ.0) RETURN
      DO 10 I=1,NPDHD
      IAD=IPDHD(I)
      IF (LPDGIV(IAD)) GOTO 10
      J=J+1
      RB(J)=RB(J)+RFPOTH(RPDH(IAD),CPDZ(IAD))-RFPOT(CPDZ(IAD))
      IF (CZI(J).EQ.CPDZ(IAD)) THEN
            CZ=CZI(J)
            RK=RFPERM(CZ)
            RTHICK=RFTOP(CZ)-RFBASE(CZ)
            RC=RPDRES(IAD)
            RSIG=RPDS(IAD)
            IF (RSIG.EQ.0.0) THEN
            RH=RPDH(IAD)-RFBASE(CZ)
            ELSE
            RH=(RFHEAD(CZ)+RPDH(IAD))*0.5-RFBASE(CZ)
            ENDIF
            RH=MIN(RH,RTHICK)
            RB(J)=RB(J)+RK*RC*RH*RSIG
      ENDIF      
  10  CONTINUE
      RETURN
      END
c
c ---------------------------------------------------------------------------
c
      SUBROUTINE PDSUB (RB,J)
c
c ---------------------------------------------------------------------------
c
      INTEGER(4) J,I,IAD
      REAL(8) RB
      INCLUDE 'pdcom.inc'
      INCLUDE 'tracom.inc'
      DIMENSION RB(1)
      IF (NPDHD.EQ.0) RETURN
      DO 10 I=1,NPDHD
      IAD=IPDHD(I)
      IF (LPDGIV(IAD)) GOTO 10
      J=J+1
      RPDS(IAD)=RPDS(IAD)+RB(J)
  10  CONTINUE
      RETURN
      END

