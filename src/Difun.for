C     Last change:  HMH  18 Dec 2006    2:20 pm
c     This file contains the following routines and functions
c
c     REAL FUNCTION RFDIPT     contribution to disch. pot. of all (3D) discs
c     REAL FUNCTION RFDIGP     same as RFDIPT, but for filling grids
c     SUBROUTINE DISCQ2D       generates discharge vector due to all (3D) discs (in 2D zone)
c          ENTRY DISCQ3D       generates discharge vector in 3D zone
c     REAL FUNCTION RFDIPC     coefficient function for disc + images
c     COMPLEX FUNCTION CDIOM   complex potential outside 3D zone
c     COMPLEX FUNCTION CDICOMC coefficient functions for complex potential
c     REAL FUNCTION RFDISC     coefficient function for disch. pot. for single disc
c     REAL FUNCTION RFDIOM     solid angle due to a disc
c     SUBROUTINE DIQIC         coefficient function for discharge vector due to a single disc
c     SUBROUTINE DIQTZ         returns radial and vertical discharge components
c
c
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDIPT (CZ)
c
c --------------------------------------------------------------------------
c
C
C     Potential contribution of all sink discs.
C     NOTE: Function is not meant for filling grids; use RFDIGP.
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RDUM,RFDIPC
      COMPLEX(8) CZ
      INCLUDE 'DICOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      RFDIPT=0.0
      IF (NDIS.EQ.0) RETURN
      DO 10 I=1,NDIS
      rdum=rfdipc(cz,i)
      RFDIPT=RFDIPT+RDIS(I)*rdum
c      write (iluout,1001) i,rdis(i),rdum
c 1001 format (' rfdipt1: i,rdis,rdum ',i3,1x,2g14.7)      
  10  CONTINUE
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDIGP (CZ)
c
c --------------------------------------------------------------------------
c
C
C     Potential contribution to disc number IEL,
C     which is passed through common.
C
      IMPLICIT NONE
      REAL(8) RFDIPC
      COMPLEX(8) CZ
      INCLUDE 'DICOM.INC'
      INCLUDE 'TRACOM.INC'
      RFDIGP=RDIS(IEL)*RFDIPC(CZ,IEL)
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      SUBROUTINE DISCQ2D (CZ,RQI)
c
c --------------------------------------------------------------------------
c
C
C     Routine adds the discharge vector due to all discs
C     to the vector <RQI>.
C
      IMPLICIT NONE
      INTEGER(4) I,J
      LOGICAL L2D
      REAL(8) RQI,RQCI,RS
      COMPLEX(8) CZ
      INCLUDE 'DICOM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION RQI(3),RQCI(3)
      L2D=.TRUE.
      GOTO 1
C     ----------------------
      ENTRY DISCQ3D (CZ,RQI)
C     ----------------------
      L2D=.FALSE.
  1   IF (NDIS.EQ.0) RETURN
      DO 20 I=1,NDIS
      CALL DIQIC (CZ,RQCI,I,L2D)
      RS=RDIS(I)
      DO 10 J=1,3
      RQI(J)=RQI(J)+RS*RQCI(J)
  10  CONTINUE
  20  CONTINUE
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDIPC (CZ,ID)
c
c --------------------------------------------------------------------------
c
C
C     Pot. coeff. for SINK disc ID and its images.
C     Note: it is efficient to call this function
C     with varying R3DZ and constant CZ.
C     Potential definition: k*phi*H
C
      IMPLICIT NONE
      INTEGER(4) ID,J
      LOGICAL LDITOP,LDIBOT
      REAL(8) RZTOP,RZBOT,RTAUI,RFTOP,RFBASE,RX,RY,RZ,RHH,
     &        RDISX,RDISY,RDISZ,RDISA,RSCO,RFSP3D,RZP,RZM,
     &        RFDISC,RF3DSL,RL,RGVREFDIST,RZP1,RZM1
      COMPLEX(8) CZ,CZ0
      INCLUDE 'DICOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION RTAUI(3)
      DATA RTAUI /3*0.0/
      IF (ID.GT.NDIS) THEN
      WRITE (ILUME,1000) ID,NDIS
      RFDIPC=0.0
      RETURN
      ENDIF
      RZTOP=RFTOP(CZ)
      RZBOT=RFBASE(CZ)
      RX=REAL(CZ)
      RY=AIMAG(CZ)
      RZ=R3DZ
      RHH=2.0*R3DH
      RDISX=RDIZ(1,ID)
      RDISY=RDIZ(2,ID)
      RDISZ=RDIZ(3,ID)
      RDISR=RDIR(ID)
      RDISA=RDIA(ID)
      RSCO=RDISA/R3DH
      RTAUI(1)=RX-RDISX
      RTAUI(2)=RY-RDISY
      RTAU2=RFSP3D(RTAUI,RTAUI)
      RTAU=SQRT(RTAU2)
      RTARA=2.0*RTAU*RDISR
C
C     far field expression (2D function)
C
      IF (RTAU.GT.RDISR+NDIRH*R3DH) THEN
      RFDIPC=0.0
C      RFDIPC=RDISA*RF3DLS(RTAU2) ! blocked and replaced by CDIOMEG
      RETURN
      ENDIF
C
      RDISR2=RDISR*RDISR
      RTMR=RTAU-RDISR
      RTPR=RTAU+RDISR
      RTDN=RTMR/RTPR
  6   RZP=RDISZ
      RZM=RDISZ
C
C     disc at upper or lower aquifer boundary
C
      LDITOP=ABS(RDISZ-RZTOP).LT.REPS
      LDIBOT=ABS(RDISZ-RZBOT).LT.REPS
      IF (LDIBOT.OR.LDITOP.and.ndimag.gt.0) THEN
      RFDIPC=RFDISC(RDISZ)
      DO J=1,NDIMAG
      RZP=RZP+RHH
      RFDIPC=RFDIPC+RFDISC(RZP)
      RZM=RZM-RHH
      RFDIPC=RFDIPC+RFDISC(RZM)
      END do
      RZP=RZP+R3DH
      RZM=RZM-R3DH
      RFDIPC=2.0*RFDIPC/RPI4
      RFDIPC=RFDIPC+RSCO*RF3DSL(RZ,RZP,RZM,RTAU2)
      RFDIPC=RFDIPC*R3DH
      cz0=CMPLX(rdisx,rdisy)
      rl=rgvrefdist(cz0)
      rfdipc=rfdipc+rdisa/rpi2*LOG(r3dh/rl)
      RETURN
      ENDIF
C
C     disc inside the aquifer (not at a boundary)
C
      RZP1=2.0*RZBOT+RHH-RDISZ
      RZM1=2.0*RZBOT-RDISZ
      RFDIPC=RFDISC(RDISZ)
      if (ndimag.eq.0) goto 21
      RFDIPC=RFDIPC+RFDISC(RZP1)
      RFDIPC=RFDIPC+RFDISC(RZM1)
      DO J=1,NDIMAG
      RZP=RZP+RHH
      RFDIPC=RFDIPC+RFDISC(RZP)
      RZP1=RZP1+RHH
      RFDIPC=RFDIPC+RFDISC(RZP1)
      RZM=RZM-RHH
      RFDIPC=RFDIPC+RFDISC(RZM)
      RZM1=RZM1-RHH
      RFDIPC=RFDIPC+RFDISC(RZM1)
      END do
  21  RZP=RZBOT+(NDIMAG+1)*RHH
      RZM=RZBOT+R3DH-(NDIMAG+1)*RHH
      RFDIPC=RFDIPC/RPI4
      if (ndimag.gt.0) RFDIPC=RFDIPC+RSCO*RF3DSL(RZ,RZP,RZM,RTAU2)
      RFDIPC=RFDIPC*R3DH
      cz0=CMPLX(rdisx,rdisy)
      rl=rgvrefdist(cz0)
      rfdipc=rfdipc+rdisa/rpi2*LOG(r3dh/rl)
      RETURN
C
 1000 FORMAT (' ***ERROR in RFDIPC: non existing disc. ID=',I5,
     .' NDIS=',I5,/)
      END
c
c --------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CDIOM (CZ)
c
c --------------------------------------------------------------------------
c
C
C     Routine returns the complex potential outside the 3D zone
C     of the sink discs
C
      IMPLICIT NONE
      INTEGER(4) I
      COMPLEX(8) CZ,CDICOMC
      INCLUDE 'DICOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
C
      CDIOM=(0.0,0.0)
      IF (NDIS.EQ.0) RETURN
      DO 10 I=1,NDIS
      CDIOM=CDIOM+RDIS(I)*CDICOMC(CZ,I)
  10  CONTINUE
      RETURN
      END      
c
c --------------------------------------------------------------------------
c
      COMPLEX(8) FUNCTION CDICOMC (CZ,I)
c
c --------------------------------------------------------------------------
c
C
C     Routine returns the complex potential coefficient outside the 3D zone
C     of the sink disc number I
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) rl,rgvrefdist
      COMPLEX(8) CZ,CZ0,CMPLX
      INCLUDE 'DICOM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
C
      CDICOMC=(0.0,0.0)
      CZ0=CMPLX(RDIZ(1,I),RDIZ(2,I))
      RTAU=ABS(CZ-CZ0)
      IF (RTAU.GT.RDIR(I)+NDIRH*R3DH) THEN
      rl=rgvrefdist(cz0)
      CDICOMC=CDICOMC+RDIA(I)*LOG((CZ-CZ0)/rl)/RPI2      !rl used to be r3dh
      ENDIF
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDISC (RZ)
c
c --------------------------------------------------------------------------
c
C
C     Pot. coeff. for a single SINK disc (called from RFDIPC).
C     The factor 4*pi is not included!
C
      IMPLICIT NONE
      REAL(8) RETA,RETA2,ROOT2,RZ,RFE,RFK,RFDIOM
      INCLUDE 'DICOM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      RETA=R3DZ-RZ
      RETA2=RETA*RETA
      ROOT2=RTAU2+RETA2+RDISR2+RTARA
      ROOT=SQRT(ROOT2)
      RK2=2.0*RTARA/ROOT2
      RK=SQRT(RK2)
      IF (ABS(RETA).LT.REPS) THEN
      IF (ABS(RTMR).LT.10.0*REPS) THEN
      RFDISC=-4.0*RDISR
      RETURN
      ENDIF
      RFDISC=2.0*(RTMR*RFK(RK)-RTPR*RFE(RK))
      RETURN
      ENDIF
      RFDISC=2.0*(RTAU2+RETA2-RDISR2)/ROOT*RFK(RK)
      RFDISC=RFDISC-2.0*ROOT*RFE(RK)+RETA*RFDIOM(RETA)
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDIOM (RETA)
c
c --------------------------------------------------------------------------
c
C
C     The SOLID ANGLE due to a disc (called from RFDISC)
C
      IMPLICIT NONE
      REAL(8) RETA,RT2,RT,RT3,RSIG,RSIG2,RFK,
     &        RALPH2,RMA2,RMA,RKP2,RFE,RFPI, D_ONE
      INCLUDE 'DICOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      IF (ABS(RETA).GT.10.0*RDISR) THEN
      RT2=RETA*RETA+RTAU2
      RT=SQRT(RT2)
      RT3=RT*RT2
      RFDIOM=RETA*RPI*RDISR2/RT3
      RETURN
      ENDIF
      D_ONE = 1.0
      RSIG=SIGN(D_ONE,RETA)
      RSIG2=SIGN(D_ONE,RTMR)
      IF (ABS(RTMR).GT.REPS) GOTO 10
      RFDIOM=RPI*RSIG
      IF (ABS(RETA).LT.REPS) GOTO 30
      RFDIOM=-RSIG*(2.0*SQRT(1.0-RK2)*RFK(RK)-RSIG2*RPI)
      GOTO 20
  10  RALPH2=2.0*RTARA/(RTAU2+RDISR2+RTARA)
      RMA2=1.0-RALPH2
      IF (RMA2.LT.REPS) THEN
      RMA=SQRT(RMA2)
      RKP2=1.0-RK2
      RFDIOM=RMA*(RFK(RK)-RFE(RK)/RKP2)
      RFDIOM=RFDIOM+RPI*0.25*(2.0-RK2*(1.0+RALPH2))/(RKP2*SQRT(RKP2))
      RFDIOM=2.0*RETA/ROOT*(RSIG2*(RFDIOM+RMA*RMA2)-RFK(RK))
      GOTO 20
      ENDIF
      RFDIOM=2.0*RETA/ROOT*(RTDN*RFPI(RALPH2,RK)-RFK(RK))
  20  IF (RTAU.GE.RDISR)  GOTO 30
      RFDIOM=RFDIOM+RSIG*RPI2
   30 CONTINUE
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      SUBROUTINE DIQIC (CZ,RQCI,ID,L2D)
c
c --------------------------------------------------------------------------
c
C
C     Coefficient function for the discharge vector
C     due to disc ID.
C
      IMPLICIT NONE
      INTEGER(4) ID,J
      LOGICAL L2D,LDITOP,LDIBOT
      REAL(8) RQCI,RTAUI,RZBOT,RFBASE,RZTOP,RFTOP,RX,RY,RZ,RHH,
     &        RDISX,RDISY,RDISZ,RDISA,RSCO,RFSP3D,RQTAU,RQZ,RZP,
     &        RZM,R,RZP1,RZM1
      COMPLEX(8) CZ
      INCLUDE 'TRACOM.INC'
      INCLUDE 'DICOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'COM3D.INC'
      DIMENSION RTAUI(3),RQCI(3)
      DATA RTAUI /3*0.0/
  1   RQCI(1)=0.0
      RQCI(2)=0.0
      RQCI(3)=0.0
      IF (ID.GT.NDIS) THEN
      WRITE (ILUME,1000) ID,NDIS
      RETURN
      ENDIF
      RZBOT=RFBASE(CZ)
      RZTOP=RFTOP(CZ)
      RX=REAL(CZ)
      RY=AIMAG(CZ)
      RZ=R3DZ
      RHH=2.0*R3DH
      RDISX=RDIZ(1,ID)
      RDISY=RDIZ(2,ID)
      RDISZ=RDIZ(3,ID)
      RDISR=RDIR(ID)
      RDISA=RDIA(ID)
      RSCO=-RDISA/R3DH
      RTAUI(1)=RX-RDISX
      RTAUI(2)=RY-RDISY
      RTAU2=RFSP3D(RTAUI,RTAUI)
      RTAU=SQRT(RTAU2)
      RTARA=2.0*RTAU*RDISR
C
C     far field expression (2D function)
C
  5   IF (RTAU.GT.RDISR+NDIRH*R3DH) THEN
      RQTAU=-RDISA/RTAU2/RPI2 ! divided by extra rtau, see after label 30
      RQZ=0.0
      GOTO 30
      ENDIF

C
      RDISR2=RDISR*RDISR
      RTMR=RTAU-RDISR
      RTPR=RTAU+RDISR
      RTDN=RTMR/RTPR
  6   RZP=RDISZ
      RZM=RDISZ
      RQTAU=0.0
      RQZ=0.0
C
C     Disc at upper or lower boundary
C
      LDIBOT=ABS(RDISZ-RZBOT).LT.REPS
      LDITOP=ABS(RDISZ-RZTOP).LT.REPS
      IF (LDIBOT.OR.LDITOP.and.ndimag.gt.0) THEN
      IF (LDIBOT) THEN
      CALL DIQTZ (RZBOT,RQTAU,RQZ)
      ELSE
      R=RZTOP+0.01*REPS
      CALL DIQTZ (R,RQTAU,RQZ)
      ENDIF
      DO 10 J=1,NDIMAG
      RZP=RZP+RHH
      CALL DIQTZ (RZP,RQTAU,RQZ)
      RZM=RZM-RHH
      CALL DIQTZ (RZM,RQTAU,RQZ)
  10  CONTINUE
      RZP=RZP+R3DH
      RZM=RZM-R3DH
      RQTAU=2.0*RQTAU
      RQZ=2.0*RQZ
      R=MAX(RTAU,REPS)
      RQTAU=RQTAU/RPI4
      RQZ=RQZ/RPI4
      CALL SEMLSQ (RZ,RZP,RZM,RTAU,RQTAU,RQZ,RSCO)      
      RQTAU=RQTAU/R
      GOTO 30
      ENDIF
C
C     Disc inside the aquifer (not at a boundary).
C
      RZP1=2.0*RZBOT+RHH-RDISZ
      RZM1=2.0*RZBOT-RDISZ
      CALL DIQTZ (RDISZ,RQTAU,RQZ)
      if (ndimag.eq.0) goto 21
      CALL DIQTZ (RZP1,RQTAU,RQZ)
      CALL DIQTZ (RZM1,RQTAU,RQZ)
      DO 20 J=1,NDIMAG
      RZP=RZP+RHH
      CALL DIQTZ (RZP,RQTAU,RQZ)
      RZP1=RZP1+RHH
      CALL DIQTZ (RZP1,RQTAU,RQZ)
      RZM=RZM-RHH
      CALL DIQTZ (RZM,RQTAU,RQZ)
      RZM1=RZM1-RHH
      CALL DIQTZ (RZM1,RQTAU,RQZ)
  20  CONTINUE
  21  R=MAX(RTAU,REPS)
      RQTAU=RQTAU/RPI4
      RQZ=RQZ/RPI4
      RZP=RZBOT+(NDIMAG+1)*RHH
      RZM=RZBOT+R3DH-(NDIMAG+1)*RHH
      if (ndimag.gt.0) CALL SEMLSQ (RZ,RZP,RZM,RTAU,RQTAU,RQZ,RSCO)
      RQTAU=RQTAU/R
C
C     Two dimensional discharge function, treat as 2D sink disc for RQTAU (overwrite 3D RQTAU)
C     However, RQZ is calculated using the 3D functions.
C
      IF (L2D) THEN
        IF (RTAU.LE.RDISR) THEN
          RQTAU=-0.5                ! divided by an extra RTAU, see label 30
        ELSE
          RQTAU=-0.5*RDISR2/RTAU2  ! divided by an extra RTAU, see label 30
        ENDIF
      ENDIF
  30  RQCI(1)=RQTAU*RTAUI(1)
      RQCI(2)=RQTAU*RTAUI(2)
      RQCI(3)=RQZ
      RETURN
C
 1000 FORMAT (' ***ERROR in DIQIC: non existing disc. ID=',I5,
     .' NDIS=',I5,/)
      END
c
c --------------------------------------------------------------------------
c
      SUBROUTINE DIQTZ (RZ,RQTAU,RQZ)
c
c --------------------------------------------------------------------------
c
C
C     Routine returns radial and vertical specific discharge
C     components of a sink disc (strength and factor 4pi are
C     not included)
C
      IMPLICIT NONE
      REAL(8) RZ,RQTAU,RQZ,RETA,RETA1,RETA2,ROOT2,RFDIOM,REMK,RFE,RFK
      INCLUDE 'DICOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      RETA=R3DZ-RZ
      RETA2=RETA*RETA
      ROOT2=RTAU2+RETA2+RDISR2+RTARA
      ROOT=SQRT(ROOT2)
      RK2=2.0*RTARA/ROOT2
      RK=SQRT(RK2)
      IF ((1.0-RK).LT.0.01*REPS) THEN
      RQTAU=-1.E10
      RQZ=RQZ-RFDIOM(RETA)
      RETURN
      ENDIF
      REMK=-0.25*RPI
      IF (RK.GT.REPS) THEN
      REMK=(-RFK(RK)+RFE(RK))/RK2
      ENDIF
      RQTAU=RQTAU+8.0*RDISR*(REMK+0.5*RFK(RK))/ROOT
      RQZ=RQZ-RFDIOM(RETA)
      RETURN
      END
