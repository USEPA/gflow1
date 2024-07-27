C     Last change:  HMH  19 Dec 2006   10:52 am
C     This file contains the following routines and functions
C
C     REAL FUNCTION RFW3PT        Specific discharge potential contribution of all part. penetrating wells
C     COMPLEX FUNCTION CFW3OM     Calculates ppwell contribution to potential in 2D zone
C     COMPLEX FUNCTION CFW3OMCOF  Coefficient function for ppwell in 2D zone
C     REAL FUNCTION RFLSF         coefficient function F
C     REAL FUNCTION RFLSG         coefficient function G
C     REAL FUNCTION RFLSE         coefficient function E
C     SUBROUTINE W3GENC           generation of common constants of F,G,E pot. coef. fcns.
C     REAL FUNCTION RFW3CP        coefficient function  (potential) for node INODE of line sink string (well) IW
C     REAL FUNCTION RFW3IP        coefficient due to the F, G, or E functions
C     SUBROUTINE W3QCO            coefficient function for specific discharge of the well.
C     SUBROUTINE W3QIM            contributions to QTAU and QZ due to the F, G, and E functions
C     SUBROUTINE W3QI2D           specific discharge vector due to all ppwells
C          ENTRY W3QI3D
C     SUBROUTINE W3QTZF           returns the radial and vertical specific discharge component for F function
C     SUBROUTINE W3QTZG           returns the radial and vertical specific discharge component for G function
C     SUBROUTINE W3QTZE           returns the radial and vertical specific discharge component for E function
c     REAL FUNCTION RF3DLS        specific discharge potential for an infinite vertical linesink
c     REAL FUNCTION RF3DSL        specific discharge potential for a pair of semi-infinite linesinks
c     SUBROUTINE SEMLSQ           returns specific discharge vector due to semi-infinite linesinks
C
C
C
C --------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFW3PT (CZ)
C
C --------------------------------------------------------------------------------------
C
C
C     Specific discharge potential contribution of all 
C     part. penetrating wells.
C     (k*phi*H)
C
      IMPLICIT NONE
      INTEGER(4) IW,INODE,IST,IEN,I
      LOGICAL L3DH,L3D,L2D
      REAL(8) RHOI,RXI,RHO2,RF3DSP,RADTST,
     &        RFBASE,RFTOP,RADTST2,RDUM1,RFW3CP,RDUM
      COMPLEX(8) CZ
      INCLUDE 'w3com.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      COMMON /FAR3D/ L2D,L3D,L3DH
      DIMENSION RHOI(3),RXI(3)
      RFW3PT=0.0
      IF(NW3.EQ.0)RETURN
      DO 20 IW=1,NW3
      INODE=IPNT(IW)
      RXI(1)=REAL(CZ)
      RXI(2)=AIMAG(CZ)
      RHOI(1)=RXI(1)-RW3ST(1,INODE)
      RHOI(2)=RXI(2)-RW3ST(2,INODE)
      RHOI(3)=0.0
      RHO2=RF3DSP(RHOI,RHOI)
      RADTST=RZONEMLTP*(RFTOP(CZ)-RFBASE(CZ))
      RADTST2=RADTST*RADTST
      IF (RHO2.GT.RADTST2) GOTO 20
      IST=IPNT(IW)
      IEN=IPNT(IW+1)-1
      DO 10 I=IST,IEN
      RDUM1=RFW3CP(I,IW,CZ)
      RDUM=RW3S(I)*RDUM1
      RFW3PT=RFW3PT+RDUM
  10  CONTINUE
  20  CONTINUE
      IF (ABS(RFW3PT).GT.1.0E10) THEN
      WRITE (ILUER,1000) CZ,R3DZ,RFW3PT
 1000 FORMAT (' ***WARNING: RFW3PT larger than 1.0E10 :'/
     &        ' CZ,R3DZ,RFW3PT ',4G14.7)      
      ENDIF
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFW3OM(CZ)
C
C --------------------------------------------------------------------------------------
C
C
C     Calculates ppwell contribution to potential in 2D zone.
C
      IMPLICIT NONE
      INTEGER(4) IW,INODE
      REAL(8) RHOI,RXI,RTAU,RTAU2,RPI8,RPI4,RPI2,RPI1,R2H,RH,RH2,RT,
     &        RHO2,RZB,RZT,RHGHT,RFBASE,RFTOP,RF3DSP,RADTST,
     &        RADTST2,RL,RGVREFDIST
      COMPLEX(8) CZ,CZ0,CFW3OMCOF
      INCLUDE 'w3com.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      COMMON /W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &               RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      INCLUDE 'tracom.inc'
      DIMENSION RHOI(3),RXI(3)
      CFW3OM=(0.0,0.0)
C            
      IF (NW3.EQ.0) RETURN
      RHGHT=RFTOP(CZ)-RFBASE(CZ)
      DO 10 IW=1,NW3
      INODE=IPNT(IW)
      RXI(1)=REAL(CZ)
      RXI(2)=AIMAG(CZ)
      RHOI(1)=RXI(1)-RW3ST(1,INODE)
      RHOI(2)=RXI(2)-RW3ST(2,INODE)
      RHOI(3)=0.0
      RHO2=RF3DSP(RHOI,RHOI)
      RADTST=RZONEMLTP*RHGHT
      RADTST2=RADTST*RADTST
      IF (RHO2.LE.RADTST2) GOTO 10
      cz0=CMPLX(RW3ST(1,INODE),RW3ST(2,INODE))
      rl=rgvrefdist(cz0)
      CFW3OM=CFW3OM+RW3Q(IW)*CFW3OMCOF(RHO2,rl) ! rl used to be rhght
   10 CONTINUE
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C      
      COMPLEX(8) FUNCTION CFW3OMCOF(RHO2,RHGHT)
C
C --------------------------------------------------------------------------------------
C
C
C     Coefficient function for ppwell in 2D zone
C
      IMPLICIT NONE
      REAL(8) RHO2,RHGHT,RHGHT2
      INCLUDE 'w3com.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      CFW3OMCOF=CMPLX(0.0,0.0)
      RHGHT2=RHGHT*RHGHT
      CFW3OMCOF=CFW3OMCOF+LOG(RHO2/RHGHT2)/(12.56637062)
      RETURN
      END
C
C     Last change:  HH   19 Oct 1999   10:38 am
      REAL(8) FUNCTION RFLSF(RZ)
C
C     F function, formula 7 of
C     `A new analytic function to model partially penetrating
C     wells', Haitjema, Water Resour. Res. 1988
C
      IMPLICIT NONE
      LOGICAL LDRT
      REAL(8) RZ,RPI8,RA,RB,RC,RD,RK,GENOLD
      COMMON/GENC/RA,RB,RC,RD,RK,GENOLD(4),LDRT
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DATA RPI8 /25.13274124/
      LDRT=.FALSE.
      CALL W3GENC(RZ)
      RFLSF=((1.0-RA)*LOG(RB)+RC)/RPI8
      IF (ABS(RFLSF).GT.1.0E10) THEN
      WRITE (ILUOUT,1000) RA,RB,RC,RFLSF
 1000 FORMAT (' ***WARNING large value in RFLSF'/
     &        'RA,RB,RC,RFLSF ',3G14.7)      
      ENDIF
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFLSG(RZ)
C
C --------------------------------------------------------------------------------------
C
C
C     G function, formula 7 of
C     `A new analytic function to model partially penetrating
C     wells', Haitjema, Water Resour. Res. 1988
C
      IMPLICIT NONE
      LOGICAL LDRT
      REAL(8) RZ,RPI8,RA,RB,RC,RD,RK,GENOLD
      COMMON/GENC/RA,RB,RC,RD,RK,GENOLD(4),LDRT
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DATA RPI8 /25.13274124/
      LDRT=.FALSE.
      CALL W3GENC(RZ)
      RFLSG=((1.+RA)*LOG(RB)-RC)/RPI8
      IF (ABS(RFLSG).GT.1.0E10) THEN
      WRITE (ILUOUT,1000) RA,RB,RC,RFLSG
 1000 FORMAT (' ***WARNING large value in RFLSG'/
     &        'RA,RB,RC,RFLSG ',3G14.7)      
      ENDIF      
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFLSE(RZ)
C
C --------------------------------------------------------------------------------------
C
C
C     E function, formula 9 of
C     `A new analytic function to model partially penetrating
C     wells', Haitjema, Water Resour. Res. 1988
C
      IMPLICIT NONE
      LOGICAL LDRT
      REAL(8) RZ,RA,RB,RC,RD,RK,GENOLD,RFK
      COMMON/GENC/RA,RB,RC,RD,RK,GENOLD(4),LDRT
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      LDRT=.TRUE.
      CALL W3GENC(RZ)
      RFLSE=RD*RFK(RK)
      IF (ABS(RFLSE).GT.1.0E10) THEN
      WRITE (ILUOUT,1000) RD,RK,RFLSE
 1000 FORMAT (' ***WARNING large value in RFLSE'/
     &        'RD,RK,RFLSE ',3G14.7)      
      ENDIF
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE W3GENC(RZ)
C
C --------------------------------------------------------------------------------------
C
C
C generation of common constants of F,G,E pot. coef. fcns.
C
      IMPLICIT NONE
      LOGICAL LDRT
      REAL(8) RZ,RA,RB,RC,RD,RK,GENOLD,
     &        RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RU2,RV2,RU,RV,SQROOT,RUV,R,RK2
      COMMON/W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &              RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      COMMON/GENC/RA,RB,RC,RD,RK,GENOLD(4),LDRT
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      SAVE
      RU2=RHO2+(RZ-RZT)**2
      RV2=RHO2+(RZ-RZB)**2
      RU=SQROOT(RU2)
      RV=SQROOT(RV2)
      RUV=RU*RV
      RA=(RZ-(RZB+RZT)*0.5)/RH
      RB=(RU+RV-2.0*RH)/(RU+RV+2.0*RH)
      RC=(RU-RV)/RH
      IF (LDRT) THEN
        R=1.0-((RU-RV)/(2.0*RH))**2
        R=MAX(R,0.0)
        RK2=R*RH2/RUV
        RK=SQROOT(RK2)
        RD=-RH/(RPI2*SQROOT(RUV))
      ENDIF
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
       REAL(8) FUNCTION RFW3CP (INODE,IW,CZ)
C
C --------------------------------------------------------------------------------------
C
C
C     Coefficient function  (potential) for node INODE of
C     line sink string (well) IW
C     Potential definition: k*phi*H
C     Two horizontal no-flow boundaries using images.
C     Note: LOGICALS LEQPT and LEQPB 
C     control alternative top or bottom equipotentials.
C
      IMPLICIT NONE
      INTEGER(4) INODE,IW,NIMAG,IST,IEN,i,j
      LOGICAL LEQPT,LEQPB,LOPOS
      REAL(8) RFLSF,RFLSG,RFLSE,RHOI,RXI,
     &        RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &              RTOP,RFTOP,RFW3CF,RFW3CG,RFW3CE,
     &              RAD,RF3DSP,RHO,SQROOT,RL,RGVREFDIST,
     &              RZ,RIM2H,RZP,RZM,RSCO,RFW3IP,RDUM,RF3DSL,R1
      COMPLEX(8) CZ,CZ0
      EXTERNAL RFLSF,RFLSG,RFLSE
      INCLUDE 'w3com.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'tracom.inc'
      COMMON/W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &              RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      COMMON/IMAG/ NIMAG,LEQPT,LEQPB
      DIMENSION RHOI(3),RXI(3)
      SAVE
      RTOP=RFTOP(CZ)
      RFW3CF=0.0
      RFW3CG=0.0
      RFW3CE=0.0
      IST=IPNT(IW)
      IEN=IPNT(IW+1)-1
      IF (LW3DRT(IW)) IEN=IEN-1
c      write(ilume,1001)lw3drt(iw),ist,ien,((i,j,rw3st(i,j),i=1,3),j=1,4)
c 1001 format ('w3fun1: lw3drt,ist,ien ',l3,2x,2(i3,2x),/,
c     &         12('rw3st(',i2,',',i2,')=',d14.7))
      RXI(1)=REAL(CZ)
      RXI(2)=AIMAG(CZ)
      RXI(3)=R3DZ
      RAD=RW3RAD(IW)
      LOPOS=RW3ST(3,IST)-RAD.LT.R3DZ.AND.RW3ST(3,IEN)+RAD.GT.R3DZ
      RHOI(1)=RXI(1)-RW3ST(1,INODE)
      RHOI(2)=RXI(2)-RW3ST(2,INODE)
      RHOI(3)=0.0
      RHO2=RF3DSP(RHOI,RHOI)
      RHO=SQROOT(RHO2)
      cz0=CMPLX(RW3ST(1,INODE),RW3ST(2,INODE))
      rl=rgvrefdist(cz0)
C
C --------- protect against point within well perimeter
C
      IF (LOPOS.AND.RHO.LT.RAD) THEN
      RHO2=RAD*RAD
      IF (RHO.LT.1.E-20) THEN
      RHOI(1)=RAD
      RHOI(2)=0.0
      ELSE
      RHOI(1)=RHOI(1)/RHO*RAD
      RHOI(2)=RHOI(2)/RHO*RAD
      ENDIF
      RHO=RAD
      ENDIF
C
      RZ=R3DZ-RTOP
      RIM2H=2.0*R3DH
      RZP=NIMAG*RIM2H+R3DH
      RZM=-RZP
      RSCO=RW3L(INODE)/R3DH
      IF (LW3DRT(IW).AND.INODE.EQ.IEN+1) THEN
C
C     use function E only
C
      IF (LW3BOT(IW)) THEN
C -------------------------well starts at aquifer bottom.
            IF (LEQPB) THEN
            WRITE (ILUER,1000)
            RETURN
            ENDIF
            IF (NIMAG.EQ.0) THEN
            WRITE (ILUER,2000)
            RETURN
            ENDIF
            RZT=RW3ST(3,IEN)-RTOP
            RH=RZT-RW3ST(3,IST)+RTOP
            RZB=RZT-2.0*RH
            RH2=RH*RH
            RZP=RZP-R3DH
            RZM=RZM+R3DH
            RFW3CE=RFW3IP(RZ,RIM2H,NIMAG,RFLSE,-1)
            GOTO 100
      ENDIF
      IF (LW3TOP(IW)) THEN
C -------------------------well ends at aquifer top.
            RZB=RW3ST(3,IST)-RTOP
            RH=RW3ST(3,IEN)-RZB-RTOP
            RZT=RZB+2.0*RH
            RH2=RH*RH
            RFW3CE=RFW3IP(RZ,RIM2H,NIMAG,RFLSE,1)
            GOTO 100
      ENDIF
C ------------------------well does not touch either aquifer boundary.
      RZB=RW3ST(3,IST)-RTOP
      RZT=RW3ST(3,IEN)-RTOP
      RH=(RZT-RZB)/2.0
      RH2=RH*RH
      RFW3CE=RFW3IP(RZ,RIM2H,NIMAG,RFLSE,0)
      GOTO 100
      ENDIF
      IF (INODE.EQ.IST) THEN
C
C     use function F only
C
      RZB=RW3ST(3,INODE)-RTOP
      RZT=RW3ST(3,INODE+1)-RTOP
      RH=0.5*(RZT-RZB)
      RH2=RH*RH
      RFW3CF=RFW3IP (RZ,RIM2H,NIMAG,RFLSF,0)
      GOTO 100
      ENDIF
      IF (INODE.EQ.IEN) THEN
C
C     use function G only
C
      RZB=RW3ST(3,INODE-1)-RTOP
      RZT=RW3ST(3,INODE)-RTOP
      RH=0.5*(RZT-RZB)
      RH2=RH*RH
      RFW3CG=RFW3IP (RZ,RIM2H,NIMAG,RFLSG,0)
      GOTO 100
      ENDIF
C
C     use both F and G
C
      RZB=RW3ST(3,INODE-1)-RTOP
      RZT=RW3ST(3,INODE)-RTOP
      RH=0.5*(RZT-RZB)
      RH2=RH*RH
      RFW3CG=RFW3IP (RZ,RIM2H,NIMAG,RFLSG,0)
      RZB=RW3ST(3,INODE)-RTOP
      RZT=RW3ST(3,INODE+1)-RTOP
      RH=0.5*(RZT-RZB)
      RH2=RH*RH
      RFW3CF=RFW3IP (RZ,RIM2H,NIMAG,RFLSF,0)
C      
  100 RFW3CP=RFW3CF+RFW3CG+RFW3CE
      IF (NIMAG.GT.0.AND..NOT.LEQPT.AND..NOT.LEQPB) THEN
      RDUM=RSCO*RF3DSL(RZ,RZP,RZM,RHO2)
      RFW3CP=RFW3CP+RDUM
      ENDIF
      RFW3CP=R3DH*RFW3CP
      r1=rw3l(inode)*LOG((r3dh*r3dh)/(rl*rl))/12.566370
      rfw3cp=rfw3cp+r1
      RETURN
 1000 FORMAT (' ***ERROR: aquifer bottom cannot be an equipotential '
     & 'when well touches the bottom.')
 2000 FORMAT (' ***ERROR: number of images cannot be zero when well '
     & 'touches the aquifer bottom.')
      END
C
C --------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFW3IP (RZ,RIM2H,NIMG,RFUNC,ICODE)
C
C --------------------------------------------------------------------------------------
C
C
C     Routine generates the coefficient due to the F, G, or E
C     function, including images
C
C
C     ICODE is used in conjunction with "E"-function call (double root elemnt)
C     ICODE = 0  well does not touch aquifer top or aquifer bottom
C     ICODE = 1  well touches aquifer top
C     ICODE = -1 well touches aquifer bottom
C
      IMPLICIT NONE
      INTEGER(4) NIMG,NIMAG,ICODE,IM
      LOGICAL LEQPT,LEQPB
      REAL(8) RZ,RIM2H,RFUNC,RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RSIGN,RZ1,RZ2,RZ1P,RZ1M,RZ2P,RZ2M,
     &        RSUM,RFAC,RZP,RZM
      COMMON/W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &              RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      COMMON/IMAG/ NIMAG,LEQPT,LEQPB
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      EXTERNAL RFUNC
C
C     Default is aquifer top and aquifer bottom are no-flow boundaries.
C     If LEQPT aquifer top is an equipotential boundary and aquifer bottom
C     is a no-flow boundary.
C     If LEQPB aquifer bottom is an equipotential boundary and aquifer top
C     is a no-flow boundary.
      IF (LEQPT) RSIGN=1.0
      IF (LEQPB) RSIGN=-1.0
C      
      IF (ICODE.EQ.0) THEN    ! well does not touch aquifer top or bottom
      RZ1=RZ
      RZ2=-RZ      
      RFW3IP=RFUNC(RZ1)
      IF (LEQPT) THEN
      RFW3IP=RFW3IP-RFUNC(RZ2)
      ELSE
      RFW3IP=RFW3IP+RFUNC(RZ2)
      ENDIF
      IF (NIMG.EQ.0) RETURN
      RZ1P=RZ1
      RZ1M=RZ1
      RZ2P=RZ2
      RZ2M=RZ2
      DO 10 IM=1,NIMG      
      RZ1P=RZ1P+RIM2H
      RZ1M=RZ1M-RIM2H
      RZ2P=RZ2P+RIM2H
      RZ2M=RZ2M-RIM2H
      RSUM=RFUNC(RZ1P)+RFUNC(RZ1M)+RFUNC(RZ2P)+RFUNC(RZ2M)
      IF (LEQPT.OR.LEQPB) THEN
      RSUM=RSIGN*RSUM
      RSIGN=-RSIGN
      ENDIF
      RFW3IP=RFW3IP+RSUM
  10  CONTINUE
      ENDIF
      IF (ICODE.EQ.-1) THEN         ! well at aquifer bottom
      IF (NIMG.EQ.0) THEN
      WRITE (ILUER,1000)
      RETURN
      ENDIF
      RFAC=0.0
      DO 20 IM=1,NIMG
      RZP=RZ+RFAC*RIM2H
      RFAC=RFAC+1.0
      RZM=RZ-RFAC*RIM2H            
      RSUM=RFUNC(RZP)+RFUNC(RZM)
      IF (LEQPT.OR.LEQPB) THEN
      RSUM=RSIGN*RSUM
      RSIGN=-RSIGN
      ENDIF
      RFW3IP=RFW3IP+RSUM
  20  CONTINUE
      ENDIF
      IF (ICODE.EQ.1) THEN          ! well at aquifer top
      RFW3IP=RFUNC(RZ)
      IF (NIMG.EQ.0) RETURN
      RZP=RZ
      RZM=RZ
      DO 30 IM=1,NIMG
      RZP=RZP+RIM2H
      RZM=RZM-RIM2H      
      RSUM=RFUNC(RZP)+RFUNC(RZM)
      IF (LEQPT.OR.LEQPB) THEN
      RSUM=RSIGN*RSUM
      RSIGN=-RSIGN
      ENDIF
      RFW3IP=RFW3IP+RSUM
  30  CONTINUE
      ENDIF      
      RETURN
 1000 FORMAT (' ***ERROR: number of images should be larger than 0 ',
     & 'when well touches aquifer bottom.')
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE W3QCO (INODE,IW,CZ,RQW3I)
C
C --------------------------------------------------------------------------------------
C
C
C     Coefficient function for specific discharge of the well.
C
C     Two horizontal no-flow boundaries using images.
C
C     Note: logicals LEQPT and LEQPB
C     control alternative top or bottom equipotentials.
C
      IMPLICIT NONE
      INTEGER(4) INODE,IW,NIMAG,IST,IEN
      LOGICAL LEQPT,LEQPB
      REAL(8) RQW3I,W3QTZF,W3QTZG,W3QTZE,
     &        RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RTAUI,RXI,RTOP,RFTOP,RQTAUF,RQTAUG,RQTAUE,
     &        RQT,RQZF,RQZG,RQZE,RQZ,RF3DSP,SQROOT,RZ,
     &        RIM2H,RZP,RZM,RSCO
      COMPLEX(8) CZ
      EXTERNAL W3QTZF,W3QTZG,W3QTZE
      COMMON/W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &              RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      INCLUDE 'w3com.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      COMMON/IMAG/ NIMAG,LEQPT,LEQPB
      DIMENSION RXI(3),RTAUI(3),RQW3I(3)
      SAVE
      RTOP=RFTOP(CZ)
      RQTAUF=0.0
      RQTAUG=0.0
      RQTAUE=0.0
      RQT=0.0
      RQZF=0.0
      RQZG=0.0
      RQZE=0.0
      RQZ=0.0
      IST=IPNT(IW)
      IEN=IPNT(IW+1)-1
      IF (LW3DRT(IW)) IEN=IEN-1
      RXI(1)=REAL(CZ)
      RXI(2)=AIMAG(CZ)
      RXI(3)=R3DZ
      RTAUI(1)=RXI(1)-RW3ST(1,INODE)
      RTAUI(2)=RXI(2)-RW3ST(2,INODE)
      RTAUI(3)=0.0
      RTAU2=RF3DSP(RTAUI,RTAUI)
      RTAU=SQROOT(RTAU2)      
      R2H=2*R3DH
      RZ=R3DZ-RTOP
      RIM2H=2.0*R3DH      
      RZP=NIMAG*RIM2H+R3DH
      RZM=-RZP
      RSCO=RW3L(INODE)/R3DH
      IF (LW3DRT(IW).AND.INODE.EQ.IEN+1) THEN
C
C     use function E only
C
      IF (LW3BOT(IW)) THEN
C -------------------------well starts at aquifer bottom.
            IF (LEQPB) THEN
            WRITE (ILUER,1000)
            RETURN
            ENDIF
            IF (NIMAG.EQ.0) THEN
            WRITE (ILUER,2000)
            RETURN
            ENDIF
            RZT=RW3ST(3,IEN)-RTOP
            RH=RZT-RW3ST(3,IST)+RTOP
            RZB=RZT-2.0*RH
            RH2=RH*RH
            RZP=RZP-R3DH
            RZM=RZM+R3DH
            CALL W3QIM (RZ,RQTAUE,RQZE,RIM2H,NIMAG,W3QTZE,-1)
            GOTO 100
      ENDIF
      IF (LW3TOP(IW)) THEN
C -------------------------well ends at aquifer top.
            RZB=RW3ST(3,IST)-RTOP
            RH=RW3ST(3,IEN)-RZB-RTOP
            RZT=RZB+2.0*RH
            RH2=RH*RH
            CALL W3QIM (RZ,RQTAUE,RQZE,RIM2H,NIMAG,W3QTZE,1)
            GOTO 100
      ENDIF
C ------------------------well does not touch either aquifer boundary.
      RZB=RW3ST(3,IST)-RTOP
      RZT=RW3ST(3,IEN)-RTOP
      RH=(RZT-RZB)/2.0
      RH2=RH*RH
      CALL W3QIM (RZ,RQTAUE,RQZE,RIM2H,NIMAG,W3QTZE,0)
      GOTO 100
      ENDIF
      IF (INODE.EQ.IST) THEN
C
C     use only F function
C
      RZB=RW3ST(3,INODE)-RTOP
      RZT=RW3ST(3,INODE+1)-RTOP
      RH=0.5*(RZT-RZB)
      RH2=RH*RH
      CALL W3QIM (RZ,RQTAUF,RQZF,RIM2H,NIMAG,W3QTZF,0)
      GOTO 100
      ENDIF
      IF (INODE.EQ.IEN) THEN
C
C     use function G only
C
      RZB=RW3ST(3,INODE-1)-RTOP
      RZT=RW3ST(3,INODE)-RTOP
      RH=0.5*(RZT-RZB)
      RH2=RH*RH
      CALL W3QIM (RZ,RQTAUG,RQZG,RIM2H,NIMAG,W3QTZG,0)
      GOTO 100
      ENDIF
C
C     use both F and G
C
      RZB=RW3ST(3,INODE-1)-RTOP
      RZT=RW3ST(3,INODE)-RTOP
      RH=0.5*(RZT-RZB)
      RH2=RH*RH
      CALL W3QIM (RZ,RQTAUG,RQZG,RIM2H,NIMAG,W3QTZG,0)
      RZB=RW3ST(3,INODE)-RTOP
      RZT=RW3ST(3,INODE+1)-RTOP
      RH=0.5*(RZT-RZB)
      RH2=RH*RH
      CALL W3QIM (RZ,RQTAUF,RQZF,RIM2H,NIMAG,W3QTZF,0)
C      
  100 IF (NIMAG.GT.0.AND..NOT.LEQPT.AND..NOT.LEQPB) 
     &CALL SEMLSQ (RZ,RZP,RZM,RTAU,RQT,RQZ,RSCO)
      IF (RTAU.GT.1.0E-21) THEN
      RQW3I(1)=-(RQTAUF+RQTAUG+RQTAUE+RQT)*RTAUI(1)/RTAU
      RQW3I(2)=-(RQTAUF+RQTAUG+RQTAUE+RQT)*RTAUI(2)/RTAU
      ENDIF
      RQW3I(3)=-(RQZF+RQZG+RQZE+RQZ)
      RETURN
 1000 FORMAT (' ***ERROR: aquifer bottom cannot be an equipotential ',
     &        'when well touches the bottom.')
 2000 FORMAT (' ***ERROR: number of images cannot be zero when well ',
     &        'touches the aquifer bottom.')
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE W3QIM  (RZ,RQT,RQZ,RIM2H,NIMG,W3SUB,ICODE)
C
C --------------------------------------------------------------------------------------
C
C
C     Routine generates the contributions to QTAU and QZ due to the
C     F, G, and E functions.
C
C     ICODE is used in conjunction with "E"-function call (double root elemnt)
C     ICODE = 0  well does not touch aquifer top or aquifer bottom
C     ICODE = 1  well touches aquifer top
C     ICODE = -1 well touches aquifer bottom
C
      IMPLICIT NONE
      INTEGER(4) NIMG,ICODE,NIMAG,IM
      LOGICAL LEQPT,LEQPB
      REAL(8) RZ,RQT,RQZ,RIM2H,
     &        RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RSIGN,RZ1,RZ2,RQTT,RQZT,RZ1P,RZ1M,RZ2P,RZ2M,
     &        RFAC,RZP,RZM
      COMMON/W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &              RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU 
      COMMON/IMAG/ NIMAG,LEQPT,LEQPB
      INCLUDE 'tracom.inc'      
      INCLUDE 'lusys.inc'
C
C     Default is aquifer top and aquifer bottom are no-flow boundaries.
C     If LEQPT aquifer top is an equipotential boundary and aquifer bottom
C     is a no-flow boundary.
C     If LEQPB aquifer bottom is an equipotential boundary and aquifer top
C     is a no-flow boundary.
      IF (LEQPT) RSIGN=1.0
      IF (LEQPB) RSIGN=-1.0
C     
      IF (ICODE.EQ.0) THEN    ! well does not touch aquifer top or bottom
      RZ1=RZ
      RZ2=-RZ
      CALL W3SUB (RZ1,RQT,RQZ,1.0)
      RQTT=0.0
      RQZT=0.0
      CALL W3SUB (RZ2,RQTT,RQZT,-1.0)
      IF (LEQPT) THEN
      RQTT=-RQTT
      RQZT=-RQZT
      ENDIF
      RQT=RQT+RQTT
      RQZ=RQZ+RQZT
      IF (NIMG.EQ.0) RETURN
      RZ1P=RZ1
      RZ1M=RZ1
      RZ2P=RZ2
      RZ2M=RZ2
      DO 10 IM=1,NIMG
      RZ1P=RZ1P+RIM2H
      RZ1M=RZ1M-RIM2H
      RZ2P=RZ2P+RIM2H
      RZ2M=RZ2M-RIM2H
      RQTT=0.0
      RQZT=0.0
      CALL W3SUB (RZ1P,RQTT,RQZT,1.0)
      CALL W3SUB (RZ1M,RQTT,RQZT,1.0)
      CALL W3SUB (RZ2P,RQTT,RQZT,-1.0)      
      CALL W3SUB (RZ2M,RQTT,RQZT,-1.0)
      IF (LEQPT.OR.LEQPB) THEN
      RQTT=RSIGN*RQTT
      RQZT=RSIGN*RQZT
      RSIGN=-RSIGN
      ENDIF
      RQT=RQT+RQTT
      RQZ=RQZ+RQZT
  10  CONTINUE
      ENDIF
      IF (ICODE.EQ.-1) THEN         ! well at aquifer bottom
      IF (NIMG.EQ.0) THEN
      WRITE (ILUER,1000)
      RETURN
      ENDIF
      RFAC=0.0
      DO 20 IM=1,NIMG
      RZP=RZ+RFAC*RIM2H
      RFAC=RFAC+1.0
      RZM=RZ-RFAC*RIM2H
      RQTT=0.0
      RQZT=0.0
      CALL W3SUB (RZP,RQTT,RQZT,1.0)
      CALL W3SUB (RZM,RQTT,RQZT,1.0)
      IF (LEQPT.OR.LEQPB) THEN
      RQTT=RSIGN*RQTT
      RQZT=RSIGN*RQZT
      RSIGN=-RSIGN
      ENDIF
      RQT=RQT+RQTT
      RQZ=RQZ+RQZT
  20  CONTINUE      
      ENDIF
      IF (ICODE.EQ.1) THEN          ! well at aquifer top
      RZ1=RZ
      CALL W3SUB (RZ1,RQT,RQZ,1.0)
      IF (NIMG.EQ.0) RETURN
      RZ1P=RZ1
      RZ1M=RZ1
      DO 30 IM=1,NIMG
      RZ1P=RZ1P+RIM2H
      RZ1M=RZ1M-RIM2H
      RQTT=0.0
      RQZT=0.0
      CALL W3SUB (RZ1P,RQTT,RQZT,1.0)
      CALL W3SUB (RZ1M,RQTT,RQZT,1.0)
      IF (LEQPT.OR.LEQPB) THEN
      RQTT=RSIGN*RQTT
      RQZT=RSIGN*RQZT
      RSIGN=-RSIGN
      ENDIF
      RQT=RQT+RQTT
      RQZ=RQZ+RQZT
  30  CONTINUE   
      ENDIF
      RETURN
 1000 FORMAT (' ***ERROR: number of images should be larger than 0 when
     & well touches aquifer bottom.')
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE W3QI2D (CZ,RQI)
C
C --------------------------------------------------------------------------------------
C
C
C     Routine adds the specific discharge vector (in 3-D zone) 
C     of all ppwells to the vector <RQI>.
C
      IMPLICIT NONE
      INTEGER(4) IW,IST,IEN,I,J
      LOGICAL LD2,LOPOS
      REAL(8) RQI,RQW3I,RXI,RHOI,RQI0,
     &        RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RAQTHICK,RFBASE,RFTOP,RF3DSP,RADTST,RADTST2,
     &        RQ1,RQ2,RHO,SQROOT,RAD,RDUM
      COMPLEX(8) CZ,CW3Z,CW3RAD,CQW3
      INCLUDE 'w3com.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      COMMON /W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &               RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      DIMENSION RQI(3),RQW3I(3),RXI(3),RHOI(3),RQI0(3)
      DATA RQI0 /0.0,0.0,0.0/
      LD2=.TRUE.              ! returns DISCHARGE
      GOTO 1
C     ---------------------
      ENTRY W3QI3D (CZ,RQI)   
C     ---------------------
      LD2=.FALSE.             ! returns SPECIFIC discharge
  1   RQW3I(1)=0.0
      RQW3I(2)=0.0
      RQW3I(3)=0.0
      IF (NW3.EQ.0) RETURN
      RAQTHICK=RFTOP(CZ)-RFBASE(CZ)
      DO 30 IW=1,NW3
      IST=IPNT(IW)
      IEN=IPNT(IW+1)-1
      RXI(1)=REAL(CZ)
      RXI(2)=AIMAG(CZ)
      RXI(3)=R3DZ
      RHOI(1)=RXI(1)-RW3ST(1,IST)
      RHOI(2)=RXI(2)-RW3ST(2,IST)
      RHOI(3)=0.0
      RHO2=RF3DSP(RHOI,RHOI)
      RADTST=RZONEMLTP*RAQTHICK
      RADTST2=RADTST*RADTST
C
C --------- protect against point within well perimeter
C
      RHO=SQROOT(RHO2)
      RAD=RW3RAD(IW)
      LOPOS=RW3ST(3,IST)-RAD.LT.R3DZ.AND.RW3ST(3,IEN-1)+RAD.GT.R3DZ
      IF (LOPOS.AND.RHO.LT.RAD) GOTO 30 ! do not calculate Qi
c
C
C     2D calculations
C
      IF (RHO2.GT.RADTST2.OR.LD2) THEN
        CW3Z=CMPLX(RW3ST(1,IST),RW3ST(2,IST))
        CW3RAD=CZ-CW3Z
        CQW3=RW3Q(IW)/RPI2/CW3RAD
        RQ1=-REAL(CQW3)
        RQ2=AIMAG(CQW3)
        IF (.NOT.LD2) THEN
          RQ1=RQ1/RAQTHICK
          RQ2=RQ2/RAQTHICK
        ENDIF
        RQI(1)=RQI(1)+RQ1
        RQI(2)=RQI(2)+RQ2
        GOTO 30
      ENDIF
c
      DO 20 I=IST,IEN  ! calculate 3D contribution
      CALL W3QCO (I,IW,CZ,RQW3I)
      DO 10 J=1,3
      RDUM=RW3S(I)*RQW3I(J)
      RQI(J)=RQI(J)+RDUM
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE W3QTZF (RZ,RQTAUF,RQZF,RSIGN)
C
C --------------------------------------------------------------------------------------
C
C
C     Returns the radial and vertical specific discharge component 
C     COEFFICIENTS: RQTAUF and RQZF
C     RSIGN is negative when RZ is imaged above the axis z=0.
C     The variables RZT, RZB, RTAU, RTAU2,and RH must be calculated prior to
C     calling this routine.
C     See "A New Analytic Function for Modeling Partially Penetrating
C     Wells", H.M. Haitjema and S.R. Kraemer, Water Resourc. Res.,
C     Vol.24, No.5, pp.683-690, May 1988.
C     The pertinent formulas are 16 and 17.
C
      IMPLICIT NONE
      REAL(8) RZ,RQTAUF,RQZF,RSIGN,
     &        RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RA,RB,RU2,RV2,RU,RV,SQROOT,RHH,RUV,
     &        RA1,RA2,RF,RD,RE,RM,RUVH
      INCLUDE 'lusys.inc'
      COMMON/W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &              RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      INCLUDE 'tracom.inc'
      RA=RZ-RZT
      RB=RZ-RZB
      RU2=RTAU2+RA**2
      RV2=RTAU2+RB**2
      RU=SQROOT(RU2)
      RV=SQROOT(RV2)
      RHH=2.0*RH
      RUV=RU*RV
      RUVH=RUV*RH
      RA1=RTAU*(RV+RU)/RUVH
      RA2=RTAU*(RV-RU)/RUVH
      RF=1.0/((((RU+RV)/RHH)**2)-1.0)
      RD=(RA*RV+RB*RU)/RUVH
      RE=RA*RV/RUVH
      RE=RE-RB*RU/RUVH
      RM=LOG((RU+RV-RHH)/(RU+RV+RHH))/RH
      RQTAUF=RQTAUF-(RA/RH*RA1*RF-RA2)/RPI8
      RQZF=RQZF+(-RA/RH*RD*RF+RE-RM)*RSIGN/RPI8
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE W3QTZG (RZ,RQTAUG,RQZG,RSIGN)
C
C --------------------------------------------------------------------------------------
C
C
C     Returns the radial and vertical specific discharge component 
C     COEFFICIENTS: RQTAUG and RQZG
C     RSIGN is negative when RZ is imaged above the axis z=0.
C     The variables RZT, RZB, RTAU, RTAU2,and RH must be calculated prior to
C     calling this routine.
C     See "A New Analytic Function for Modeling Partially Penetrating
C     Wells", H.M. Haitjema and S.R. Kraemer, Water Resourc. Res.,
C     Vol.24, No.5, pp.683-690, May 1988.
C     The pertinent formulas are 16 and 17.
C
      IMPLICIT NONE
      REAL(8) RZ,RQTAUG,RQZG,RSIGN,
     &        RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RA,RB,RU2,RV2,RU,RV,RUV,RUVH,SQROOT,
     &        RHH,RA1,RA2,RF,RD,RE,RM
      INCLUDE 'lusys.inc'
      COMMON/W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &              RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      INCLUDE 'tracom.inc'
      RA=RZ-RZT
      RB=RZ-RZB
      RU2=RTAU2+RA**2
      RV2=RTAU2+RB**2
      RU=SQROOT(RU2)
      RV=SQROOT(RV2)
      RHH=2.0*RH
      RUV=RU*RV
      RUVH=RUV*RH
      RA1=RTAU*(RV+RU)/RUVH
      RA2=RTAU*(RV-RU)/RUVH
      RF=1.0/((((RU+RV)/RHH)**2)-1.0)
      RD=(RA*RV+RB*RU)/RUVH
      RE=RA*RV/RUVH
      RE=RE-RB*RU/RUVH
      RM=LOG((RU+RV-RHH)/(RU+RV+RHH))/RH
      RQTAUG=RQTAUG+(RB/RH*RA1*RF-RA2)/RPI8
      RQZG=RQZG+(RB/RH*RD*RF-RE+RM)*RSIGN/RPI8
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE W3QTZE (RZ,RQTAUE,RQZE,RSIGN)
C
C --------------------------------------------------------------------------------------
C
C
C     Returns radial and vertical specific discharge component 
C     COEFFICIENTS: RQTAUE and RQZE.
C     RSIGN is negative when RZ is imaged above the axis z=0.
C     The variables RZT, RZB, RTAU, RTAU2, RH, and RH2 must be 
C     calculated prior to calling this routine.
C     See "A New Analytic Function for Modeling Partially Penetrating
C     Wells", H.M. Haitjema and S.R. Kraemer, Water Resourc. Res.,
C     Vol.24, No.5, pp.683-690, May 1988.
C     The pertinent formulas are 16, 17 and 18.
C     Note: sign of Gamma derivatives are wrong; add minus sign.
C
      IMPLICIT NONE
      REAL(8) RZ,RQTAUE,RQZE,RSIGN,RPI1D4,
     &        RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RDELTA,RA,RB,RU2,RV2,RU,RV,RUV,RUVH,
     &        SQROOT,RHH,R,RK2,RK,RKP2,RKP12,RFKRK,
     &        RFK,RC1,RH1,RJ,RFERK,RFE,RL,RG
      INCLUDE 'lusys.inc'
      COMMON/W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &              RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      DATA RPI1D4 /0.785398163/
      INCLUDE 'tracom.inc'
      RDELTA=0.03
      RA=RZ-RZT
      RB=RZ-RZB
      RU2=RTAU2+RA**2
      RV2=RTAU2+RB**2
      RU=SQROOT(RU2)
      RV=SQROOT(RV2)
      RHH=2.0*RH
      RUV=RU*RV
      R=1.0-((RU-RV)/RHH)**2
      R=MAX(R,0.0)
      RK2=R*RH2/RUV
      RK=SQROOT(RK2)
      RKP2=1.0-RK2
      RKP12=-1.0/RKP2
      RFKRK=RFK(RK)
      RC1=RV2*RA/RUV
      RC1=RC1+RU2*RB/RUV
      RH1=RH/(RUV*SQROOT(RUV))
      RJ=(RU2+RV2)/RUV
      RFERK=RFE(RK)
      RL=RFERK*RKP12
      RG=RPI1D4
      IF(RK.GE.RDELTA) THEN
      RG=RFERK/RK2
      RG=RG-RKP2*RFKRK/RK2
      ENDIF
      RQTAUE=RQTAUE-RH1*RTAU*((1.0-(0.5*RJ))*RG*RKP12+RL*RJ)/RPI4
      RQZE=RQZE-RH1*((RA+RB-RC1)*RG*RKP12+2.0*RL*RC1)*RSIGN/RPI8
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RF3DLS (RTAU2)
C
C --------------------------------------------------------------------------------------
C
C
C     Pot. coeff. function for an infinite vertical line sink.
C
      IMPLICIT NONE
      REAL(8) RTAU2,RTOLD,RLSOLD,RPI4,R
      INCLUDE 'com3d.inc'
      SAVE
      DATA RTOLD,RLSOLD /0.0,0.0/
      DATA RPI4 /12.56637062/
      RF3DLS=RLSOLD
      IF (RTAU2.EQ.RTOLD) RETURN
      R=RTAU2/(R3DH*R3DH)
      RF3DLS=LOG(R)/RPI4
      RTOLD=RTAU2
      RLSOLD=RF3DLS
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RF3DSL (RZ,RZP,RZM,RTAU2)
C
C --------------------------------------------------------------------------------------
C
C
C     Pot. coeff. function for a pair of semi-infinite line sinks.
C     See `Comparing a three-dimensional and a Dupuit-Forchheimer
C     solution for a circular recharge area in a confined aquifer',
C     Haitjema, Journal of Hydrology, 91 (1987) 83-101.
C     See formula B6 without sigma. (Note there is a sign error
C     from B5 to B6, consequently the reference to source in the
C     paper should be read as sink.)
C
      IMPLICIT NONE
      REAL(8) RZ,RZP,RZM,RTAU2,RPI4,R1,R2,RU2,RV2,RU,RV,
     &        SQROOT
      INCLUDE 'com3d.inc'
      DATA RPI4 /12.56637062/
      R1=RZP-RZ
      R2=RZ-RZM
      RU2=R1*R1+RTAU2
      RV2=R2*R2+RTAU2
      RU=SQROOT(RU2)
      RV=SQROOT(RV2)
      RF3DSL=LOG((R1+RU)*(R2+RV)/(R3DH*R3DH))/RPI4
      RETURN
      END
C
C --------------------------------------------------------------------------------------
C
      SUBROUTINE SEMLSQ (RZ,RZP,RZM,RTAU,RQTAU,RQZ,RSCO)
C
C --------------------------------------------------------------------------------------
C
C
C     Routine returns radial and vertical specific discharge component
C     coefficient functions due to semi-infinite line sinks.
C     See `Comparing a three-dimensional and a Dupuit-Forchheimer
C     solution for a circular recharge area in a confined aquifer',
C     Haitjema, Journal of Hydrology, 91 (1987) 83-101.
C     Added to RQTAU and RQZ are the derivatives of formula B6 with respect
C     to TAU and Z, respectively, multiplied by RSCO which contains
C     the factor h/H.
C     See formulas C17 and C18 without sigma. (Note reference to source
C     in the paper should be read as sink, see also routine RF3DlS.)
C
      IMPLICIT NONE
      REAL(8) RZ,RZP,RZM,RTAU,RQTAU,RQZ,RSCO,RPI4,
     &        R1,R2,RTAU2,RU2,RV2,RU,RV,SQROOT,RDUMT,RDUMZ
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DATA RPI4 /12.56637062/
      R1=RZP-RZ
      R2=RZ-RZM
      RTAU2=RTAU*RTAU
      RU2=R1*R1+RTAU2
      RV2=R2*R2+RTAU2
      RU=SQROOT(RU2)
      RV=SQROOT(RV2)
      RDUMT=(RTAU/(RU*(R1+RU))+RTAU/(RV*(R2+RV)))*RSCO/RPI4
      RQTAU=RQTAU+RDUMT
      RDUMZ=(1.0/RV-1.0/RU)*RSCO/RPI4
      RQZ=RQZ+RDUMZ
      RETURN
      END




