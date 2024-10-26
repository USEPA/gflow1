C     Last change:  HMH   4 Apr 2007    7:35 pm
c     This file contains the following routines or functions
c
c     COMPLEX FUNCTION COMEGA       all functions that contribute to the complex potential
c     COMPLEX FUNCTION CFBIGZ       returns "big" Z for complex line integrals (maps z onto Z-plane)
c     COMPLEX FUNCTION CFSMALLZ     returns the inverse of "big" Z (maps Z back to physical plane)
c     REAL FUNCTION RFPOT           returns the real part of the complex potential: the discharge potential
c     REAL FUNCTION RFPSI           returns the imaginary part of the complex potential: the stream function
c     REAL FUNCTION RFHEAD          returns the head measured with respect to the datum mean sea level (msl)
c     REAL FUNCTION RFDSCH          returns the absolute value of the discharge vector
c     SUBROUTINE DISCH              returns the components if the discharge vector and the vertical specific discharge
c          ENTRY DISCHSPEC          returns specific discharge vector due to 2D functions (called by SDISCH)
c     SUBROUTINE SDISCH             returns the complete specific discharge vector (incl. approx. vertical discharge)
c     SUBROUTINE VELOC              returns the velocity vector (incl. approx. vertical velocity)
c     REAL FUNCTION RFK             complete elliptic integral of the first kind K(k)
c     REAL FUNCTION RFE             complete elliptic integral of the second kind E(k)
c     REAL FUNCTION RFPI            complete elliptic integral of the third kind Pi(k,alpha)
c     REAL FUNCTION ERFC            complementary error function Erfc(x)
c     SUBROUTINE LEGEND             Legendre functions: Pn(0), Pn(x) and P'n(x)
c     REAL FUNCTION RE1             Exponential integral function E1(x)
c     FUNCTION RF3DSP               scalar product function (3D vectors)
c     REAL FUNCTION SQROOT          shell for intrinsic square root function
C     SUBROUTINE VDIF2D             difference between vectors RV1 and RV2 (two dimensional).
c     SUBROUTINE VDIF3D             difference between vectors RV1 and RV2 (three dimensional).
c     SUBROUTINE VSUM2D             sum of two vectors (two dimensional)
c     SUBROUTINE VSUM3D             sum of two vectors (three dimensional)
c     REAL FUNCTION VABS2D          returns absolute value of vector (two dimensional)
c     REAL FUNCTION VABS3D          returns absolute value of vector (three dimensional)
c     REAL FUNCTION RFSP2D          scalar product of two 2D vectors
c     REAL FUNCTION RFSCALAR        returns the scalar product of 2D vectors (in complex numbers)
c     FUNCTION RFSP3D               scalar product function (three dimensional)
C
C
C ----------------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION COMEGA (CZ)
C
C ----------------------------------------------------------------------------------
C
C     Combines all functions that contribute to the complex potential.
C     If a function contributes to the complex potential only outside
C     a particular domain, COMEGA is set zero inside that domain.
c
c     On July 21, 2005 the function cflk_omega has been added. It contributes
c     to the complex potential, but also adds a real Poisson solution to the real part
c     of omega when CZ is inside a leakage element. This design deviates from the
c     design for all other functions.
C
      IMPLICIT NONE
      COMPLEX(8) CZ,CFUFOM,CFW3OM,CFWEOM,CFPDOM,CFLSU,
     +           CFLSG,CDBOM,CDIOM,cflk_omega,cflk_subomega
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      COMEGA=CDBOM(CZ)
c     write (iluer,1001) comega
 1001 format (' comega1: DB',2(d14.7))
      COMEGA=COMEGA+CFUFOM(CZ)
c     write (iluer,1002) comega
 1002 format (' comega2: UF',2(d14.7))
      COMEGA=COMEGA+CFW3OM(CZ)
c     write (iluer,1003) comega
 1003 format (' comega3: W3',2(d14.7))
      COMEGA=COMEGA+CFWEOM(CZ)
c     write (iluer,1004) comega
 1004 format (' comega4: WE',2(d14.7))
      COMEGA=COMEGA+CFPDOM(CZ)
c     write (iluer,1005) comega
 1005 format (' comega5: PD',2(d14.7))
      COMEGA=COMEGA+CFLSU(CZ)
c     write (iluer,1006) comega
 1006 format (' comega6:LSU',2(d14.7))
      COMEGA=COMEGA+CFLSG(CZ)
c     write (iluer,1007) comega
 1007 format (' comega7:LSG',2(d14.7))
      COMEGA=COMEGA+CDIOM(CZ)
c     write (iluer,1008) comega
 1008 format (' comega8: DI',2(d14.7))
      COMEGA=COMEGA+cflk_omega(cz)
c     write (iluer,1009) comega
 1009 format (' comega9: MG',2(d14.7))
      COMEGA=COMEGA+ cflk_subomega(cz)
c     write (iluer,1010) comega
 1010 format (' comega10:SG',2(d14.7))
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFBIGZ (CZ,CZZ1,CZZ2)
C
C ----------------------------------------------------------------------------------
C
C
C     Routine calculates "BIG Z" for complex line integrals
C
      IMPLICIT NONE
      COMPLEX(8) CZ,CZZ1,CZZ2
      INCLUDE 'lusys.inc'
C
      CFBIGZ=(2.0*CZ-(CZZ1+CZZ2))/(CZZ2-CZZ1)
C
      RETURN
      END
C ----------------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CFSMALLZ (CZ,CZ1,CZ2)
C
C ----------------------------------------------------------------------------------
C
C
C     Routine calculates "SMALL Z" for complex line integrals (is inverse of "BIG Z")
C
      IMPLICIT NONE
      COMPLEX(8) CZ,CZ1,CZ2
      INCLUDE 'lusys.inc'
C
      CFSMALLZ=0.5*(CZ*(CZ2-CZ1)+CZ1+CZ2)
C
      RETURN
      END

C
C ----------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFPOT(CZ)
C
C ----------------------------------------------------------------------------------
C
C     Function returns the POTENTIAL at CZ, R3DZ.
C     Combines all functions that generate real potentials and adds
C     the real part of the complex potential. In domains where a function
C     contributes to the complex potential the real function is set to zero.
C      
      IMPLICIT NONE
      REAL(8) RZ0,RFPOT0,RTW0,RTIME,RTWTIM,RFW3PT,
     &        RFDIPT,RFTWPT,RFDBPT,RFCONPT,RFPDPT
      COMPLEX(8) CZ,CZ0,COMEGA
      INCLUDE 'com3d.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      DATA CZ0,RZ0,RFPOT0,RTW0 /(1.0D21,1.0D21),3*1.0D21/
      RFPOT=RFPOT0
      IF (CZ.EQ.(1.0D21,1.0D21)) THEN
      CZ0=CZ
      RZ0=1.0D21
      RFPOT0=1.0D21
      RTW0=1.0D21
      RETURN
      ENDIF
      RTIME=RTWTIM()
      IF (CZ.EQ.CZ0.AND.R3DZ.EQ.RZ0.AND.RTIME.EQ.RTW0) THEN
      RETURN
      ENDIF
      RFPOT=REAL(COMEGA(CZ))
c     write (iluer,1001) rfpot
 1001 format ('rfpot1: OME: ',d14.7)
      RFPOT=RFPOT+RFTWPT(CZ)
c     write (iluer,1002) rfpot
 1002 format ('rfpot2: +TW: ',d14.7)
      RFPOT=RFPOT+RFDIPT(CZ)
c     write (iluer,1003) rfpot
 1003 format ('rfpot3: +DI: ',d14.7)
      RFPOT=RFPOT+RFW3PT(CZ)
c     write (iluer,1004) rfpot
 1004 format ('rfpot4: +W3: ',d14.7)
      RFPOT=RFPOT+RFPDPT(CZ)
c     write (iluer,1005) rfpot
 1005 format ('rfpot5: +PD: ',d14.7)
      RFPOT=RFPOT+RFDBPT(CZ)
c     write (iluer,1006) rfpot
 1006 format ('rfpot6: +DB: ',d14.7)
      RFPOT=RFPOT+RFCONPT(CZ)
c     write (iluer,1007) rfpot
 1007 format ('rfpot7: +CS: ',d14.7)
      CZ0=CZ
      RZ0=R3DZ
      RTW0=RTIME
      RFPOT0=RFPOT
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFPSI (CZ)
C
C ----------------------------------------------------------------------------------
C
C
C     Function returns the imaginary part of the complex potential
C      
      IMPLICIT NONE
      COMPLEX(8) CZ,COMEGA
      INCLUDE 'tracom.inc'
      EXTERNAL COMEGA
      RFPSI=AIMAG(COMEGA(CZ))
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFHEAD (CZ)
C
C ----------------------------------------------------------------------------------
C
C
C     Function returns the HEAD at CZ, R3DZ.
C      
      IMPLICIT NONE
      REAL(8) RDUM,RFPOT,RFHEDP
      COMPLEX(8) CZ
      INCLUDE 'tracom.inc'
      RDUM=RFPOT(CZ)
      RFHEAD=RFHEDP(RDUM,CZ)
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFDSCH (CZ)
C
C ----------------------------------------------------------------------------------
C
C
C     Function returns the absolute value of the discharge vector
C      
      IMPLICIT NONE
      REAL(8) RQI,SQROOT
      COMPLEX(8) CZ
      INCLUDE 'tracom.inc'
      DIMENSION RQI(3)
      RQI(1)=0.0
      RQI(2)=0.0
      RQI(3)=0.0
      CALL DISCH (CZ,RQI)
      RFDSCH=SQROOT(RQI(1)*RQI(1)+RQI(2)*RQI(2))
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      SUBROUTINE DISCH (CZ,RQI)
C
C ----------------------------------------------------------------------------------
C
C
C     Routine returns the DISCHARGE components in RQI.
C     -----
C     NOTE:
C     -----
C     RQI(3) is set equal to the vertical specific discharge,
C     whether approximate (Dupuit-Forchheimer) or exact (3D).
C
      IMPLICIT NONE
      REAL(8) RQI
      COMPLEX(8) CZ
      INCLUDE 'lusys.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'tracom.inc'
      DIMENSION RQI(3)
      RQI(1)=0.0
      RQI(2)=0.0
      RQI(3)=0.0
C
C     3D functions
C      
      CALL DISCQ2D (CZ,RQI)
      CALL W3QI2D (CZ,RQI)
      GOTO 10
C     ------------------------
      ENTRY DISCHSPEC (CZ,RQI)      ! used when spec. disch. is requested
C     ------------------------      ! what is returned are all discharge vector components
C                                     without the 3D specific discharges.
C     2D functions                    These must be added in the calling function.
C
      RQI(1)=0.0
      RQI(2)=0.0
      RQI(3)=0.0  
  10  CONTINUE
      CALL DBQI (CZ,RQI)
      CALL PDQI (CZ,RQI)
      CALL GVQI (CZ,RQI)
      CALL LSQI (CZ,RQI)
      CALL WELQI(CZ,RQI)
      CALL TWQI (CZ,RQI)
      CALL LKQI (CZ,RQI)
      call lksubqi (cz,rqi)
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C      
      SUBROUTINE SDISCH (CZ,RQI)
C
C ----------------------------------------------------------------------------------
C
C
C     Routine returns the SPECIFIC DISCHARGE in RQI.
C
      IMPLICIT NONE
      LOGICAL LREV,LEND,lfinterface
      REAL(8) RQI,RHGHT,RFHGHT,RSTEP,RDS0,RPOT,RFPOT,RQ30,
     &        RHED,RFHEDP,RTOP,RFTOP,RQS2,RFBASE,rbloc,rfperm,
     &        rfinterface,rintelv,rhss,rsgs,rsgf,rfac1,rfac2
      COMPLEX(8) CZ
      INCLUDE 'lusys.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'dwmn.inc'
      INCLUDE 'tracom.inc'
      DIMENSION RQI(3)
      RHGHT=RFHGHT(CZ)
      IF (RHGHT.LT.0.001) THEN
C           CALL TONE              ! NOT AVAILABLE IN BATCH MODEE
C            IF (LOPEN) CALL PLSETCUR (7,1,1)   ! NOT AVAILABLE IN BATCH MODE
            WRITE (ILUER,1000)
            CALL GETSTEP (RDS0,RSTEP,LEND,LREV)
            LEND=.TRUE.
            CALL SETSTEP (RSTEP,LEND)
            RETURN
      ENDIF
      CALL DISCHSPEC (CZ,RQI)
      RPOT=RFPOT(CZ)
      RHED=RFHEDP(RPOT,CZ)
      RTOP=RFTOP(CZ)
      rbloc=rfbase(cz)
      RINTELV=RFINTERFACE(CZ)
c      write (ilume,1001) rintelv,rbloc
c 1001 format (' sdisch1: rintelv,rbloc ',2(d14.7))
c      write (ilume,1003) rbloc,rintelv,r3dz,rhed
c 1003 format (' sdisch3: rbloc,rintelv,r3dz,rhed ',4(d14.7))
      IF (LFINTERFACE().AND.RINTELV.GT.RBLOC) THEN ! interface is present
         CALL interfacedata (rhss,rsgs,rsgf,rfac1,rfac2)
         RPOT=RFPERM(CZ)*RHGHT*RHGHT
         RQS2=RQI(1)*RQI(1)+RQI(2)*RQI(2)
         RQ30=RFAC2/RFAC1*RQS2/RPOT
c      write (ilume,1002) r3dz,rpot,rhed,rtop,rbloc,rintelv,rqs2,rq30
c 1002 format (' sdisch2: r3dz,rpot,rhed,rtop,rbloc,rintelv,rqs2,rq30',/,
c     &        8(d14.7))
         IF (RHED.GE.RTOP) THEN  ! ----------------- confined interface flow
           RQI(3)=RQI(3)-RQ30*(R3DZ-RINTELV)/RHGHT+RQ30
         ELSE  ! ----------------------------------- unconfined interface flow
           RQI(3)=RQI(3)-RQ30*RFAC1*(R3DZ-RINTELV)/RHGHT+RQ30
         END IF
      ELSE
        IF (RHED.LT.RTOP) THEN  ! no interface present, us original logic
C -------------------------------------add term to q3 if unconfined
        RQS2=RQI(1)*RQI(1)+RQI(2)*RQI(2)
        RQI(3)=RQI(3)-0.5*RQS2/RPOT*(R3DZ-RBLOC)/RHGHT
      ENDIF
      ENDIF
      RQI(1)=RQI(1)/RHGHT
      RQI(2)=RQI(2)/RHGHT
C
C     Add 3D functions
C
      CALL DISCQ3D (CZ,RQI)
      CALL W3QI3D (CZ,RQI)  ! only valid under confined flow conditions.
C                             The vertical flow due to recharge and leakage is calculated above,
      RETURN                ! while the vertical flow due to the 3D functions are added here.
 1000 FORMAT (' ***ERROR calculating qi: saturated aquifer thickness ',
     &'= 0.')
      END
C
C ----------------------------------------------------------------------------------
C    
      SUBROUTINE VELOC (CZ,RVI)
C
C ----------------------------------------------------------------------------------
C     
C
C     Routine returns the VELOCITY in RVI.
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RVI,RQI,RPOR,RFPOR
      COMPLEX(8) CZ
      INCLUDE 'com3d.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      DIMENSION RQI(3),RVI(3)
      RPOR=RFPOR (CZ)
      RQI(1:3)=0.0d0
      CALL SDISCH (CZ,RQI)
      DO 10 I=1,3
      RVI(I)=RQI(I)/RPOR
  10  CONTINUE
      RETURN
      END
C ---------------------------------------------------------------------------------
c 
C     COMPLETE ELLIPTIC INTEGRALS OF THE FIRST, SECOND AND THIRD KIND.
C
C     References:
C
C     "On the direct numerical calculation of elliptic functions
C     and integrals", L.V. King, Cambridge University Press,
C     London, 1924
C
C     "Handbook of elliptic integrals for engineers and scientists",
C     P.F.Byrd and M.D. Friedman, second edition, Springer-Verlag,
C     New York, 1971.
C
C     "Table of Integrals, Series and Products", I.S. Gradshteyn
C     and I.M. Ryzhik, Academic Press, New York, 1980.
C
C
C ----------------------------------------------------------------------------------
C 
      REAL(8) FUNCTION RFK(RK)
C
C ----------------------------------------------------------------------------------
C 
C
C     Complete elliptic integral of the first kind
C
C     King, 1924 formula before 27 on page 8
C     see also Byrd and Friedman page 298 and 299
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RK,RA,RB,RC,RFKOLD,RKOLD,RKP2,RKP,SQROOT,RDUM
      COMMON /ELINT/ RA(8),RB(8),RC(8)
      INCLUDE 'lusys.inc'
      SAVE
      DATA RKOLD,RFKOLD /2.0,1.E25/
      IF (RK.LT.0.0.OR.RK.GT.1.0) THEN
      RFK=0.0
      WRITE (ILUER,1000) RK
      GOTO 20
      ENDIF
      RFK=RFKOLD
      IF (RK.EQ.RKOLD) GOTO 20
      RKOLD=RK      
      IF (RK.EQ.1.0) THEN
      RFK=1.E25
      RFKOLD=RFK
      GOTO 20
      ENDIF
      IF (RK.EQ.0.0) THEN
      RFK=1.57079633
      RFKOLD=RFK
      GOTO 20
      ENDIF
      RA(1)=1
      RKP2=1.0-RK*RK
      RKP=SQROOT(RKP2)
      RB(1)=RKP
      RC(1)=RK
      DO 10 I=1,7
      RA(I+1)=0.5*(RA(I)+RB(I))
      RB(I+1)=SQROOT(RA(I)*RB(I))
      RC(I+1)=0.5*(RA(I)-RB(I))
  10  CONTINUE
      RFK=1.57079633/RA(8)
      IF (RKP.LT.0.01) THEN
      RDUM=1.38629436-LOG(RKP)
      RFK=RDUM+0.25*(RDUM-1.0)*RKP2
      RFK=RFK+0.5625*(RDUM-1.166666666666)*RKP2*RKP2
      ENDIF
  20  RFKOLD=RFK
      RETURN
 1000 FORMAT (' ***ERROR IN RFK: MODULUS OUT OF RANGE. K=',G11.4)
      END
C
C ----------------------------------------------------------------------------------
C 
      REAL(8) FUNCTION RFE(RK)
C
C ----------------------------------------------------------------------------------
C 
C
C     Complete elliptic integral of the second kind.
C
C     References:
C
C     King, 1924 formula 27 on page 8
C     Byrd and Friedman, 1971 formula 900.10 on page 299
C
C
      IMPLICIT NONE
      INTEGER(4) I,J
      REAL(8) RK,RKOLD,RFEOLD,RA,RB,RC,RK0,RFK,RDUM,RKP2,RKP,SQROOT
      COMMON /ELINT/ RA(8),RB(8),RC(8)
      INCLUDE 'lusys.inc'
      SAVE      
      DATA RKOLD,RFEOLD /2.0,1.E25/
      IF (RK.LT.0.0.OR.RK.GT.1.0) THEN
      RFE=0.0
      WRITE (ILUER,1000) RK
      GOTO 15
      ENDIF
      RFE=RFEOLD
      IF (RK.EQ.RKOLD) GOTO 15
      RKOLD=RK
      IF (RK.EQ.0.0) THEN
      RFE=1.570796327
      GOTO 15
      ENDIF
      IF (RK.EQ.1.0) THEN
      RFE=1.0
      GOTO 15
      ENDIF
      RK0=RFK(RK)
      RDUM=0.0
      J=1
      DO 10 I=1,8
      RDUM=RDUM+J*RC(I)*RC(I)
      J=J*2
  10  CONTINUE
      RFE=RK0*(1.0-(0.5*RDUM))
      IF (RK.GT.0.99) THEN
      RFE=1.0
      RKP2=1.0-RK*RK
      IF (RKP2.GT.0.0) THEN
      RKP=SQROOT(RKP2)
      ELSE
      RKP=0.0
      ENDIF
      IF (RKP.LT.1.E-10) GOTO 15
      RDUM=LOG(4.0/RKP)
      RFE=RFE+0.5*(RDUM-0.5)*RKP2
      RFE=RFE+0.1875*(RDUM-1.0833333333333)*RKP2*RKP2
      ENDIF
  15  RFEOLD=RFE
      RETURN
 1000 FORMAT (' ***ERROR IN RFE: MODULUS OUT OF RANGE. K=',G11.4)
      END
C
C ----------------------------------------------------------------------------------
C 
      REAL(8) FUNCTION RFPI (RALPH2,RK)
C
C ----------------------------------------------------------------------------------
C 
C
C     Complete elliptic integral of the third kind.
C     ( only circular case: k<alpha<1 , and the case k=alpha )
C
C     References:
C
C     Byrd and Friedman, 1971 definition 110.08 on page 10, and second form.
C                        of 412.01 on page 227.
C     Circular case k<alpha<1 (case III in King, 1924 on page 19)
C     King, 1924 formula 75 on page 15 and formula 123 on page 22
C
      IMPLICIT NONE
      INTEGER(4) I,J
      REAL(8) RALPH2,RK,RA,RB,RC,RSN,RCN,RPSI,
     &        RKOLD,RALOLD,RPIOLD,REPS,RK2,RKP2,
     &        RFE,RK0,RFK,RFAC,TETA,SQROOT,RDUM,RCOS,RZ
      COMMON /ELINT/ RA(8),RB(8),RC(8)
      INCLUDE 'lusys.inc'
      DIMENSION RSN(7),RCN(7),RPSI(8)
      SAVE
      DATA RKOLD,RALOLD,RPIOLD,REPS /0.0,0.0,1.57079633,0.000001/
      RFPI=RPIOLD
      IF (RK.EQ.RKOLD.AND.RALPH2.EQ.RALOLD) GOTO 25
      RKOLD=RK
      RALOLD=RALPH2
      RK2=RK*RK
      RFPI=0.0
      RKP2=1.0-RK2
      IF (ABS(RALPH2-RK2).LT.REPS) THEN
      RFPI=RFE(RK)/RKP2
      GOTO 25
      ENDIF
      IF (RK2.GT.RALPH2) THEN
      WRITE (ILUER,1000) RALPH2,RK
      GOTO 25
      ENDIF
      IF (RK.LT.0.0.OR.RALPH2.GT.1.0) THEN
      WRITE (ILUER,2000) RALPH2,RK
      GOTO 25
      ENDIF
      RK0=RFK(RK)
      RFAC=(RALPH2-RK2)/RKP2/RALPH2
      TETA=ASIN(SQROOT(RFAC))
      RDUM=RB(1)*RB(1)/RA(1)/RA(1)
      RCOS=SQROOT(1-RDUM*RFAC)
      RPSI(1)=TETA
      DO 10 I=1,7
      RSN(I)=RB(I)*SIN(RPSI(I))/RA(I)
      RCN(I)=SQROOT(1.0-RSN(I)*RSN(I))
      RPSI(I+1)=0.5*(RPSI(I)+ASIN(RSN(I)))
  10  CONTINUE
      RZ=0.0
      J=2
      DO 20 I=1,5
      J=J*2
      RZ=RZ+(J-2)*RC(I+2)*RSN(I+1)/RCN(I+1)/RCN(I+2)
  20  CONTINUE
      RFPI=2.0*RK0*RCOS*(RA(8)*SIN(RPSI(8))+RZ)/(RKP2*SIN(2.0*RPSI(1)))
  25  RPIOLD=RFPI
      RETURN
 1000 FORMAT (' ***ERROR IN RFPI: K**2 LARGER THAN ALPHA**2',/,
     &        ' ALPHA**2=',G11.4,'  K=',G11.4)
 2000 FORMAT (' ***ERROR IN RFPI: ARGUMENTS OUT OF RANGE.'/,
     &        ' ALPHA**2=',G11.4,'  K=',G11.4)
      END
C
C ----------------------------------------------------------------------------------
C 
      REAL(8) FUNCTION ERFC (X)
C
C ----------------------------------------------------------------------------------
C 
C
C     COMPLEMENTARY ERROR FUNCTION.
C     Rational approximation for 0<=X<infinite.
C     Reference: Abramowitz & Stegun formula 7.1.26 (page 299).
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) X,RA,RX,RX2,RT,RT0,RP
      DIMENSION RA(5)
      DATA RA /0.254829592,-0.284496736,1.421413741,
     &         -1.453152027,1.061405429/
      DATA RP /0.3275911/
C
      ERFC=2.0
      IF (X.LT.-8.8) RETURN
      ERFC=0.0
      IF (X.GT.8.8) RETURN
      RX=ABS(X)
      RX2=RX*RX
      RT=1.0/(1.0+RP*RX)
      RT0=RT
      DO 10 I=1,5
      ERFC=ERFC+RT*RA(I)
      RT=RT*RT0
  10  CONTINUE
      ERFC=ERFC*EXP(-RX2)
      IF (X.LT.0.0) ERFC=2.0-ERFC
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C 
      SUBROUTINE LEGEND (RX,IORDER)
C
C ----------------------------------------------------------------------------------
C 
C
C     Legendre functions: Pn(0), Pn(x) and P'n(x).
C     Ref. Murry Spiegel, "Mathematical Handbook",
C     Schaum's outline series, McGraw Hill Book Company,
C     1968, page 147.
C     At first call RP0(I) is filled with the Legendre functions
C     with argument zero, see formula 25.32.
C     RP0(1)=0, RP0(1)=-1/2, RP0(3)=0, RP0(4)=1*3/(2*4),...
C     The array RP(I) is filled with Legendre functions with
C     argument RX, see formula 25.20.
C     The array RP1(I) is filled with the derivatives of the
C     Legendre functions, see formula 25.23.
C
      IMPLICIT NONE
      INTEGER(4) IORDER,I
      LOGICAL LFIRST
      REAL(8) RX,RP,RP0,RP1
      COMMON /LEG/ RP(100),RP0(100),RP1(100)
      INCLUDE 'lusys.inc'
      SAVE
      DATA LFIRST /.TRUE./
      IF (IORDER.GT.100) THEN
      WRITE (ILUER,1000) IORDER
 1000 FORMAT ('***ERROR in LEGEND: IORDER=',I5,' (max. 100)')
      RETURN
      ENDIF
      IF (LFIRST) THEN
      RP0(1)=0.0
      RP0(2)=-0.5
      DO 5 I=4,100,2
      RP0(I-1)=0.0
      RP0(I)=-RP0(I-2)*(I-1)/I
  5   CONTINUE
      LFIRST=.FALSE.
      ENDIF
      RP(1)=RX
      RP(2)=0.5*(3.0*RX*RX-1.0)
      RP1(1)=1
      RP1(2)=3.0*RX
      DO 10 I=3,IORDER
      RP(I)=((2*(I-1)+1)*RX*RP(I-1)-(I-1)*RP(I-2))/I
      RP1(I)=RP1(I-2)+(2*I-1)*RP(I-1)
  10  CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RE1 (R)
C
C ----------------------------------------------------------------------------------
C
C
C     EXPONENTIAL INTEGRAL FUNCTION E1.
C     REFERENCE ABRAMOWITZ & STEGUN, 5.1.53 AND 5.1.56
C     NOTE: REST TERM SMALLER THAN THE MACHINE ACCURACY.
C
      IMPLICIT NONE
      REAL(8) R,R0,R1,R2,R3,R4,R5
      INCLUDE 'lusys.inc'
      IF (R.EQ.0.0) THEN
        RE1=1.0E+31
        WRITE(ILUER,*) 'ARG = 0.0 IN ROUTINE RE1'
        RETURN
        ENDIF
      R2=R*R
      R3=R*R2
      R4=R*R3
      IF (R.GT.1.0) GOTO 10
      R5=R*R4
      RE1=-0.57721566+0.99999193*R-0.24991055*R2+0.05519968*R3
      RE1=RE1-0.00976004*R4+0.00107857*R5-LOG(R)
      RETURN
  10  RE1=R4+8.5733287401*R3+18.059016973*R2+8.6347608925*R+0.2677737343
      R1=R4+9.5733223454*R3+25.6329561486*R2+21.0996530827*R+3.958496228
      R0=0.0
      IF (R.LT.150.0) R0=EXP(-R)
      RE1=RE1*R0/(R*R1)
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RF3DSP(RUI,RVI)
C
C ----------------------------------------------------------------------------------
C
C
C   scalar product function (3D)
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RUI,RVI
      INCLUDE 'tracom.inc'
      DIMENSION RUI(3),RVI(3)
      RF3DSP=0.0
      DO 10 I=1,3
      RF3DSP=RF3DSP+RUI(I)*RVI(I)
   10 CONTINUE
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      REAL(8) FUNCTION SQROOT(RX)
C
C ----------------------------------------------------------------------------------
C
      IMPLICIT NONE
      REAL(8) RX
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      SQROOT=0.0
      IF (RX.LT.1.0E-30.AND.RX.GE.0.0) SQROOT=0.0
      IF (RX.GE.1.0E-30) SQROOT=SQRT(RX)
      IF (RX.LT.0.0) WRITE (ILUER,1000) RX
 1000 FORMAT (' ***ERROR, negative argument in SQROOT:',G11.4)
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFSCALAR(CZ1,CZ2)
C
C ----------------------------------------------------------------------------------
C
C
C     real function returns the scalar product of two vectors,
C     which x and y components are the real and imaginary part of
C     CZ1 and CZ2, respectively.
C
      IMPLICIT NONE
      COMPLEX(8) CZ1,CZ2
      RFSCALAR=REAL(CZ1)*REAL(CZ2)+AIMAG(CZ1)*AIMAG(CZ2)
      RETURN
      END
C
C ----------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFSP3D (RUI,RVI)
C
C ----------------------------------------------------------------------------------
C
C
C   scalar product function (three dimensional).
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RUI,RVI
      SAVE
      DIMENSION RUI(3),RVI(3)
      RFSP3D=0.0
      DO 10 I=1,3
      RFSP3D=RFSP3D+RUI(I)*RVI(I)
   10 CONTINUE
      RETURN
      END

