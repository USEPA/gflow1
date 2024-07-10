C     Last change:  HMH   4 Apr 2007    4:40 pm
C    This file contains the following routines and functions:
C
c    subroutine w3matsize     updates the matrix array dimensions
C    SUBROUTINE W3CZC         generates control point vector for ppwells
C    SUBROUTINE W3MAT         generates matrix coefficients for ppwells
C    SUBROUTINE W3KNO         generates known vector for ppwells
C    SUBROUTINE W3SUB         substitutes solution vector into strengths parameters for ppwells
c    subroutine w3update      updates the potentials at collocation points
C    REAL FUNCTION RFW3MATC   coefficient for ppwell strength parameters   
C
c -------------------------------------------------------------------
c
      subroutine w3matsize (M,N)
c
c -------------------------------------------------------------------
c
c     Routine updates the matrix array dimensions M and N
c
      implicit none
      INTEGER i,iw,M,N,ist,ien,ilast,N0
      INCLUDE 'w3com.inc'
c
      nw3strengths=0
      if (nw3.gt.0) then
        N0=N
        do iw=1,nw3
          IST=IPNT(IW)
          ILAST=IPNT(IW+1)-1
          IEN=ILAST
          IF (LW3Q(IW)) IEN=ILAST-1
          DO i=IST,IEN
             M=M+1
             N=N+1
          end do
        end do
        nw3strengths=N-N0
      end if
      return
c
      end subroutine
c
c -------------------------------------------------------------------------------------
c
      subroutine w3keep(drb,j,m)
c
c --------------------------------------------------------------------------------------
c
c
c     add number of 3D well strength parameters to array counter and set corresponding
c     values of drb array equal to 0. The drb array will contain any strength differences between the
c     solution and the actual stored strength parameters, see lskeepsigma and lkkeepstrength routines.
c
      INTEGER i,j,m
      REAL(8) drb
      DIMENSION drb(m)
      include 'w3com.inc'
      if (nw3strengths.gt.0) then  ! we have 3D wells in the matrix
       do i=1,nw3strengths
         j=j+1
         drb(j)=0.0d0
       end do
      endif
      return
      end subroutine

C
C --------------------------------------------------------------------
C
      SUBROUTINE W3CZC (CZI,M,N,DRFAC,CALPH,ITYPE)
C
C --------------------------------------------------------------------
C
C
C     Control points
C
      IMPLICIT NONE
      INTEGER(4) M,N,ITYPE,IW,IST,ILAST,IEN,I
      REAL(8) DRFAC
      COMPLEX(8) CZI,CALPH
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION CZI(*),DRFAC(4,*),CALPH(*),ITYPE(*)
      IF (NW3.EQ.0) RETURN
      DO 20 IW=1,NW3
      IST=IPNT(IW)
      ILAST=IPNT(IW+1)-1
      IEN=ILAST
      IF (LW3Q(IW)) IEN=ILAST-1
      DO 10 I=IST,IEN
      M=M+1 ! add one equation for each well collocation point
      N=N+1 ! add one strength parameter for each collocation point
      CZI(M)=CMPLX(RW3CP(1,I),RW3CP(2,I))
      R3DCPZ(M)=RW3CP(3,I)
      ITYPE(M)=1
      DRFAC(1,M)=1.0D00
      DRFAC(4,M)=1.0D0
      CALPH(M)=CMPLX(0.0D0,0.0D0)
      IF (LW3Q(IW)) THEN
            ITYPE(M)=-1
            CALPH(M)=CMPLX(RW3CP(1,I+1),RW3CP(2,I+1))
            R3DZA(M)=RW3CP(3,I+1)
      ENDIF
  10  CONTINUE
  20  CONTINUE
      RETURN
      END
C
C --------------------------------------------------------------------
C
      SUBROUTINE W3GENMAT (DRA,CZI,M,N,J,CALPH,ITYPE)   ! TO BE DONE
C
C --------------------------------------------------------------------
C
C
C     Matrix generation without any corrections, hence equations make Phi or normal flows
C
      IMPLICIT NONE
      INTEGER(4) M,N,J,ITYPE,IW,IST,ILAST,IEN,I,IEQS,IEQ,J1,I1,II
      LOGICAL LNEG
      REAL(8) DRA,RQW3I,RQW3E,RFW3MATC,RAD,RFNFW3CO
      COMPLEX(8) CZI,CALPH,CZ,CZA,CZ0
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION RQW3I(3),RQW3E(3)
      DIMENSION DRA(M,*),CZI(*),CALPH(*),ITYPE(*)
      IF(NW3.EQ.0) RETURN
      DO 50 IW=1,NW3
      IST=IPNT(IW)
      ILAST=IPNT(IW+1)-1
      IEN=ILAST
      IF (LW3Q(IW)) IEN=ILAST-1   ! one equation less than the number of control points when Q specified
      DO 40 I=1,M
      CZ=CZI(I)
      CZA=CALPH(I)
      R3DZ=R3DCPZ(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between CZ & CZA, positive to the left when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request for zero matrix coefficient
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      J1=J
      GOTO (1,2,3,4,5,1),IEQ
    1 DO II=IST,IEN       ! provide potential at CZ
      J1=J1+1
      R3DZ=R3DCPZ(I)
      DRA(I,J1)=RFW3MATC(II,IW,CZ)
      IF (LW3Q(IW)) DRA(I,J1)=DRA(I,J1)-   ! necessary to get last strength parameter included in \Phi
     &              RW3L(II)/RW3L(ILAST)*RFW3MATC(ILAST,IW,CZ) ! done to avoid a continuity equation
      END DO
      GOTO 9
    2 CONTINUE
      IF (LNEG) THEN
        CZ0=CMPLX(RW3ST(1,IST),RW3ST(2,IST))
        RAD=RW3RAD(IW)
        DO II=IST,IEN
        J1=J1+1
        IF (.NOT.LW3Q(IW)) THEN
          DRA(I,J1)=RW3L(II)*RFNFW3CO(CZ0,RAD,CZ,CZA) ! add normal flow across CZ and CZA
        END IF
        END DO
        GOTO 9
      ELSE  ! when discharge specified flow across CZ and CZA is part of "given" functions in RFNFW3
C ----------------provide PSI: not yet implemented
      GOTO 9
      END IF
    3 CONTINUE
C ----------------discharge normal and parallel unit vector CZA: not yet implemented
      GOTO 9
    4 CONTINUE     ! then what is this ????? henk 7/9/03
      IF (LW3Q(IW)) CALL W3QCO (ILAST,IW,CZ,RQW3E)
      DO II=IST,IEN
      CALL W3QCO (II,IW,CZ,RQW3I)
      J1=J1+1
      DRA(I,J1)=DRA(I,J1)+RQW3I(3)
      IF (LW3Q(IW)) DRA(I,J1)=DRA(I,J1)-RW3L(II)/RW3L(ILAST)*RQW3E(3)
      END DO
      GOTO 9
    5 CONTINUE     ! provide total discharge for continuity equation
      IF (LW3Q(IW)) GOTO 9
      DO II=IST,IEN
      J1=J1+1
      DRA(I,J1)=DRA(I,J1)+RW3L(II)
      END DO
      GOTO 9
   9  CONTINUE
  40  CONTINUE
      J=J1
!      WRITE (ILUME,1000) N,J1
  50  CONTINUE
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)      
      END
C
C --------------------------------------------------------------------
C
      SUBROUTINE W3MATCORRECT (DRA,CZI,M,N,J,CALPH,ITYPE)     !  TO be done
C
C --------------------------------------------------------------------
C
C
C     Matrix
C
      IMPLICIT NONE
      INTEGER(4) M,N,J,ITYPE,IW,IST,ILAST,IEN,I,IEQS,IEQ,J1,I1,II
      LOGICAL LNEG
      REAL(8) DRA,RQW3I,RQW3E,RFW3MATC,RAD,RFNFW3CO
      COMPLEX(8) CZI,CALPH,CZ,CZA,CZ0
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION RQW3I(3),RQW3E(3)
      DIMENSION DRA(M,*),CZI(*),CALPH(*),ITYPE(*)
      IF(NW3.EQ.0) RETURN
      DO 50 IW=1,NW3
      IST=IPNT(IW)
      ILAST=IPNT(IW+1)-1
      IEN=ILAST
      IF (LW3Q(IW)) IEN=ILAST-1  ! one equation less then number of control points when Q specified
      DO 40 I=1,M
      CZ=CZI(I)
      CZA=CALPH(I)
      R3DZ=R3DCPZ(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across line between Cz & CZA, positive to the left when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request for zero matrix coefficient
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      J1=J
      GOTO (1,9,9,9,9,6),IEQ
    1 DO II=IST,IEN       ! provide potential at CZ
      J1=J1+1
      IF (LNEG) THEN         ! subtract potential at CZA
       R3DZ=R3DZA(I)
       DRA(I,J1)=DRA(I,J1)-RFW3MATC(II,IW,CZA)
       IF (LW3Q(IW)) DRA(I,J1)=DRA(I,J1)+
     &               RW3L(II)/RW3L(ILAST)*RFW3MATC(ILAST,IW,CZA)
      ENDIF
      END DO
      GOTO 9
   6  DO II=IST,IEN
      J1=J1+1
      DRA(I,J1)=0.0D+0
      END DO
      GOTO 9
   9  CONTINUE
  40  CONTINUE
      J=J1
!      WRITE (ILUME,1000) N,J1
  50  CONTINUE
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)      
      END

C
C --------------------------------------------------------------------
C      
      SUBROUTINE W3KNO (DRB,J,CZI,CALPH,ITYPE)
C
C --------------------------------------------------------------------
C
C
C     Known vector
C      
      IMPLICIT NONE
      INTEGER(4) J,ITYPE,IW,IST,IEN,I,IEQS,IEQ,II
      LOGICAL LNEG
      REAL(8) DRB,RFPOT,RHEAD,RPOT,RFPOTH,rdum1,rdum2
      COMPLEX(8) CZI,CALPH,CZ,CZA
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION DRB(*),CZI(*),CALPH(*),ITYPE(*)
      IF (NW3.EQ.0) RETURN
      DO 20 IW=1,NW3
      IST=IPNT(IW)
      IEN=IPNT(IW+1)-1
      IF (LW3Q(IW)) IEN=IEN-1 ! one equation less then number of control points when Q specified
      DO 10 I=IST,IEN
      J=J+1
      CZ=CZI(J)
      CZA=CALPH(J)
      R3DZ=R3DCPZ(J)
      IEQS=ITYPE(J)
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      GOTO (1,9,9,9,9,9),IEQ
   1  IF (LNEG) THEN
c      rdum1=rfpot(cz)
      rdum2=dw3pot(i)
      DRB(J)=DRB(J)-rdum2
      R3DZ=R3DZA(J)
c      rdum1=rfpot(cza)
      if (i.ne.ien) then
       rdum2=dw3pot(i+1)
      else
       rdum2=rfpot(cza)  ! no \Phi equation for last control point when Q specified
      end if
      DRB(J)=DRB(J)+rdum2
      R3DZ=R3DCPZ(J)
      ELSE   
      RHEAD=RW3HED(IW)
      RPOT=RFPOTH(RHEAD,CZ)
c      rdum1=rfpot(cz)
      rdum2=dw3pot(i)
      DRB(J)=DRB(J)+RPOT-rdum2
      ENDIF
c      write (ilume,1001) lneg,iw,i,rdum1,rdum2
c 1001 format ('w3kno1: lneg,iw,i,rdum1,rdum2 ',l3,2(i3,1x),2(d14.7,1x))
      GOTO 10
   9  CONTINUE
      write (ilume,1000) itype(j)
  10  CONTINUE
  20  CONTINUE
      RETURN
 1000 FORMAT (' ***ERROR: Illegal equation type in W3KNO: ITYPE=',I5)
      END
C
C --------------------------------------------------------------------
C
      SUBROUTINE W3SUB (DRB,J)
C
C --------------------------------------------------------------------
C
C
C     Substitution of strength parameters
C
      IMPLICIT NONE
      INTEGER(4) J,IW,IST,ILAST,IEN,I
      REAL(8) DRB,RS
      INCLUDE 'W3COM.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION DRB(*)
      IF(NW3.EQ.0) RETURN
      DO 20 IW=1,NW3
      IST=IPNT(IW)
      ILAST=IPNT(IW+1)-1
      IEN=ILAST
      IF (LW3Q(IW)) IEN=ILAST-1
      RS=RW3Q(IW)
      DO 10 I=IST,IEN
      J=J+1
      RW3S(I)=RW3S(I)+DRB(J)
      IF (LW3Q(IW)) RS=RS-RW3S(I)*RW3L(I)
  10  CONTINUE
      IF (LW3Q(IW)) THEN
            RS=RS-RW3S(ILAST)*RW3L(ILAST)
            RW3S(ILAST)=RW3S(ILAST)+RS/RW3L(ILAST)
      ELSE
            RW3Q(IW)=0.0
            DO 15 I=IST,IEN
            RW3Q(IW)=RW3Q(IW)+RW3S(I)*RW3L(I)
  15        CONTINUE              
      ENDIF      
  20  CONTINUE
      RETURN
      END
C
C --------------------------------------------------------------------
C      
      SUBROUTINE W3UPDATE_check (J,M,N,CZI,ITYPE)
C
C --------------------------------------------------------------------
C
C
C     Known vector
C      
      IMPLICIT NONE
      INTEGER(4) J,ITYPE,IW,IST,IEN,I,IEQS,IEQ,II,M,N
      LOGICAL LNEG
      REAL(8) RFPOT,rdum1,rdum2
      COMPLEX(8) CZI,CZ
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION CZI(m),ITYPE(m)
      IF (NW3.EQ.0) RETURN
      DO 20 IW=1,NW3
      IST=IPNT(IW)
      IEN=IPNT(IW+1)-1
      IF (LW3Q(IW)) IEN=IEN-1 ! one equation less then number of control points when Q specified
      DO 10 I=IST,IEN
      J=J+1
      CZ=CZI(J)
      R3DZ=R3DCPZ(J)
      IEQS=ITYPE(J)
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      GOTO (1,9,9,9,9,9),IEQ
   1  rdum1=rfpot(cz)
      rdum2=dw3pot(i)
      write (ilume,1000) j,rdum1,rdum2
 1000 format ('w3update_check: j,rfpot,dw3pot(i) ',/,
     &         i4,2(d14.7,2x))
      GOTO 10
   9  CONTINUE
      WRITE (ILUER,2000) ITYPE(J)
  10  CONTINUE
  20  CONTINUE
      RETURN
 2000 FORMAT('***ERROR: Illegal equation type in W3UPDATE_check: ITYPE='
     &       ,I5)
      END
C
C --------------------------------------------------------------------
C      
      SUBROUTINE W3UPDATE (DRSCR,J,M,N,CZI,CALPH,ITYPE,lincludesub,
     &                     isubsign)
C
C --------------------------------------------------------------------
C
C
C     Known vector
C      
      IMPLICIT NONE
      INTEGER(4) J,ITYPE,IW,IST,IEN,I,IEQS,IEQ,II,M,N,
     &           nsolOut,isubsign
      LOGICAL LNEG,lincludesub,laddsubcells,
     &           lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut
      REAL(8) DRSCR,RFPOT,RHEAD,RPOT,RFPOTH,rsubsign
      COMPLEX(8) CZI,CALPH,CZ,CZA,cflk_subomega
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION DRSCR(*),CZI(*),CALPH(*),ITYPE(*)
      IF (NW3.EQ.0) RETURN
      laddsubcells=lincludesub
      rsubsign=REAL(isubsign)
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      DO 20 IW=1,NW3
      IST=IPNT(IW)
      IEN=IPNT(IW+1)-1
      IF (LW3Q(IW)) IEN=IEN-1 ! one equation less then number of control points when Q specified
      DO 10 I=IST,IEN
      J=J+1
      CZ=CZI(J)
      CZA=CALPH(J)
      R3DZ=R3DCPZ(J)
      IEQS=ITYPE(J)
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      GOTO (1,9,9,9,9,9),IEQ
   1  if (nsolOut.eq.0) then
        dw3pot(i)=rfpot(cz)
      else
        dw3pot(i)=dw3pot(i)+drscr(j)
        if (laddsubcells) then
          dw3pot(i)=dw3pot(i)+rsubsign*REAL(cflk_subomega(cz)) ! sub-cell contributions not in drscr(j)
        end if
      endif
      GOTO 10
   9  CONTINUE
      WRITE (ILUER,1000) ITYPE(J)
  10  CONTINUE
  20  CONTINUE
      RETURN
 1000 FORMAT (' ***ERROR: Illegal equation type in W3UPDATE: ITYPE=',I5)
      END

C
C -------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFW3MATC (INODE,IW,CZ)
C
C     Generates coefficient for ppwell strength parameters 
C     Used in W3MAT only
C      
      IMPLICIT NONE
      INTEGER(4) IW,INODE
      REAL(8) RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &        RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU,
     &        RHOI,RXI,RHGHT,RFBASE,RFTOP,RF3DSP,
     &        RADTST,RADTST2,RFW3CP,RL,RGVREFDIST
      COMPLEX(8) CZ,CFW3OMCOF,CZ0
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      COMMON /W3PASS/RZT,RZB,RHO2,RT,RH,RH2,R2H,
     &               RPI1,RPI2,RPI4,RPI8,RTAU2,RTAU
      DIMENSION RHOI(3),RXI(3)
      RHGHT=RFTOP(CZ)-RFBASE(CZ)
      RXI(1)=REAL(CZ)
      RXI(2)=AIMAG(CZ)
      RHOI(1)=RXI(1)-RW3ST(1,INODE)
      RHOI(2)=RXI(2)-RW3ST(2,INODE)
      RHOI(3)=0.0
      RHO2=RF3DSP(RHOI,RHOI)
      RADTST=RZONEMLTP*RHGHT
      RADTST2=RADTST*RADTST
      IF (RHO2.LE.RADTST2) THEN     ! 3D zone
      RFW3MATC=RFW3CP (INODE,IW,CZ)
      ELSE
      cz0=CMPLX(RW3ST(1,INODE),RW3ST(2,INODE))
      rl=rgvrefdist(cz0)                          ! 2D zone
      RFW3MATC=RW3L(INODE)*REAL(CFW3OMCOF(RHO2,rl)) ! rl used to be rhght
      ENDIF
      RETURN
      END

