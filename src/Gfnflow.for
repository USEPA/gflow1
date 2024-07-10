C     Last change:  HMH  16 Nov 2008    2:46 am
c     This file contains the following routines and functions
c
c     REAL function rfnormalflow      driver for calculating the flow across a straight line element
c     REAL FUNCTION RFBRANCH          returns jump in PSI across line for a unit withdrawal at a point
c     SUBROUTINE BRANCHCUT            returns jump in PSI across line for linesink or linedipole
c     REAL FUNCTION RFAREATRIANGLE    returns area for triangle
c     REAL FUNCTION RFNUMNF           returns a numerical approximation of the flow across a straight line element
c     LOGICAL FUNCTION LINSECTCIRCLE  true if line intersects a circle
c     LOGICAL FUNCTION LINSECTLINE    true if two lines intersect
c     REAL*8 FUNCTION DFA             returns signed area of triangle 
c
c
c
c-------------------------------------------------------------------------------------
c
      REAL(8) function rfnormalflow(cz1,cz2)
c
c-------------------------------------------------------------------------------------
c
c
c      Function returns the total flow across a line between cz1 and cz2.
c      The flow is positive from left to right across the line when viewed from cz1.
c      The function is called by subroutine "extract".
c      NOTE: the transient wells are not represented in this function.
c
      implicit none
      REAL(8) rfnfgv,rfnfpd,rfnfdb,rfnfls,rfnfwl,rfnflk,rfnflksub, ! 2D functions
     &     rfnfdi,rfnfw3                       ! 3D functions
      COMPLEX(8)  cz1,cz2
      INCLUDE 'MAIN.INC'
      INCLUDE 'lusys.inc'
c
      RFNORMALFLOW=0.0
      RFNORMALFLOW=RFNORMALFLOW+RFNFGV(CZ1,CZ2) ! uniform flow
      RFNORMALFLOW=RFNORMALFLOW+RFNFPD(CZ1,CZ2) ! sink discs (2D)
      RFNORMALFLOW=RFNORMALFLOW+RFNFWL(CZ1,CZ2) ! wells (steady state, 2D)
      RFNORMALFLOW=RFNORMALFLOW+RFNFLS(CZ1,CZ2) ! line sinks
      RFNORMALFLOW=RFNORMALFLOW+RFNFDB(CZ1,CZ2) ! line doublets/dipoles
      RFNORMALFLOW=RFNORMALFLOW+RFNFW3(CZ1,CZ2) ! partially penetrating wells
      RFNORMALFLOW=RFNORMALFLOW+RFNFDI(CZ1,CZ2) ! partially penetrating wells
      RFNORMALFLOW=RFNORMALFLOW+RFNFLK(CZ1,CZ2) ! leakage and recharge grid
      RFNORMALFLOW=RFNORMALFLOW+RFNFLKSUB(CZ1,CZ2) ! leakage in subgrids
      return
      end
c
c-------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFBRANCH(CBRANCH,CZ1,CZ2)
c
c-------------------------------------------------------------------------------------
c
C
C     Returns branch cut in PSI for a unit withdrawal (well with pumping rate 1.0) at CBRANCH.
C
      IMPLICIT NONE
      REAL(8) RXW,RYW,RXINT,RX1,RY1,RX2,RY2
      COMPLEX(8) CBRANCH,CZ1,CZ2
      include 'lusys.inc'
C
      RXW=REAL(CBRANCH)
      RYW=AIMAG(CBRANCH)
      RFBRANCH=0.0
      RXINT=0.0
      RX1=REAL(CZ1)
      RY1=AIMAG(CZ1)
      RX2=REAL(CZ2)
      RY2=AIMAG(CZ2)
      IF ((RY1.GE.RYW).AND.(RY2.LT.RYW)) THEN  ! CZ1 above and CZ2 below branch cut line
        RXINT=RX1-(RY1-RYW)*(RX2-RX1)/(RY2-RY1)
        IF (RXINT .LT. RXW) RFBRANCH=-1.0      ! intersection point on branch cut (left of CBRANCH)
      ENDIF
      IF ((RY1.LT.RYW).AND.(RY2.GE.RYW)) THEN  ! CZ1 below and CZ2 above branch cut line
        RXINT=RX1-(RY1-RYW)*(RX2-RX1)/(RY2-RY1)
        IF (RXINT .LT. RXW) RFBRANCH=1.0       ! intersection point on branch cut (left of CBRANCH)
      ENDIF
!      WRITE (ILUME,1001) CBRANCH,CZ1,CZ2,RFBRANCH
! 1001 FORMAT (' leaving RFBRANCH: cbranch,cz1,cz2,rfbranch ',
!     &        3(2(e14.7),2x),2x,e14.7)
      RETURN
      END
c
c-------------------------------------------------------------------------------------
c
      SUBROUTINE BRANCHCUT(CZ1,CZ2,CZS,CZE,CZ0,RBRANCH,LBRANCH)
c
c-------------------------------------------------------------------------------------
c
C
C     LBRANCH=.TRUE. if line CZ1 - CZ2 intersects the line CZS - CZE at CZ0,
c     where CZ0 is between CZ1 and CZ2 and between CZS and CZE
c     RBRANCH contains the flow in the branch cut for a constant strenght line sink,
c     which is on the line element CZS and CZE the same as for a linear strength line dipole.
c     Functions is called by RFBCLSU (in LSFUN.FOR) and DBINTERSECT (in DBFUN.FOR).
c     Input:
c            CZ1        starting point of line segment
C            CZ2        end point of line segment
C            CZS        starting point of line sink or line dipole
C            CZE        end point of line sink or line dipole
C     Output:
C            CZ0        intersection point on line sink or line dipole
C            RBRANCH    flow in branch cut for unit strength line sink
C            LBRANCH    TRUE when CZ0 between CZS and CZE (on line sink or line dipole)
C
      IMPLICIT NONE
      LOGICAL LBRANCH
      REAL(8) RBRANCH,RBX1,RBY1,RBX2,RBY2,RXINT
      COMPLEX(8) CZ1,CZ2,CZS,CZE,CZ0,CBZ1,CBZ2,CFBIGZ
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
      LBRANCH=.FALSE.
      RBRANCH=0.0
      CBZ1=CFBIGZ(CZ1,CZS,CZE) ! map CZ1 onto reference plane
      CBZ2=CFBIGZ(CZ2,CZS,CZE) ! map CZ2 onto reference plane
      RBX1=REAL(CBZ1)
      RBY1=AIMAG(CBZ1)
      RBX2=REAL(CBZ2)
      RBY2=AIMAG(CBZ2)
      IF ( RBY1 .GE. 0.0 .AND. RBY2 .LT. 0.0 .OR.
     .     RBY2 .GE. 0.0 .AND. RBY1 .LT. 0.0 ) THEN  ! straddling the X-axis
        IF ( RBX2 .NE. RBX1 ) THEN  ! calculate intersection with X-axis
          RXINT=RBX1-RBY1*(RBX2-RBX1)/(RBY2-RBY1)
        ELSE
          RXINT=RBX1
        ENDIF      
C       NOTE: The branch cut for a unit linesink varies from 2.0 at Z=(-1,0) to zero at Z=(1,0)
        IF (RXINT .LE. -1.0) THEN
          RBRANCH=2.0
        ENDIF
        IF (RXINT .GT. -1.0 .AND. RXINT .LE. 1.0 ) THEN ! between CZS and CZE
          RBRANCH=(1.0-RXINT)
          LBRANCH=.TRUE.
          CZ0=0.5*(CMPLX(RXINT,0.0)*(CZE-CZS)+CZS+CZE) ! map intersection back to physical plane
        ENDIF
        IF ( RBY1 .LT. 0.0 ) THEN ! note: was .GT. in SLWL routine
          RBRANCH=RBRANCH*ABS(CZE-CZS)*0.5 ! map branch cut jump back onto physical plane
        ELSE
          RBRANCH=-RBRANCH*ABS(CZE-CZS)*0.5 ! map branch cut jump back onto physical plane
        ENDIF
      ENDIF
      RETURN
      END
c
c-------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFAREATRIANGLE(CZ1,CZ2,CZ3)
c
c-------------------------------------------------------------------------------------
c
C
C     Returns area of a triangle spanned by the three coordinate pairs CZ1, CZ2 and CZ3
c     The area is signed and used to calculate the flow across CZ1 and CZ2 due to a sink disc
c     of constant sink density 1, centered at CZ3
c     Function is used for Sink Discs (2D & 3D) and recharge inhomogeneities.
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RA0,DR1,DR2,DS1,DS2,DW3
      COMPLEX(8) CZ1,CZ2,CZ3,CR,CS
      INCLUDE 'TRACOM.INC'
c     Vector CR is the vector from CZ1 to the center of the recharge circle
      CR = CZ3-CZ1      ! To center of circle
      DR1 = REAL(CR)
      DR2 = AIMAG(CR)
c     Vector CS is the vector from CZ1 to CZ2
      CS = CZ2-CZ1
      DS1 = REAL(CS)
      DS2 = AIMAG(CS)
c     DW3 is the z-component of the curl of R and S
      DW3 = DR1*DS2 - DR2*DS1
c     The signed area bounded by R and S is 1/2 of W3
      RFAREATRIANGLE = 0.5*DW3
      END
c
c-------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFNUMNF(CZ1,CZ2)
c
c-------------------------------------------------------------------------------------
c
C
C     Returns a numerical approximation of the total flow normal to the line CZ1 - CZ2
c     Used only for debugging: comparing results to that of "rfnormalflow".
C     Integration of Qn is done by Simpson's rule
C
      implicit none
      INTEGER(4) NSECT,I
      REAL(8) RDIS,RQ1,RQ2,RQ3,RQI
      COMPLEX(8) CZ1,CZ2,CZ,CZDIS,CZNORM
      DIMENSION RQI(3)

C
      NSECT=200   ! number of steps along line segment (must be even).
      CZDIS=(CZ2-CZ1)/NSECT
      RDIS=ABS(CZDIS)  ! delta s along line segment CZ1 - CZ2
      CZNORM=-(0.0,1.0)*CZDIS/RDIS ! unit normal vector pointing from left to right when at CZ1
      CZ=CZ1
      CALL DISCH(CZ,RQI)
      RQ3=REAL(CZNORM)*RQI(1)+AIMAG(CZNORM)*RQI(2) ! scalar product of unit normal with Qi
      RFNUMNF=0.0
C
      DO I=1,NSECT,2
      RQ1=RQ3
      CZ=CZ+CZDIS
      CALL DISCH(CZ,RQI)
      RQ2=REAL(CZNORM)*RQI(1)+AIMAG(CZNORM)*RQI(2) ! scalar product of unit normal with Qi
      CZ=CZ+CZDIS
      CALL DISCH(CZ,RQI)
      RQ3=REAL(CZNORM)*RQI(1)+AIMAG(CZNORM)*RQI(2) ! scalar product of unit normal with Qi
      RFNUMNF=RFNUMNF+RQ1+4*RQ2+RQ3
      END DO
      RFNUMNF=RDIS*RFNUMNF/3.0
      RETURN
      END
c
c-------------------------------------------------------------------------------------
c
	LOGICAL FUNCTION LINSECTCIRCLE(CZ0,RAD,CZ1,CZ2,CZ3,CZ4)
c
c-------------------------------------------------------------------------------------
c
C
C     Function is TRUE when the line through CZ1 and CZ2 intersects the circle with
C     radius RAD and center CZ0 and if the intersection points CZ3 and CZ4 are inbetween
C     the points CZ1 and CZ2. Function is called by PD, WL, W3 and DI modules.
C
      IMPLICIT NONE
      REAL(8) RAD,RM,RN,RA,RB,RC,RD,RY,RX0,RY0,RX1,RY1,RX2,RY2,
     &        RX3,RY3,RX4,RY4,RS1,RS2,RTEST
      COMPLEX(8) CZ0,CZ1,CZ2,CZ3,CZ4
      include 'lusys.inc'
C
c      write (ilume,1001) cz0,rad,cz1,cz2
c 1001 format (' Entering LINSECTCIRCLE: cz0,rad,cz1,cz2 ',/,
c     & 2(e14.7),2x,e14.7,2x,2(2(e14.7),2x))
	RX1=REAL(CZ1-CZ0)
      RY1=AIMAG(CZ1-CZ0)
      RX2=REAL(CZ2-CZ0)
      RY2=AIMAG(CZ2-CZ0)
C
      IF (ABS(RX2-RX1).LT.1.0E-10) THEN ! vertical line
       RD=(RX1+RX2)/2.0
       IF (ABS(RD).LT.RAD) THEN ! line intersects disc boundary
        RS1=SIGN(1.0,RY1)
        RS2=SIGN(1.0,RY2)
        IF (RS1*RS2.GT.0.0) THEN ! intersection not on segment CZ1 and CZ2
          LINSECTCIRCLE=.FALSE.
          RETURN
        END IF
c        Intersection on segment between CZ1 and CZ2
        RY=SQRT(RAD*RAD-RD*RD)
        CZ3=CMPLX(RD,RS1*RY)+CZ0
        CZ4=CMPLX(RD,RS2*RY)+CZ0
        LINSECTCIRCLE=.TRUE.
       ELSE
        LINSECTCIRCLE=.FALSE.
        RETURN
       END IF
      ELSE       ! sloping line
       RM=(RY2-RY1)/(RX2-RX1) ! slope
       RN=RY1-RM*RX1
       RA=RM*RM+1.0
       RB=2.0*RM*RN
       RC=RN*RN-RAD*RAD
       RD=RB*RB-4.0*RA*RC
       IF (RD.GT.0.0) THEN ! line intersects disc boundary
         IF (RX1.LT.RX2) THEN
           RX3=(-RB-SQRT(RD))/(2.0*RA)
           RX4=(-RB+SQRT(RD))/(2.0*RA)
         ELSE
           RX3=(-RB+SQRT(RD))/(2.0*RA)
           RX4=(-RB-SQRT(RD))/(2.0*RA)
         END IF
         RY3=RM*RX3+RN
         RY4=RM*RX4+RN
         CZ3=CMPLX(RX3,RY3)+CZ0
         CZ4=CMPLX(RX4,RY4)+CZ0
         RTEST=ABS(CZ3-CZ1)+ABS(CZ3-CZ2)-ABS(CZ2-CZ1)
         IF (ABS(RTEST).LT.0.01) THEN ! CZ3 between CZ1 and CZ2
           LINSECTCIRCLE=.TRUE.
           RETURN
         END IF
         LINSECTCIRCLE=.FALSE.
       ELSE
        LINSECTCIRCLE=.FALSE.
       END IF
      ENDIF
      RETURN
      END
c
c-------------------------------------------------------------------------------------
c
      LOGICAL FUNCTION LINSECTLINE(CZS,CZE,CZ1,CZ2,CZ0)
c
c-------------------------------------------------------------------------------------
c
C
C     Function is TRUE if the intersection point CZ0 is inbetween CZ1 & CZ2 and CZS & CZE
c     Note: CZ0 will only be set if intersection between CZ1 and CZ2
C
      IMPLICIT NONE
      REAL(8) RBX1,RBY1,RBX2,RBY2,RXINT
      COMPLEX(8) CZS,CZE,CZ1,CZ2,CZ0,CBZ1,CBZ2,CFBIGZ
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      CBZ1=CFBIGZ(CZ1,CZS,CZE) ! map CZ1 onto reference plane for CZS-CZE
      CBZ2=CFBIGZ(CZ2,CZS,CZE) ! map CZ2 onto reference plane for CZS-CZE
      RBX1=REAL(CBZ1)
      RBY1=AIMAG(CBZ1)
      RBX2=REAL(CBZ2)
      RBY2=AIMAG(CBZ2)
      LINSECTLINE=.FALSE.
      IF ( RBY1 .GE. 0.0 .AND. RBY2 .LT. 0.0 .OR. ! x-axis is now the line through CZS and CZE
     .     RBY2 .GE. 0.0 .AND. RBY1 .LT. 0.0 ) THEN  ! intersection is between CZ1 and CZ2
        IF ( RBX2 .NE. RBX1 ) THEN  ! calculate intersection with X-axis
          RXINT=RBX1-RBY1*(RBX2-RBX1)/(RBY2-RBY1)
        ELSE
          RXINT=RBX1
        ENDIF
        CZ0=CMPLX(RXINT,0.0)*0.5*(CZE-CZS)+0.5*(CZS+CZE) ! map intersection back to physical plane
        IF (RXINT.GE.-1.0.AND.RXINT.LE.1.0) THEN  ! intersection is also between CZS and CZE
        LINSECTLINE=.TRUE.
        END IF
      ENDIF
      RETURN
      END
c
c-------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION DFA(CZZ1,CZZ2)
c
c-------------------------------------------------------------------------------------
c
C
C     Function returns the signed area of the triangle spanned by the vectors stored on the
C     complex numbers CZZ1 and CZZ2. The sign of the area is determined using the right-hand
C     rule when going from vector CZZ1 to CZZ2.
C
      IMPLICIT NONE
      REAL(8) DX1,DY1,DX2,DY2
      COMPLEX(8) CZZ1,CZZ2
      DX1=REAL(CZZ1)
      DY1=AIMAG(CZZ1)
      DX2=REAL(CZZ2)
      DY2=AIMAG(CZZ2)
      DFA=0.5*(DX1*DY2-DX2*DY1)
      RETURN
      END
      


