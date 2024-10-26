C     Last change:  HMH  18 Dec 2006    2:59 pm
c     This file contains the following routines and functions
c
c     REAL FUNCTION RFNFPD      returns flow across a straight line element
c     REAL FUNCTION RFNFPDCO    normal flow coefficient function
c     REAL FUNCTION RFDPPDI     normal flow coefficient function when inside disc
c     REAL FUNCTION RFDPPDO     normal flow coefficient function when outside disc
c
c
c
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFNFPD(CZ1,CZ2)
C
C ----------------------------------------------------------------------------
C
c
c     Return all LAPLACE and POISSON contributions to the "Delta-PSI".
c     Includes branch cut correction.
c
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RFNFPDCO
      COMPLEX(8) CZ1,CZ2
      INCLUDE 'pdcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
c
      RFNFPD = 0.0
      IF (NPD.EQ.0) RETURN
      DO I=1,NPD
      RFNFPD=RFNFPD+RPDS(I)*RFNFPDCO(I,CZ1,CZ2)
      END DO
      RETURN
      END
c
c --------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFNFPDCO(I,CZ1,CZ2)
C
C ----------------------------------------------------------------------------
C
c
c     "Delta PSI" coefficient function for sink disc (pond) number i - branch cut.
c     Adapted from Vic Kelson (test SLWL)
c
      IMPLICIT none
      INTEGER(4) I
      LOGICAL LINSECTCIRCLE
      REAL(8) RAD,RR1,RR2,RA,RB,RC,RD,CABS,RFSCALAR,RFDPPDO,RFDPPDI
      COMPLEX(8) CZ1,CZ2,CZ3,CZ4,CV12,CV21,CR1,CR2
      INCLUDE 'pdcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
c
      RFNFPDCO = 0.0
      RAD=RPDR(I)
      CR1=CZ1-CPDZ(I)
      RR1=ABS(CR1)
      CR2=CZ2-CPDZ(I)
      RR2=ABS(CR2)
      IF ( RR1.LE.RAD .AND. RR2.LE.RAD ) THEN ! CZ1 and CZ2 inside
          RFNFPDCO=RFDPPDI(I,CZ1,CZ2)
          RETURN
      ELSEIF ( RR1.GT.RAD .AND. RR2.LT.RAD) THEN  ! CZ1 outside, CZ2 inside
          CV12=CZ2-CZ1 
          RA=ABS(RFSCALAR(CR1,CV12)/ABS(CV12))
          RB=ABS(CR1)
          RC=SQRT(RB*RB-RA*RA)
          RD=SQRT(RAD*RAD-RC*RC)
          CZ3=(RA-RD)*CV12/ABS(CV12)+CZ1 ! intersection of CZ1 - CZ2 with boundary of disc
          RFNFPDCO=RFDPPDO(I,CZ1,CZ3)+RFDPPDI(I,CZ3,CZ2)
          RETURN
      ELSEIF ( RR2.GT.RAD .AND. RR1.LT.RAD) THEN  ! CZ1 inside, CZ2 outside
          CV21=CZ1-CZ2
          RA=ABS(RFSCALAR(CR2,CV21)/ABS(CV21))
          RB=ABS(CR2)
          RC=SQRT(RB*RB-RA*RA)
          RD=SQRT(RAD*RAD-RC*RC)
          CZ3=(RA-RD)*CV21/ABS(CV21)+CZ2 ! intersection of CZ1 - CZ2 with boundary of disc
          RFNFPDCO=RFDPPDI(I,CZ1,CZ3)+RFDPPDO(I,CZ3,CZ2)
          RETURN
      ELSEIF ( RR1.GT.RAD .AND. RR2.GT.RAD ) THEN ! CZ1 and CZ2 outside
          IF (LINSECTCIRCLE(CPDZ(I),RAD,CZ1,CZ2,CZ3,CZ4)) THEN ! line segment intersect disc boundary
            RFNFPDCO=
     &      RFDPPDO(I,CZ1,CZ3)+RFDPPDI(I,CZ3,CZ4)+RFDPPDO(I,CZ4,CZ2)
          ELSE  ! line segment entirely outside disc
            RFNFPDCO=RFDPPDO(I,CZ1,CZ2)
          END IF
      ENDIF
      RETURN
      END
c
c ------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDPPDI(I,CZ1,CZ2)
C
C ----------------------------------------------------------------------------
C
C
C     "Delta-PSI" equivalent for a sink disc with both CZ1 and CZ2 inside the disc.
c     Function returns the "delta-PSI" for disc i with unit sink density.
c     NOTE: there is no test that CZ! and CZ2 are actually inside the disc!!
c     Adapted from Vic Kelson (test SLWL)
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RA0,RFAREATRIANGLE,DR1,DR2,DS1,DS2,DW3
      COMPLEX(8) CZ1,CZ2,CR,CS
      INCLUDE 'pdcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      RFDPPDI=RFAREATRIANGLE (CZ1,CZ2,CPDZ(i))
      RETURN
      END

C
C ------------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFDPPDO(I,CZ1,CZ2)
C
C ----------------------------------------------------------------------------
C
C
c     real function returns the Delta-PSI for a sink disc - branch cut, both CZ1 and CZ2 are outside the disc
c     Delta-PSI is calculated for disc i with unit sink density.
c     NOTE: there is no test for CZ1 and CZ2 to be actually outside the disc!!
C
      IMPLICIT NONE
      INTEGER(4) I
      REAL(8) RFBRANCH,RPI,RAD
      COMPLEX(8) CZ1,CZ2,CR1,CR2
      INCLUDE 'pdcom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      DATA RPI /3.1415926/
c     Compute the delta-psi for the sink disc being on the outside.
C     NOTE: The coefficient function CFPDCO is not used to avoid a point just inside
C     the disc (numerical inaccuracies) to miss a PSI allocation.
      RAD=RPDR(I)
      CR1=CZ1-CPDZ(I)
      CR2=CZ2-CPDZ(I)
      RFDPPDO=0.5*RAD*RAD*AIMAG(LOG(CR1/RAD)-LOG(CR2/RAD))
C     Add flow in branch cut of PSI, in case CZ1 - CZ2 crosses branch cut.
      RFDPPDO=RFDPPDO+RPI*RPDR(I)*RPDR(I)*RFBRANCH(CPDZ(I),CZ1,CZ2)
      RETURN
      END

