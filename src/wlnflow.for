C     Last change:  HMH  18 Dec 2006    3:06 pm
C     This file contains the following routines and functions
C
C     REAL FUNCTION RFNFWL    returns the flow across a line due to all wells
C     REAL FUNCTION RFNFWLCO  coefficient function for the flow across a line due to well i.
C     
C
C ------------------------------------------------------------------------------
C
      REAL FUNCTION RFNFWL(CZ1,CZ2)
C
C ------------------------------------------------------------------------------
C
C     Returns the flow across a line between CZ1 and CZ2 due to all wells.
C     The flow is positive if from left to right when seen from CZ1
C
      IMPLICIT NONE
      INCLUDE 'WLCOM.INC'
      INCLUDE 'TRACOM.INC'
      INTEGER(4) I
      REAL(8) RFNFWLCO
      COMPLEX(8) CZ1,CZ2
C
      RFNFWL=0.0
      IF (NWL.EQ.0) RETURN
      DO I=1,NWL
      RFNFWL=RFNFWL+RWLQ(I)*RFNFWLCO(I,CZ1,CZ2)
      END DO
      RETURN
      END
c
c ------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFNFWLCO(I,CZ1,CZ2)
C
C ------------------------------------------------------------------------------
C
c     Coefficient function for the flow between CZ1 and CZ2 due to well I.
c     Adapted from Vic Kelson (test SLWL)
c
      IMPLICIT NONE
      INTEGER(4) I
      LOGICAL LINSECTCIRCLE
      REAL(8) RFBRANCH,RDIS,CABS
      COMPLEX(8) CZ1,CZ2,CFWLOMC,CDUM1,CDUM2
      INCLUDE 'WLCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
c
      RFNFWLCO=0.0
c     test for point inside well and return zero flow if true
      RDIS=ABS(CZ1-CWLZ(I))
      IF (RDIS.LT.RWLR(I)) RETURN ! cz1 is inside, segment is radial for well; return 0
      RDIS=ABS(CZ2-CWLZ(I))
      IF (RDIS.LT.RWLR(I)) RETURN ! cz2 is inside, segment is radial for well; return 0
      IF (LINSECTCIRCLE(CWLZ(I),RWLR(I),CZ1,CZ2,CDUM1,CDUM2))  RETURN ! if line intersects well: return 0
c     calculate flow across line through CZ1 and CZ2
      RFNFWLCO=AIMAG(cfwlomc(cz1,i))-AIMAG(cfwlomc(cz2,i)) +
     &        rfbranch(cwlz(i),cz1,cz2)  ! add flow in branch cut
      end
