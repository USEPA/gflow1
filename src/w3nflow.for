C     Last change:  HMH  18 Dec 2006    3:03 pm
C     This file contains the following routines and functions
C
C     REAL FUNCTION RFNFW3    returns the flow across a line due to all wells
C     REAL FUNCTION RFNFW3CO  coefficient function for the flow across a line due to well i.
C     
C
C ------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFNFW3(CZ1,CZ2)
C
C ------------------------------------------------------------------------------------------------
C
C     Function returns the total flow, integrated over the aquifer height, across the
C     line between CZ1 and CZ2 due to a partially penetrating well.
C
      IMPLICIT NONE
      INTEGER(4) I,INODE
      LOGICAL LINSECTCIRCLE
      REAL(8) RFBRANCH,RFNFW3CO
      COMPLEX(8) CZ,CZ1,CZ2,C1,C2,CDUM1,CDUM2
      INCLUDE 'W3COM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      RFNFW3=0.0
      IF (NW3.EQ.0) RETURN
      DO I=1,NW3
      INODE=IPNT(I)
      CZ=CMPLX(RW3ST(1,INODE),RW3ST(2,INODE))
      RFNFW3=RFNFW3+RW3Q(I)*RFNFW3CO(CZ,RW3RAD(I),CZ1,CZ2)
      END DO
      RETURN
      END
C
C ------------------------------------------------------------------------------------------------
C
      REAL(8) FUNCTION RFNFW3CO(CZ,RAD,CZ1,CZ2)
C
C ------------------------------------------------------------------------------------------------
C
C     Function returns flow across the line between CZ1 & CZ2 due to a ppwell at CZ with
C     radius RAD and pumping rate 1.
C
      IMPLICIT NONE
      INTEGER(8) I,INODE
      LOGICAL LINSECTCIRCLE
      REAL(8) RFBRANCH,RAD
      COMPLEX(8) CZ,CZ1,CZ2,C1,C2,CDUM1,CDUM2
      INCLUDE 'W3COM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      RFNFW3CO=0.0
      C1=CZ1-CZ
      IF (ABS(C1).LT.RAD) RETURN ! if CZ1 inside well do not add flow for this line segment
      C2=CZ2-CZ
      IF (ABS(C2).LT.RAD) RETURN ! if CZ2 inside well do not add flow for this line segment
      IF (LINSECTCIRCLE(CZ,RAD,CZ1,CZ2,CDUM1,CDUM2)) RETURN ! if line intersects well radius do not add flow
      RFNFW3CO=AIMAG(LOG(C1)-LOG(C2))/6.2831853 ! add Delta PSI
      RFNFW3CO=RFNFW3CO+RFBRANCH(CZ,CZ1,CZ2)  ! add flow in branch cuts
      RETURN
      END
