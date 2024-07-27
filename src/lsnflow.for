C     Last change:  HMH  18 Dec 2006    2:57 pm
c       This file contains the following routines and functions:
c
C	REAL FUNCTION RFNFLS      calculates the flow between two points due to all line-sinks
C	REAL FUNCTION RFNFLSCO    coefficient functions for one line-sink with sink density 1
C	REAL FUNCTION RFBCLSU     calculates the jump in PSI for a branch cut due to line-sink I
c
C
c --------------------------------------------------------------------------------------------------
C
      REAL FUNCTION RFNFLS(CZ1,CZ2)
C
c --------------------------------------------------------------------------------------------------
C
c     Return flow across line between CZ1 and CZ2 for all line sinks
c
      IMPLICIT NONE
      INTEGER I
      REAL RFNFLSCO,RNF
      COMPLEX CZ1,CZ2
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
      RFNFLS = 0.0
      IF (NLS.EQ.0) RETURN
      DO  I=1,NLS
      RNF=RLSIG(I)*RFNFLSCO(I,CZ1,CZ2)
      RFNFLS=RFNFLS+RNF
c      write (ilume,1001) i,rlsig(i),rnf,rfnfls
c 1001 format (' RFNFLS: i,rlsig,rnf,rfnfls ',i4,2x,e14.7,2x,2(e16.9,2x))
      END DO
      RETURN
      END
C
c --------------------------------------------------------------------------------------------------
C
      REAL FUNCTION RFNFLSCO(I,CZ1,CZ2)
C
c --------------------------------------------------------------------------------------------------
C
c     Coefficient function for flow across line between CZ1 and CZ2
c     due to line sink I with sink density 1
c
      IMPLICIT NONE
      INTEGER I
      REAL RFBCLSU
      COMPLEX CZ1,CZ2,COMLS
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
      RFNFLSCO = AIMAG(COMLS(CZ1,CLSZS(I),CLSZE(I))) -
     .           AIMAG(COMLS(CZ2,CLSZS(I),CLSZE(I))) ! calculate Delta Psi
C      Note: the constant CLSCONST(I) will cancel out, hence is not added to COMLS
c      write (ilume,1001) i,cz1,cz2,rfnflsco
c 1001 format (' RFNFLSCO: i,cz1,cz2, delta psi ',/,
c     . i4,2x,2(2(e14.7),2x),e16.9)
      RFNFLSCO=RFNFLSCO+RFBCLSU(I,CZ1,CZ2) ! add flow in branch cut
c      write (ilume,1002) i,cz1,cz2,rfnflsco
c 1002 format (' RFNFLSCO: i,cz1,cz2, delta psi and flow in branch cut ',/,
c     . i4,2x,2(2(e14.7),2x),e16.9)
      RETURN
      END
C
c --------------------------------------------------------------------------------------------------
C
      REAL FUNCTION RFBCLSU(I,CZ1,CZ2)
C
c --------------------------------------------------------------------------------------------------
C
c     Real function, returns the branch cut in the stream function from
c     the first point to the second for a single linesink of unit strength.
c
      IMPLICIT COMPLEX(C), LOGICAL(L)
      INCLUDE 'LSCOM.INC'
      INCLUDE 'TRACOM.INC'
      RFBCLSU=0.0
      call BRANCHCUT(CZ1,CZ2,CLSZS(I),CLSZE(I),CZ0,RBRANCH,LBRANCH)
      RFBCLSU=RBRANCH
      RETURN
      END
