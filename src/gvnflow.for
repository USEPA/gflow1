C     Last change:  HMH   2 Jan 2001   10:17 pm
c       This file contains the following routines and functions:
c
c       RFNFGV   calculates flow across line from CZ1 to CZ2 due to uniform flow
C
C --------------------------------------------------------------------------------------------
C
	REAL(8) FUNCTION RFNFGV(CZ1,CZ2)
C
C --------------------------------------------------------------------------------------------
C
C     Return flow across line from CZ1 to CZ2 due to uniform flow component
C
      IMPLICIT NONE
      COMPLEX(8) CZ1,CZ2,CFUFOM
      RFNFGV=AIMAG(CFUFOM(CZ1))-AIMAG(CFUFOM(CZ2))
      END
