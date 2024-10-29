C     Last change:  HMH  14 Mar 2003   12:13 pm
C     This file contains the following routines and functions
C
C     Logical FUNCTION lwlinfo    routine returns information for the 2D well
C
C ------------------------------------------------------------------------------
C
      Logical FUNCTION lwlinfo (CZ,CZWELL,RWELL,RQWELL,AWELL)
C
C ------------------------------------------------------------------------------
C
C     Routine returns information for the 2D well at CZ
C     Routine is called by "lwellinfo" in gftrace.for
C
      implicit none
      INTEGER(4) iw
      REAL(8) rwell,rqwell,rdis
      COMPLEX(8) cz,czwell
      CHARACTER(16) AWELL
      INCLUDE 'wlcom.inc'
      INCLUDE 'lusys.inc'
      lwlinfo=.false.
      IF (NWL.EQ.0) RETURN
      DO IW=1,NWL
      RDIS=ABS(CZ-CWLZ(IW))
      IF (RDIS.LT.rwlr(iw)) THEN ! within the radius, close enough
      CZWELL=CWLZ(IW)
      RWELL=RWLR(IW)
      RQWELL=RWLQ(IW)
      AWELL=AWLAB(IW)
      lwlinfo=.true.
      return
      ENDIF
      END do
      RETURN
      END
