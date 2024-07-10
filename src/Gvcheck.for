C     Last change:  HMH  23 Jun 2005    5:46 pm
c     This file contains the following routines and functions:
c
c     GVERROR   reports the error in head at the reference point
C
C--------------------------------------------------------------------------------------------------
C
      SUBROUTINE GVERROR (RERMAX)
C
C--------------------------------------------------------------------------------------------------
C
C     Routine reports % error in head at the reference point.
C     (called in SOLUT)
C 
      IMPLICIT NONE
      REAL(8) RDUM,RERMAX,RFBASE,RERRH,RFHEAD
      INCLUDE 'GVCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      rdum=RHEAD0-rfbase(crefz)
      IF (rdum.NE.0.0) THEN
       RERRH=ABS(RHEAD0-RFHEAD(CREFZ))/(ABS(rdum))*100
       IF (.NOT.LUCON) WRITE (ILUME,1000) RERRH
       WRITE (*,1000) RERRH
      ELSE
       RERRH=ABS(RHEAD0-RFHEAD(CREFZ))
       IF (.NOT.LUCON) WRITE (ILUME,2000) RERRH
       WRITE (*,2000) RERRH
      ENDIF
      if (rerrh.ge.rconverge_reference) lquit=.false. ! do not abort iterations yet
      RERMAX=MAX(RERMAX,RERRH)
      RETURN
 1000 FORMAT ('    error in head at the reference point:           ',
     & E11.4,' %')
 2000 FORMAT ('    error in head at the reference point: (abs.)    ',
     & E11.4)
      END

