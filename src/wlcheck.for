C     Last change:  HMH  19 Aug 2004    9:39 pm
C     This file contains the following routines and functions
C  
C     SUBROUTINE WLERROR   maximum error in the boundary conditions specified at wells
C
C
C ----------------------------------------------------------------------------
C
      SUBROUTINE WLERROR (RERMAX)
C
C ----------------------------------------------------------------------------
C
C     Routine calculates and reports the maximum error in the
C     boundary conditions specified at wells.
C
      IMPLICIT NONE
      INTEGER(4) IW,IHED
      REAL(8) RERMAX,RERRH,RBAS,RS,RF,RFFF,RDUM,RFBASE,RFHEAD
      COMPLEX(8) CZ
      INCLUDE 'wlcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      IF (NWL.EQ.0) RETURN
      RERRH=0.0
      IHED=0
      DO 10 IW=1,NWL
      IF (LWLH(IW)) THEN
            IHED=IHED+1
            CZ=CWLZ(IW)+CMPLX(0.0d0,RWLR(IW)) ! (make sure is same as in WLCZCZ)
            rbas=rfbase(cz)
            RS=RWLH(IW)-rbas
            RF=RFHEAD(CZ)-rbas
            RFFF=MAX(ABS(RF+RS),1.0E-10)      ! avoid division by zero
            RDUM=ABS(RF-RS)/(RFFF/2.0)*100
            RWLERRH(IW)=RDUM    ! save for extract routine
            RERRH=MAX(RDUM,RERRH)
      ELSE                                ! check for wells pumped dry
            CZ=CWLZ(IW)+CMPLX(0.0d0,RWLR(IW)) ! (make sure is same as in WLCZCZ)
            RBAS=RFBASE(CZ)
            RF=RFHEAD(CZ)
            RWLH(IW)=RF
            IF (ABS(RF-RBAS).LT.0.001) THEN
              WRITE (ILUER,2000) AWLAB(IW)
              if (.not.lucon) write (ilume,2000) awlab(iw)
              write (*,2000) awlab(iw)
              RWLH(IW)=-9999.0  ! flag dry wells for GUI
            ENDIF
      ENDIF
  10  CONTINUE
      IF (IHED.GT.0) then
      if (.not.lucon) WRITE (ILUME,1000) IHED,RERRH
      write (*,1000) ihed,rerrh
      if (rerrh.ge.rconverge_well_2D) lquit=.FALSE. ! do not yet abort iterations
      endif
      RERMAX=MAX(RERMAX,RERRH)
      RETURN
 1000 FORMAT (' ',I3,' head specified wells:               max. error=',
     &        E11.4,' %')
 2000 FORMAT (' *** WARNING: well ',A16,' pumped dry!')
      END

