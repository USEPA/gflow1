C     Last change:  HMH  19 Aug 2004    9:39 pm
C     This file contains the following routines and functions
C  
C     SUBROUTINE W3ERROR   maximum error in the boundary conditions specified at ppwells
C
C
C --------------------------------------------------------------------
C
      SUBROUTINE W3ERROR (RERMAX)
C
C --------------------------------------------------------------------
C
C
C     Reporting maximum error
C
      IMPLICIT NONE
      INTEGER(4) IW,NWQ,NWH,IEN,IST,I,NPT
      REAL(8) RERMAX,RHDI,RW3ERRQ,RW3ERRH,RQ,
     &        RHDSUM,RBAS,RFBASE,RFHEAD,RERR,RHDAV,RFHGHT,RHEDS
      COMPLEX(8) CZ
      INCLUDE 'W3COM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION RHDI(100)
      RW3ERRQ=0.0
      RW3ERRH=0.0
      IF (NW3.EQ.0) RETURN
      NWQ=0.0
      NWH=0.0
      DO 50 IW=1,NW3
      IEN=IPNT(IW+1)-1
      IST=IPNT(IW)
      RQ=0.0
      RHDSUM=0.0
      NPT=0
      cz=CMPLX(rw3cp(1,ist),rw3cp(2,ist))
      rbas=rfbase(cz)
      DO 10 I=IEN,IST,-1      ! calculate heads and discharge
      CZ=CMPLX(RW3CP(1,I),RW3CP(2,I))
      R3DZ=RW3CP(3,I)
      RQ=RQ+RW3S(I)*RW3L(I)
      NPT=NPT+1
      RHDI(NPT)=rfhead(cz)-rbas   ! head measured with respect to aquifer base
      RHDSUM=RHDSUM+RHDI(NPT)
  10  CONTINUE
      IF (LW3Q(IW)) THEN            ! discharge specified well
      NWQ=NWQ+1
      RERR=100.0
      IF (RW3Q(IW).NE.0.0) RERR=ABS((RW3Q(IW)-RQ)/RW3Q(IW))*100.0
      RW3ERRQ=MAX(RW3ERRQ,RERR)
      RHDAV=RHDSUM/NPT          ! average head along the well bore (with respect to aquifer base)
      IF (RHDAV.LE.RFHGHT(CZ)) THEN    !  ppwell in unconfined aquifer zone
        IF (.NOT.LUCON) WRITE (ILUME,4000) AW3LAB(IW)
        WRITE (*,4000) AW3LAB(IW)
        WRITE (ILUER,4000) AW3LAB(IW)
      END IF
      IF (RHDAV.LE.0.0) THEN    !  ppwell pumped dry
        IF (.NOT.LUCON) WRITE (ILUME,3000) AW3LAB(IW)
        WRITE (*,3000) AW3LAB(IW)
        WRITE (ILUER,3000) AW3LAB(IW)
      END IF
      IF (RHDAV.NE.0.0) THEN  ! skip this if well is pumped dry
        NPT=0.0
        DO 20 I=IEN,IST,-1
        NPT=NPT+1
        RERR=100.0
        RERR=ABS((RHDI(NPT)-RHDAV)/RHDAV)*100.0
        RW3ERRQ=MAX(RW3ERRQ,RERR)
  20    CONTINUE
      ENDIF
      ELSE                          ! head specified well
      NWH=NWH+1
      RHEDS=RW3HED(IW)-rbas
      NPT=0
      DO 30 I=IEN,IST,-1
      NPT=NPT+1
      RERR=100.0
      IF (RHEDS.NE.0.0) RERR=ABS((RHEDS-RHDI(NPT))/RHEDS)*100.0
      RW3ERRH=MAX(RW3ERRH,RERR)
  30  CONTINUE
      ENDIF
  50  CONTINUE
      IF (NWQ.NE.0) then
      if (.not.lucon) WRITE (ILUME,1000) NWQ,RW3ERRQ
      write (*,1000) nwq,rw3errq
      if (rw3errq.ge.rconverge_well_3D) lquit=.FALSE. ! do not yet abort iterations
      endif
      IF (NWH.NE.0) then
      if (.not.lucon) WRITE (ILUME,2000) NWH,RW3ERRH
      write (*,2000) nwh,rw3errh
      if (rw3errh.ge.rconverge_well_3D) lquit=.FALSE. ! do not yet abort iterations
      endif
      RERR=MAX(RW3ERRQ,RW3ERRH)
      RERMAX=MAX(RERMAX,RERR)
      RETURN
 1000 FORMAT (' ',I3,' disch. spec. part. penetr. wells:   max. error=',
     &        E11.4,' %')
 2000 FORMAT (' ',I3,' head spec.   part. penetr. wells:   max. error=',
     &        E11.4,' %')
 3000 FORMAT (' WARNING: Part. pen. well ',a16,' is pumped dry!')
 4000 FORMAT (' WARNING: Part. pen. well ',a16,' in unconfined zone!')
      END