C     Last change:  HMH  18 Dec 2006    2:17 pm
c --------------------------------------------------------------------------------
c
c     This file contains the following routines or functions
c
c     DBNEAR    handles pathline trace near a inhomogeneity or slurry wall
c       
c
c --------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBNEAR (CZ,CZNEW,RZ0,RZNEW,RDTWALL,CZCROSS,LREDO)
c
c --------------------------------------------------------------------------------------------------------
c
C
C     Routine adjust the end point of a path line segment to stop just across
C     a slurry wall if it crossed it, not if it moves underneath it.
C     For that case the extra residence time through the wall is given to RDTWALL.
C
C     Input:     CZ       current point of the pathline
C                CZNEW    projected end point of the pathline segment
C                RZ0      elevation of the pathline at CZ
C                RZNEW    elevation of the pathline at the point CZNEW
C                LREDO    some module wants to redo the calculations of CZNEW & RZNEW
C
C     Output:    CZNEW    position of end point of the pathline segment
C                RNEW     elevation of the pathline at the new CZNEW
C                RDTWALL  added residence time due to passage through a slurry wall.
C                CZCROSS  intersection of trace with linedoublet
C
C     NOTE: pathline segment may only intersect 1 line doublet at the time!!! Is not checked.
C
C
      IMPLICIT NONE
      INTEGER(4) INOD,INOD1,INODL,INODP1,ISTR
      LOGICAL LINSECTLINE,l3dend,l3drev,lredo
      REAL(8) RZ0,RZNEW,RDTWALL,RB,RKC,RC,RHEAD,RHEAD1,RHEAD2,RFPOT,
     &     RBIGX,RS,RBASE,RPOR,RP,RW,RFHEAD,RFBASE,RFHGHT,RFPERM,RFHEDP,
     &     RFPOR,RFDBCONDUCTANCE,RSTEP,RH1,RH2,RHGHT,RPOT,RPOT1,RPOT2,
     &     rk,rh3,rzdiv,rsd0,RFTOP
      COMPLEX(8) CZ,CZNEW,CZS,CZE,CZ0,CZCROSS
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
C
      RDTWALL=0.0
      IF (NDB.EQ.0.or.lredo) RETURN
      call getstep (rsd0,rstep,l3dend,l3drev)
      if (l3dend) RETURN ! streamline ended elsewhere
      DO ISTR=1,NDBSTR
      IF (IDOMAINTYPE(ISTR).EQ.2.OR.IDOMAINTYPE(ISTR).EQ.3) THEN ! slurry wall, look for intersection
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      DO INOD=INOD1,INODL
      IF (INOD.EQ.INODL.AND.IDOMAINTYPE(ISTR).EQ.3) THEN
        INODP1=INOD1
      ELSE
        INODP1=INOD+1
      END IF
      CZS=CDBZ(INOD)
      CZE=CDBZ(INODP1)
      IF (LINSECTLINE(CZS,CZE,CZ,CZNEW,CZ0)) THEN ! found intersection
      CZCROSS=CZ0
      CALL DBPREP (CZ0)
      RBIGX=REAL(CDBZEL(INOD))              ! X of CZ0 on line doublet
      RS=-0.5*DBSTR(INOD)*(RBIGX-1.0)      ! S due to starting node
      RS=RS+0.5*DBSTR(INODP1)*(RBIGX+1.0)    ! S due to end node
      RS=RS-DBQSTR(INOD)*(RBIGX*RBIGX-1.0) ! S due to parabol. strength
      CZ0=CZ0+0.001*(CZ0-CZ)  ! ensure CZ0 is just across from line doublet
      RC=RFDBCONDUCTANCE(ISTR,CZ0)
      IF (RC.EQ.0.0) THEN ! crossed no-flow boundary
        WRITE (ILUER,1000) ISTR,CZ0
        RDTWALL=1.0E+21  ! set infinite residence time
        L3DEND=.TRUE.    ! end pathline trace
        RSTEP=ABS(CZ0-CZ)
        CALL SETSTEP (RSTEP,L3DEND) ! adjust step
      END IF
      RBASE=RFBASE(CZ0)
      RPOT=RFPOT(CZ0)
      RHEAD=RFHEDP(RPOT,CZ0)-RBASE ! head at CZ0 (relative to aquifer base)
      RHGHT=RFTOP(CZ0)-RFBASE(CZ0)  ! aquifer thickness (wet  or dry)
      RZNEW=RZ0+ABS(CZ-CZ0)/ABS(CZ-CZNEW)*(RZNEW-RZ0) ! elevation estimate just before CZ0
      RKC=RDBK(ISTR)
      RW=RDBW(ISTR)
      RB=RDBB(ISTR)-RBASE     ! distance between aquifer base and slurry wall bottom
      RB=MAX(RB,0.0)          ! distance is zero when slurry wall ends below aquifer base
      RP=RDBP(ISTR)
      RPOR=RFPOR(CZ0)
      IF (L3DREV) THEN ! pathline jumps up, CZ0 is on upgradient side of wall
        RHEAD1=RHEAD
        RPOT2=RPOT-ABS(RS)
        RHEAD2=RFHEDP(RPOT2,CZ0)-RBASE
        RH1=MIN(RHEAD1,RHGHT)
        RH2=MIN(RHEAD2,RHGHT)
        RZNEW=(RH1/RH2)*(RZNEW-RBASE)+RBASE
      ELSE             ! pathline jumps down, CZ0 is on downgradient side of wall
        RHEAD2=RHEAD
        RPOT1=RPOT+ABS(RS)
        RHEAD1=RFHEDP(RPOT1,CZ0)-RBASE
        RH1=MIN(RHEAD1,RHGHT)
        RH2=MIN(RHEAD2,RHGHT)
        RZNEW=(RH2/RH1)*(RZNEW-RBASE)+RBASE
      END IF
      rk=rfperm(cz0)
      RHEAD=0.5*(RHEAD1+RHEAD2)
      RHEAD=MIN(RHEAD,RHGHT)
      rh3=rh2
      if (l3drev) rh3=rh1
      rzdiv=rb*rk*rh3/((rhead-rb)*rkc+rb*rk)+rbase ! dividing streamline
      IF (RZNEW.GT.rzdiv) THEN ! flow through slurry wall
          RDTWALL=RW*RW*RP*RHEAD*RFPERM(CZ0)/(RKC*ABS(RS))
      END IF
      RETURN
      ENDIF  ! linsectline
      END DO
      END IF ! idomaintype
      END DO
 1000 FORMAT (' ***ERROR in DBNEAR: pathline crossed no-flow boundary,',
     & /,' string ',I3,' trace aborted at ',2(E14.7,1x))
c
      END SUBROUTINE

