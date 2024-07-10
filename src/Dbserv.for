C     Last change:  HMH  29 Apr 2015    4:37 pm
c --------------------------------------------------------------------------------
c
c     This file contains the following routines or functions
c
c     DBPREP    driver for assigning data to common block for line doublet functions
c     DBFILL    fills common blocks for line doublets on either side of a single node
c     LINDOM    true if inside a inhomogeneity, passes relevant data for that inhomogeneity
c     LFDBIN    true is inside inhomogeneity string
c     RFDBGM    returns exfiltration rate at point due to all inhomogeneities
c     LDBOUTSIDE    used to determine if point is in the inactive zone for a closed slurry wall (horizontal barrier)
c     DBCHECKBOTTOM stops program if bottom elevation of inhomogeneity is at or above aquifer top
c     setldbmatrix_true  called in SOLUT to avoid point shifts (points can be on vertices)
c     setldbmatrix_false called in SOLUT upon exit to allow for points to be shifted away from vertices
C     dbcollocation_prep  routine flags common boundaries for transmissivity inhomogeneity domains and
c                         determines k+, k-, b+, and b- for each line-doublet
c     ldbsharednode    logical function set true when a node is shared by one or more inhomogeneity domains.
c     set_idbcode      routine sorts out the 7 cases for wich node may be shared
c     SET_DBPROPERTIES  Routine assigns properties to k or b at a node.
c     db_cancellakerecharge creates a recharge only domain to cancel recharge inside a lake.
c
c ----------------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBPREP (CZ)
c
c ----------------------------------------------------------------------------------------------------------------
c
C
C     Routine fills common /DOUBL/ with element related data used by
C     potential functions RFDBF, RFDBS, and RFDBG, and by the derivative
C     functions RFDBFD, RFDBSD, RFDBGD.
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD1,INODL,INOD
      LOGICAL LCZ,LFDBIN
      COMPLEX(8) CZ,CZ0
      INCLUDE 'DBCOM.INC'
      INCLUDE 'MAIN.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      COMMON /temp/ lcz  ! debugging, remove
      DATA CZ0 /(1.0E21,1.0E21)/
      IF (NDB.EQ.0) RETURN
      IF (NDBSTR.EQ.0) THEN
      WRITE (ILUER,1000) NDB
      RETURN
      ENDIF
      IF (CZ.EQ.(1.0E21,1.0E21)) THEN
C -----------this call is from <SOLUT>, <DBKNO> or <MAIN> (Load command) to reset CZ0 after a new solution.
C            CZ0 normally contains CZ from previous call.
        CZ0=CZ
        RETURN
      ENDIF
C      IF (lsol.and.CZ.EQ.CZ0) RETURN !  Test for redundant execution of this routine (lsol because CZ may not differ but data may in interactive version)
c
c     !!!!!!!!!!!!!!!! replacing the statement above by the one below is only OK for BATCH use of GFLOW1.EXE 
c
      IF (CZ.EQ.CZ0) RETURN   ! OK for batch version (speeds up solutions)  --------------- 06/22/20000
      do inod=1,ndb  ! initialize
      cddpad(inod)=(0.0d0,0.0d0)
      end do
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN  ! open string
      DO INOD=INOD1,INODL
      CALL DBFILL (INOD1,INOD,0,CZ) ! 0 for last point flags open string, see also DBMAT
      END DO
      ELSE                        ! closed string
      DO INOD=INOD1,INODL
      CALL DBFILL (INOD1,INOD,INODL,CZ)
      END DO
      ENDIF
  20  CONTINUE
      DO 30 ISTR=1,NDBSTR  ! Set logicals for CZ inside a string
      LDBINS(ISTR)=LFDBIN(ISTR)
  30  CONTINUE      
      CZ0=CZ   ! update CZ0
      RETURN
 1000 FORMAT (' ***ERROR in DBPREP: ',I4,' doublet nodes, but no ',
     &        'strings encountered!')
      END
c
c ----------------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBFILL (INOD1,INOD,INODL,CZ)
c
c ----------------------------------------------------------------------------------------------------------------
c
C
C     Routine fills common block to calculate functions at node INOD
C     If CZ at node INOD ln((z-z3)/(z-z1)) is given to CDDPAD
C     to be given to CFDBF(INOD)
C     Routine is called by DBPREP and DBMAT
C     NOTE: INODL=0 flags an open string.
C
      IMPLICIT NONE
      INTEGER(4) INOD1,INOD,INODL,INODM1,INODP1,I
      LOGICAL LCZ,L2
      REAL(8) DUM
      COMPLEX(8) CDDUM,CDD1,CDD2,CDDBZ,CZ,CFBIGZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      COMMON /temp/ lcz  ! debugging, remove
      IF (INOD.EQ.INOD1) THEN
      INODM1=INODL         ! will be zero for open string
      ELSE
      INODM1=INOD-1
      ENDIF
      IF (INOD.EQ.INODL) THEN ! will not happen for open string, INODL=0 passed
      INODP1=INOD1
      ELSE
      INODP1=INOD+1
      ENDIF
C
C Note: this routine is generating coefficient arrays (in common) used to calculate complex potentials
C       and derivative functions. The coefficients are calculated for the line doublet that follows
C       node INOD, exept for the LOG term, which may be calculated for the previous and current line doublet
C       at once for the case that CZ coincides with the starting node of the current line doublet.
C
      CDDBZ=CFBIGZ(CZ,CDBZ(INOD),CDBZ(INODP1))
      CDBZEL(INOD)=CDDBZ
      LDBFAR(INOD)=ABS(CDDBZ).GT.2.0D0
      IF (LDBFAR(INOD)) THEN
C -----------calculate terms for farfield expansion
       CDDUM=1.0D0
       DO I=1,NDBTRM
       CDDUM=CDDUM/CDDBZ
       CDBZIV(I,INOD)=CDDUM
       END DO
      ENDIF
C ------------calculate complex logarithm  NOTE: this logic must remain intact in order to accomodate the various
c                                                calls to LINSIDE etc., which uses vertices!!!
       LDBNOD(INOD)=ABS(cz-cdbz(inod)).LT.DBEPS ! too close to starting node
       L2=ABS(cz-cdbz(inodp1)).LT.DBEPS           ! too close to ending node
       IF (LDBNOD(INOD).OR.L2) THEN
        CDDBLN(INOD)=(0.0D0,0.0D0)
        IF (LDBNOD(INOD).AND.INODM1.NE.0) THEN   ! skip if CZ on first point of open string
          CDDUM=CMPLX((CZ-CDBZ(INODP1))/(CZ-CDBZ(INODM1)))
          CDDPAD(inod)=LOG(CDDUM)
          IF (AIMAG(CDDPAD(inod)).LT.0.0D0)
     &      CDDPAD(inod)=CDDPAD(inod)+CDI*D2PI
        ENDIF
        IF (LDBNOD(INOD).AND.INODM1.EQ.0) THEN   ! set CDDPAD zero when at first point of open string
           CDDPAD(inod)=(0.0D0,0.0D0)
        END IF
       ELSE
        CDD1=CDDBZ-1.0D0
        CDD2=CDDBZ+1.0D0
        CDDUM=CDD1/CDD2
        CDDBLN(INOD)=LOG(CDDUM)
        IF (AIMAG(CDDUM).EQ.0.0D0) THEN
         DUM=CDDUM
         IF (DUM.LT.0.0D0) THEN
          DUM=CDDBLN(INOD)
          CDDBLN(INOD)=CMPLX(DUM,DPI)
         ENDIF
        ENDIF
       ENDIF
C --------------calculate z2-z1 for derivative functions
      CDBZ21(INOD)=CDBZ(INODP1)-CDBZ(INOD)
c
c     test use of ldbnod etc.
c
c      write (ilume,1001) inod,ldbnod(inod),cddpad(inod)
c 1001 format ('dbfill1: inod,ldbnod,cddpad ',i4,1x,l3,1x,2(d14.7))
      RETURN
      END
c
c
c ----------------------------------------------------------------------------------------------------------------
c
      LOGICAL FUNCTION Linside_inhom (CZ,RKI,RP,RBI,RT)
c
c ----------------------------------------------------------------------------------------------------------------
c
C     This function is a reduced verion of the original function LINDOM,
c     to be used only if properties inside a inhomogeneity domain are needed.
c     It is called in GFSERV, where conductivity, base, porosity, and aquifer top
c     are being set.  Created on 9/3/04.
C     Function is true when CZ is inside an inhomogeneity domain.
C
c
!       Added 9/2/99: If inside conductivity is negative, the previous (non-negative)
!       conductivity is given to the inside.
!       Also the previous aquifer base and porosity values are returned.
!       This modification makes recharge only inhomogeneities transparent with respect
!       to transmissivity and porosity. Consequently, recharge only inhomogeneities are
!       allowed to overlap other inhomogeneity domains.
!
!       Added 7/4/2000: Recharge only inhomogeneities are redefined as domains for which both the
!       conductivity and the base elevation are not set (flagged by -9999.0).
!       Warning: if k=-9999.0 and b<>-9999.0 the outside k value is given to rdbk(istr), effectively
!       redefining the inhomogeneity domain as one with both k and b set.
!       This precludes interactive use of GFLOW1.EXE.!! Best to redesign this logic!!
!
!       Modiefied on 1/15/2013 to always set the inside porosity (logic was missing porosity when T was not modified)
!
C     The variable RT is given the global aquifer top (not allowed to vary).
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD
      REAL(8) RKI,RP,RBI,RT
      COMPLEX(8) CZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'GVCOM.INC'   ! change this to a GCPAR call
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
C
      Linside_inhom=.FALSE.
      RKI=RK
      RP=RPOR
      RBI=RBASE
      RT=RBASE+RH   ! RT is always the same in GFLOW version 3.0
      IF (NDB.EQ.0) RETURN
      CALL DBPREP (CZ)
      DO 10 ISTR=1,NDBSTR
c      write (iluer,1001) cz,istr,idomaintype(istr)
 1001 format (' Linside_inhom1: cz,istr,idomaintype '2(E14.7),2(i3))
      IF (IDOMAINTYPE(ISTR).NE.1) GOTO 10 ! skip if a slurry wall (type = 2 or 3)
      IF (LDBINS(ISTR)) THEN
      RP=RDBP(ISTR) ! Always set inside pororosity (modified on 1/15/2013)
c      write (iluer,1002) rdbk(istr),rdbb(istr)
 1002 format (' Linside_inhom2: rdbk(istr),rdbb(istr) ',2(E14.7))
        if (RDBK(istr).NE.-9999.0.AND.RDBB(istr).NE.-9999.0) then ! not a recharge only domain       PROBLEM >> Missing porosity only inhomogeneity
          RKI=RDBK(ISTR)
          RBI=RDBB(ISTR)
        endif   ! keep previous inside properties if recharge inhomogeneity only
         Linside_inhom=.TRUE.  ! are in one or more general or recharge only inhomogeneities
      ENDIF
c      write (iluer,1003) rbi,Linside_inhom
 1003 format (' Linside_inhom3 ',E14.7,l3)
  10  CONTINUE

      RETURN
      END      

c
c ----------------------------------------------------------------------------------------------------------------
c
      LOGICAL FUNCTION LINDOM (ISTR0,RKI,RKO,RP,RBI,RBO,RT)
c
c ----------------------------------------------------------------------------------------------------------------
c
C     This function has been redesigned to allow for common inhomogeneity
c     boundaries. It's use is now limited to solution routines in DBMAT.FOR where
c     when both inside and outside properties are needed and to DBCHECK.FOR.
c     Note: the domain will be seen as inside itself. Hence, the inside and outside properties
c     are returned when the domain is only by itself.
c
C     Function is true when string ISTR0 is inside an inhomogeneity domain.
C     The inside properties for inhomogeniety ISTR0 are given to RKI, RBI, RP, and RT,
c     while the properties of the inhomogeneity in which it occurs are given to
c     RKO, AND RBO.
c     Changes made on 9/3/04.
c
c     *********** REDESIGN THIS LOGIC ************
c
c          KEEP K=B=-9999.0 TEST REMOVE THE REST, SEE DBFAC IN DBMAT.FOR
c

!       Added 9/2/99: If inside conductivity is negative, the previous (non-negative)
!       conductivity are given to the inside and outside k.
!       Also the previous aquifer base and porosity values are returned.
!       This modification makes recharge only inhomogeneities transparent with respect
!       to transmissivity and porosity. Consequently, recharge only inhomogeneities are
!       allowed to overlap other inhomogeneity domains.
!
!       Added 7/4/2000: Recharge only inhomogeneities are redefined as domains for which both the
!       conductivity and the base elevation are not set (flagged by -9999.0).
!       Warning: if k=-9999.0 and b<>-9999.0 the outside k value is given to rdbk(istr), effectively
!       redefining the inhomogeneity domain as one with both k and b set.
!       This precludes interactive use if GFLOW1.EXE.!! Best to redsign this logic!!
!
C     The variable RT is given the global aquifer top (not allowed to vary).
C
      IMPLICIT NONE
      LOGICAL LDOMINSIDEDOM
      INTEGER(4) ISTR,ISTR0,INOD
      REAL(8) RKO,RKI,RP,RBO,RBI,RT
      INCLUDE 'DBCOM.INC'
      INCLUDE 'GVCOM.INC'   ! change this to a GCPAR call
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
C
      LINDOM=.FALSE.
      RKO=RK
      RKI=RK
      RP=RPOR
      RBO=RBASE
      RBI=RBASE
      RT=RBASE+RH   ! RT is always the same in GFLOW version 3.0
      IF (NDB.EQ.0) RETURN
      DO 10 ISTR=1,NDBSTR
      IF (IDOMAINTYPE(ISTR).NE.1) GOTO 10 ! skip if not a transmissivity inhomogeneity
      IF (LDOMINSIDEDOM(ISTR0,ISTR)) THEN
        if (RDBK(istr).NE.-9999.0.AND.RDBB(istr).NE.-9999.0) then ! not a recharge only domain
          RKO=RKI
          RBO=RBI
          RKI=RDBK(ISTR)
          RBI=RDBB(ISTR)
          RP=RDBP(ISTR)
        endif   ! keep previous inside/outside transmissivities and porosity if recharge inhomogeneity only
         LINDOM=.TRUE.  ! are in one or more general or recharge only inhomogeneities
      ENDIF
  10  CONTINUE
      RETURN
      END      
c
c ----------------------------------------------------------------------------------------------------------------
c
      LOGICAL FUNCTION LFDBIN (ISTR)
c
c ----------------------------------------------------------------------------------------------------------------
c
C
C     Function value is .TRUE. if CZ is inside string ISTR.
C     Prior to execution of this routine DBPREP (CZ) should
C     have been called. NOTE: this function is called at the
C     end of the routine DBPREP and should not need further execution.
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD1,INOD,NOD
      REAL(8) RDUM
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
C
      IF (ISTR.LT.1.OR.ISTR.GT.NDBSTR) THEN
      WRITE (ILUER,1000) ISTR,NDBSTR
      LFDBIN=.FALSE.
      RETURN
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN ! cannot be inside open string
      LFDBIN=.FALSE.
      RETURN
      ENDIF
      INOD1=IDBSTA(ISTR)
      NOD=INOD1+NDBSTI(ISTR)-1
      RDUM=0.0
      DO 10 INOD=INOD1,NOD
      RDUM=RDUM+AIMAG(CDDBLN(INOD))
  10  CONTINUE
      LFDBIN=RDUM.GT.1.0
C --------if inside RDUM~2pi, else RDUM~0, test for >1!        
      RETURN
 1000 FORMAT (' ***ERROR in LFDBIN: ISTR=',I3,', but should be ',
     &'between 0 and ',I3)
      END
c
c ----------------------------------------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFDBGM (CZ)
c
c ----------------------------------------------------------------------------------------------------------------
c
C
C    Function returns the total exfiltration rate at CZ due to
C    inhomogeneities.
C
      IMPLICIT NONE
      INTEGER(4) ISTR
      LOGICAL LFDBIN
      COMPLEX(8) CZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      RFDBGM=0.0
      IF (NDBSTR.LE.0) RETURN
      CALL DBPREP (CZ)
      DO 10 ISTR=1,NDBSTR
      IF (LFDBIN(ISTR)) RFDBGM=RFDBGM+RDBGAM(ISTR)
  10  CONTINUE
      RETURN
      END
c
c ----------------------------------------------------------------------------------------------------------------
c
      LOGICAL FUNCTION LDBOUTSIDE (RX,RY)
c
c ----------------------------------------------------------------------------------------------------------------
c
C
C     Function is true when rx,ry is INSIDE a slurry wall, which is declared an
C     "outside" feature; where the model domain is outside the slurry wall.
C     Function is true when rx,ry is OUTSIDE a slurry wall, which is declared an
C     "inside" feature; where the model domain is inside the slurry wall.
C
      IMPLICIT NONE
      INTEGER(4) ISTR
      LOGICAL LFDBIN,LIN,L1,L2
      REAL(8) RX,RY
      COMPLEX(8) CZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      LDBOUTSIDE=.FALSE.
      IF (NDB.EQ.0) RETURN
      DO ISTR=1,NDBSTR
      IF (IDOMAINTYPE(ISTR).EQ.3.AND.ABS(RDBT(ISTR)).EQ.9999.0) THEN ! closed slurry wall with only one model domain
        CZ=CMPLX(RX,RY)
        CALL DBPREP(CZ)
        LIN=LFDBIN(ISTR)
        L1=(LIN.AND.RDBT(ISTR).EQ.+9999.0)      ! inside a slurry wall with outside model domain
        L2=(.NOT.LIN.AND.RDBT(ISTR).EQ.-9999.0) ! outside a slurry wall with inside model domain
        IF (L1.OR.L2) THEN
          LDBOUTSIDE=.TRUE.
          RETURN
        END IF
      END IF
      END DO
      RETURN
      END
c
c -----------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBCHECKBOTTOM ()
c
c -----------------------------------------------------------------------------------------------------
c
C
C     Routine scans for bottom elevations that are equal or higher than the aquifer top
C     at the control points of inhomogeneities or slurry walls.
C     If found, program execution is halted.
C
      IMPLICIT none
      INTEGER(4) INOD,INOD1,INODL,ISTR
      REAL(8) RBOT,RTOP1,RTOP2,RFTOP
      COMPLEX(8) CZ1,CZ2,CZZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      IF (NDB.EQ.0) RETURN
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1  ! last node of closed string, one but last node of open string
      RBOT=RDBB(ISTR)             ! get bottom elevation for the string
      IF (RBOT.EQ.-9999.0) GOTO 20 ! skip if the bottom is to be assigned a local elevation
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
      IF (RDBTFACNOD(INOD1).EQ.1.0E21) GOTO 20 ! k inside = k outside or "recharge only", skip domain
c
c      Hydraulic conductivity inhomogeneity
c
        DO INOD=INOD1,INODL
        CZ1=CDBZ(INOD)  ! first collocation point is line doublet starting point
        IF (INOD.EQ.INODL) THEN
          CZZ=CDBZ(INOD)+CDBZ(INOD1)
        ELSE
          CZZ=CDBZ(INOD)+CDBZ(INOD+1)
        ENDIF
        CZ2=0.5*CZZ    ! second collocation point is line doublet mid point
        RTOP1=RFTOP(CZ1)  ! get aquifer top elevation at the first collocation point
        RTOP2=RFTOP(CZ2)  ! get aquifer top elevation at the second collocation point
C
        IF (RBOT.GE.RTOP1.OR.RBOT.GE.RTOP2) THEN
          AMESS(1)='Inhomogeneity has bottom at or above aquifer top.'
          AMESS(2)='Correct problem and solve again.'
          CALL HALT(2)  ! stop program execution
        END IF
C
        END DO
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
c
c      Open slurry wall  (open horizontal barrier)
c
        DO INOD=INOD1,INODL
        IF (INOD.EQ.INOD1) THEN ! first interval on center of line doublet
          CZ1=CDBZ(INOD)+0.001*(CDBZ(INOD+1)-CDBZ(INOD)) ! first point near start of string
        ELSE
          CZ1=CZ2
        END IF
        IF (INOD.EQ.INODL) THEN ! last point near end of string
          CZ2=CDBZ(INOD)+0.999*(CDBZ(INOD+1)-CDBZ(INOD))
        ELSE                    ! second interval straddles the line doublet end point
          CZ2=CDBZ(INOD)+0.75*(CDBZ(INOD+1)-CDBZ(INOD))
          CZ1=CZ2
          CZ2=CDBZ(INOD+1)+0.25*(CDBZ(INOD+2)-CDBZ(INOD+1))
        END IF
        RTOP1=RFTOP(CZ1)  ! get aquifer top elevation at the first collocation point
        RTOP2=RFTOP(CZ2)  ! get aquifer top elevation at the second collocation point
C
        IF (RBOT.GE.RTOP1.OR.RBOT.GE.RTOP2) THEN
          AMESS(1)='Open barrier has bottom at or above aquifer top.'
          AMESS(2)='Correct problem and solve again.'
          CALL HALT(2)  ! stop program execution
        END IF
C
        END DO
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
c
c      Closed slurry wall     (rdboffsetn=0.5 and rdboffsetc=0.05
c
        CZ1=CDBZ(INODL)+(1.0-RDBOFFSETN)*(CDBZ(INOD1)-CDBZ(INODL)) ! first interval straddles first node
        CZ2=CDBZ(INOD1)+RDBOFFSETN*(CDBZ(INOD1+1)-CDBZ(INOD1))
        RTOP1=RFTOP(CZ1)  ! get aquifer top elevation at the first collocation point
        RTOP2=RFTOP(CZ2)  ! get aquifer top elevation at the second collocation point
C
        IF (RBOT.GE.RTOP1.OR.RBOT.GE.RTOP2) THEN
          AMESS(1)='Closed barrier has bottom at or above aquifer top.'
          AMESS(2)='Correct problem and solve again.'
          CALL HALT(2)  ! stop program execution
        END IF
C
        DO INOD=INOD1,INODL-2
        CZ1=CDBZ(INOD)+RDBOFFSETC*(CDBZ(INOD+1)-CDBZ(INOD))  ! next interval on center of line doublet
        CZ2=CDBZ(INOD)+(1.0-RDBOFFSETC)*(CDBZ(INOD+1)-CDBZ(INOD))
        RTOP1=RFTOP(CZ1)  ! get aquifer top elevation at the first collocation point
        RTOP2=RFTOP(CZ2)  ! get aquifer top elevation at the second collocation point
C
        IF (RBOT.GE.RTOP1.OR.RBOT.GE.RTOP2) THEN
          AMESS(1)='Closed barrier has bottom at or above aquifer top.'
          AMESS(2)='Correct problem and solve again.'
          CALL HALT(2)  ! stop program execution
        END IF
C
        CZ1=CDBZ(INOD)+(1.0-RDBOFFSETN)*(CDBZ(INOD+1)-CDBZ(INOD))
        CZ2=CDBZ(INOD+1)+RDBOFFSETN*(CDBZ(INOD+2)-CDBZ(INOD+1)) ! next interval straddles end node of line doublet
        RTOP1=RFTOP(CZ1)  ! get aquifer top elevation at the first collocation point
        RTOP2=RFTOP(CZ2)  ! get aquifer top elevation at the second collocation point
C
        IF (RBOT.GE.RTOP1.OR.RBOT.GE.RTOP2) THEN
          AMESS(1)='Closed barrier has bottom at or above aquifer top.'
          AMESS(2)='Correct problem and solve again.'
          CALL HALT(2)  ! stop program execution
        END IF
C
        END DO
        CZ1=CDBZ(INODL-1)+RDBOFFSETC*(CDBZ(INODL)-CDBZ(INODL-1))  ! next interval on center of one but last line doublet
        CZ2=CDBZ(INODL-1)+(1.0-RDBOFFSETC)*
     &           (CDBZ(INODL)-CDBZ(INODL-1))
        RTOP1=RFTOP(CZ1)  ! get aquifer top elevation at the first collocation point
        RTOP2=RFTOP(CZ2)  ! get aquifer top elevation at the second collocation point
C
        IF (RBOT.GE.RTOP1.OR.RBOT.GE.RTOP2) THEN
          AMESS(1)='Closed barrier has bottom at or above aquifer top.'
          AMESS(2)='Correct problem and solve again.'
          CALL HALT(2)  ! stop program execution
        END IF
C
        CZ1=CDBZ(INODL-1)+(1.0-RDBOFFSETN)*
     &           (CDBZ(INODL)-CDBZ(INODL-1))  ! next interval straddles last node (INODL)
        CZ2=CDBZ(INODL)+RDBOFFSETN*(CDBZ(INOD1)-CDBZ(INODL))
        RTOP1=RFTOP(CZ1)  ! get aquifer top elevation at the first collocation point
        RTOP2=RFTOP(CZ2)  ! get aquifer top elevation at the second collocation point
C
        IF (RBOT.GE.RTOP1.OR.RBOT.GE.RTOP2) THEN
          AMESS(1)='Closed barrier has bottom at or above aquifer top.'
          AMESS(2)='Correct problem and solve again.'
          CALL HALT(2)  ! stop program execution
        END IF
C
        CZ1=CDBZ(INODL)+RDBOFFSETC*(CDBZ(INOD1)-CDBZ(INODL))  ! last interval on center of last line doublet
        CZ2=CDBZ(INODL)+(1.0-RDBOFFSETC)*(CDBZ(INOD1)-CDBZ(INODL))
        RTOP1=RFTOP(CZ1)  ! get aquifer top elevation at the first collocation point
        RTOP2=RFTOP(CZ2)  ! get aquifer top elevation at the second collocation point
C
        IF (RBOT.GE.RTOP1.OR.RBOT.GE.RTOP2) THEN
          AMESS(1)='Closed barrier has bottom at or above aquifer top.'
          AMESS(2)='Correct problem and solve again.'
          CALL HALT(2)  ! stop program execution
        END IF
C
      END IF
  20  CONTINUE
      RETURN
      end
c
c -----------------------------------------------------------------------------------------------------
c
      SUBROUTINE setldbmatrix_true
c
c -----------------------------------------------------------------------------------------------------
c
C
C    routine is called in SOLUT to set ldmatrix true,
c    which signals that no point shifts will occur to avoid vertices when evaluating PHI.
c    During matrix solution procedures we have collocation points on vertices, which requires
c    removing singularities by considering both line-doublets that share the vertex.
c    After the matrix solution procedure, this option is disabled and line-doublets are evaluated
c    one by one, see subroutine setldbmatrix_false()
C
      IMPLICIT none
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      IF (NDB.EQ.0) RETURN
      ldbmatrix=.true.
      return
      END
c -----------------------------------------------------------------------------------------------------
c
      SUBROUTINE setldbmatrix_false
c
c -----------------------------------------------------------------------------------------------------
c
C
C    routine is called in SOLUT to set ldmatrix false,
c    which signals that from here on point shifts will occur to avoid vertices.
c    This is necessary to allow line-doublets to be evaluated one by one, without looking at
c    nearby line-doublets. These point shifts occur in CDBOM (and also in DBQI).
C
      IMPLICIT none
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      IF (NDB.EQ.0) RETURN
      ldbmatrix=.false.
      return
      END
c ----------------------------------------------------------------------------------------------------------------
c
      SUBROUTINE dbcollocation_prep
c
c ----------------------------------------------------------------------------------------------------------------
c
C
C     Routine sorts out which line-doublets in inhomogeneity domains are in common and
c     how to organize the collocation points. Duplicate line-doublets will be flagged to be skipped.
c     The routine also finds and stores inside and outside
c     properties for each line-doublet in the matrix.
c
c     idbcode=  -1  error, node has not been visited
c
c                0  node is shared by another earlier visited node and on a "common boundary;" not to be included in matrix
c
c                1  node is not on a "common boundary" and to be included in the matrix OR
c                   node is on a "common boundary" between 2 or more domains and to be included in the matrix
c                   (coinciding nodes on other domains are flagged 0, hence not included in the matrix)
c
c                2  node is shared by one or two domains, but not on a "common boundary."
c
c     Note: A node is on a "COMMON BOUNDARY" iff the previous and the next node are shared by that same boundary,
c           in other words if both adjacent line doublets are coinciding.
c
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD1,INODL,INOD,IDOMAIN1,IDOMAIN2,NDOM,
     &           INODP1,INODM1,istr0,idom
      LOGICAL LDBSHAREDNODE,LINDOM
      REAL(8) RKI,RKO,RP,RBI,RBO,RT,RKK0,RHH0,RHED0,RBB0,RPP0,
     &        RHIT,RH1,RFTOP
      COMPLEX(8) CZ,CZ1,CZ2,CZ3,CZZ0
      DIMENSION IDOMAIN1(20),IDOMAIN2(20)
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      IF (NDB.EQ.0) RETURN
      IF (NDBSTR.EQ.0) THEN
      WRITE (ILUER,1000) NDB
      RETURN
      ENDIF
C
      DO INOD=1,NDB   ! debugging; all must be changed, see end of routine
      IDBCODE(INOD)=-1
      END DO
C
C     k=-9999.0 and b=-0.9999 recharge inhomogeneity only, leave -9999.0 in place
C     k=-9999.0 and b=real value: base jump, set k to local hydraulic conductivity
C     k=real value and b=-0.9999: cond. jump, set b to local base elevation
C
        DO ISTR=1,NDBSTR     ! sweep through all line-doublet strings  ?? NOT SURE WHY, WE ARE SETTING PROPERTIES IN A
C                                                                         SPECIAL ROUTINE BELOW (AFTER SETTING IDBCODE).
          IF (IDOMAINTYPE(ISTR).EQ.1) THEN  ! act on all inhomogeneity domains
            IF (RDBK(ISTR).EQ.-9999.0D0.AND.RDBB(ISTR).NE.-9999.0D0)THEN  ! b jump only
                IF (LINDOM (ISTR,RKI,RKO,RP,RBI,RBO,RT)) THEN ! inside an inhomogeneity find local k
                  RDBK(ISTR)=RKI
                ELSE
                  CALL GVPAR (RKK0,RHH0,RHED0,RBB0,RPP0,CZZ0) ! not inside any inhomogeneity set global k
                  RDBK(ISTR)=RKK0
                END IF
      ELSE IF (RDBK(ISTR).NE.-9999.0D0.AND.RDBB(ISTR).EQ.-9999.0D0) THEN  ! k jump only
                IF (LINDOM (ISTR,RKI,RKO,RP,RBI,RBO,RT)) THEN ! inside an inhomogeneity find local b
                  RDBB(ISTR)=RBI
                ELSE
                  CALL GVPAR (RKK0,RHH0,RHED0,RBB0,RPP0,CZZ0) ! not inside any inhomogeneity set global b
                  RDBB(ISTR)=RBB0
                END IF
             END IF
          END IF
        END DO
C
      DO  ISTR=1,NDBSTR  ! Once again sweep through all line-doublet strings
       ldbrechargeonly(istr)=
     &  rdbk(istr).EQ.-9999.0.and.rdbb(istr).EQ.-9999.0 ! no k or b changed in GUI
       IF (IDOMAINTYPE(ISTR).EQ.1) THEN  ! proceed, this is an inhomogeneity domain.
        INOD1=IDBSTA(ISTR)
        INODL=INOD1+NDBSTI(ISTR)-1
c        write (iluer,1001) istr,RDBK(ISTR),RDBB(ISTR)
c 1001 format (' dbcollocation_prep1: istr=',i3,' RDBK,RDBB ',2(D14.7))
        DO INOD=INOD1,INODL  ! sweep through all nodes (vertices) of current domain (string)
        if (ldbrechargeonly(istr)) then
         idbcode(inod)=0  ! keep out of matrix
        else
         NDOM=0
         IF (LDBSHAREDNODE(ISTR,INOD,IDOMAIN1,NDOM)) THEN ! node is shared by at least 1 domain
          DO IDOM=1,NDOM ! check for common boundaries with each of these domains.
           INODM1=INOD-1
           IF (INOD.EQ.INOD1) INODM1=INODL     ! Hmm - these 6 lines can be outside this loop, I think
           INODP1=INOD+1
           IF (INOD.EQ.INODL) INODP1=INOD1
           cz1=cdbz(inodm1)
           cz2=cdbz(inod)
           cz3=cdbz(inodp1)
           istr0=idomain1(idom) ! string numbers of domains on which node occurs are stored in array "idomain."
c
           call set_idbcode (INOD,cz1,cz2,cz3,istr0) ! This is the only place where the idbcode is set!!
c
          END do
         ELSE
          if (idbcode(inod).EQ.-1) IDBCODE(INOD)=1   ! node has not been visited and is not shared by another domain
         ENDIF
        endif
        END DO
       ENDIF
      END DO
c
c    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> TEST OF MISTAKEN IDBCODE ASSIGNMENT!  R E M O V E WHEN DONE!
c
C       idbcode(10)=2
C
C
      DO ISTR=1,NDBSTR ! Once again sweep through all strings to check for missed nodes (= debugging).
       IF (IDOMAINTYPE(ISTR).EQ.1) THEN
       INOD1=IDBSTA(ISTR)
       INODL=INOD1+NDBSTI(ISTR)-1
       DO INOD=INOD1,INODL
c        write (iluer,1003) istr,inod,idbcode(inod),cdbz(inod)
c 1003   format (' dbcollocation_prep3: ISTR=',I3,' IDBCODE(',I3,')=',I3,
c     &                                      ' CZ=',2(d14.7))
        IF (IDBCODE(INOD).EQ.-1) THEN
         WRITE (ILUER,2000) INOD
        END IF
       END DO
       END IF
      END DO
C
c     set properties
c
      DO ISTR=1,NDBSTR
       IF (IDOMAINTYPE(ISTR).EQ.1) THEN
       INOD1=IDBSTA(ISTR)
       INODL=INOD1+NDBSTI(ISTR)-1
       DO INOD=INOD1,INODL
       CALL SET_DBPROPERTIES (ISTR,INOD) ! store inside and outside k and b values
c       write (iluer,1002) inod,idbcode(inod),cdbz(inod),
c     &                                       rdbki(inod),rdbko(inod),
c     &                                       rdbbi(inod),rdbbo(inod)
c 1002  format (' dbcollocation_prep2:',
c     &         ' inod=',i3,' idbcode=',i3,' cz1=',2(d14.7),\,
c     &         ' rdbki,rdbko,rdbbi,rdbbo ',4(d14.7))
c       write (iluer,1004)
c 1004  format (' new line')
      END DO
       END IF
      END DO

c
c     set initial values for the head at collocation points
c
      DO ISTR=1,NDBSTR   ! Once again sweep through all strings, now to set heads (for base jump cases only)
       IF (IDOMAINTYPE(ISTR).EQ.1) THEN
         INOD1=IDBSTA(ISTR)
         INODL=INOD1+NDBSTI(ISTR)-1
         DO INOD=INOD1,INODL
         if (idbcode(inod).gt.0) then ! line-doublet will be in matrix
          if (rdbbi(inod).ne.rdbbo(inod)) then ! we have a base jump
           CZ=CDBZ(INOD1)
           RHIT=RFTOP(CZ)
           IF (RDBH(ISTR).EQ.-9999.0D0) THEN           ! flag for no head estimate given
             CALL GVPAR (RKK0,RHH0,RHED0,RBB0,RPP0,CZZ0)    ! use head at the reference point
             RH1=RHED0
           ELSE
             RH1=RDBH(ISTR)        ! RDBH is an estimate of the average head inside the domain
           END IF
           RH1=MIN(RH1,RHIT)     ! select the initial heads
           rdbhnod(inod)=rh1
           rdbhctr(inod)=rh1
           rdbhend(inod)=rh1
          end if
         end if
       END DO
       END IF
      END DO
c
      RETURN
 1000 FORMAT (' ***ERROR in DBCOLLOCATION_PREP: ',I4,
     &        ' doublet nodes, but no strings encountered!')
 2000 FORMAT (' ***ERROR in DBCOLLOCATION_PREP: IDBCODE(',
     &        I5,' not set.')
      END
c ----------------------------------------------------------------------------------------------------------------
c
      LOGICAL function ldbsharednode(istr0,inod0,idomain,ndom)
c
c ----------------------------------------------------------------------------------------------------------------
c
C
C     Function is true is the node at INOD0 is also occuring on one or more successive domain boundaries
c     If true, the domain (string) numbers are stored in IDOMAIN
C
      IMPLICIT NONE
      INTEGER(4) ISTR,istr0,inod0,INOD1,INODL,INOD,IDOMAIN,ndom
      DIMENSION IDOMAIN(20)
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      idomain(1:20)=0
      ldbsharednode=.false.
      IF (NDB.EQ.0) RETURN
      IF (NDBSTR.EQ.0) THEN
      WRITE (ILUER,1000) NDB
      RETURN
      ENDIF
      if (istr0.eq.ndbstr) return ! we are already in the last string
      ndom=0
      DO 20 ISTR=istr0+1,NDBSTR  ! don't look at previous strings, those have been visited already.
       IF (IDOMAINTYPE(ISTR).EQ.1) THEN  ! proceed, is an inhomogeneity domain.
        INOD1=IDBSTA(ISTR)
        INODL=INOD1+NDBSTI(ISTR)-1
        DO INOD=INOD1,INODL
        if (ABS(cdbz(inod)-cdbz(inod0)).lt.1.0d-4) then  ! found coinciding node
         cdbz(inod)=cdbz(inod0) ! force nodes to be identical
         ndom=ndom+1
         idomain(ndom)=istr
         ldbsharednode=.true.
         GOTO 20        ! go on to next domain, see if there are more conciding nodes
        end if
        END DO
       ENDIF
  20  continue
      RETURN
 1000 FORMAT (' ***ERROR in LDBSHAREDNODE: ',I4,
     &        ' doublet nodes, but no strings encountered!')
      END

c ----------------------------------------------------------------------------------------------------------------
c
      subroutine set_idbcode (INOD0,CZ1,CZ2,CZ3,ISTR)
c
c ----------------------------------------------------------------------------------------------------------------
c
C   Routine sorts out the 7 cases for wich the node CZ2 may be shared with a node on string ISTR.
c
c   Description of the cases:
c
c   case 0: domains only share the point CZ2, while orientation may be the same or opposite
c   case 1: both the previous and current line-doublet are shared and the orientation is the same
c   case 2: only curent line-doublet is shared and the orientation is the same
c   case 3: only previous line-doublet is shared and the orientation is the same
c   case 4: like case 1, but orientation is opposite
c   case 5: like case 2, but orientation is opposite
c   case 6: like case 3, but orientation is opposite
C
c   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> THIS ROUTINE SEEMS INCOMPLETE OR IN ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C                  it works for two abutting domains, but not when there is a third domain inside one of the two
c
C   input parameters:
c   INOD0       node for which idbcode is to be determined
c   CZ1,CZ2,CZ3 are three successive points on the curent string, whereby CZ2 is the node in question.
c   ISTR        string on which to look for a matching node with INOD0
c
c   output parameters: none   (we are assigning values to IDBCODE, which is in a common block in DBCOM.INC
c
c
C
      IMPLICIT NONE
      LOGICAL l11,l13,l31,l33,lnot0,l0not0,l0not2
      INTEGER(4) ISTR,INOD0,INOD1,INODL,INOD,INODP1,INODM1,icase
      COMPLEX(8) CZ1,CZ2,CZ3,C1,C2,C3
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
        INOD1=IDBSTA(ISTR)
        INODL=INOD1+NDBSTI(ISTR)-1
c        write (iluer,1001) inod0,cz1,cz2,cz3,istr
c 1001 format (' set_idbcode1: inod0,cz1,cz2,cz3,istr0 ',/,
c     &        i3,1x,3(2(d14.7),1x),i3)
        DO INOD=INOD1,INODL
        IF (CDBZ(INOD).EQ.CZ2) THEN ! found the matching node, now see if the surrounding ones match also.
          INODM1=INOD-1
          IF (INOD.EQ.INOD1) INODM1=INODL
          INODP1=INOD+1
          IF (INOD.EQ.INODL) INODP1=INOD1
          C1=CDBZ(INODM1)
          C2=CDBZ(INOD)      ! So C2 coincides with CZ2
          C3=CDBZ(INODP1)
c          write (iluer,1002) cz1,cz2,cz3,c1,c2,c3
c 1002 format (' set_idbcode2:',/,
c     &       ' cz1,cz2,cz3 '6(d14.7),/,
c     &       ' c1, c2, c3  '6(d14.7))
          l11=ABS(c1-cz1).lt.1.0d-4   ! previous node on other string coincides with previous node on current string (same orientation)
          l13=ABS(c1-cz3).lt.1.0d-4   ! previous node on other string coincides with next node on current string (opposite orientation)
          l31=ABS(c3-cz1).lt.1.0d-4   ! next node on other string coincides with previous node on current string (opposite orientation)
          l33=ABS(c3-cz3).lt.1.0d-4   ! next node on other string coincides with next node on current string (same orientation)
          lnot0=idbcode(inod).ne.0    ! matching node on other string not already excluded from the matrix
          l0not0=idbcode(inod0).ne.0  ! matching node on current string not already excluded from the matrix
          l0not2=idbcode(inod0).ne.2
c                               ! case 0
            IF (L0NOT0) IDBCODE(INOD0)=2
            IF (LNOT0) IDBCODE(INOD)=2    ! start by assuming case 0 and set codes for both matching nodes to 2 - this may get changed below.
            ICASE=0
          IF (l11.and.l33) then ! case 1   C1.EQ.CZ1.AND.C3.EQ.CZ3
            IF (L0NOT0.and.l0not2) IDBCODE(INOD0)=1
            IF (LNOT0) IDBCODE(INOD)=0
            icase=1
          END IF
          IF (.not.l11.and.l33) then ! case 2    C1.NE.CZ1.AND.C3.EQ.CZ3
            IF (L0NOT0) IDBCODE(INOD0)=2
            IF (LNOT0) IDBCODE(INOD)=0
            icase=2
          END IF
          IF (l11.and..not.l33) THEN ! case 3    C1.EQ.CZ1.AND.C3.NE.CZ3
            IF (L0NOT0) IDBCODE(INOD0)=2
            IF (LNOT0) IDBCODE(INOD)=2
            icase=3
          END IF
          IF (l13.and.l31) THEN ! case 4   C1.EQ.CZ3.AND.C3.EQ.CZ1
            IF (L0NOT0.and.l0not2) IDBCODE(INOD0)=1
            IF (LNOT0) IDBCODE(INOD)=0
            icase=4
          END IF
          IF (l13.and..not.l31) THEN ! case 5   C1.EQ.CZ3.AND.C3.NE.CZ1
            IF (L0NOT0) IDBCODE(INOD0)=2
            IF (LNOT0) IDBCODE(INOD)=2
            icase=5
          END IF
          IF (.not.l13.and.l31) THEN ! case 6  C1.NE.CZ3.AND.C3.EQ.CZ1
            IF (L0NOT0) IDBCODE(INOD0)=2
            IF (LNOT0) IDBCODE(INOD)=0
            icase=6
          END IF
c          write (iluer,1003) icase
c 1003 format (' set_idbcode3: icase=',i3)
        END IF
        END DO
      RETURN
      END subroutine
C
C -------------------------------------------------------------------------------------------------
C
      SUBROUTINE SET_DBPROPERTIES (ISTR,INOD)
C
C -------------------------------------------------------------------------------------------------
C
c     Routine assigns properties to k or b at a node.
c
c     It finds the properties by looking just to the left (inside) and to the right (outside) of the center of each line-doublet
c     This is a robust way to handle multiple domains including domains with comon boundaries.
c     NOTE: INOD is the address of the first node of the line-doublet.
C
      IMPLICIT NONE
      LOGICAL Linside_inhom
      INTEGER(4) ISTR,INOD1,INODL,INOD,INODP1
      REAL RK0,RKI,RB0,RBI,RP,RT,RHED0,RP0,RH0
      COMPLEX(8) CZ1,CZ2,CZ,CFSMALLZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
c
      IF (IDBCODE(INOD).GT.0) THEN ! is in matrix, set properties
       CALL GVPAR (RK0,RH0,RHED0,RB0,RP0,CZ) ! get background properties
c       write (iluer,1001) RK0,RB0
c 1001  format (' set_dbproperties1: rk0,rb0 ',2(d14.7))
       INOD1=IDBSTA(ISTR)
       INODL=INOD1+NDBSTI(ISTR)-1
       INODP1=INOD+1
       IF (INOD.EQ.INODL) INODP1=INOD1
       CZ1=CDBZ(INOD)
       CZ2=CDBZ(INODP1)    ! defines a line-doublet
c
       CZ=CMPLX(0.0,1.0d-6) ! point just on the inside opposite the line doublet center
       CZ=CFSMALLZ (CZ,CZ1,CZ2)
c       write (iluer,1002) cz,cz1,cz2
c 1002  format (' set_dbproperties2: cz,cz1,cz2 ',6(d14.7))
       IF (Linside_inhom (CZ,RKI,RP,RBI,RT)) THEN
        RDBKI(INOD)=RKI   ! set to k in last domain in which the point falls
        RDBBI(INOD)=RBI   ! set to b in last domain in which the point falls
       ELSE
        RDBKI(INOD)=RK0   ! not in a domain, set background k
        RDBBI(INOD)=RB0   ! not in a domain, set background b
       ENDIF
c
       CZ=CMPLX(0.0,-0.01) ! point just on the outside opposite the line doublet center (? much larger distance then before)
       CZ=CFSMALLZ (CZ,CZ1,CZ2)
c       write (iluer,1003) cz,cz1,cz2
c 1003  format (' set_dbproperties3: cz,cz1,cz2 ',6(d14.7))
       IF (Linside_inhom (CZ,RKI,RP,RBI,RT)) THEN
        RDBKO(INOD)=RKI
        RDBBO(INOD)=RBI
       ELSE
        RDBKO(INOD)=RK0
        RDBBO(INOD)=RB0
       ENDIF
      ENDIF
      RETURN
C
      END SUBROUTINE
c
c --------------------------------------------------------------------------
c
      subroutine db_cancellakerecharge (nodes,cztemp,alabtemp)
c
c --------------------------------------------------------------------------
c
c     Routine creates a recharge only inhomogeneity with vertices stored in the
c     array cztemp and labels made out of the line-sink labels stored in alabtemp.
c     The recharge will be set to -1.0 * the current recharge.
c
c
      implicit none
      INTEGER nodes,i,istart,iend,inod
      CHARACTER(16) alabtemp,alabend
      REAL(8) rfpor,rfdbgm,rftop
      COMPLEX(8) cztemp,czex,czplus,cz
      DIMENSION cztemp(*),alabtemp(*)
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
c
      if (nodes.gt.0) then ! OK, we have data
       if (ndbstr+1.le.ndbsmx.and.ndb+nodes.le.ndbzmx) then ! OK, we have space
        ndbstr=ndbstr+1
        ndbsti(ndbstr)=nodes
        idbsta(ndbstr)=ndb+1
        idomaintype(ndbstr)=1
        ldbrechargeonly(ndbstr)=.true.
        rdbk(ndbstr)=-9999.0
        rdbb(ndbstr)=-9999.0
        do i=1,nodes
         ndb=ndb+1
         cdbz(ndb)=cztemp(i)
         iend=LEN(alabtemp(i))
         alabend=alabtemp(i)(3:iend)
         adblab(ndb)="IN"//TRIM(alabend) ! use line-sink labels with LS replaced by IN
         idbcode(ndb)=0
        end do
        call dbdoublpt(ndbstr,inod,czex)
        if (inod.gt.0) then
         write (iluer,3000) czex
        end if
        call dborien() ! ensure that domain is to the left when traversing the boundary.
        istart=idbsta(ndbstr)
        czplus=cdi*(cdbz(istart+1)-cdbz(istart))*0.01
        cz=0.5*(cdbz(istart)+cdbz(istart+1))+czplus ! let's hope this is inside the lake
        rdbp(ndbstr)=rfpor(cz)
        rdbt(ndbstr)=rftop(cz)
        rdbgam(ndbstr)=-rfdbgm(cz) ! assign minus the current recharge to nullify it.
        call dbrech(ndbstr) ! generate the recharge coefficients
        call dbsort() ! sort the domain so they are ordered from outside to inside
       else
        write (iluer,2000) ndbsmx,ndbzmx
       endif
      else
        write (iluer,1000) nodes
      end if
      return
c
 1000 format (' ***ERROR in db_cancellakerecharge: ',
     &'illegal number of vertices, nodes=',i5)
 2000 format (' ***ERROR in db_cancellakerecharge: ',
     &'not enough space left in inhomogeneiety domains.',/,
     &'maximum strings=',i5,' maximum line-doublets=',i5 )
 3000 format (' ***ERROR in dbcancellakerecharge: ',
     &'zero length line-doublet encountered, removed ',
     &'vertex at ',2(d14.7))
      end subroutine

