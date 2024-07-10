C     Last change:  HMH   1 Sep 2017    5:26 pm
c --------------------------------------------------------------------------------
c
c     This file contains the following routines or functions
c
c     DBERROR    driver for reporting errors at control points
c     RFDBER     returns the error at control points for a string
c     rfdbslocal returns the strength at point x on the line-doublet
c
c -------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBERROR (RERMAX)
c
c -------------------------------------------------------------------------------------------------------
c
C
C     Routine calculates the error in the boundary condition at all doublet
C     controll points. The maximum relative error is reported.
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD1,INODL,ICOUNT,INOD,INODS,INODE,INODM1,INODP1
      LOGICAL LREPOR,LINDOM,LCONDUCTANCE,lshift1,lshift2,ls,lc,le
      REAL(8) RERMAX,RERROR,RERINH,RT,RBO,RBI,RP,RKO,RKI,RFPOT,
     &        RDUM,RFDBER,RERSLUR,RTOT1,RTOT2,R1,R2,RTOT3,RBZ,
     &        RFDBNFCONDITION,RFNORMALFLOW,R3,rfhead,rhloc,rftop,rfhfp
      COMPLEX(8) CZ,CZA,CZB,cfsmallz
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      IF (NDB.EQ.0) RETURN
      RERROR=0.0
      LREPOR=.FALSE.
      DO 20 ISTR=1,NDBSTR
      INOD1=IDBSTA(ISTR)
      LREPOR=.TRUE.
      INODL=INOD1+NDBSTI(ISTR)-1
      IF (IDOMAINTYPE(ISTR).EQ.1) THEN
C
C     Hydraulic conductivity inhomogeneities
C
      IF (ldbrechargeonly(istr)) GOTO 20 ! "recharge only", skip domain
      RERINH=0.0
      DO INOD=INOD1,INODL            ! now calculate and store the errors in BC at node and center
      if (idbcode(inod).gt.0) then   ! line-doublet is indeed in the matrix
       INODS=INOD
       INODE=INOD+1
       IF (INOD.EQ.INODL) INODE=INOD1
       lshift1=idbcode(inods).eq.2
       lshift2=idbcode(inode).eq.2.or.idbcode(inode).eq.0
       ls=rdbtfacnod(inods).ne.1.0d21
       lc=rdbtfacctr(inods).ne.1.0d21
       le=rdbtfacend(inods).ne.1.0d21
       if (ls) then
        CZ=CDBZ(INOD)        ! first collocation point
        RBZ=-1.0D0
        IF (lshift1) THEN
         rbz=rbz+rdboffsetinhom
         CZ=CMPLX(RBZ,0.0D0)
         CZ=CFSMALLZ(CZ,CDBZ(INODS),CDBZ(INODE))
        END IF
        RDUM=RFDBER (CZ,rbz,ISTR,INODS,INODE) ! error at or near first collocation point
        rdberrstrt(inods)=rdum  ! save for use in inhomextract
        RERINH=MAX(RDUM,RERINH)
       endif
       if (lc) then
        CZ=0.5*(CDBZ(INODS)+CDBZ(INODE)) ! center collocation point
        RDUM=RFDBER (CZ,0.0,ISTR,INODS,INODE)
        rdberrcntr(inods)=rdum  ! save for use in inhomextract
        RERINH=MAX(RDUM,RERINH)
       endif
c
       IF (lshift2.and.le) THEN ! there is also a collocation point near the end of the line-doublet
        RBZ=1.0D0-rdboffsetinhom
        CZ=CMPLX(RBZ,0.0D0)
        CZ=CFSMALLZ(CZ,CDBZ(INODS),CDBZ(INODE))
        RDUM=RFDBER (CZ,rbz,ISTR,INODS,INODE)
        rdberrend(inods)=rdum  ! save for use in inhomextract (not currently used)
        RERINH=MAX(RDUM,RERINH)
       END IF
      endif
      END DO
      if (.not.lucon) WRITE (ILUME,1000) adblab(inod1)(1:10),RERINH
      write (*,1000) adblab(inod1)(1:10),rerinh
      if (rerinh.ge.rconverge_inhomogeneities) lquit=.false. ! do not abort iterations yet
      RERROR=MAX(RERINH,RERROR)
      ENDIF
      IF (IDOMAINTYPE(ISTR).EQ.2) THEN
c
c     Open slurry walls
c
      RERSLUR=0.0
      lconductance=.false.
      rtot1=0.0
      rtot2=0.0
      icount=0
      DO INOD=INOD1,INODL
      IF (INOD.EQ.INOD1) THEN
      CZA=CDBZ(INOD1)+0.001*(CDBZ(INOD+1)-CDBZ(INOD))
      ELSE
      CZA=CDBZ(INOD)+0.25*(CDBZ(INOD+1)-CDBZ(INOD))
      if (rdbtfacnod(inod).gt.0.0) lconductance=.true.
      END IF
      if (rdbtfacctr(inod).gt.0.0) lconductance=.true.
      IF (INOD.EQ.INODL) THEN
      CZB=CDBZ(INODL)+0.999*(CDBZ(INODL+1)-CDBZ(INODL))
      ELSE
      CZB=CDBZ(INOD)+0.75*(CDBZ(INOD+1)-CDBZ(INOD))
      END IF
      RDUM=RFDBNFCONDITION(ISTR,CZA,CZB)
      r1=RFNORMALFLOW(CZA,CZB)
      r2=0.5*(rdum+r1)
      rtot1=rtot1+r1
      rtot2=rtot2+ABS(r2)
      icount=icount+1
      rdum=rdum-r1
      rdberrcntr(inod)=rdum  ! save for use in inhomextract
      RERSLUR=MAX(ABS(RDUM),RERSLUR) ! select maximum error
      IF (INOD.LT.INODL) THEN
      CZA=CZB
      CZB=CDBZ(INOD+1)+0.25*(CDBZ(INOD+2)-CDBZ(INOD+1))
      RDUM=RFDBNFCONDITION(ISTR,CZA,CZB)
      r1=RFNORMALFLOW(CZA,CDBZ(INOD+1))+RFNORMALFLOW(CDBZ(INOD+1),CZB)
      r2=0.5*(rdum+r1)
      rtot1=rtot1+r1
      rtot2=rtot2+ABS(r2)
      icount=icount+1
      rdum=rdum-r1
      rdberrstrt(inod+1)=rdum  ! save for use in inhomextract
      RERSLUR=MAX(ABS(RDUM),RERSLUR) ! select maximum error
      END IF
      END DO
      if (lconductance) then
      rerslur=rerslur/rtot2*icount ! calculate percentage of average flow
      WRITE (ILUME,2010) adblab(inod1)(1:10),RERSLUR
      write (ILUME,2020) rtot1
      WRITE (*,2010) adblab(inod1)(1:10),RERSLUR
      write (*,2020) rtot1
      if (rerslur.ge.rconverge_barriers_resistance) lquit=.false. ! do not abort iterations yet
      else
      WRITE (ILUME,2000) adblab(inod1)(1:10),RERSLUR
      write (ILUME,2020) rtot1
      WRITE (*,2000) adblab(inod1)(1:10),RERSLUR
      write (*,2020) rtot1
      if (rerslur.ge.rconverge_barriers_noflow) lquit=.false. ! do not abort iterations yet
      endif
      RERROR=MAX(RERSLUR,RERROR)
      END IF
      IF (IDOMAINTYPE(ISTR).EQ.3) THEN
c
c     Closed slurry walls
c
      RERSLUR=0.0
      lconductance=.false.
      rtot1=0.0
      rtot2=0.0
      rtot3=0.0
      icount=0
      INODM1=INODL
      DO INOD=INOD1,INODL
      IF (INOD.EQ.INODL) THEN
      INODP1=INOD1
      ELSE
      INODP1=INOD+1
      END IF
      if (rdbtfacnod(inod).gt.0.0) lconductance=.true.
      if (rdbtfacctr(inod).gt.0.0) lconductance=.true.
      CZA=CDBZ(INODM1)+(1.0-RDBOFFSETN)*(CDBZ(INOD)-CDBZ(INODM1))
      CZB=CDBZ(INOD)+RDBOFFSETN*(CDBZ(INODP1)-CDBZ(INOD))
      RDUM=RFDBNFCONDITION(ISTR,CZA,CZB)
      r1=RFNORMALFLOW(CZA,CDBZ(INOD))+RFNORMALFLOW(CDBZ(INOD),CZB)
      r3=0.5*(rdum+r1)
      rtot1=rtot1+r1
      rtot3=rtot3+ABS(r3)
      icount=icount+1
      rdum=rdum-r1
      rdberrstrt(inod)=rdum  ! save for use in inhomextract
      RERSLUR=MAX(ABS(RDUM),RERSLUR) ! select maximum error
      CZA=CDBZ(INOD)+RDBOFFSETC*(CDBZ(INODP1)-CDBZ(INOD))
      CZB=CDBZ(INOD)+(1.0-RDBOFFSETC)*(CDBZ(INODP1)-CDBZ(INOD))
      RDUM=RFDBNFCONDITION(ISTR,CZA,CZB)
      r2=RFNORMALFLOW(CZA,CZB)
      r3=0.5*(rdum+r2)
      rtot2=rtot2+r2
      rtot3=rtot3+ABS(r3)
      icount=icount+1
      rdum=rdum-r2
      rdberrcntr(inod)=rdum  ! save for use in inhomextract
      RERSLUR=MAX(ABS(RDUM),RERSLUR) ! select maximum error
      INODM1=INOD
      END DO
      if (lconductance) then
      rerslur=rerslur/rtot3*icount  ! calculate percentage of average flow
      WRITE (ILUME,3010) adblab(inod1)(1:10),RERSLUR
      write (ilume,3020) rtot1
      WRITE (*,3010) adblab(inod1)(1:10),RERSLUR
      write (*,3020) rtot1
      if (rerslur.ge.rconverge_barriers_resistance) lquit=.false. ! do not abort iterations yet
      else
      WRITE (ILUME,3000) adblab(inod1)(1:10),RERSLUR
      write (ilume,3020) rtot1
      WRITE (*,3000) adblab(inod1)(1:10),RERSLUR
      write (*,3020) rtot1
      if (rerslur.ge.rconverge_barriers_noflow) lquit=.false. ! do not abort iterations yet
      endif
      RERROR=MAX(RERSLUR,RERROR)
      END IF
  20  CONTINUE
      RERMAX=MAX(RERMAX,RERROR)
      RETURN
 1000 FORMAT ('  inhomogeneity domain ',A10,' has as max. error=',
     &        E11.4,' %')
 2000 FORMAT ('    open slurry wall   ',A10,' has as max. error=',
     &        E11.4)
 2010 FORMAT ('    open slurry wall   ',A10,' has as max. error=',
     &        E11.4,' %')
 2020 FORMAT (' Total flow across slurry wall = ',
     &        E11.4)
 3000 FORMAT ('  closed slurry wall   ',A10,' has as max. error=',
     &       E11.4)
 3010 FORMAT ('  closed slurry wall   ',A10,' has as max. error=',
     &       E11.4,' %')
 3020 FORMAT (' Total flow across slurry wall = ',
     &       E11.4)
 4000 FORMAT (' ',I3,' doublets in walls and domains:     max. error=',
     &       E11.4,' %')
      END
C
c
c -------------------------------------------------------------------------------------------------------
c
      REAL FUNCTION RFDBER(CZ,RBIGX,ISTR,INODS,INODE)
c
c -------------------------------------------------------------------------------------------------------
c
C
C     Function returns the percent ERROR in the HEAD at boundary point CZ
C     of domain number ISTR.  Transmissivity inhomogeneities only!
C     RBIGX is the location of CZ on the element (-1 <= RBIGX <= +1).
C     INODS and INODE are the starting and endpoint of the element, repectively.
C     Potential at control points stored in error arrays. Average potential
C     at control points passed as rdbpotaverage
c
c     Note: "error in potential" logic changed to "error in head" logic in June 2002
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INODS,INODE
      LOGICAL LOW
      REAL(8) RPOT,RFPOT,RHEAD,RFHEDP,RKI,RKO,RBI,RBO,RHI,RHO,RP,RFBASE,
     &     RTI,RTO,RDPHI,RFHGHT,RTOP,RFTOP,RBIGS,RBIGX,RERROR,RT,
     &     RPOTFAC,rhtemp,rfhfp,rtt,rhloc,rht,rfdbslocal,
     &     rh1,rh2,rha
      COMPLEX(8) CZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
C
      RFDBER=0.0
      IF (ldbrechargeonly(istr)) THEN ! ERROR =  "recharge only" domain
        write (iluer,2200)
        RETURN
      ENDIF
        R3DZ=RFBASE(CZ)+RFHGHT(CZ) ! default z at saturated aquifer top
        IF (RBIGX.EQ.-1.0D0.OR.RBIGX.EQ.-1.0D0+rdboffsetinhom) THEN
          RPOT=RDBERRSTRT(INODS) ! CZ is at or near start of line-doublet, get potential stored here in DBUPDATE
        ELSE IF (RBIGX.EQ.0.0D0) THEN
          RPOT=RDBERRCNTR(INODS) ! CZ is at center, get potential stored here in DBUPDATE
        ELSE IF (RBIGX.EQ.1.0D0-rdboffsetinhom) THEN
          RPOT=RDBERREND(INODS)  ! CZ is near end of line-doublet, get potential stored here in DBUPDATE
        END IF
        IF (RPOT.LT.0.0) THEN             ! negative potential: skip
          WRITE (ILUER,1100) CZ,RPOT
          rpot=0.0
       END IF
c                 At this point both Ti and To should have proper values.
      rbigs=rfdbslocal (rbigx,dbstr(inods),dbqstr(inods),dbsen(inods))
      rbi=rdbbi(inods)
      rbo=rdbbo(inods)
      rki=rdbki(inods)
      rko=rdbko(inods)
      rt=rftop(cz)
      IF (RBI.LE.RBO) THEN  ! select inside head to test accuracy
        RHLOC=RT-RBI
        RHEAD=RFHFP (RPOT,RKI,RHLOC,RBI)! inside potential generated in DBERROR
        rh1=RHEAD ! inside head relative to inside base elevation
        RHT=RT-RBO  ! now calculate outside head
        RPOT=RPOT-RBIGS   ! make outside potential from inside potential
        RHO=RFHFP (RPOT,RKO,RHT,RBO)
        rh2=RHO+rbo-rbi ! outside head relative to inside base elevation
        if (rho.lt.1.0d-3.and.rh1.lt.rh2) then   !NOTE: may need tolerance instead of zero
          rh2=rh1 ! inside head below outside aquifer base and outside head at aquifer base, force 0 error
        end if
        rpot=rdbpotaverage(istr)
        rha=RFHFP (rpot,rki,rhloc,rbi) ! inside average head relative to inside base elevation
      ELSE                  ! select outside head to test accuracy
        RHLOC=RT-RBO
        RHEAD=RFHFP (RPOT,RKO,RHLOC,RBO)! outside potential generated in DBERROR
        rh1=RHEAD ! outside head relative to outside base elevation
        RHT=RT-RBI  ! now calculate inside head
        RPOT=RPOT+RBIGS   ! make inside potential from outside potential
        RHI=RFHFP (RPOT,RKI,RHT,RBI)
        rh2=RHI+rbi-RBO ! inside head relative to outside base elevation
        if (rhi.lt.1.0d-3.and.rh1.lt.rh2) then    !NOTE: may need tolerance instead of zero
           rh2=rh1 ! outside head below inside aquifer base and inside head at aquifer base, force 0 error
        end if
        rpot=rdbpotaverage(istr)-RBIGS ! make outside average potential
        rha=RFHFP (rpot,rko,rhloc,rbo) ! outside average head relative to outside base elevation
      END IF
      rha=RFHFP (rdbpotaverage(istr),rki,rhloc,rbi) ! average inside head (surrogate)
      if (rha.ne.0.0) then           ! protection added against division by zero on 9/1/17
        rfdber=ABS(rh1-rh2)/rha*100  ! new test procedure (% error in head)
      else
        rfdber=0.0
      end if
c      write (ilume,1002) cz,rh1,rh2,rha,rfdber
c 1002 format ('rfdbr2: cz,rh1,rh2,rha,rfdber (new) ',
c     &          2(d14.7),1x,3(d14.7),1x,d14.7)
      RETURN
 1100 FORMAT (' ***WARNING in RFDBR called from DBERROR:',
     &/, 'Potential at ',2F11.1,' is ',G11.4,)
 2200 format (' ***ERROR in RFDBER call (from DBERROR):',/,
     &/,' Recharge only inhomogeneity, point is skipped!')
      END
c
c-----------------------------------------------------------------------------------------
c
      REAL(8) function rfdbslocal (rbigx,rs1,rmu,rs2)
c
c ----------------------------------------------------------------------------------------
c
c     function returns the strength at point RBIGX on the line-doublet
c
      implicit none
      include 'lusys.inc'
      REAL(8) rbigx,rs1,rmu,rs2
c
      if (rbigx.GE.-1.0d0.and.rbigx.le.1.0d0) then
        rfdbslocal=0.5d0*rs1*(1.0d0-RBIGX)+0.5d0*rs2*(RBIGX+1.0d0)
     &           +rmu*(1.0d0-RBIGX*RBIGX)
      else
        write (iluer,1000) rbigx
      endif
c
      return
 1000 format (' ***ERROR in RFDBSLOCAL: x-value out of range.',/,
     & ' x=',d14.7,' (-1<=x<=+1). Strength is set to 0.')
      end

