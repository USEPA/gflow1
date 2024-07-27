C     Last change:  HMH  17 Jan 2018    1:56 am
c       This file contains the following routines and functions:
c
c       subroutine lskeepsigma   store current strength parameters to compare with those from surface water routines
c       subroutine lsmatsize     updates the matrix array dimensions
C	SUBROUTINE LSCZC         generates control point vector for line-sinks
C	SUBROUTINE LSGENMAT      generates matrix coefficients for line-sinks (no corrections, no delta phi)
C       SUBROUTINE LSMATCORECT   correct for resistance and subtract equations for galleries
C       SUBROUTINE ls_remove_equations  "remove" equations by setting all coefficients zero, except diagonal element
C	SUBROUTINE LSKNO         generates known vector for line-sinks
C	SUBROUTINE LSSUB         substitutes solution vector into strength parameters for line-sinks
C	SUBROUTINE LSUPDATE      stores the calculated heads at the centers of line-sinks
C       REAL FUNCTION rflstable   returns the outlet flow rate using a linear interpolation between stage levels and outlet flow rates
c                                stored in the tables RLSTAGE(i,itable) and rlstable(i,itable)
C       SUBROUTINE LSREADTABLES  reads the lake stage versus outflow rate tables for
c                                      stream sections (strings) that are lake outlet streams
c       LOGICAL function lflslakeiterations: controls iteration process for lakes
c       subroutine set_lakedone() initializes the logical LAKEDONE (tracom.inc)
c
c-------------------------------------------------------------------------------------------------------
c
      subroutine lskeepsigma (drb,j,m)
c
c-------------------------------------------------------------------------------------------------------
c
c     Store current line-sink strengths in DRB and subtract those found at the end of LSCZC.
c     This will create an array DRB with all line-sink strength changes that occur in LSSTREAM and LSCZC.
c     A set of "update" routine calls will adjust all known vector arrays to reflect these changes.
c     Set j=0 before calling the routine, at least when LSCZC is first in the list of control point routines.
c     This routine must be called at the end of the update calls just after the solution substitution.
c
      implicit none
      INTEGER j,m,i,iad
      REAL(8) drb
      DIMENSION drb(m)
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
c
      if (nlsh.eq.0) return
      do i=1,nlsh
      iad=klspth(i)
      j=j+1
      drb(j)=rlsig(iad) ! new sigmas will be subtracted in LSCZC
      end do
c
      end subroutine
c
c --------------------------------------------------------------------------------------------------------
c
      subroutine lsmatsize (M,N)
c
c --------------------------------------------------------------------------------------------------------
c
c     Routine updates the array dimensions M and N for the number of equations and strength parameters
c
      implicit none
      INTEGER i,j,M,N
      INCLUDE 'lscom.inc'
      if (nlsh.gt.0) then
       do i=1,nlsh ! loop over all head specified line-sinks (includes drains and galleries)
         M=M+1
         N=N+1
       end do
      end if
      return
c
      end subroutine
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSCZC(CZI,M,N,DRB,DRFAC,CALPH,ITYPE,LSURFWAT,
     &             lDirectFromDisk,lskip,ltimer)
C
C---------------------------------------------------------------------------------------------------------
C
c     Routine generates the solution arrays: czi, itype, and drfac.
c     Removes or reintroduces line-sinks depending on flow conditions.
c     Strength parameters for line-sinks that are removed may be changed.
c     The solution vector DRB() is given the strength changes, for use in lsupdate_head
c     to be called after the collocation point vector is constructed
c
c     Note: use of stored matrix and solution. Make sure LSMAT is first in matrix, so that
c           the ith head specified line-sink corresponds to the ith matrix equation
c
c
      IMPLICIT NONE
      INTEGER(4) M,N,ITYPE,IAD,I,ICPT,nsolOut
      LOGICAL LNET,LSURFWAT,lDirectFromDisk,
     &        lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut,lskip,ltimer
      REAL(8) DRFAC,RDIS,RC,RHEAD,RFHEDP,RSPERC,RS,RFPERM,RTOL,DRB,RFPOT
      COMPLEX(8) CZI,CALPH,CZ0,CZOLD,CZ0FIRST
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION CZI(*),DRFAC(4,*),CALPH(*),ITYPE(*),DRB(*),lskip(*)
      IF (NLSH.EQ.0) RETURN
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      lsfirst=nsolOut.eq.0
      LSURFWAT=.FALSE.
      IF (LSFIRST) CALL LSREADTABLES ! must read and store the stage tables for lakes
      DO 10 I=1,NLSH
      IAD=KLSPTH(I)
      M=M+1 ! add one equation per head specified line-sink
      N=N+1 ! add one strength parameter (sink density) for each head specified line-sink
      RDIS=RLSDEP(IAD)
      RC=RLSRES(IAD)
      IF (RC.GT.0.0D0) THEN
        IF (RLSWID(IAD).GT.0.0D0) LSURFWAT=.TRUE.
      END IF
      LNET=KLSDN(IAD).NE.IAD              ! line sink part of network
      CZ0=0.5*(CLSZS(IAD)+CLSZE(IAD))
      if (lsfirst) then
        rlspot(iad)=rfpot(CZ0) ! necessary, since not yet set at this point
      endif
      if (lsinlet(iad)) lsgiv(iad)=.true.
      IF (LSDRAIN(IAD).and..not.lsfirst) THEN     ! drain feature
        RHEAD=RFHEDP(RLSPOT(IAD),CZ0)
        IF (RHEAD.GT.RLSH(IAD)) THEN   ! head above drain elevation, add to solution procedure
          LSGIV(IAD)=.FALSE.
        ELSE                           ! head below drain elevation, remove from solution procedure
          LSGIV(IAD)=.TRUE.
          rlsig(iad)=0.0d+0
        END IF
      ELSE                       ! no drain feature, head specified or part of stream network?
      IF ((.NOT.LSBASE.OR..NOT.LNET).AND..NOT.LSINLET(IAD).AND.
     &     LSGIV(IAD)) THEN ! no reset in lsbase
c        write (iluout,1001) lsbase,lnet,lsgiv(iad),iad
c 1001   format (' lsczc1: lsbase,lnet,lsgiv,iad ',3l4,i5)
c        write (iluout,1002) rc,alslab(iad)
        IF (RC.GT.0.0) THEN             ! check for head above res. layer base
c          write (iluout,1002) rc,alslab(iad)        
c 1002     format (' lsczc2: rc,alslab ',g14.7,a16)          
          RHEAD=RFHEDP(RLSPOT(IAD),CZ0)
c          write (iluout,1003) cz0,rhead,rlsh(iad),rdis
          IF (RHEAD.GT.RLSH(IAD)-RDIS) THEN                 ! head is above res. layer base
c            write (iluout,1003) cz0,rhead,rlsh(iad),rdis          
c 1003       format (' lsczc3: cz0,rhead,rlsh,rdis ',2g11.4,3g14.7)
            RLSIG(IAD)=RLSWID(IAD)*(RHEAD-RLSH(IAD))/RC  ! set current strength
            LSGIV(IAD)=.FALSE.  ! bring back into the solution procedure
c            write (iluout,1004) rlswid(iad),rlsig(iad),lsgiv(iad)
c 1004       format (' lsczc4: rlswid,rlsig,lsgiv ',2g14.7,l5)            
          ELSE                           ! head is NOT above resistance layer base
            rsperc=-RLSWID(IAD)*RDIS/RC
            IF (.NOT.LSBASE) then
              RLSIG(IAD)=rsperc ! head below res. layer
            endif
          ENDIF
        ELSE
c          write (iluout,1005) rhead,rlsh(iad)
          IF (RHEAD.GT.RLSH(IAD)) THEN          ! head is above stream level
c            write (ilume,1005) rhead,rlsh(iad)
c 1005       format (' lsczc5: rhead,rlsh ',2g14.7)
            RLSIG(IAD)=0.0
            LSGIV(IAD)=.FALSE.    ! bring back into the solution procedure
          ENDIF
        ENDIF  
      ENDIF
      IF (.NOT.LSFIRST.AND..NOT.LSGIV(IAD).AND..NOT.LSGALLERY(IAD)) THEN
C                              If not a given strength (set in GENSTREAMFLOW),
C                              check for unsaturated zone underneath stream.
C                              Except if at first iteration.
        IF (RDIS*RC.NE.0.0) THEN
          RS=-RDIS/RC*RLSWID(IAD)      ! resistance line sinks
        ELSE
          RS=-RFPERM(CZ0)*RLSWID(IAD)  ! no resistance
          IF (RLSWID(IAD).EQ.0.0) RS=-1.0E20 ! is flag, do not limit infiltration
        ENDIF
        IF (RLSIG(IAD).LE.RS) THEN
          RLSIG(IAD)=rs
          LSGIV(IAD)=.TRUE.
        ENDIF
      ENDIF
      END IF
c      write (iluout,1006) iad,alslab(iad),rlsig(iad),lsgiv(iad)
c 1006 format (' lsczc6: iad,alslab,rlsig,lsgiv ',/,i5,a16,g14.7,l5)      
      DRFAC(1,M)=1.0D00 ! not to be used
      if (LSGIV(IAD)) then
        lskip(M)=.true.
      else
        lskip(M)=.false.
      endif
      IF (LSGALLERY(IAD)) THEN   ! gallery line-sink
        IF (RLSQ(IAD).EQ.-9999.0D0) THEN  ! successive line-sinks in gallery
           ITYPE(M)=-1
           CZI(M)=CZ0
           CALPH(M)=CZOLD
           CZOLD=CZ0
           DRFAC(2,M)=REAL(CZ0FIRST)  ! Mark each successive line-sink with the control point of the first line-sink
           DRFAC(3,M)=AIMAG(CZ0FIRST)
        ELSE ! first line-sink in gallery
           ITYPE(M)=6
           CZI(M)=CZ0
           CALPH(M)=(0.0D0,0.0D0)
           CZOLD=CZ0
           CZ0FIRST=CZ0
           DRFAC(2,M)=REAL(CZ0FIRST)
           DRFAC(3,M)=AIMAG(CZ0FIRST)
        END IF
      ELSE  ! head specified line-sink or drain or stream
      ITYPE(M)=1
      CZI(M)=CZ0
      CALPH(M)=(0.0D0,0.0D0)
      END IF
      if (lsfirst) then
C --------------------------------  Check for nearby control points
        RTOL=0.01*ABS(CLSZS(IAD)-CLSZE(IAD))  ! 1% of current line sink length
        DO 5 ICPT=1,M-1
        IF (ABS(CZI(M)-CZI(ICPT)).LT.RTOL)
     &   WRITE (ILUER,1000) IAD,ALSLAB(IAD),ICPT,CZI(ICPT)
   5    CONTINUE
      endif
c
c     Here we are creating an array of differences between the original line-sink strengths that followed
c     from a matrix solution process and the modified values in first GENSTREAMFLOW and subsequently LSCZC.
c
      drb(M)=rlsig(iad)-drb(M) ! store strength update in solution vector for use in "lsupdate_head"
  10  CONTINUE
      RETURN
 1000 FORMAT (' ***WARNING: line sink ',I3,' with label ',A16,/
     &' may be too close to control point # ',I3,' =',2G11.4)
      END
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSGENMAT (DRA,CZI,M,N,J,DRFAC,CALPH,ITYPE,
     &                  lDirectFromDisk,ltimer)
C
C---------------------------------------------------------------------------------------------------------
c
c     Generate initial matrix coefficients (without corrections)
c
c     ra(i,j) i=row number (equation number) j=element number (column number)
c
c     Note: use of stored matrix and solution. Make sure LSGENMAT is first in matrix, so that
c           the ith head specified line-sink corresponds to the ith matrix equation
c
c     IMPORTANT
c     The idea is to let every line-sink generate a potential coefficient at its own collocation point
C     so the matrix equations can be used in "lsupdate_head" for fast recalculation of all potentials
C     at line-sink collocation points. The actual matrix coefficients necessary for the matrix solution
c     will be created by modifying these initial matrix coefficients in "lsmatcorrect" (after "lsupdate_head"
c     has been called). 
c
      IMPLICIT NONE
      INTEGER(4) M,N,J,ITYPE,II,IAD,I,IEQS,IEQ
      LOGICAL LNEG,lDirectFromDisk,ltimer
      REAL(8) DRA,DRFAC,RFNFLSCO
      COMPLEX(8) CZI,CALPH,CZ1,CZ2,CZ0,CZ,CZA,COMLS,CDLS
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRA(M,1),CZI(1),CALPH(1),ITYPE(1),DRFAC(4,1)
      IF (NLSH.EQ.0) RETURN
      DO 20 II=1,NLSH
      IAD=KLSPTH(II)
      CZ1=CLSZS(IAD)
      CZ2=CLSZE(IAD)
      CZ0=0.5D0*(CZ1+CZ2)
      J=J+1
!      WRITE (ILUME,1000) N,J
      DO 10 I=1,M
      CZ=CZI(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across the line between CZ & CZA, positive from left to right when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request for zero matrix element
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      GOTO (1,2,3,4,5,1),IEQ  ! note itype 6 is interpreted here as itype=1, to get desired potential
c                               coefficient at the line-sink collocation point.
  1   DRA(I,J)=DRA(I,J)+REAL(COMLS(CZ,CZ1,CZ2)+clsconst(iad)) ! provide potential at CZ
c                       note: handle difference in PHI in lsmatcorrect
c      IF (LNEG) THEN
c      DRA(I,J)=DRA(I,J)-REAL(COMLS(CZA,CZ1,CZ2)+clsconst(iad)) ! subtract potential at CZA
c      ENDIF
      GOTO 9
  2   IF (LNEG) THEN
      DRA(I,J)=DRA(I,J)+RFNFLSCO(IAD,CZ,CZA) ! provide flow across CZ & CZA (sink density 1)
      ELSE
      DRA(I,J)=DRA(I,J)+AIMAG(COMLS(CZ,CZ1,CZ2)+clsconst(iad))  ! provide PSI at CZ
      END IF
      GOTO 9
  3   DRA(I,J)=DRA(I,J)+AIMAG(CDLS(CZ,CZ1,CZ2)*CONJG(CZA)) ! provide Q normal to unit vector CZA
      GOTO 9
  4   DRA(I,J)=DRA(I,J)+REAL(CDLS(CZ,CZ1,CZ2)*CONJG(CZA))  ! provide Q parallel to unit vector CZA
      GOTO 9
  5   DRA(I,J)=DRA(I,J)+ABS(CZ2-CZ1)    ! provide total discharge for continuity equation
C     if itype=6 first store potential coefficient (to be used in lsupdate_head) proper matrix elements
C     created in lsmatcorrect
  9   continue
  10  CONTINUE
  20  CONTINUE
      RETURN
 1000 FORMAT ('+Generating',I4,' equations, doing equation #: ',I4)
      END
c
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSMATCORRECT(DRA,CZI,DRB,DRSCR,M,N,J,DRFAC,CALPH,ITYPE,
     &                  lDirectFromDisk,ltimer)
C
C---------------------------------------------------------------------------------------------------------
c
c     Correct matrix coefficients when resistance linesinks
c
c     ra(i,j) i=row number (equation number) j=element number (column number)
c
c     Note: use of stored matrix and solution. Make sure LSMAT is first in matrix, so that
c           the ith head specified line-sink corresponds to the ith matrix equation
c
      IMPLICIT NONE
      INTEGER(4) M,N,J,ITYPE,II,IAD,I,IEQS,IEQ,ikeep
      LOGICAL LNEG,lresistance,lDirectFromDisk,ltimer,
     &        LFINTERFACE,LINTERFACEFLOW
      REAL(8) DRA,DRFAC,RK,RFPERM,RTHICK,RFBASE,RFTOP,RC,RW,RH,rh0,
     &        RHEAD,RFHEDP,RFNFLSCO,RALPHA,RBETA,RHSS,RGSS,RGFF,RF1,RF2,
     &        RFPOTH_confined,RFPOTH_unconfined,rbase,rtop,DRB,DRSCR
      COMPLEX(8) CZI,CALPH,CZ1,CZ2,CZ0,CZ,CZA,COMLS,CDLS,CZ0FIRST,
     &           cz0start
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRA(M,*),CZI(*),CALPH(*),ITYPE(*),DRFAC(4,*),DRB(*),
     &          DRSCR(*)
      LINTERFACEFLOW=LFINTERFACE()
      IF (NLSH.EQ.0) RETURN
C
      cz0start=CMPLX(1.0d21,1.0d21)
      DO 20 II=1,NLSH
      IAD=KLSPTH(II)
      CZ1=CLSZS(IAD)
      CZ2=CLSZE(IAD)
      CZ0=0.5D0*(CZ1+CZ2)
      RC=RLSRES(IAD)
      RW=RLSWID(IAD)
      if (lsgallery(iad).and.rlsq(iad).NE.-9999.0d0) then ! first line-sink in gallery
       cz0start=cz0
      end if
      lresistance=rc*rw.ne.0.0D0
      J=J+1
C      if (lsgiv(iad)) GOTO 20 ! do not correct if equation is to be removed
C
C      BLOCKED OFF THIS STATEMENT TO ENSURE CORRECTION FOR DRAINS BE GIVEN TO DECOMPOSED MATRIX.
C
!      WRITE (ILUME,1000) N,J
      DO 10 I=1,M
      CZ=CZI(I)
      CZA=CALPH(I)
      IEQS=ITYPE(I)
C
C ITYPE=+1 potential specified at CZ
C       -1 difference in potential specified: PHI(CZ)-PHI(CZA)
C       +2 stream function specified at CZ
C       -2 flow across the line between CZ & CZA, positive from left to right when at CZ
C       +3 discharge component normal to the unit vector CZA (rotated to the left)
C       +4 discharge component parallel to the unit vector CZA
C        5 continuity equation: provide total discharge
C        6 request for zero matrix element
C
      LNEG=IEQS.LT.0
      IEQ=IABS(IEQS)
      if (ieq.eq.1) then
      IF (lresistance.AND.CZ.EQ.CZ0) THEN       ! provide correction term when at own collocation point
C
C  At line sink there is resistance between stream and aquifer:
C  phi(cz)-phi(specified)/c=qz (qz=specific discharge into stream)
C  qz*w=s (s=strength of line sink and w=width of stream)
C
        RK=RFPERM(CZ0)
        RHEAD=RFHEDP(RLSPOT(IAD),CZ0) ! actual head at line-sink center
        rbase=rfbase(cz0)
        rtop=rftop(cz0)
        rthick=rtop-rbase
        ralpha=1.0d+0
        rbeta=0.0d+0
        IF (LINTERFACEFLOW) THEN  ! interface flow
          CALL INTERFACEDATA (RHSS,RGSS,RGFF,RF1,RF2)
          RHSS=RHSS-RBASE   ! sea level with respect to aquifer base
c           write (iluer,1002) rhss,rgss,rgff,rf1,rf2,rhead
c 1002  FORMAT(' LSMAT2: rhss,rgss,rgff,rf1,rf2,rhead ',/,6(D14.7))
          if (rhead-rbase.lt.rf1*rhss) then ! interface present
            IF (RHEAD.GE.rtop) THEN ! confined
              RALPHA=RF2/RF1
              RBETA=-RF2*RHSS+RTHICK
            ELSE                         ! unconfined
              RALPHA=RF2
              RBETA=-RF2*RHSS
            END IF
          endif
        END IF
c            write (iluer,1005) lsfirst,lsgallery(iad)
c 1005 format ('lsfirst,lsgallery ',2L5)
            IF (LSFIRST) THEN
              RH=RLSH(IAD)-rbase  ! specified head with respect to the aquifer base at line-sink center
              rh=MIN(rh,rthick)
              IF (LSGALLERY(IAD)) RH=RLSHMIN(IAD)-rbase
            ELSE
                RH=(RHEAD+RLSH(IAD))*0.5D0-rbase   ! average head with respect to aquifer base
                IF (RHEAD.GE.rtop.AND..NOT.LINTERFACEFLOW) RH=RTHICK
              if (LSGALLERY(IAD)) then
                rh0=rlshmin(iad)-rbase ! added to try and improve
                RH=MAX(RH,rh0) ! NOTE: rh0 used to be 0.0D+0
              endif
            ENDIF
            RH=RALPHA*RH+RBETA
            DRA(I,J)=DRA(I,J)-RK*RC*RH/RW
c            WRITE (ILUER,1003)I,J,CZ,CZA,CZ0,RK,RC,RH,RW,RLSH(IAD),
c     &                        RHEAD,RLSHMIN(IAD),RALPHA,RBETA
c 1003 FORMAT (' LSMAT3: I,J,CZ,CZA,CZ0 ',2(I3),6(D14.7),/,
c     &        ' RK,RC,RH,RW,RLSH,RHEAD ',6(D14.7),/,' RLSHMIN ',D14.7,/,
c     &        'RALPHA,RBETA ',2(D14.7))
C             write (iluer,1003) i,j,rk,rc,rh,rw,iad,dra(i,j)
C 1003 format (' lsmat3: i,j,rk,rc,rh,rw,dra(i,j) ',
C     &          2(i4),4(d14.7),i4,d14.7)
      ENDIF
      IF (LNEG) THEN
      DRA(I,J)=DRA(I,J)-REAL(COMLS(CZA,CZ1,CZ2)+clsconst(iad)) ! subtract potential at CZA
      IF (lresistance.AND.CZA.EQ.CZ0) THEN       ! provide correction term when at own collocation point
C
C  At line sink there is resistance between stream and aquifer:
C  phi(cz)-phi(specified)/c=qz (qz=specific discharge into stream)
C  qz*w=s (s=strength of line sink and w=width of stream)
C
        RK=RFPERM(CZA)
        RHEAD=RFHEDP(RLSPOT(IAD),CZA)
        rbase=rfbase(czA)
        rtop=rftop(czA)
        rthick=rtop-rbase
        ralpha=1.0d+0
        rbeta=0.0d+0
        IF (LINTERFACEFLOW) THEN  ! interface flow
          CALL INTERFACEDATA (RHSS,RGSS,RGFF,RF1,RF2)
          RHSS=RHSS-RBASE   ! sea level with respect to aquifer base
          if (rhead-rbase.lt.rf1*rhss) then ! interface present
            IF (RHEAD.GE.rtop) THEN ! confined
              RALPHA=RF2/RF1
              RBETA=-RF2*RHSS+RTHICK
            ELSE                         ! unconfined
              RALPHA=RF2
              RBETA=-RF2*RHSS
            END IF
          endif
        END IF
            IF (LSFIRST) THEN
              RH=RLSH(IAD)-rbase
              RH=MIN(RH,RTHICK)
              IF (LSGALLERY(IAD)) RH=RLSHMIN(IAD)-rbase
            ELSE
               RH=(RHEAD+RLSH(IAD))*0.5D0-rbase
               IF (RHEAD.GE.rtop.AND..NOT.LINTERFACEFLOW)RH=RTHICK
              if (LSGALLERY(IAD)) then
                rh0=rlshmin(iad)-rbase
                RH=MAX(RH,rh0)    ! NOTE: rh0 used to be 0.0D+0
              endif
            ENDIF
            RH=RALPHA*RH+RBETA
            DRA(I,J)=DRA(I,J)+RK*RC*RH/RW
C            WRITE (ILUER,1004)I,J,CZ,CZA,CZ0,RK,RC,RH,RW,RLSH(IAD),RHEAD
C     &                       ,RALPHA,RBETA
C 1004 FORMAT (' LSMAT4: I,J,CZ,CZA,CZ0 ',2(I3),6(D14.7),/,
C     &        ' RK,RC,RH,RW,RLSH,RHEAD ',6(D14.7),
C     &        ' RALPHA,RBETA ',2(d14.7))
      ENDIF
      ENDIF
      endif
      if (ieq.eq.6) THEN ! make zero coefficients
      IF (LSGALLERY(IAD)) THEN
        CZ0FIRST=CMPLX(DRFAC(2,I),DRFAC(3,I))
c        write (ilume,1001) cz0start,cz0first
c 1001 format (' lsmat1: cz0start,cz0first ',2(e11.4),2x,2(e11.4))
        IF (ABS(CZ0FIRST-CZ0START).LT.0.00001) THEN
          DRA(I,J)=ABS(CZ2-CZ1)  ! continuity equation for gallery or lake
        ELSE
          DRA(I,J)=0.0D0
        ENDIF
      ELSE
          DRA(I,J)=0.0D0
      END IF
      end if
  10  CONTINUE
  20  CONTINUE
c
      RETURN
 1000 FORMAT ('+Correcting',I4,' equations, doing equation #: ',I4)
      END
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE ls_remove_equations (DRA,M,N,ltimer)
C
C---------------------------------------------------------------------------------------------------------
c
c     remove matrix equations for given line-sinks
c
c     ra(i,j) i=row number (equation number) j=element number (column number)
c
c     Note: use of stored matrix and solution. Make sure LSMAT is first in matrix, so that
c           the ith head specified line-sink corresponds to the ith matrix equation
c
      IMPLICIT NONE
      INTEGER(4) M,N,II,IAD,I,iticks,iticks1,iticks2,iremove
      LOGICAL ltimer
      REAL(8) DRA
      DIMENSION DRA(M,N)
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
c
      IF (NLSH.EQ.0) RETURN
c     if (ltimer) call timer  (iticks1)
      i=0
      iremove=0
      DO II=1,NLSH ! check all head specified line-sinks (with or without resistance, galleries, drains, lakes)
      i=i+1
      IAD=KLSPTH(II)
      if (lsgiv(iad)) THEN ! remove equation
        dra(i,1:n)=0.0d0
        dra(i,i)=1.0d+01
        iremove=iremove+1
      end if
      ENDDO
      write (ilume,1000) iremove
      if (ltimer) then
        !call timer (iticks2)
        iticks=iticks2-iticks1
      if (ltimer) write (ilume,1001) iticks
 1001 format(' ls_remove_equations execution time=                ',i10,
     &       ' E-2 seconds.')
      end if
      RETURN
 1000 format (' Number of equations removed = ',i10)
      END
C
C ---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSKNO (DRB,J,CZI,ITYPE,lDirectFromDisk)
c
c ---------------------------------------------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER(4) J,IAD,I,J0,ITYPE,iadfirst
      LOGICAL LRESISTANCE,lDirectFromDisk,
     &        LFINTERFACE,LINTERFACEFLOW
      REAL(8) DRB,RSIG,RC,RW,RHEAD,RFHEDP,RH,RFBASE,rbase,RK,rh0,rhed,
     &      RFPERM,RTHICK,RFTOP,rtop,RPOTB,RFPOTH,RCORRECTION,
     &      RCORRECTIONOLD,RALPHA,RBETA,RHSS,RGSS,RGFF,RF1,RF2,RHS,
     &      RFPOTH_confined,RFPOTH_unconfined,rcunconfined,rcconfined,
     &      rcond,rdum,rgvfac1,rgvfac2,rhedtip
      COMPLEX(8) CZI,CZ0,CZ
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      DIMENSION DRB(*),CZI(*),ITYPE(*)
      LINTERFACEFLOW=LFINTERFACE()
      IF (NLSH.EQ.0) RETURN
      iadfirst=0
      DO 10 I=1,NLSH
      IAD=KLSPTH(I)
      J=J+1
      CZ0=0.5D0*(CLSZS(IAD)+CLSZE(IAD))
      IF (LSGIV(IAD)) then
          drb(j)=0.0d0    !  "skip" line-sinks with a given sink density; force delta sigma equal to zero
      else
      RSIG=RLSIG(IAD)
      CZ=CZI(J)
      RC=RLSRES(IAD)
      RW=RLSWID(IAD)
      LRESISTANCE=RC*RW.NE.0.0D0
      RK=RFPERM(CZ)
      rbase=rfbase(cz)
      rtop=rftop(cz)
      RTHICK=RFTOP(CZ)-rbase
      RHEAD=RFHEDP(RLSPOT(IAD),CZ)   ! actual head at the line-sink center
      IF (LRESISTANCE) THEN  ! resistance line-sink, calculate correction term
        ralpha=1.0d+0
        rbeta=0.0d+0
        IF (LINTERFACEFLOW) THEN  ! interface flow
          CALL INTERFACEDATA (RHSS,RGSS,RGFF,RF1,RF2)
          RHSS=RHSS-RBASE   ! sea level with respect to aquifer base
          if (rhead-rbase.lt.rf1*rhss) then ! interface present
            IF (RHEAD.GE.RFTOP(CZ)) THEN ! confined
              RALPHA=RF2/RF1
              RBETA=-RF2*RHSS+RTHICK
            ELSE                         ! unconfined
              RALPHA=RF2
              RBETA=-RF2*RHSS
            END IF
          endif
        END IF
        IF (lsfirst) THEN
           RH=RLSH(IAD)-RBASE
           IF (LSGALLERY(IAD)) RH=RLSHMIN(IAD)-RBASE
        ELSE
           RH=(RHEAD+RLSH(IAD))*0.5D0-RBASE  ! "effective" saturated thickness for unconfined conditions
           rh0=rlshmin(iad)-rbase
           IF (RHEAD.GE.RTOP.AND..NOT.LINTERFACEFLOW) RH=RTHICK
           IF (LSGALLERY(IAD)) RH=MAX(RH,rh0)        ! NOTE: rh0 used to be 0.0D+0
        ENDIF
        RH=RALPHA*RH+RBETA
        RCORRECTION=RK*RC*RH*RSIG/RW
c            WRITE (ILUER,1004) RK,RC,RH,RSIG,RW,RLSH(IAD),RALPHA,RBETA
c 1004 FORMAT (' LSKNO4: RK,RC,RH,RSIG,RW,RLSH,RALPHA,RBETA ',8G11.4)
      END IF
c
c     -----------------------------
c
      IF (ITYPE(J).EQ.1) THEN   ! head specified line-sink or stream or drain
c    ----------------------------------------------------------------------- make potential for specified head PHI_s
      IF (RSIG.NE.0.0D0.AND.LRESISTANCE) THEN   ! near-field line-sink
c
c       convert the stream water level (rlsh(iad)) into a potential using the equation
c       for unconfined or confined flow based on the actual head (rhead) in the aquifer.   HMMM-- Does not include interface flow or not on the actual head!!!
c
        RHED=RLSH(IAD)-rbase  ! specified head at line-sink center with repsect to aquifer base
       IF (RHEAD.LT.RTOP) THEN   ! unconfined conditions
c          RPOTB=RFPOTH_unconfined (RLSH(IAD),CZ)   ! replaced below by code from this sub, but with RHEAD as the interface test head
         IF (LINTERFACEFLOW) THEN ! --------- interface flow -------------
         CALL INTERFACEDATA (RHSS,RGSS,RGFF,RGVFAC1,RGVFAC2)
         RHS=RHSS-RBASE   ! sea level with respect to aquifer base
         RHEDTIP=RGVFAC1*RHS
c      write (ilume,1001) rhed,rhs,rgvfac1,rgvfac2,rhedtip
c 1001 FORMAT(' rfpoth1: rhed,rhs,rgvfac1,rgvfac2,rhedtip ',5(d14.7))
          IF (RHEAD-rbase.GE.RHEDTIP) THEN ! no interface
            RPOTB=0.5D0*RK*RHED*RHED
c      write (ilume,1004) rfpoth
c 1004 format (' rfpoth4 unconfined, no interface: rfpoth=',d14.7)
          ELSE ! interface present
            RALPHA=RGVFAC2
            RBETA=-RGVFAC2*RHS
            RCUNCONFINED=0.5D0*RK*RGVFAC1*RHS*RHS
            RDUM=(RHED+RBETA/RALPHA)
            RPOTB=0.5D0*RK*RALPHA*RDUM*RDUM+RCUNCONFINED
c      write (ilume,1005) ralpha,rbeta,rcunconfined,rdum,rfpoth
c 1005 format(' rfpoth5 unconfined+interface: ralpha,rbeta,rcunconfined,'
c     &        ,' rdum,rfpoth=',/,5(d14.7))
          END IF
         ELSE  ! --------------- original logic without interface -------------
          RPOTB=0.5D0*RK*RHED*RHED
          IF (RHED.LE.0.0D0) RPOTB=0.0D0
         ENDIF
       ELSE                           ! confined conditions
c          RPOTB=RFPOTH_confined (RLSH(IAD),CZ) ! replaced by code from this sub, but with RHEAD as the interface test head
        IF (LINTERFACEFLOW) THEN ! --------- interface flow -------------
          CALL INTERFACEDATA (RHSS,RGSS,RGFF,RGVFAC1,RGVFAC2)
          RHS=RHSS-RBASE   ! sea level with respect to aquifer base
          RHEDTIP=RGVFAC1*RHS
c      write (ilume,1001) rhead,rhs,rgvfac1,rgvfac2,rhedtip
c 1001 FORMAT(' rfpoth1: rhead,rhs,rgvfac1,rgvfac2,rhedtip ',5(d14.7))
           IF (RHEAD-rbase.GE.RHEDTIP) THEN ! no interface   NOTE: Based on the actual head, not the specified head to be comnverted into a potential
            RCOND=0.5D0*RK*RTHICK*RTHICK
            RPOTB=RK*RTHICK*RHED-RCOND
c      write (ilume,1002) rcond,rfpoth
c 1002 format (' rfpoth2 confined, no interface: rcond,rfpoth=',2(d14.7))
           ELSE ! interface present
            RALPHA=RGVFAC2/RGVFAC1
            RBETA=-RGVFAC2*RHS+RTHICK
            RCCONFINED=RK*RGVFAC1*RTHICK*RHS
     &                 -0.5D0*RK*RGVFAC1*RTHICK*RTHICK
            RDUM=(RHED+RBETA/RALPHA)
            RPOTB=0.5D0*RK*RALPHA*RDUM*RDUM+RCCONFINED
c      write (ilume,1003) ralpha,rbeta,rcconfined,rdum,rfpoth
c 1003 format (' rfpoth3 confined+interface: ralpha,rbeta,rcconfined,',
c     &        ' rdum,rfpoth=',/,5(d14.7))
           END IF
        ELSE  ! --------------- original logic without interface -------------
          RCOND=0.5D0*RK*RTHICK*RTHICK
          RPOTB=RK*RTHICK*RHED-RCOND
        ENDIF
       ENDIF
      ELSE                                     ! far-field line-sink
        RPOTB=RFPOTH(RLSH(IAD),CZ)
      ENDIF
c      write (iluer,1005) rhead,rtop,rlsh(iad),rpotb
c 1005 format (' lskno5: rhead,rtop,rlsh(iad),rpotb ',4(D14.7))
      DRB(J)=DRB(J)+RPOTB-RLSPOT(IAD)
C      IF (CZI(J).EQ.CZ0) THEN        !  are at control point
c      IF (CZI(J).NE.CZ0) THEN
c        WRITE (ILUER,1001) J,CZI(J),CZ0
c 1001 FORMAT(' ***ERROR in LSKNO: equation ',i3,/,' collocation point ',
c     &2(D14.7),' does not coincide with line-sink center ',2(D14.7))
c      END IF
       IF (LRESISTANCE)  DRB(J)=DRB(J)+RCORRECTION  ! correction term
c       if (lresistance) write (iluer,1006) rpotb,rlspot(iad),rcorrection
c     &                                    ,drb(j)
c 1006 format (' lskno6: rpotb,rlspot(iad),rcorrection,drb(j) ',4(D14.7))
      ENDIF
C      END IF
c
c     -------------------------------------------
c
      IF (ITYPE(J).EQ.-1) THEN  ! successive line-sinks in gallery
c        WRITE (ILUER,1002)
c 1002 FORMAT(' ***WARNING: in LSKNO it is assumed that galleries do ',
c     &       'not cross inhomogeneity boundaries')
        DRB(J)=DRB(J)+RLSPOT(iadfirst)-RLSPOT(IAD)
        DRB(J)=DRB(J)+RCORRECTION-RCORRECTIONOLD
c        WRITE (ILUER,1003) RCORRECTION,RCORRECTIONOLD
c 1003 FORMAT (' LSKNO3: RCORRECTION,RCORRECTIONOLD ',2(D14.7))
        DRB(J0)=DRB(J0)-ABS(CLSZS(IAD)-CLSZE(IAD))*RSIG
      END IF
c
c     --------------------------------------
c
      IF (ITYPE(J).EQ.6) THEN  ! first line-sink in gallery
            DRB(J)=DRB(J)+RLSQ(IAD)
            iadfirst=iad
            J0=J
            DRB(J0)=DRB(J0)-ABS(CLSZS(IAD)-CLSZE(IAD))*RSIG
            RCORRECTIONOLD=RCORRECTION
      END IF
      endif
  10  CONTINUE
      RETURN
      END
C
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSSUB (DRB,J)
C
C---------------------------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) J,I,IAD,j0
      REAL(8) DRB,RSIGAV,RTEST
      INCLUDE 'lscom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      DIMENSION DRB(*)
      IF (NLSH.EQ.0) RETURN
      j0=j
      DO 10 I=1,NLSH
      IAD=KLSPTH(I)
      J=J+1
      RLSIG(IAD)=RLSIG(IAD)+DRB(J) ! solution vector contains strength increment
  10  CONTINUE
C
C     Filter out very small strength parameters
C  
      RSIGAV=0.0
      DO 20 I=1,NLSH    ! calculate average absolute strength
      IAD=KLSPTH(I)
      RSIGAV=RSIGAV+ABS(RLSIG(IAD))
  20  CONTINUE
      RSIGAV=RSIGAV/NLSH
      RTEST=1.0E-06*RSIGAV    ! test for strengths which are E-6 times average
      j=j0
      DO 30 I=1,NLSH
      IAD=KLSPTH(I)
      j=j+1
      IF (ABS(RLSIG(IAD)).LE.RTEST) THEN  ! set small strengths equal to zero
      RLSIG(IAD)=0.0
      LSGIV(IAD)=.TRUE.
      ENDIF
c
c      IF (LSDRAIN(IAD).AND.RLSIG(IAD).LT.0.0) RLSIG(IAD)=0.0 ! prevent recharging drains
! this statement is blocked off to preserve water balance in the groundwater flow solution.
  30  CONTINUE
      lsfirst=.false.
      RETURN
      END
c
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSUPDATE (DRSCR,j,M,N,lsubinclude,isubsign)
C
C---------------------------------------------------------------------------------------------------------
C
C     Routine updates head specified line-sink data between iterations:
C     1) Update the array RLSPOT, which contains the calculated
C     potential at the centers of the line sinks.
C     2) Update the head in the array RLSH for galleries with resistance
c
c
c     NOTE: make sure that DRA contains uncorrected matrix coefficients
c           make sure that DRB contains the latest solution vector
c
c
C
      IMPLICIT NONE
      INTEGER(4) I,j,M,N,iad,nsolOut,isubsign
      LOGICAL LRESISTANCE,lflksolving,
     &        lsubinclude,laddsubcells,
     &        lsolOut,loadsolOut,linalreadyOut,
     &        lErrorReportOut,lDirectfromDiskOut
      REAL(8) RFPOT,RHEAD,RC,RW,RFHEDP,DRSCR,rdum,rsubsign
      COMPLEX(8) CZ,cflk_subomega
      CHARACTER(8) aBasenameOut
      CHARACTER(16)aDateTimeOut
      DIMENSION DRSCR(*)
      INCLUDE 'lscom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
C
      IF (NLSH.EQ.0) RETURN
      laddsubcells=lsubinclude
      rsubsign=REAL(isubsign)
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      lsfirst=nsolOut.eq.0
      DO I=1,NLSH
      j=j+1
      IAD=KLSPTH(I)
      CZ=(CLSZS(IAD)+CLSZE(IAD))*0.5
      IF (LSFIRST) THEN   
        RLSPOT(IAD)=RFPOT(CZ) ! will give bad result, but necessary to build on (use in LSKNO)
      ELSE
        RLSPOT(IAD)=RLSPOT(IAD)+drscr(j) ! ADD NEW CONTRIBUTION
        if (laddsubcells) then
          RLSPOT(IAD)=RLSPOT(IAD)+rsubsign*REAL(cflk_subomega(cz)) ! add/subtract sub-cell contributions
        end if
      ENDIF
c      if (i.eq.1) then
c        rdum=rfpot(cz)
c        write (iluer,1001) i,rlspot(iad),rdum
c 1001 format (' lsupdate1: i,rlspot,rfpot -->',i3,2x,d14.7,2x,d14.7)
c      end if
      IF (LSGALLERY(IAD)) THEN  ! update RLSH on all line-sinks in galleries
       RC=RLSRES(IAD)
       RW=RLSWID(IAD)
       if (lsfirst) then
         rhead=rlsh(iad)
       else
         RHEAD=RFHEDP(RLSPOT(IAD),CZ)
       end if
       lresistance=rc*rw.ne.0.0D0
       IF (LRESISTANCE) THEN
         RHEAD=RHEAD-RLSIG(IAD)*RC/RW
       END IF
       RLSH(IAD)=RHEAD
      END IF
      END do
      END
c
c---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSUPDATE_check (j,M,N)
C
C---------------------------------------------------------------------------------------------------------
C
C     Routine is for debugging purposes only
c     A comparison is made between RLSPOT and RFPOT
c
c
C
      IMPLICIT NONE
      INTEGER(4) I,iad,j,M,N
      REAL(8) RFPOT,rdum
      COMPLEX(8) CZ
      INCLUDE 'lscom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
C
      IF (NLSH.EQ.0) RETURN
      DO I=1,NLSH
      j=j+1
      IAD=KLSPTH(I)
      CZ=(CLSZS(IAD)+CLSZE(IAD))*0.5
      rdum=rfpot(cz)
      write (ilume,1001) j,rlspot(iad),rdum
 1001 format (' LSUPDATE_check1: j,rlspot(iad),rfpot(cz) ',
     &          i4,2x,2(d14.7))
      enddo
      return
      END
c
c
C---------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSOUTLETUPDATE
C
C---------------------------------------------------------------------------------------------------------
C
C     Routine updates restances for first line-sink in a lake outlet stream
C     to force the proper outflow.
C
C     **************   CURRENTLY NOT USED *********************
C
C     The procedure is not as robust as forcing the first line sink in an
C     outlet string to extract the outflow. This is done in GENSTREAMFLOW.
C
      IMPLICIT NONE
      INTEGER(4) I,ISTR
      LOGICAL LRESISTANCE
      REAL(8) RHEAD,RC,RW,RFHEDP,RLENGTH,RQ,
     &        rflstable
      COMPLEX(8) CZ
      INCLUDE 'lscom.inc'
      INCLUDE 'tracom.inc'
      INCLUDE 'lusys.inc'
      IF (.NOT.LSFIRST) THEN  ! do not force flow yet in first iteration
      DO ISTR=1,NLSTRING     !  Update resistance for outlet line-sink
        I=KLSTRING(ISTR)
        IF (LSOUTLET(I)) THEN  ! found outlet string, now calculate resistance
          CZ=0.5D+0*(CLSZS(I)+CLSZE(I))
          RHEAD=RFHEDP(RLSPOT(I),CZ)
          RLENGTH=ABS(CLSZE(I)-CLSZS(I))
          RQ=rflstable(RHEAD,ISTR)
          IF (RQ.EQ.0.0D0) THEN
          WRITE (ILUER,1000) ALSLAB(I)
            RLSRES(I)=1.0D5
          ELSE
            RLSRES(I)=(RHEAD-RLSH(I))*RLSWID(I)*RLENGTH/RQ
          END IF
          IF (RLSRES(I).LT.0.0D+0) THEN  ! lake stage (RHEAD) below outlet head (RLSH(IAD))
            RLSRES(I)=1.0D5
            WRITE (ILUER,2000) ALSLAB(I)
          END IF
            write (iluer,1001) rhead,rlsh(i),rq,rlsres(i)
 1001 format (' LSOUTLETUPDATE1:  rhead,rlsh,rq,rlsres ',4(D14.7))
        END IF
      END DO
      END IF
      RETURN
 1000 FORMAT (' ***WARNING: Calculated outflow rate is zero'
     &        ' for line-sink ',a16,/,' large resistance set, but'
     &        ' some water will flow back into the lake.')
 2000 FORMAT (' ***WARNING: Lake stage dropped below outlet head'
     &        ' for line-sink ',a16,/,' large resistance set, but'
     &        ' some water will flow back into the lake.')
      END
c
c -----------------------------------------------------------------------------------
c
      REAL(8) function rflstable(rhead,istr)
c
c       The function returns the linear interpolation of RLSTABLE(i,itable)
c       in the *.STG file using the argument RHEAD and RLSTAGE(i,itable).
c       The array RLSTABLE(i,itable) may be the outflow rate of an outlet stream
c       or the area of a lake, depending on where this function is called, see
c       the routine GENSTREAMFLOW in the file LSSTREAM.FOR
c
c
      implicit none
      INTEGER istr,itable,ilength,i
      REAL(8) rhead,rdh,rdq
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
c
      itable=nlstab(istr)
      ilength=nlstablength(istr)
      if (itable.le.0.or.ilength.le.0) then ! illegal data
        rflstable=0.0D+0  ! return arbitrarily a zero flow rate
        write (iluer,1000) istr,itable,ilength
        return
      end if
      if (rhead.le.rlstage(1,itable)) then ! rhead below lowest stage
        rflstable=rlstable(1,itable)
      else if (rhead.ge.rlstage(ilength,itable)) then ! rhead above highest stage
        rflstable=rlstable(ilength,itable)
      else   ! rhead between stages
        do i=2,ilength
        if
     &(rhead.le.rlstage(i,itable).and.rhead.ge.rlstage(i-1,itable)) then
          rdh=rlstage(i,itable)-rlstage(i-1,itable)
          rdq=rlstable(i,itable)-rlstable(i-1,itable)
      rflstable=rlstable(i-1,itable)+(rhead-rlstage(i-1,itable))*rdq/rdh
        end if
        end do
      end if
      return
 1000 format (' ***ERROR in rflstable: incorrect stage and discharge',
     & ' table data:',/,'istr,itable,ilength=',3I5)
      end
C
C ---------------------------------------------------------------------------------
C
      subroutine lsreadtables
c
c     Routine reads the lake stage versus outflow rate tables for
c     stream sections (strings) that are lake outlet streams or lakes.
c
      implicit none
      INTEGER istr,iad,itable,itablelength,ilu,itemp,i
      LOGICAL lret
      REAL(8) rdum1,rdum2,rvar
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
c
      if (lsfirst.and.nlstring.gt.0) THEN ! see is there are tables to be read
      itable=0
      itablelength=0
      ilu=2
      do istr=1,nlstring
      iad=klstring(istr)
      if (lsoutlet(iad).or.lslake(iad)) THEN ! found outlet or lake string, now read table
        itable=itable+1
        afile=alstblfilename(istr)
        call filech ('.STG')   ! check or add extension
        call opafil(ilu,-2,lret)
        if (lret) then  ! failed to open file
           write (iluer,1000) afile
         amess(1)='Could not open table file for outlet stream or lake.'
         amess(2)='Check existence of table file.'
           call halt(2) ! stop program execution for batch version
        end if
        itemp=iluin ! temporarily store iluin and replace for this read
        iluin=ilu
  10    call inline
        rdum1=rvar(1)
        rdum2=rvar(2)
        if (.not.lerror) then
          itablelength=itablelength+1
          rlstage(itablelength,itable)=rdum1
          rlstable(itablelength,itable)=rdum2
          GOTO 10
        end if
        lerror=.FALSE. ! assume end of file
        close (iluin)
        iluin=itemp ! retore iluin
        if (itablelength.eq.0) THEN ! empty file, abort program execution
          write (iluer,1000) alslab(iad),afile
          AMESS(1)='Empty table file encountered.'
          AMESS(2)='Check table files for outlet streams.'
          CALL HALT(2) ! stop program execution for batch version
        end if
        nlstab(istr)=itable
        nlstablength(istr)=itablelength
C
c      write (iluer,1001) istr,itable,itablelength
 1001 format (' rflstable1: istr,itable,ilength ',3I5)
c      write (iluer,1002) (i,rlstage(i,itable),rlstable(i,itable),
c     &i=1,itablelength)
 1002 format (' rflstable2:  i      stage     outflow or area'
     &         ,/,4(12x,I3,1x,2(D14.7),/))
C
        itablelength=0
      end if
      end do
      end if
 1000 format (' ***ERROR in LSREADTABLES: Could not open stage table '
     &        'file ',A16,/,' Make sure file is in project path '
     &        '(see Project Settings in GUI).')
 2000 format(' ***ERROR in LSREADTABLES: Line-sink string starting with'
     &        ' line-sink ',A16,/,' has no or illegal data in table '
     &        ' file ',A16)
      end subroutine
c
c ----------------------------------------------------------------------------------------
c
      subroutine set_lakedone()
c
c     Routine checks for the presence of lake features.
c     If present the logical LAKEDONE is set FALSE as there is no solution yet.
c     If not present the logical LAKEDONE is set TRUE, since no lake convergence criteriua is to be met.
c
      implicit none
      INTEGER istr,iad
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
c
      lakedone=.false.
      do istr=1,nlstring
      iad=klstring(istr)
      if (lslake(iad)) THEN ! found lake
        lakedone=.false.
      endif
      enddo
      end subroutine
c
c
c ----------------------------------------------------------------------------------------
c
      LOGICAL function lflslakeiterations(niter)
c
c     Routine generates a new estimate for lake levels during solution procedure
c
c     Input:
c            niter                   number of outer loop iterations to find lake level
c
c     Output:
c            lflslakeiterations      .TRUE. if lakes found, niter not exceeded and heads adjusted
c
c     This function should be called AFTER the first set of groundwater and
c     surface water iterations (inner loop)
c
      implicit none
      INTEGER i,niter,istr,iad,iend,inext,istart
      LOGICAL lakeiterations
      REAL(8) rh,rh1,rh2,rq1,rq2,rq,rqf,rflstable
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      DATA lakeiterations /.false./
c
      lflslakeiterations=.false.
      if (nlakeiterations.lt.niter) then ! not done yet
c      write (iluer,1004) nlakeiterations,niter
c 1004 format (' lflslakeiterations4: nlakeiterations,niter ',2(i3,2x))
       do istr=1,nlstring
       iad=klstring(istr)
       if (lslake(iad)) THEN ! found lake
        lakeiterations=.true.
        rq=rlsq(iad)         ! first calculate water balance
        rh=rlsh(iad)
  5     inext=klsdn(iad)
        if (inext.gt.0) then ! next line-sink
         rq=rq+rlsq(inext)
         rqf=rlsbf(inext)+rlsof(inext) ! get "streamflow" in case this is the last line-sink
c         write (iluer,1002) ISTR,INEXT,RQ,RQF
c 1002 FORMAT (' lflslakeiterations2: istr,inext,rq,rqf ',
c     &       2(i3,1x),2(d14.7,1x))
         iad=inext
         GOTO 5
        else                 ! at end of string
         rq=rq+rqf           ! add "streamflow" in lake
         rq=rq-rflstable(rh,istr)*rlsevap(istr) ! subtract evapotranspiration
         rlsh1(istr)=rlsh2(istr)    ! shift h & q
         rlsq1(istr)=rlsq2(istr)
         rlsh2(istr)=rh             ! store h & q
         rlsq2(istr)=rq
        end if
c         write (iluer,1111) rflstable(rh,istr),rlsevap(istr)
c 1111 format (' lflslakeiterations3:rflstable(rh,istr),rlsevap(istr)',
c     &2(d14.7))
c         write (iluer,1003) istr,rlsh1(istr),rlsh2(istr),
c     &                      rlsq1(istr),rlsq2(istr)
c 1003 format (' lflslakeiterations3: istr,rlsh1,rlsh2,rlsq1,rlsq2 ',
c     &I3,2x,4(d14.7,1x))
        if (nlakeiterations.eq.0) then ! use second estimated lake level (shifted into rlsh1)
         rh=rlsh1(istr)
        else                  ! calculate new lake level
         rh1=rlsh1(istr)
         rh2=rlsh2(istr)
         rq1=rlsq1(istr)
         rq2=rlsq2(istr)
         rh=rh1-rq1*(rh2-rh1)/(rq2-rq1)   ! secant method for estimating a new lake stage
        end if
c         write (iluer,1001) istr,nlakeiterations,rh
c 1001 format (' lflslakeiterations1: istr,nlakeiterations,rh ',
c     &         i3,1x,i3,1x,d14.7)
         istart=klstring(istr)
         iend=klstrend(istr)
         do i=istart,iend
         rlsh(i)=rh
         end do
       end if
       end do
       if (lakeiterations) then
        nlakeiterations=nlakeiterations+1
        lflslakeiterations=.true.
       end if
      end if
      return
      end

