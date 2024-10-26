C     Last change:  HMH   8 Aug 2005    9:57 pm
c       This file contains the following subroutines and functions:
c
C	SUBROUTINE LSNETWORK       generates a stream network
C	SUBROUTINE GENSTREAMFLOW   generates streamflow and eliminates some line-sinks from solution
c       SUBROUTINE LSOPENEND       checks for stream reaches that are not connected in the network
c       SUBROUTINE lsread_relaxation_file  reads a table with relaxation factors versus iterations for
c                                          streamflow calcutations.
c
C
C----------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSNETWORK
C
C----------------------------------------------------------------------------------------------------
C
C
C     Routine connects the streamsections entered with the stream command.
C     The downstream end of a string is linked to the upstream end of
C     a line sink in another string, unless the downstream end is a 
C     "stream end"!
C     The orientation of a string is determined by inspecting the
C     "down stream address" of the first line sink in the string.
C     Stream flow orientation is determined in RELSPT called in LSIN.
C     If a string consists of only one line sink its orientation is 
C     undetermined and the networking routine tries to link either end to 
C     a line sink in an other string. Single line sink strings are not prefered
C     and should be limited to model (short) branches in stream net works.
C     The network routine will only try to link a line sink string to another line sink
C     string if the latter line sink is down stream (has a lower head).
C     This limitation does not hold when linking a string to a "Lake string".
C     NOTE: upstream pointer of a line sink just downstream of a branch
C     will point to the last connected string. Consequently, when moving
C     upstream through the network, you will branch off along the last
C     entered strings. Baseflow routing, however, is done DOWNSTREAM!!
C     For baseflow routing, move downstream from all headwaters
C     (start addresses of strings with KLSUP pointing to itself) and
C     add groundwater inflow to RLSBF for each line sink.
C
      IMPLICIT NONE
      INTEGER(4) ISTR,IAD1,IADLAST,ISTR1,INEXT1,IAD0,ISTR2,IAD2,
     &           INEXT2
      LOGICAL LDN1,LUP1,LSIBR1,LDN2,LUP2,LSIBR2,LSKIP
      REAL(8) RDIS,RMIN
      COMPLEX(8) CZS,CZE
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
c      character(1) adum
C
C      write (iluout,1000) (i,klsup(i),klsdn(i),i=1,nls)
C 1000 format (' lsnetwork0: entering LSNETWORK.',/,
C     &' i,klsup,klsdn ',/,11(3i5,/))      
      IF (NLSTRING.LE.1) RETURN
      DO 2 ISTR=1,NLSTRING      ! unlink strings
      IAD1=KLSTRING(ISTR)
      IADLAST=KLSTREND(ISTR)
      KLSUP(IAD1)=IAD1
      IF (KLSDN(IADLAST).NE.0) KLSDN(IADLAST)=IADLAST
  2   CONTINUE      
C      write (iluout,1001) (i,klsup(i),klsdn(i),i=1,nls)
C 1001 format (' lsnetwork1: i,klsup,klsdn (links undone)',/,
C     &11(3(i5,2x),/))
      DO 20 ISTR1=1,NLSTRING    ! start linking strings
C      write (ilume,1007) istr1
c      read (iluin,1008) adum      
c      if (adum.eq.'R'.or.adum.eq.'r') return      
C 1007 format (' Press <Enter> to link string ',I3,' else "R"')
c 1008 format (a1)      
      IAD1=KLSTRING(ISTR1)
      INEXT1=KLSDN(IAD1)
      lsibr1=.false.
      IF (INEXT1.EQ.0) GOTO 20  ! single line-sink in end stream, don't try to link
      LDN1=INEXT1.GT.IAD1       ! linesinks oriented downstream
      LUP1=INEXT1.LT.IAD1       ! linesinks oriented upstream
      IF (INEXT1.EQ.IAD1) THEN  ! single line sink in string, try link both ends
      LSIBR1=.TRUE.
      LDN1=.TRUE.           ! start linking end node
      ENDIF
c      write (iluout,1002) istr1,iad1,inext1,ldn1,lup1,lsibr1
c 1002 format (' lsnetwork2: istr1,iad1,inext1,ldn1,lup1,lsibr1 ',3i5,
c     &3l3,/,' start of new string to be linked to other strings')
   5  IF (INEXT1.EQ.IAD1) THEN  ! open string end, search for nearest node
      RMIN=1.0E21
      IAD0=0
   6  IF (LDN1) CZE=CLSZE(IAD1) ! line sinks oriented downstream
      IF (LUP1) CZE=CLSZS(IAD1) ! line sinks oriented upstream
c      write (iluout,1003) inext1,klsup(inext1),klsdn(inext1),lup1,ldn1,
c     &cze
c 1003 format (' lsnetwork3: inext1,klsup(inext1),klsdn(inext1),lup1, ',
c     &'ldn1,cze ',/,3i5,2l5,2G11.4,/,
c     &' found open end, try to connect CZE to other string')
      DO 10 ISTR2=1,NLSTRING
c      write (ilume,1009) istr1,istr2
c      read (iluin,1008) adum            
c      if (adum.eq.'R'.or.adum.eq.'r') return
c 1009 format (' Press <Enter> to link string ',I3,' to string ',I3,
c     &' else "R".')
      IF (ISTR2.EQ.ISTR1) GOTO 10  ! do not try to link string to itself
      IAD2=KLSTRING(ISTR2)
      INEXT2=KLSDN(IAD2)
      LDN2=INEXT2.GT.IAD2       ! determine string orientation based on 
      LUP2=INEXT2.LT.IAD2       ! first two line sinks in the string
      LSIBR2=.FALSE.
      IF (INEXT2.EQ.IAD2) THEN  ! single line sink in string, try both ends
      LSIBR2=.TRUE.
      LDN2=.TRUE.                  ! start linking to end node
      ENDIF
  7   IF (LDN2) CZS=CLSZS(IAD2) ! line sinks oriented downstream
      IF (LUP2) CZS=CLSZE(IAD2) ! line sinks oriented upstream
      RDIS=ABS(CZE-CZS)
      LSKIP=RLSH(IAD2).GE.RLSH(IAD1) ! do not try to connect to upstream ls
      if (lslake(iad2)) lskip=.FALSE. ! overrule when linking to lake string
      LSKIP=LSKIP.OR.IAD1.EQ.IAD2    ! do not try to connect ls to itself
c      WRITE (ILUOUT,1004) ISTR1,ISTR2,IAD1,IAD2,CZS,CZE,lskip
c 1004 FORMAT (' ISTR1,ISTR2,IAD1,IAD2,CZS,CZE,lskip ',4I3,4G11.4,l5)      
      IF (RDIS.LT.RMIN.AND..NOT.LSKIP) THEN     ! may be correct node, save
      IAD0=IAD2
      RMIN=RDIS
      ENDIF
      IF (LSIBR2) THEN          ! single line sink in string
      LUP2=.TRUE.               ! try also to link to start node
      LSIBR2=.FALSE.
      GOTO 7
      ENDIF
      INEXT2=KLSDN(IAD2)        ! move down the string
c      write (iluout,1005) iad0,rmin,iad2,inext2
c 1005 format (' lsnetwork5: iad0,rmin,iad2,inext2 ',i5,G11.4,2i5,/,
c     &' if not open string and not end string, continue search on same '
c     &'string (setting iad2=inext2)')      
      IF (INEXT2.EQ.IAD2) GOTO 10      ! open string, go on to next string
      IF (INEXT2.EQ.0) GOTO 10         ! end of stream, go on to next string
      IAD2=INEXT2
      GOTO 7      ! continue to check for node proximity on same string
  10  CONTINUE
      IF (LSIBR1) THEN     ! single line sink in string
      LUP1=.TRUE.          ! try to link start node also
      LSIBR1=.FALSE.
      GOTO 6
      ENDIF 
      IF (IAD0.GT.0) THEN  ! link string ISTR1 to node IAD0 in string ISTR2
        KLSDN(IAD1)=IAD0     
        KLSUP(IAD0)=IAD1
      ENDIF
      ELSE                 ! not yet at string end, proceed down stream
      IAD1=INEXT1 
      INEXT1=KLSDN(IAD1)        ! search downstream for open end
      IF (INEXT1.EQ.0) GOTO 20  ! end of stream start at new string
      GOTO 5
      ENDIF
  20  CONTINUE             ! continue to link up next string
c      write (iluout,1006) (i,klsup(i),klsdn(i),i=1,nls)
c 1006 format (' lsnetwork6: leaving LSNETWORK.',/,
c     &' i,klsup,klsdn ',/,11(3i5,/))        
C ---- check for open ends (LSOPENEND) deferred to LSCUR call or GENSTREAM
C ---- since more streams may be added at this point.
      RETURN
      END
C
C----------------------------------------------------------------------------------------------------
C
      SUBROUTINE GENSTREAMFLOW (icode,lErrorReport,nsol,
     &                          ilstablelength,niterarray,relaxarray)
C
C----------------------------------------------------------------------------------------------------
C
C     Routine generates baseflow and overland flow in stream networks,
C     using the line sink discharges and overland inflow rates, respectively.
C     If ICODE = 0 negative streamflow is not corrected.
C     If ICODE = 1 go through the stream network to adjust line sink discharges
C                  in order to avoid negative stream flows and overinfiltration.
C     lErrorReport=true report streamflow error
C     Streamflow is defined as the sum of overland flow and baseflow
C
C
      IMPLICIT NONE
      INTEGER(4) ICODE,I,ISTR,IAD,INEXT,ILS,ILAKE,
     &           ilstablelength,niterarray,nsol
      LOGICAL LCHANGE,LGIVEN,LSIGNIFICANT,lreintroduce,lErrorReport
      REAL(8) RLSRELAX,RLENGTH,RBFSUM,ROFSUM,RBFLOCAL,ROFLOCAL,
     &        RBFCORRECT,RBFCORTEMP,RFLOW,RQS,RHEAD,RFHEDP,RHT,
     &        RQMAX,RFPERM,RS,RDQ,RDSIG,rflstable,RQ,RDIST,RL,RQF,
     &        relaxarray
      COMPLEX(8) CZ
      DIMENSION niterarray(*),relaxarray(*)
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      lreintroduce=.true.  ! may be used to suppress reintroduction of stream sections
c                            during some initial iterations. Does not work well.
c      write (ilume,1001) nsol,lreintroduce
c 1001 format (' genstreamflow1: ncall, lreintroduce ',I4,1x,l4)
      do i=1,ilstablelength
c      if (nsol.ge.niterarray(i)) rlsrelax=relaxarray(i) ! temporary for when extra solve in solut is bypassed
      if (nsol-1.ge.niterarray(i)) rlsrelax=relaxarray(i)
      end do
      write (ilume,1000) rlsrelax
      IF (NLSTRING.EQ.0) THEN
      WRITE (ILUER,3000)
      RETURN
      ENDIF
      CALL LSOPENEND                      ! check and correct open ends
      DO I=1,NLS
      IF (LSLAKE(I)) RLSQ(I)=0.0D0 ! reset any lake outflows (via outlet streams)
      RLSBF(I)=0.123456E+21  ! flag: not yet included in base flow calculations
      RLSOF(I)=0.0           ! zero out overland flow
      END DO
C -------------------------   initial baseflow and overland flow calculations.
      DO 20 ISTR=1,NLSTRING
      IAD=KLSTRING(ISTR)
      IF (KLSUP(IAD).EQ.IAD) THEN   ! (1) head water, start calculations
      RLENGTH=ABS(CLSZE(IAD)-CLSZS(IAD))
      IF (LSOUTLET(IAD)) THEN ! (2) first line sink in outlet string, force sink density
       ILAKE=0                ! check if outlet stream is near a line-sink lake feature
       RL=RLENGTH
       DO I=1,NLS
       IF (I.NE.IAD) THEN
        RDIST=ABS(CLSZS(I)-CLSZS(IAD))
c        write (iluer,1004) iad,i,rdist,rl,lslake(i)
c 1004 format (' genstreamflow4: iad,i,rdist,rl,lslake ',
c     &          2(I4),2(d14.7),l4)
        IF (RDIST.LT.RL.AND.LSLAKE(I)) THEN ! we are close to a lake feature
         ILAKE=I
         RL=RDIST
        END IF
       END IF
       END DO
       IF (ILAKE.GT.0) THEN ! outlet stream connects to a line-sink lake feature
        LSCONNECT(IAD)=.TRUE.
        RHEAD=RLSH(ILAKE)
        RQ=RFLSTABLE(RHEAD,ISTR) ! use table for outlet stream to get outflow rate
        RLSOVHW(ISTR)=RQ         ! give outflow rate to stream as "end inflow"
        RLSQ(ILAKE)=-RQ              ! take outflow rate out of lake line-sink
       ELSE                 ! outlet stream not near line-sink lake, assume high-k lake
        LSCONNECT(IAD)=.FALSE.
        CZ=0.5D+0*(CLSZS(IAD)+CLSZE(IAD))
        RHEAD=RFHEDP(RLSPOT(IAD),CZ)
        RQ=RFLSTABLE(RHEAD,ISTR)
        RLSIG(IAD)=RQ/RLENGTH    ! give outflow rate to first line-sink and make "given"
        LSGIV(IAD)=.TRUE.
       ENDIF
c       write (iluer,1001) alslab(iad),rhead,rq
c 1001 format (' GENSTREAMFLOW1: alslab,rhead,rq ',a16,2(D14.7))
      ENDIF ! (2)
      RLSBF(IAD)=RLSIG(IAD)*RLENGTH+RLSOVHW(ISTR) ! baseflow + headwaterinflow = baseflow
      RBFSUM=RLSBF(IAD)
C add overland inflow to first element
      RLSOF(IAD)=RLSOFSIG(IAD)*RLENGTH ! overlandflow
      ROFSUM=RLSOF(IAD)
  5   INEXT=KLSDN(IAD)
      IF (INEXT.EQ.0)  GOTO 20      ! end of stream, next string
      RLENGTH=ABS(CLSZE(INEXT)-CLSZS(INEXT))   ! we are now moving downgradient along the stream
      RBFLOCAL=RLSIG(INEXT)*RLENGTH
      ROFLOCAL=RLSOFSIG(INEXT)*RLENGTH
      IF (RLSBF(INEXT).EQ.0.123456E+21) THEN ! add line sink strength to base flow
      RBFSUM=RBFSUM+RBFLOCAL
      ROFSUM=ROFSUM+ROFLOCAL
c
c this commented out code may not be correct at this stage. Perhaps negative baseflow is to be prevented lateron.
c
c      if (rbfsum.lt.0.0) then ! do not allow for negative baseflow, reduce overland flow instead
c        rofsum=rofsum+rbfsum
c        rbfsum=0.0
c      end if
      RLSBF(INEXT)=RBFSUM
      RLSOF(INEXT)=ROFSUM
      ELSE                          ! add new flow components only
      RLSBF(INEXT)=RLSBF(INEXT)+RBFSUM
      RLSOF(INEXT)=RLSOF(INEXT)+ROFSUM
      ENDIF
c      write (iluer,1002)inext,alslab(inext),lsgiv(inext),lsinlet(inext),
c     &      rlength,rlsig(inext),rlsbf(inext),rlsof(inext),rbfsum,rofsum
c 1002 format (' genstreamflow2: inext,alslab,lsgiv,lsinlet ', i3,a16,2l3
c     &     ,/,' rlength,rlsig,rlsbf,rlsof,rbfsum,rofsum ',/,6(d14.7))
      IAD=INEXT
      GOTO 5
      ENDIF ! (1)
  20  CONTINUE
      DO 21 I=1,NLS                 ! remove flags
      IF (RLSBF(I).EQ.0.123456E+21) RLSBF(I)=0.0
  21  CONTINUE      
      IF (icode.EQ.0) GOTO 45 ! no corrections, calculate max. flows and return
C ----------------------------- ELSE go through the network to eliminate
C                               negative stream flows and overinfiltration
      LSBASE=.TRUE.
      LCHANGE=.FALSE.
      DO 30 ISTR=1,NLSTRING    ! ********* START OF STREAMFLOW CORRECTIONS **************
      IAD=KLSTRING(ISTR)
      IF (LSLAKE(IAD)) GOTO 30 ! do not "correct" flow in lake features
      IF (KLSUP(IAD).EQ.IAD) THEN ! head water, update baseflow/strength---(0)
      RBFCORRECT=0.0
  25  continue
      IF (LSLAKE(IAD)) GOTO 30 ! do not "correct" flow in lake features
c      write (iluer,1001) iad,lsinlet(iad),alslab(iad)
c 1001 format (' genstreamflow1: iad,lsinlet,alslab ',i3,l3,1x,a16)
      IF (LSINLET(IAD)) THEN ! @ last element of inlet stream, infiltrate all streamflow
         RLENGTH=ABS(CLSZE(IAD)-CLSZS(IAD))
         RQS=(RLSBF(IAD-1)+RLSOF(IAD-1))/RLENGTH
         RLSIG(IAD)=-RQS
         RLSBF(IAD)=0.0D+0
         RLSOF(IAD)=0.0D+0
c         write (iluer,1011) iad,rlsig(iad),rlsbf(iad),rlsof(iad),
c     &    alslab(iad)
c 1011 format (' genstream11: iad,rlsig,rlsbf,rlsof,alslab',/,
c     & i3,3(d14.7),a16)
      GOTO 30  !  this is the last line-sink of an end stream, go to next head water
      ENDIF
      IF (LSOUTLET(IAD).AND..NOT.LSCONNECT(IAD)) THEN ! first element of string with FORCED strength, goto next element
C                                              this is a outlet stream connected to a high-k lake, not a line-sink lake
      INEXT=KLSDN(IAD)
      IF (INEXT.EQ.0) GOTO 30 ! end of stream, next string
      IAD=INEXT
      GOTO 25
      ENDIF
      RLSBF(IAD)=RLSBF(IAD)+RBFCORRECT ! next element: correct baseflow
      RBFCORTEMP=0.0 ! correction term based on this element only
      RFLOW=RLSBF(IAD)+RLSOF(IAD) ! get streamflow
      RLENGTH=ABS(CLSZE(IAD)-CLSZS(IAD))  
      RQS=RFLOW/RLENGTH
      CZ=0.5*(CLSZS(IAD)+CLSZE(IAD))
      RHEAD=RFHEDP(RLSPOT(IAD),CZ)
      IF (RLSRES(IAD).EQ.0.0) THEN ! --------------------------------------(1)
        RHT=RLSH(IAD)
        RQMAX=RFPERM(CZ)*RLSWID(IAD)
        IF (RLSWID(IAD).EQ.0.0) RQMAX=1.0E20 ! do not limit infiltration
      ELSE ! --------------------------------------------------------------(1)
        RHT=RLSH(IAD)-RLSDEP(IAD)
        RQMAX=RLSDEP(IAD)/RLSRES(IAD)*RLSWID(IAD)
      ENDIF ! -------------------------------------------------------------(1)
      IF (RQS.LT.-1.0E-6*RQMAX) THEN !!!! negative stream flow, correct ---(2)   will go wrong if c=w=0 (RQMAX=1.E20)
        if (lErrorReport) write (ILUME,5000) RFLOW,ALSLAB(IAD)
        if (lErrorReport) write (*,5000) RFLOW,ALSLAB(IAD)
        LCHANGE=.TRUE.
        IF (RLSIG(IAD).LT.RQS) THEN ! -------------------------------------(3)
          RBFCORTEMP=-RFLOW*rlsrelax
          RLSIG(IAD)=RLSIG(IAD)-RQS*rlsrelax
          LSGIV(IAD)=.TRUE.
        ELSE ! ------------------------------------------------------------(3)
          if (ABS(rlsig(iad)).lt.1.0d-06) then                                   ! new logic, idea is that neg. streamflow is remnant of some headinflow
            RBFCORTEMP=-rflow*rlsrelax                                           ! change due to, for instance, lake outflow changes
          else
            RBFCORTEMP=-RLSIG(IAD)*RLENGTH*rlsrelax
            RLSIG(IAD)=rlsig(iad)-rlsig(iad)*rlsrelax
          endif
          LSGIV(IAD)=.TRUE.
        ENDIF ! -----------------------------------------------------------(3)
      ELSE   !!!! positive stream flow; check proper infiltration ---------(2)
        LGIVEN=LSGIV(IAD)
        IF (LGIVEN) THEN !!!  Proper inf. rate? ---------------------------(4)
          IF (RHEAD.GT.RHT) THEN !! may be overinfiltrating----------------(5)
            IF (RLSRES(IAD).EQ.0.0.and.lreintroduce) THEN ! overinfiltrating----------------(6)
              LSGIV(IAD)=.FALSE. ! reintroduce in GW solution
              RLSIG(IAD)=0.0
              if (lErrorReport) write (ILUME,4010) ALSLAB(IAD)
              if (lErrorReport) write (*,4010) ALSLAB(IAD)
              LCHANGE=.TRUE.
            ELSE !! resistance specified line sink ------------------------(6)
              RS=(RLSH(IAD)-RHEAD)/RLSRES(IAD)*RLSWID(IAD)
              RDQ=RS+RLSIG(IAD)
              IF (RDQ.GT.0.0) THEN ! underinfitrating ? -------------------(7)
                RDSIG=MIN(RDQ,RQS)
                IF (RDSIG.GT.0.0) THEN ! yes underinfiltrating ------------(7)
                  LCHANGE=.TRUE.
                  if (abs(rlsig(iad)).gt.1.0e-5) then
                    lsignificant=abs(rdsig/rlsig(iad)).gt.1.0e-5
                  else
                    lsignificant=abs(rdsig).gt.1.0e-5
                  endif
                  if (lsignificant) then
                  if (lErrorReport) write (ILUME,4501) RDSIG,ALSLAB(IAD)
                  if (lErrorReport) write (*,4501) RDSIG,ALSLAB(IAD)
                  endif
                ENDIF
              ENDIF  
              IF (RDQ.LT.0.0.and.lreintroduce) THEN ! overinfiltrating ---------------------(7)
                RDSIG=RDQ
                LCHANGE=.TRUE.
                LSGIV(IAD)=.FALSE.
                if (abs(rlsig(iad)).gt.1.0e-5) then
                    lsignificant=abs(rdsig/rlsig(iad)).gt.1.0e-5
                else
                    lsignificant=abs(rdsig).gt.1.0e-5
                endif
                if (lsignificant) then
                if (lErrorReport) write (ILUME,4004) RDSIG,ALSLAB(IAD)
                if (lErrorReport) write (*,4004) RDSIG,ALSLAB(IAD)
                endif
              ENDIF
            RLSIG(IAD)=RLSIG(IAD)-RDSIG*rlsrelax
            RBFCORTEMP=-RDSIG*RLENGTH*rlsrelax
            ENDIF ! -------------------------------------------------------(6)
          ELSE  !! head below res. layer ----------------------------------(5)
            RDQ=RQMAX+RLSIG(IAD)
            IF (RDQ.GT.0.0) THEN ! more of streamflow can be infiltrated---(8)
              RDSIG=MIN(RDQ,RQS)
              IF (RDSIG.GT.0.0) THEN
                LCHANGE=.TRUE.
                if (abs(rlsig(iad)).gt.1.0e-5) then
                    lsignificant=abs(rdsig/rlsig(iad)).gt.1.0e-5
                else
                    lsignificant=abs(rdsig).gt.1.0e-5
                endif
                if (lsignificant) then
                if (lErrorReport) write (ILUME,4500) RDSIG,ALSLAB(IAD)
                if (lErrorReport) write (*,4500) RDSIG,ALSLAB(IAD)
                endif
              ENDIF
            ENDIF
            IF (RDQ.LT.0.0) THEN ! overinfiltrating------------------------(8)
              RDSIG=RDQ
              LCHANGE=.TRUE.
              if (abs(rlsig(iad)).gt.1.0e-5) then
                lsignificant=abs(rdsig/rlsig(iad)).gt.1.0e-5
              else
                lsignificant=abs(rdsig).gt.1.0e-5
              endif
              if (lsignificant) then
              if (lErrorReport) write (ILUME,4002) RDSIG,ALSLAB(IAD)
              if (lErrorReport) write (*,4002) RDSIG,ALSLAB(IAD)
              endif
            ENDIF  ! -------------------------------------------------------(8)
            RLSIG(IAD)=RLSIG(IAD)-RDSIG*rlsrelax
            RBFCORTEMP=-RDSIG*RLENGTH*rlsrelax
          ENDIF ! ---------------------------------------------------------(5)
        ELSE !!! currently not given, check if inf. rate is possible ------(4)
          IF (RLSIG(IAD).LT.0.0) THEN !! infiltrating, else no check ------(9)
            IF (.NOT.(RLSRES(IAD).EQ.0.0.AND.RHEAD.GT.RLSH(IAD))) THEN !---(10)
C           else set 'given' in previous iteration: don't touch
              RS=MIN(RQS,RQMAX)
              IF (RLSIG(IAD).LT.-RS) THEN ! limit infiltration ------------(11)
                RDSIG=RS+RLSIG(IAD)
                if (abs(rlsig(iad)).gt.1.0e-5) then
                  lsignificant=abs(rdsig/rlsig(iad)).gt.1.0e-5
                else
                  lsignificant=abs(rdsig).gt.1.0e-5
                endif
                if (lsignificant) then
                if (lErrorReport) write (ILUME,4003) RDSIG,ALSLAB(IAD)
                if (lErrorReport) write (*,4003) RDSIG,ALSLAB(IAD)
                endif
                LCHANGE=.TRUE.
                RBFCORTEMP=-RDSIG*RLENGTH*rlsrelax
                RLSIG(IAD)=rlsig(iad)-rdsig*rlsrelax
                LSGIV(IAD)=.TRUE.
              ELSE ! ------------------------------------------------------(11)
              ENDIF ! -----------------------------------------------------(11)
            ELSE ! --------------------------------------------------------(10)
            ENDIF ! -------------------------------------------------------(10)
          ELSE ! ----------------------------------------------------------(9)
          ENDIF ! ---------------------------------------------------------(9)
        ENDIF ! -----------------------------------------------------------(4)
      ENDIF ! -------------------------------------------------------------(2)
      RLSBF(IAD)=RLSBF(IAD)+RBFCORTEMP  ! correct baseflow in current element
      RBFCORRECT=RBFCORRECT+RBFCORTEMP  ! update correction term for next elements
      INEXT=KLSDN(IAD)
      IF (INEXT.EQ.0) GOTO 30 ! end of stream, next string
      IAD=INEXT
      GOTO 25
      ENDIF ! -------------------------------------------------------------(0)
  30  CONTINUE
      IF (LCHANGE) WRITE (ILUME,2000)
  45  RLSBFMX=-1.0E21
      RLSOFMX=-1.0E21
      DO 50 ILS=1,NLS
      RLSBFMX=MAX(RLSBFMX,RLSBF(ILS))
      RLSOFMX=MAX(RLSOFMX,RLSOF(ILS))
  50  CONTINUE
      RETURN
 1000 FORMAT (' Relaxation factor in use:',D14.7)
 2000 FORMAT (' Resolve groundwater flow.')
 3000 FORMAT (' *** Found no stream networks!')
 4001 FORMAT (' Found overinfiltration    ',G11.4,' at ',A16) 
 4002 FORMAT (' Found overinfiltration    ',G11.4,' at ',A16) 
 4003 FORMAT (' Found overinfiltration    ',G11.4,' at ',A16)   
 4004 FORMAT (' Found overinfiltration    ',G11.4,' at ',A16)    
 4010 FORMAT (' Found overinfiltration at ',A16) 
 4500 FORMAT (' Found underinfiltration   ',G11.4,' at ',A16)  
 4501 FORMAT (' Found underinfiltration   ',G11.4,' at ',A16)  
 5000 FORMAT (' Found negative streamflow ',G11.4,' at ',A16)
      END
C
C-----------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSOPENEND
C
C-----------------------------------------------------------------------------------------------------
C
C     Routine check and correct open ends in stream network
C
      IMPLICIT NONE
      INTEGER(4) ISTR,IADLAST
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
C
      IF (NLSTRING.EQ.0) RETURN
      DO 30 ISTR=1,NLSTRING 
      IADLAST=KLSTREND(ISTR)
        IF (KLSDN(IADLAST).EQ.IADLAST) THEN
        WRITE (ILUER,2000) ALSLAB(IADLAST)
        KLSDN(IADLAST)=0
        ENDIF
  30  CONTINUE        
      RETURN
 2000 FORMAT (' ***ERROR in line sink module: line sink ',A16,/,
     &' is not declared a stream end, but also not connected!',/,
     &' May be due to missing END statement following STREAM command',
     &' during input,',/,
     &' or nearby stream sections all have higher water elevations.',/,
     &' Line sink has been changed into a stream end!')
      END
c
C ---------------------------------------------------------------------------------
C
      subroutine lsread_relaxation_file
c
c     Routine reads the file "relax.tab" if it exists. If the table does not exists the
c     relaxation factor in the surface water routine is set to 1.0
c     If the table does exists it passes the specified iteration numbers and associated relaxation
c     factors to the surface water solution routine GENSTREAMFLOW
c
c     Format of the file "relax.tab":
c
c     * comment statement (* in first column)
c     * iteration #     relaxation factor
c          1              0.5
c          5              0.8
c          8              1.0
c     quit
c
c     The effect of this file is as follows: During the first 5 iterations a (under)relaxation factor
c     of 0.5 is used. During the next 3 iterations a relaxation factor of 0.8 is used, while after
c     that the relaxation factor is set to 1.0 (no relaxation).
c     NOTE: The relaxation factor should NEVER be set larger than 1.0
c

      implicit none
      INTEGER ilu,itemp,i,ivar,idum1,ierr
      LOGICAL lret,lnorelaxfile
      REAL(8) rdum1,rdum2,rvar
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
      INCLUDE 'main.inc'
c
      lnorelaxfile=.false.
      ilstablelength=0
      ilu=2
      afile='relax.tab'
      itemp=iluin ! temporarily store iluin and replace for this read
c      call opafil(ilu,-2,lret)
        iluin=ilu
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',ERR=20,IOSTAT=IERR)
  10    call inline
        if (aline(1).eq.'Q'.or.aline(1).eq.'q') GOTO 20 ! end of data in file
        idum1=ivar(1)
        rdum2=rvar(2)
        if (.not.lerror) then
          ilstablelength=ilstablelength+1
          niterarray(ilstablelength)=idum1
          if (rdum2.gt.1.0) then
            write (iluer,2000) rdum2
            rdum2=1.0
          end if
          relaxarray(ilstablelength)=rdum2
          GOTO 10
        end if
        if (lerror) then
          write (iluer,3000)
          lerror=.false.
        end if
  20    close (iluin)
        iluin=itemp ! restore iluin
        if (ierr.ne.0) then  ! failed to open file, assume it is not present
          write (ilume,1000)
          lnorelaxfile=.true.
        endif
        if (ilstablelength.eq.0.or.lnorelaxfile) THEN ! set relaxation factor to default
          ilstablelength=1
          niterarray(ilstablelength)=1
          relaxarray(ilstablelength)=1.0
        end if
c
C
c      write (iluer,1001) ilstablelength
c 1001 format (' lsread_relaxation_file: ilstablelength ',I5)
c      write (iluer,1002) (i,niterarray(i),relaxarray(i),
c     &i=1,ilstablelength)
c 1002 format (' lsread_relaxation_file:  i   niterarray  relaxarray'
c     &         ,/,4(26x,I3,5x,i3,5x,D14.7,/))
C
      if (.not.lnorelaxfile) then
        write (ilume,4000)
        write (ilume,5000) (niterarray(i),relaxarray(i),
     &                      i=1,ilstablelength)
      end if
 1000 format (' No file "relax.tab" found, relaxation factor for'
     &        ' streamflow is set to 1.0')
 2000 format(' ***ERROR in the file "relax.tab"!',/,
     &       'The relaxation factor in line ',i3,' is ',D14.7,/,
     &       'The relaxation factor has been reset to 1.0')
 3000 format (' ***ERROR in lsread_relaxation_file:',/,
     &        'Illegal syntax in the file "relax.tab"')
 4000 format(' Found file "relax.tab" The following table will be used',
     &       /,' during streamflow calculations:',/,
     &         ' iteration   relaxation factor')
 5000 format (10(3x,i3,5x,d14.7,/))
      end subroutine

