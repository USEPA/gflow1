C     Last change:  HMH  29 Oct 2012    1:12 pm
c     This file contains the following routines and functions:
c
C     LSERROR  calculates errors at line-sink control points and reports maximum error
c     subroutine lslakewaterbalance calculates and writes lake water balance to disk and concole
c
C
C-------------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSERROR (RERMAX)
C
C-------------------------------------------------------------------------------------------------------------
C
C     Routine calculates and reports the maximum error in the
C     boundary conditions specified at line sinks.
C     Note: Galleries without resistance are tested on constant head condition along gallery
C           Galleries with resistance are tested on matching strength condition at control points, like
C           a head-specified line-sink with resistance.
C
      IMPLICIT NONE
      INTEGER(4) I,K,IHED,IRES,IDR,IDH,IGH,IGR
      LOGICAL lpump
      REAL(8) RERMAX,RERRH,RERRS,RERDS,RERDH,RBASE,RFBASE,RS,
     &        RF,RFHEDP,RC,RFFF,RDUM,RW,RSIG,RSIGC,RSSS,
     &        RERGH,RERGR
      COMPLEX(8) CZ
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'tracom.inc'
      IF (NLSH.EQ.0) RETURN
      RERRH=0.0D0
      RERRS=0.0D0
      RERDS=0.0D0
      RERDH=0.0D0
      RERGH=0.0D0
      RERGR=0.0D0
      IHED=0
      IRES=0
      IDR=0
      IDH=0
      IGH=0
      IGR=0
      DO I=1,NLS
      RLSERROR(I)=0.0D0 ! set all errors to zero.
      END DO
      DO 10 I=1,NLSH
      K=KLSPTH(I)
      IF (LSGIV(K)) GOTO 10
      IF (LSDRAIN(K).AND.RLSIG(K).EQ.0.0D0) THEN ! active drain with zero sink density
        WRITE (ILUER,3000) ALSLAB(K)  ! issue warning, may cause large error
      END IF
      CZ=0.5*(CLSZS(K)+CLSZE(K))
      rbase=rfbase(cz)
      RS=RLSH(K)-rbase
      RF=RFHEDP(RLSPOT(K),CZ)-rbase
      RC=RLSRES(K)
      IF (RC.EQ.0.0D0) THEN
C ---------Head specified line sinks without resistance
        RFFF=MAX(ABS(RF+RS),0.02D0)      ! avoid division by zero
        RDUM=ABS(RF-RS)/(RFFF/2.0D0)*100
        RLSERROR(K)=RDUM
        IF (LSDRAIN(K)) THEN        !  drain
          RERDH=MAX(RDUM,RERDH)
          IDH=IDH+1
        ELSE IF (LSGALLERY(K)) THEN ! gallery, tested on equal heads
          if (rlsq(k).NE.-9999.0d0) then
            lpump=rlsq(k).ge.0.0d0  ! are we pumping or injecting?
          end if
          if (rlsh(k).lt.rlshmin(k).and.lpump) then ! line-sink pumped dry
           write (*,7000) alslab(k)
           if (.not.lucon) write (ilume,7000) alslab(k)
          end if
          IF (RLSQ(K).EQ.-9999.0D0) THEN
            RF=RLSH(K-1)-RFBASE(CZ)
            RFFF=MAX(ABS(RF+RS),0.02D0)      ! avoid division by zero
            RDUM=ABS(RF-RS)/(RFFF/2.0D0)*100
            RLSERROR(K)=RDUM
            RERGH=MAX(RDUM,RERGH)
          END IF
          IGH=IGH+1
        ELSE                        ! other
          RERRH=MAX(RDUM,RERRH)
          IHED=IHED+1
        END IF
      ELSE
C ---------Head specified line sinks with resistance      
        RW=RLSWID(K)
        RSIG=RLSIG(K)
        RSIGC=(RF-RS)*RW/RC
        RSSS=MAX(ABS(RSIG+RSIGC),10.D-10) ! avoid division by zero
        RDUM=ABS((RSIG-RSIGC)/(RSSS/2.0D0))*100
        RLSERROR(K)=RDUM
        IF (LSDRAIN(K)) THEN          ! drains
          RERDS=MAX(RDUM,RERDS)
          IDR=IDR+1
        ELSE IF (LSGALLERY(K)) THEN   ! galleries, tested on matching strength
          RERGR=MAX(RDUM,RERGR)
          IGR=IGR+1
        ELSE                          ! other
          RERRS=MAX(RDUM,RERRS)
          IRES=IRES+1
        END IF
      ENDIF
  10  CONTINUE
      IF (IHED.NE.0) then
        if (.not.lucon) WRITE (ILUME,1000) IHED,RERRH
        write (*,1000) ihed,rerrh
        if (rerrh.ge.rconverge_linesinks) lquit=.FALSE. ! do not abort iterations yet
      endif
      IF (IRES.NE.0) then
        if (.not.lucon) WRITE (ILUME,2000) IRES,RERRS
        write (*,2000) ires,rerrs
        if (rerrs.ge.rconverge_linesinks_resistance) lquit=.FALSE. ! do not abort iterations yet
      endif
      IF (IDH.NE.0) THEN
        IF (.NOT.LUCON) WRITE (ILUME,4000) IDH,RERDH
        WRITE (*,4000) IDH,RERDH
        if (rerdh.ge.rconverge_linesinks) lquit=.FALSE. ! do not abort iterations yet
      END IF
      IF (IDR.NE.0) THEN
        IF (.NOT.LUCON) WRITE (ILUME,5000) IDR,RERDS
        WRITE (*,5000) IDR,RERDS
        if (rerds.ge.rconverge_linesinks_resistance) lquit=.FALSE. ! do not abort iterations yet
      END IF
      IF (IGH.NE.0) THEN
        IF (.NOT.LUCON) WRITE (ILUME,6000) IGH,RERGH
        WRITE (*,6000) IGH,RERGH
        if (rergh.ge.rconverge_linesinks) lquit=.FALSE. ! do not abort iterations yet
      END IF
      IF (IGR.NE.0) THEN
        IF (.NOT.LUCON) WRITE (ILUME,6500) IGR,RERGR
        WRITE (*,6500) IGR,RERGR
        if (rergr.ge.rconverge_linesinks_resistance) lquit=.FALSE. ! do not abort iterations yet
      END IF
      RERMAX=MAX(RERMAX,RERRH)
      RERMAX=MAX(RERMAX,RERRS)
      RERMAX=MAX(RERMAX,RERDH)
      RERMAX=MAX(RERMAX,RERDS)
      RERMAX=MAX(RERMAX,RERGH)
      RERMAX=MAX(RERMAX,RERGR)
      RETURN
 1000 FORMAT (' ',I4,' line sinks without resistance:     max. error=',
     &        E11.4,' %')
 2000 FORMAT (' ',I4,' line sinks with resistance:        max. error=',
     &        E11.4,' %')
 4000 FORMAT (' ',I4,' drain elements without resistance: max. error=',
     &        E11.4,' %')
 5000 FORMAT (' ',I4,' drain elements with resistance:    max. error=',
     &        E11.4,' %')
 6000 FORMAT (' ',I4,' gallery elements no resistance:    max. error=',
     &        E11.4,' %')
 6500 FORMAT (' ',I4,' gallery elements with resistance:  max. error=',
     &        E11.4,' %')
 7000 format (' ***WARNING: Line-sink ',a16,' in gallery pumped dry!')
 3000 FORMAT (' ***WARNING in LSERROR: drain with label ',a16,/,
     &'is active, but has a zero sink density. This may result in a',/,
     &'relatively large error at the drain (line sink) center',/)
      END
c
c ----------------------------------------------------------------------------------------
c
      subroutine lslakewaterbalance(lakereportin,lErrorReport)
c
c    Routine writes the water balance data for lake features to a file.
c    lakereport = .TRUE.  write to the file "lakes.out" and "basename.lks"
c    lErrorReport = .TRUE.  write to message.log and console
c
      implicit none
      INTEGER i,niter,istr,iad,iend,inext,ilu,ierr,nsolOut
      LOGICAL lakereport,lakereportin,lErrorReport,l1,l2,l3,l4,l5
      REAL(8) rh,rq,rqo,rqf,rflstable,rarea,rqep,rlen,
     &        rsin,rsout,rin,rout,rdeltaq,rdeltaqpercent
      CHARACTER(8) aBasenameOut
      CHARACTER(16) aDateTimeOut
      INCLUDE 'lscom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
      INCLUDE 'tracom.inc'
c

c
      lakereport=lakereportin
      if (lErrorReport.or.lakereport) then ! something to write, go ahead...
      if (.not.lakereport) GOTO 1 ! do not open file if no lakereport
      ilu=2   ! use logical unit 2
      OPEN (UNIT=ILU,ERR=50,FILE="lakes.out",STATUS='UNKNOWN',
     &      IOSTAT=IERR)                                                ! open file as "unknown"
      GOTO 1
  50  write (ilume,100) ierr                                           ! could not open file; no lake report
      write (*,100) ierr
      lakereport=.false.
      GOTO 1
c
   1   lakedone=.true.
      do istr=1,nlstring
       iad=klstring(istr)
       if (lslake(iad)) THEN ! found lake
        if (lErrorReport) write (ilume,500)    ! write header for lake water balance
        if (lErrorReport) write (*,500)
        rlen=ABS(clsze(iad)-clszs(iad))
        rsin=0.0d0
        rsout=0.0d0
        if (rlsig(iad).gt.0.0d0) rsin=rlsig(iad)*rlen    ! collect all groundwater inflow into the lake
        if (rlsig(iad).lt.0.0d0) rsout=-rlsig(iad)*rlen  ! collect all groundwater outflows from the lake
        rq=-rlsq(iad)         ! collect all flows into outlet streams
        rh=rlsh(iad)         ! get lake stage
        if (lErrorReport) write (ilume,1000) alslab(iad)              ! write first label
        if (lErrorReport) write (*,1000) alslab(iad)
        if (lakereport) write (ilu,1100)  alslab(iad)
  5     inext=klsdn(iad)
        if (inext.gt.0) then ! next line-sink
         rlen=ABS(clsze(inext)-clszs(inext))
         if (rlsig(inext).gt.0.0d0) rsin=rsin+rlsig(inext)*rlen    ! collect all groundwater inflow into the lake
         if (rlsig(inext).lt.0.0d0) rsout=rsout-rlsig(inext)*rlen  ! collect all groundwater outflows from the lake
         rq=rq-rlsq(inext)      ! collect all flows into outlet streams
         rqo=rlsof(inext)       ! get overland flow in case this is the last line-sink
         rqf=rlsbf(inext)-rsin+rsout ! get stream inflow in case this is the last line-sink
         iad=inext
         GOTO 5
        else                 ! at end of string
         if (lErrorReport) write (ilume,2000) alslab(iad)              ! write last label
         if (lErrorReport) write (*,2000) alslab(iad)
         if (lakereport) write (ilu,2100)  alslab(iad)
         if (lErrorReport) write (ilume,3000) rh                       ! write lake stage
         if (lErrorReport) write (*,3000) rh
         if (lakereport) write (ilu,3100) rh
         rarea=rflstable(rh,istr)
         if (lErrorReport) write (ilume,4000) rarea                    ! write lake area
         if (lErrorReport) write (*,4000) rarea
         if (lakereport) write (ilu,4100) rarea
         rqep=rarea*rlsevap(istr) ! (evaporation - precipitation) times lake area
         if (lErrorReport) write (ilume,4500) rqep                     ! write net evaporation (not to report)
         if (lErrorReport) write (*,4500) rqep
         if (lErrorReport) write (ilume,5000) rsin                     ! write groundwater inflow
         if (lErrorReport) write (*,5000) rsin
         if (lakereport) write (ilu,5100) rsin
         if (lErrorReport) write (ilume,6000) rsout                    ! write groundwater outflow
         if (lErrorReport) write (*,6000) rsout
         if (lakereport) write (ilu,6100) rsout
         if (lErrorReport) write (ilume,7000) rqf                      ! write stream inflow
         if (lErrorReport) write (*,7000) rqf
         if (lakereport) write (ilu,7100) rqf
         if (lErrorReport) write (ilume,8000) rq                       ! write stream outflow
         if (lErrorReport) write (*,8000) rq
         if (lakereport) write (ilu,8100) rq
         if (lErrorReport) write (ilume,9000) rqo                      ! write overland flow
         if (lErrorReport) write (*,9000) rqo
         if (lakereport) write (ilu,9100) rqo
         rin=rsin+rqf+rqo
         rout=rsout+rq+rqep
         rdeltaq=ABS(rin-rout)   ! use absolute value for water balance deficiency
         rdeltaqpercent=ABS(rdeltaq/rin*100)
         if (lErrorReport) write (ilume,9400) rdeltaq                  ! write water balance deficiency
         if (lErrorReport) write (*,9400) rdeltaq
         if (lakereport) write (ilu,9410) rdeltaq
         if (lErrorReport) write (ilume,9600) rdeltaqpercent           ! write percent water balance deficiency
         if (lErrorReport) write (*,9600) rdeltaqpercent
         if (lakereport) write (ilu,9610) rdeltaqpercent
         if (rconverge_lake_waterbalance.lt.rdeltaqpercent) then
           lakedone=.false.
         endif
c         write (ilume,1001) rconverge_lake_waterbalance,rdeltaqpercent,
c     &                      lakedone
c 1001 format (' lslakewaterbalance: rconverge_lake_waterbalance,'
c     &        'rdeltaqpercent,lakedone ',2(d14.7),l3)
        end if
       end if
       end do
c
      if (lakereport) close (ilu)              ! close lake report file
c
      if (lakereport) then
c
c
c ------------OK, now write to the file "basename.lks"
c
      call GetMainData (l1,l2,l3,l4,l5,
     &           aBasenameOut,aDateTimeOut,nsolOut)     ! get the basename
c
      ilu=2   ! use logical unit 2
      afile=TRIM(aBasenameOut) // ".lks"
      OPEN (UNIT=ILU,ERR=55,FILE=afile,STATUS='UNKNOWN',
     &      IOSTAT=IERR)                                                ! open file as "unknown"
      GOTO 2
  55  write (ilume,105) ierr                                           ! could not open file; no lake report
      write (*,105) ierr
      lakereport=.false.
      return     ! get out of this routine because cannot create/open file
c
   2  continue
      do istr=1,nlstring
       iad=klstring(istr)
       if (lslake(iad)) THEN ! found lake
        write (ilu,500)  ! write header for lake water balance
        rlen=ABS(clsze(iad)-clszs(iad))
        rsin=0.0d0
        rsout=0.0d0
        if (rlsig(iad).gt.0.0d0) rsin=rlsig(iad)*rlen    ! collect all groundwater inflow into the lake
        if (rlsig(iad).lt.0.0d0) rsout=-rlsig(iad)*rlen  ! collect all groundwater outflows from the lake
        rq=-rlsq(iad)         ! collect all flows into outlet streams
        rh=rlsh(iad)         ! get lake stage
        write (ilu,1000)  alslab(iad)
  6     inext=klsdn(iad)
        if (inext.gt.0) then ! next line-sink
         rlen=ABS(clsze(inext)-clszs(inext))
         if (rlsig(inext).gt.0.0d0) rsin=rsin+rlsig(inext)*rlen    ! collect all groundwater inflow into the lake
         if (rlsig(inext).lt.0.0d0) rsout=rsout-rlsig(inext)*rlen  ! collect all groundwater outflows from the lake
         rq=rq-rlsq(inext)      ! collect all flows into outlet streams
         rqo=rlsof(inext)       ! get overland flow in case this is the last line-sink
         rqf=rlsbf(inext)-rsin+rsout ! get stream inflow in case this is the last line-sink
         iad=inext
         GOTO 6
        else                 ! at end of string
         write (ilu,2000)  alslab(iad)
         write (ilu,3000) rh
         rarea=rflstable(rh,istr)
         write (ilu,4000) rarea
         rqep=rarea*rlsevap(istr) ! (evaporation - precipitation) times lake area
         write (ilu,4500) rqep
         write (ilu,5000) rsin
         write (ilu,6000) rsout
         write (ilu,7000) rqf
         write (ilu,8000) rq
         write (ilu,9000) rqo
         rin=rsin+rqf+rqo
         rout=rsout+rq+rqep
         rdeltaq=ABS(rin-rout)   ! use absolute value for water balance deficiency
         rdeltaqpercent=ABS(rdeltaq/rin*100)
         write (ilu,9400) rdeltaq
         write (ilu,9600) rdeltaqpercent
        end if
       end if
       end do
c
      close (ilu)              ! close lake report file
      endif
c
      endif

      return
  100 format (' ***ERROR opening file "lakes.out", I/O error: ',I3)
  105 format (' ***ERROR opening file "<basename>.lks", I/O error: ',I3)
  500 format (' ',/,'           Lake water balance:')
 1000 format ('first GFLOW label: ',A16)
 1100 format (A16)
 2000 format (' last GFLOW label: ',A16)
 2100 format (A16)
 3000 format ('       lake stage: ',D17.10)
 3100 format (D17.10)
 4000 format ('        lake area: ',D17.10)
 4100 format (D17.10)
 4500 format ('  evap. - precip.: ',D17.10) ! not in "lakes.out", but in "basename.lks"
 5000 format ('   groundwater in: ',D17.10)
 5100 format (D17.10)
 6000 format ('  groundwater out: ',D17.10)
 6100 format (D17.10)
 7000 format ('   stream flow in: ',D17.10)
 7100 format (D17.10)
 8000 format ('  stream flow out: ',D17.10)
 8100 format (D17.10)
 9000 format (' overland flow in: ',D17.10)
 9100 format (D17.10)
 9400 format ('water balance deficiency: ',D17.10)
 9410 format (D17.10)
 9600 format ('water balance deficiency: ',D17.10,' %')
 9610 format (D17.10)
      end
