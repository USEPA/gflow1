C     Last change:  HMH  21 Oct 2014    3:38 pm
c       This file contains the following routines and functions:
c
C	BLOCK DATA LSDAT    initializes line-sink common blocks
C	SUBROUTINE LSIN     handles input of line-sinks as defined in the .dat file
C	SUBROUTINE RELSPT   input service routine for LSIN
c
c -------------------
c lsmat.for contains:
c -------------------
c
C	SUBROUTINE LSCZC     generates control point vector for line-sinks
C	SUBROUTINE LSMAT     generates matrix coefficients for line-sinks
C	SUBROUTINE LSKNO     generates known vector for line-sinks
C	SUBROUTINE LSSUB     substitutes solution vector into strength parameters for line-sinks
C	SUBROUTINE LSUPDATE  stores the calculated heads at the centers of line-sinks
c
c -------------------
c lsfun.for contains:
c -------------------
c
C	SUBROUTINE LSQI         calcuates the discharge vector for all line-sinks
C	COMPLEX FUNCTION COMLS  calculates the complex potential for lins-sink between CZ1 and CZ2  c                               with sink density 1
C	COMPLEX FUNCTION CFLSG  complex potential for all line-sinks with given sink density
C	COMPLEX FUNCTION CFLSU  complex potential for all head specified line-sinks
C	COMPLEX FUNCTION CDLS   derivative of COMLS
C	COMPLEX FUNCTION CDLSG  derivative of CFLSG
C	COMPLEX FUNCTION CDLSU  derivative of CFLSU
C	REAL FUNCTION RFLSPT    real discharge potential due to all line-sinks
c
c ---------------------
c lsnflow.for contains:
c ---------------------
c
C	REAL FUNCTION RFNFLS      calculates the flow between two points due to all line-sinks
C	REAL FUNCTION RFNFLSCO    coefficient functions for one line-sink with sink density 1
C	REAL FUNCTION RFBCLSU     calculates the jump in PSI for a branch cut due to line-sink I
c
c ------------------
c lsio.for contains:
c ------------------
c
C	subroutine lsextracthead    writes data for head specified line sinks and drains to .xtr file
C	SUBROUTINE LSDATIO          writes all line sink data to a ".dat" file
C	SUBROUTINE LSIO             reads or writes contents of LSCOM.INC from or to .sol file
C	SUBROUTINE LSNKOUT          writes contents of LSCOM formatted to ILUOUT for debugging purposes
c
c ---------------------
c lscheck.for contains:
c ---------------------
c
C     LSERROR  calculates errors at line-sink control points and reports maximum error
c
c --------------------
c lsnear.for contains:
c --------------------
c
C	SUBROUTINE LSNEAR   handles particle tracing near line-sinks
c
c --------------------
c lsserv.for contains:
c --------------------
c
C	SUBROUTINE NOCONJUNCTIVE   sets flag for conjuctive solutions false
c       SUBROUTINE LSPREP          calculate width, CLSCONST, and cancel recharge on lakes
c       real function rflsnearflow returns the streamflow in the nearest line-sink
c       subroutine ls_setwidth     routine sets the representative width parameter for resistance line-sinks
c
c ----------------------
c lsstream.for contains:
c ----------------------
c
C	SUBROUTINE LSNETWORK       generates a stream network
C	SUBROUTINE GENSTREAMFLOW   generates streamflow and eliminates some line-sinks from solution
c       SUBROUTINE LSOPENEND       checks for stream reaches that are not connected in the network
c       SUBROUTINE lsread_relaxation_file  reads a table with relaxation factors versus iterations for
c                                          streamflow calcutations.c
c ---------------------
c lsextra.for contains:
c ---------------------
c
C	SUBROUTINE LSPLOT
C	SUBROUTINE PLOTLINESINK
C	SUBROUTINE PLOTLSNODES
C	SUBROUTINE PLOTFLOW
C	SUBROUTINE LSCUR
C	SUBROUTINE LSCHEK
c       REAL FUNCTION RFBCLS
c       REAL FUNCTION RFBCLSS
c       LSNETPLOT
c       DOWNSTREAMPLOT
C
c
C-------------------------------------------------------------------------------------------------------
C
	BLOCK DATA LSDAT
C
C-------------------------------------------------------------------------------------------------------
C
      INCLUDE 'LSCOM.INC'
      DATA RO2PI / .159154943D0/
      DATA NLS,RLSIG/ 0,NLSMAX*.0/
      DATA LSGIV,LSDRAIN /NLSMAX*.FALSE.,NLSMAX*.FALSE./
      DATA LSGALLERY /NLSMAX*.FALSE./
      DATA NLSIG,NLSH /0,0/
      DATA RLSW0,RLSC0,RLSD0 /0.0D0,0.0D0,0.0D0/
      DATA LSHIRE,LSHIPE /2*.FALSE./
      DATA LSEND,LSBASE,LSFIRST /2*.FALSE.,.TRUE./
      DATA LSINLET,LSOUTLET /NLSMAX*.FALSE.,NLSMAX*.FALSE./
      DATA LSCONNECT /NLSMAX*.FALSE./
      DATA LSLAKE /NLSMAX*.FALSE./
      DATA LSINLETSTRING,LSOUTLETSTRING /2*.FALSE./
      DATA LSLAKESTRING /.FALSE./
      DATA ILSPLOT /0/
      DATA RLSQ /NLSMAX*0.0D0/
      DATA NLSTAB,NLSTABLENGTH /NLSTMX*0,NLSTMX*0/
      DATA RLSH1,RLSH2 /NLSTMX*0.0D0,NLSTMX*0.0D0/
      DATA RLSQ1,RLSQ2 /NLSTMX*0.0D0,NLSTMX*0.0D0/
      DATA RLSH0,RLSEVAP /NLSTMX*0.0D0,NLSTMX*0.0D0/
      DATA NLAKEITERATIONS /1/
      DATA ilsbound0 /0/
      END
C
C-------------------------------------------------------------------------------------------------------
C
      SUBROUTINE LSIN (LSOL,RA,IRA,RSCR)
C
C-------------------------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'LSCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'GRID.INC'
      INTEGER(4) IRA,ILU,INEXT,IAD,ISPRSTR,ICODE,JUMP,I,IWMIN,IWMAX,
     &           NSPRFILES,ISTR0,idum,ivar
      INTEGER :: rvalues(8)
      LOGICAL LSOL,LINSTREAM,LINLINESINKS,LBAD,LDRAIN,LNET,LRET,LGALLERY
      REAL(8) RA,RSCR,RDUM,RVAR,RHED,RLEN,RHW,ROF,RL,RWMIN,RWMAX,RFHEAD,
     &        rqgallery,RH0,RH2
      COMPLEX(8) CZ
      CHARACTER(1) AWORD(68),APAR1(16),APAR2(21),APAR3(30)
      CHARACTER(8) ADATE
      CHARACTER(11) ATIME
      CHARACTER(5) azone
      CHARACTER(80) ASPRFILE(20)
      INCLUDE 'MATCH.INC'
      DIMENSION RA(IRA,*),RSCR(*),ISPRSTR(20)
      DATA AWORD/ 'D','I','S','C',' ',
     .            'D','R','A','I',' ',
     .            'G','A','L','L',' ',
     .            'H','E','A','D',' ',
     .            'S','T','R','E',' ',
     .            '?',' ',
     .            'Q','U','I','T',' ',
     .            'W','I','D','T',' ',
     .            'R','E','S','I',' ',
     .            'D','E','P','T',' ',
     .            'P','U','M','P',' ',
     .            'H','I','G','H',' ',
     .            'C','U','R','S',' ',
     .            'P','L','O','T',' ',
     .            ATERM/
      DATA APAR1/ 'N','O','N','E',' ',
     .            'R','E','C','H',' ',
     .            'P','E','R','C',' ',
     .            ATERM/
      DATA APAR2/ 'L','A','Y','O',' ',
     .            'B','A','S','E',' ',
     .            'O','V','E','R',' ',
     .            'S','T','R','E',' ',
     .            ATERM/
      DATA APAR3/ 'E','N','D',' ',
     .            'O','V','E','R',' ',
     .            'I','N','L','E',' ',
     .            'O','U','T','L',' ',
     .            'L','A','K','E',' ',
     .            'L','E','V','E',' ',
     .            ATERM/
      LERROR=.FALSE.
      LMISS=.FALSE.
      ICODE=0
      LINSTREAM=.FALSE.
      LINLINESINKS=.FALSE.
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
 10   IF (LERROR.OR.LMISS) WRITE(ILUER,9001) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE. 
      IF (ICODE.NE.0) THEN
      JUMP=IABS(ICODE)
      GOTO 15
      ENDIF
      if (lucon) then
      WRITE (ILUME,8000) NLSMAX,NLSTMX
      WRITE (ILUME,8001) RLSW0,RLSC0,RLSD0
      IF (.NOT.LSHIRE.AND..NOT.LSHIPE) WRITE (ILUME,8002)
      IF (LSHIRE.AND.LSHIPE) WRITE (ILUME,8003)
      IF (LSHIPE.AND..NOT.LSHIRE) WRITE (ILUME,8004)
      IF (ILSPLOT.EQ.0) WRITE (ILUME,8005)
      IF (ILSPLOT.EQ.1) WRITE (ILUME,8006)
      IF (ILSPLOT.EQ.2) WRITE (ILUME,8007)
      IF (ILSPLOT.EQ.3) WRITE (ILUME,8008)                  
      WRITE (ILUME,8100)
      endif
      CALL INLINE
 11   CALL MATCH(AWORD,1,JUMP,LBAD)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      IF(.NOT.LBAD) GOTO 15
      GOTO (10,13), JUMP
 13   WRITE(ILUER,9002) ALINE2
      LERROR=.FALSE.
      GOTO 10
 15   GOTO(120,150,160,170,190,250,300,400,500,600,650,700,800,900),JUMP
C
C   Input discharge specified linesinks
C
 120  if (lucon) WRITE(ILUME,9005)
      CALL INLINE
      RDUM=RVAR(1)
      icode=1
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 11
      ENDIF
      IF (LMISS) GOTO 10
      LSOL=.FALSE.
      CALL RELSPT(ICODE,ROF,RHW,RQGALLERY,RH0,RH2)
      LINLINESINKS=.TRUE.
      GOTO 120
C 
C   Input head specified linesinks
C
 150  LDRAIN=.TRUE. ! line sink is a drain feature
      LGALLERY=.FALSE.
      GOTO 180
 160  LGALLERY=.TRUE.  ! line-sink is a gallery feature
      LDRAIN=.FALSE.
      GOTO 180
 170  LDRAIN=.FALSE.
      LGALLERY=.FALSE.
 180  if (lucon) WRITE(ILUME,9006)
      CALL INLINE
      RDUM=RVAR(1)
      IF (LDRAIN) THEN
      ICODE=2
      ELSE IF (LGALLERY) THEN
      ICODE=3
      ELSE
      ICODE=4
      END IF
      IF (LERROR) THEN
        LERROR=.FALSE.
        LMISS=.FALSE.
        GOTO 11
      ENDIF
      IF (LMISS) GOTO 10
      LSOL=.FALSE.
      CALL RELSPT(ICODE,ROF,RHW,RQGALLERY,RH0,RH2)
      IF (LGALLERY) RQGALLERY=-9999.0D0
      LINLINESINKS=.TRUE.
      GOTO 180
C
C     Stream sections for stream networks
C
  190 IF (ICODE.EQ.5) GOTO 195
      LINSTREAM=.TRUE.
      LINLINESINKS=.TRUE.
      if (lucon) WRITE (ILUME,9040)           ! first string entry
      LSEND=.FALSE.
      LSINLETSTRING=.FALSE.
      LSOUTLETSTRING=.FALSE.
      LSLAKESTRING=.FALSE.
      ROF=0.0
      RHW=RVAR(2)      ! head water inflow OR lake evapotranspiration
      IF (LERROR) THEN
      LERROR=.FALSE.
      RHW=0.0
      ENDIF
  195 if (lucon) WRITE (ILUME,9050)
      ICODE=5
      CALL INLINE
  200 if (lucon) print*,' '
      RDUM=RVAR(1)
      IF (LERROR) THEN
      LERROR=.FALSE.
      CALL MATCH (APAR3,1,JUMP,LBAD)
      IF (LBAD) GOTO 11
      IF (JUMP.EQ.1) THEN     ! end of stream command
      LSEND=.TRUE.
      GOTO 195
      ENDIF
      IF (JUMP.EQ.2) THEN     ! overlandflow command
      ROF=RVAR(2)
      IF (LERROR) GOTO 10
      GOTO 195
      ENDIF
      IF (JUMP.EQ.3) THEN     ! end stream forms lake inlet (high-k lakes only)
      LSEND=.TRUE.
      LSINLETSTRING=.TRUE.
      GOTO 195
      END IF
      IF (JUMP.EQ.4) THEN     ! stream forms lake outlet
      CALL GETFN(2) ! read table filename into AFILE
      ATEMP=AFILE
      IF (LERROR) GOTO 10
      LSOUTLETSTRING=.TRUE.
      GOTO 195
      END IF
      IF (JUMP.EQ.5) THEN     ! stream forms lake
      CALL GETFN(2) ! read table filename into AFILE
      ATEMP=AFILE
      IF (LERROR) GOTO 10
      LSLAKESTRING=.TRUE.
      GOTO 195
      END IF
      IF (JUMP.EQ.6) THEN   ! read lake levels
      RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RH0=RDUM       ! lake bottom elevation
      RDUM=RVAR(3)
      IF (LERROR) GOTO 10
      RH2=RDUM       ! second head estimate
      GOTO 195
      END IF
      ENDIF
      IF (LMISS) GOTO 10
      LSOL=.FALSE.
      RDUM=RVAR(3)
      IF (LMISS) THEN    ! no parameters: end of string
      LERROR=.FALSE.
      LMISS=.FALSE.
      ICODE=-5
      ELSE
      ICODE=5
      ENDIF
      CALL RELSPT (ICODE,ROF,RHW,RQGALLERY,RH0,RH2)
      IF (ICODE.EQ.-5) THEN  ! already read next command in RELSPT
        DO 202 I=1,80
        ALINE(I)=ALINE2(I)
  202   CONTINUE
        CALL TIDY
        GOTO 11
      ENDIF
      IF (ICODE.EQ.0) GOTO 10   ! abort entering string.
      GOTO 200
C
C   Help
C
 250  AFILE='                '
      AFILE(1:9)='LSHLP.HLP'
C      CALL HELP                      ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
c   quit command
C
C   Return, after checking for reasonable width and after networking
C   the stream sections entered by use of the stream command.
C      
 300  lmiss=.false.
      lerror=.false.
      IF (LINLINESINKS) THEN
      DO 310 I=1,NLS
      RL=ABS(CLSZE(I)-CLSZS(I))
      RWMIN=0.0001*RL
      IWMIN=INT(10.0*RWMIN+0.5)
      RWMIN=REAL(IWMIN)/10.0
      RWMIN=MAX(RWMIN,0.1)
      RWMAX=0.3*RL
      IWMAX=INT(10.0*RWMAX-0.5)
      RWMAX=REAL(IWMAX)/10.0
      RWMAX=MAX(RWMAX,0.1)
c      IF (RLSWID(I).GT.RWMAX)
c     & WRITE (ILUER,3001) ALSLAB(I),RLSWID(I),RWMAX
      IF (RLSRES(I).EQ.0.0.AND.RLSWID(I).LT.0.0) THEN
        WRITE (ILUER,3003) ALSLAB(I),RLSWID(I)
        RLSWID(I)=0.0
        GOTO 310
      ENDIF
      IF (RLSRES(I).NE.0.0.AND.RLSWID(I).LE.0.0) THEN
        WRITE (ILUER,3004) ALSLAB(I)
        RLSRES(I)=0.0
        RLSWID(I)=0.0
        GOTO 310
      ENDIF
c      IF (RLSRES(I).NE.0.0.AND.RLSWID(I).LT.RWMIN)
c     & WRITE (ILUER,3002) ALSLAB(I),RLSWID(I),RWMIN
 310  CONTINUE
      ENDIF
      IF (LINSTREAM) CALL LSNETWORK
      RETURN
C
C     set default width
C
 400  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RLSW0=RDUM
      idum=ivar(3)
      if (lerror.or.lmiss) then
        ilsbound0=0
        lerror=.false.
        lmiss=.false.
      else
        ilsbound0=idum
      end if
      GOTO 10
C
C     set default resistance
C
 500  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RLSC0=RDUM
      GOTO 10
C 
C     set default depth (distance between water level in stream
C     and bottom of resistance layer
C
 600  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RLSD0=RDUM
      GOTO 10
C
C     set pumping rate for a gallery
C
 650  RDUM=RVAR(2)
      IF (LERROR) THEN
        WRITE (ILUER,9003)
        RQGALLERY=0.0D0
      ELSE
        RQGALLERY=RDUM
      END IF
      GOTO 10
C
C     highlight recharging or percolating line sinks
C
  700 CALL MATCH (APAR1,2,JUMP,LBAD)
      ICODE=0
      LERROR=LBAD
      IF (LERROR) GOTO 10
      IF (JUMP.EQ.1) THEN
      LSHIRE=.FALSE.
      LSHIPE=.FALSE.
      ENDIF
      IF (JUMP.EQ.2) THEN
      LSHIRE=.TRUE.
      LSHIPE=.TRUE.
      ENDIF
      IF (JUMP.EQ.3) THEN
      LSHIRE=.FALSE.
      LSHIPE=.TRUE.
      ENDIF
      GOTO 10
C
C     Cursor      
C
 800  IF (.NOT.LSOL) THEN
C        CALL TONE            ! NOT AVAILABLE IN BATCH MODE
        WRITE (ILUER,8200)
        LERROR=.TRUE.
      ENDIF
      IF (.NOT.LSBASE) THEN
        LNET=.FALSE.
        DO 810 I=1,NLS
        IF (KLSDN(I).NE.I) LNET=.TRUE.
 810    CONTINUE 
        IF (LNET) THEN
C          CALL TONE           ! NOT AVAILABLE IN BATCH MODE
          WRITE (ILUER,8210)
          LERROR=.TRUE.
        ENDIF
      ENDIF
      CALL LSOPENEND
C      CALL LSCUR (RA,IRA,RSCR,LSOL,NSPRFILES,ISPRSTR,ASPRFILE) ! NOT AVAILABLE IN BATCH MODE
      NSPRFILES=0 ! in batch mode force spreadsheet files to zero
      do i=1,20
      ISPRSTR(i)=0.0  ! no string addresses for spreadsheet files in batch mode
      ASPRFILE(i)=''  ! no filenames when in batch mode
      end do
      IF (.NOT.LSOL) LINLINESINKS=.TRUE.
      IF (NSPRFILES.GT.0) GOTO 950
      GOTO 10
C
C     Plot layout/baseflow/overlandflow/streamflow
C
 900  CALL MATCH (APAR2,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      ILSPLOT=JUMP-1
      GOTO 10
C
C     Spreadsheet file   (only accessible in interactive mode with data from LSCUR routine)
C                        may be accessed by adding a command and assigning ISPRSTR and ASPRFILE arrays
C      
 950  IF (NSPRFILES.GT.0) THEN
      ILU=2
      DO 960 I=1,NSPRFILES
      ISTR0=ISPRSTR(I)
      AFILE=ASPRFILE(I)
      CALL CRAFIL (ILU,-4,LRET,'.SPS')
      IF (LRET) GOTO 10
      call date_and_time (adate,atime,azone,rvalues)
      WRITE (ILU,2200) ADATE,ATIME,ATITLE
 2200 FORMAT (' *Date: ',A8,2X, 'Time: ',A11,/,
     &        ' *Title: ',A16,' Stream reach data written by GFLOW1.',/)
      WRITE (ILU,9760)
      IAD=KLSTRING(ISTR0)
  957 RLEN=ABS(CLSZS(IAD)-CLSZE(IAD))
      CZ=0.5*(CLSZS(IAD)+CLSZE(IAD))
      RHED=RFHEAD(CZ)
      WRITE (ILU,9761) RLSBF(IAD),RLSOF(IAD),RHED,RLSH(IAD),RLEN,
     &RLSIG(IAD),RLSRES(IAD),RLSWID(IAD),RLSDEP(IAD),ALSLAB(IAD)
      INEXT=KLSDN(IAD)
      IF (INEXT.EQ.0) GOTO 958 ! end of stream
      IF (INEXT.EQ.IAD) THEN
      WRITE (ILUER,9600) ALSLAB(IAD)
      CLOSE (ILU)
      GOTO 960
      ENDIF
      IAD=INEXT
      GOTO 957
  958 CLOSE (ILU)
      if (lucon) WRITE (ILUME,9765)  AFILE
  960 CONTINUE
      if (lucon) WRITE (ILUME,9800)
      CALL INLINE
c      CALL CLEARSCREEN      ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
      ENDIF      
C 
 3001 FORMAT (' ***WARNING in line sink module:',/,
     &' width of linesink ',A16,' is ',G11.4,/,
     &' recommended maximum width is ',G11.4,/,
     &' (For wide streams in nearfield consider parallel linesinks.)',/)
 3002 FORMAT (' ***WARNING in line sink module:',/,
     &' width of linesink ',A16,' is ',G11.4,
     &' Width below ',G11.4,' may lead to unstable solution.',/,
     &' (Width = zero is O.K. for line sinks without resistance.)',/)
 3003 FORMAT (' ***ERROR in line sink module:',/,
     &' width of linesink ',A16,' is ',G11.4,
     &' set to ZERO!',/)     
 3004 FORMAT (' ***ERROR in line sink module:',/,
     &' zero width and non-zero resistance for linesink ',A16,/,
     &' Both resistance and width are set to ZERO!',/)     
 8000 FORMAT('                   -------- LINE SINK module --------'/
     &       ' Maximum number of line sinks:',I4,/
     &       ' Maximum number of streams:',I4,/
     &       ' Available commands:',/,
     &       ' <F1> = Help',/,' HEAD',/,' DISCHARGE',/,
     &       ' STREAM [headwater inflow L^3/T]')
 8001 FORMAT (' WIDTH ',G11.4,/,' RESISTANCE ',G11.4,/,
     &' DEPTH ',G11.4,/,' CURSOR')
 8002 FORMAT (' HIGHLIGHT NONE (recharging/percolating)')
 8003 FORMAT (' HIGHLIGHT RECHARGING (percolating/none)')
 8004 FORMAT (' HIGHLIGHT PERCOLATING (recharging/none)')
 8005 FORMAT (' PLOT LAYOUT (baseflow/overlandflow/streamflow)')
 8006 FORMAT (' PLOT BASEFLOW (layout/overlandflow/streamflow)')
 8007 FORMAT (' PLOT OVERLANDFLOW (layout/baseflow/streamflow)')
 8008 FORMAT (' PLOT STREAMFLOW (layout/baseflow/overlandflow)')
 8100 FORMAT (' <Esc> or QUIT',/,' >')
 8200 FORMAT (' ***WARNING: no valid solution, data likely in ERROR !'/)
 8210 FORMAT (' ***WARNING: stream network detected, but no',
     &' conjunctive solution!',/, ' Data likely in ERROR !'/)
 9000 FORMAT(80A1)
 9001 FORMAT(' ***ILLEGAL OR MISSING PARAMETERS in line sink module:',/,
     &       ' ',80A1,/)
 9002 FORMAT(' *** ILLEGAL COMMAND in line sink module:',/,
     &       ' ',80A1,/)
 9003 FORMAT(' ***WARNING: no discharge specified for gallery, ',
     &       'discharge set to zero.')
 9005 FORMAT(' X1,Y1 X2,Y2 SIGMA (extr. rate L^2/T) [WIDTH] [LABEL]',/)
 9006 FORMAT(' X1,Y1,X2,Y2 HEAD [WIDTH] [RESISTANCE] [DEPTH] [LABEL]',/)
 9040 FORMAT(' Enter nodes, heads and optional parameters in string.',/
     &       ' Parameters apply to line sink following the node.',/ 
     &       ' When at last node enter X,Y only!',/,
     &' Include "end" command to indicate last downstream string of a',
     &' stream.',/,' Include "overlandflow (Q)" to add overland flow;',
     &' Q in L^3/T for string',/)
 9050 FORMAT(' X,Y HEAD [WIDTH] [RESISTANCE] [DEPTH] [LABEL]',/)
 9600 FORMAT (' ***ERROR: open end (label=',A16,')')
 9760 FORMAT (' baseflow  ,overlandflow,  head    , waterlevel,  length' 
     &,'   ,  strength , resistance, width     ,   depth   ,    label') 
 9761 FORMAT (9(G11.4,','),A16)      
 9765 FORMAT (' File ',A,' has been written.')
 9766 FORMAT (' ***ERROR in closing file ',A,' IOSTAT=',I5) 
 9800 FORMAT (' Press any key to continue.')
      END
C
C-------------------------------------------------------------------------------------------------------
C
      SUBROUTINE RELSPT(ICOD,ROF,RHW,RQGALLERY,RH0,RH2)
C
C-------------------------------------------------------------------------------------------------------
C
C     Input:
C     ICODE is absolute ICOD
C     ICODE = 1 discharge specified line sinks
C     ICODE = 2 drain feature
C     ICODE = 3 gallery feature
C     ICODE = 4 head specified line sinks
C     ICODE = 5 strings of head specified line sinks
C     ICOD = -5 end node for string of line sinks
C     ROF       overland inflow   (when ICODE =5)
C     RHW       head water inflow or evapotranspiration (only when ICODE =5)
C     RQGALLERY total pumping rate for a gallery
C     RH0       bottom elevation of lake
C     RH2       second head estimate of lake
C
C     Output:
C     ICOD = 0 zero line sink in string, abort string input.
C      
      IMPLICIT NONE
      INCLUDE 'LSCOM.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      INTEGER(4) ICOD,ICODE,I3,I4,I5,I6,I7,JUMP,ISTART,IOFF,I
      LOGICAL LFIRST,LBAD,LEQUALHD
      REAL(8) ROF,RHW,RTOL,RFGTOL,RDUM1,RVAR,RDUM2,RDUM3,RDUM4,
     &        RLEN,RFGRTOL,RQGALLERY,RH0,RH2,RDUM
      COMPLEX(8) CDUM1,CVAR,CDUM2,CTEMP
      CHARACTER(1) APAR (30),AWORD(68)
      DATA LFIRST /.TRUE./
      DATA APAR/  'E','N','D',' ',
     .            'O','V','E','R',' ',
     .            'I','N','L','E',' ',
     .            'O','U','T','L',' ',
     .            'L','A','K','E',' ',
     .            'L','E','V','E',' ',
     .            ATERM/
      DATA AWORD/ 'D','I','S','C',' ',     ! 1
     .            'D','R','A','I',' ',     ! 2
     .            'G','A','L','L',' ',     ! 3
     .            'H','E','A','D',' ',     ! 4
     .            'S','T','R','E',' ',     ! 5
     .            '?',' ',                 ! 6
     .            'Q','U','I','T',' ',     ! 7
     .            'W','I','D','T',' ',     ! 8
     .            'R','E','S','I',' ',     ! 9
     .            'D','E','P','T',' ',     ! 10
     .            'P','U','M','P',' ',     ! 11
     .            'H','I','G','H',' ',     ! 12
     .            'C','U','R','S',' ',     ! 13
     .            'P','L','O','T',' ',     ! 14
     .            ATERM/
      save
      ICODE=IABS(ICOD)
      RTOL=0.001*RFGRTOL()  ! reduced tolerance to avoid unecessary warnings for short line-sinks (7/7/00)
      IF (ICODE.EQ.5.AND.NLSTRING.GE.NLSTMX) THEN
      WRITE (ILUER,2000) NLSTMX   ! string arrays at capacity, abort
      RETURN
      ENDIF
      CDUM1=CVAR(1)   ! read x1,y1
      IF (ICODE.NE.5) THEN
        CDUM2=CVAR(3) ! read x2,y2 from same line (not a stream feature)
        I3=5
        I4=6
        I5=7
        I6=8
        I7=9
      ELSE  ! stream feature has 2 reals (x1,y2) less on a line
        I3=3
        I4=4
        I5=5
        I6=6
        I7=7
      ENDIF
      IF (ICOD.GT.0) RDUM1=RVAR(I3)  ! read head at line sink center
      IF (LERROR.OR.LMISS) THEN
      WRITE (ILUER,9001) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      RETURN
      ENDIF
      RDUM2=RVAR(I4)                ! try to read optional width parameters
      IF (LERROR.OR.LMISS) THEN     ! not found, set defaults for:
      RDUM2=RLSW0                   ! width
      RDUM3=RLSC0                   ! resistance
      RDUM4=RLSD0                   ! depth
      CALL GETFN (I4)     ! read optional label (following head)
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 5
      ENDIF
      RDUM3=RVAR(I5)      ! try to read optional resistance parameter
      IF (LERROR.OR.LMISS) THEN     ! not found, set defaults for:
      RDUM3=RLSC0                   ! resistance
      RDUM4=RLSD0                   ! depth
      CALL GETFN (I5)     ! read optional label (following width parameter)
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 5
      ENDIF
      RDUM4=RVAR(I6)                ! try to read optional depth parameter
      IF (LERROR.OR.LMISS) THEN     ! not found, set default for:
      RDUM4=RLSD0                   ! depth
      CALL GETFN (I6)     ! read optional label (following resistance parameter)
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 5
      ENDIF
      CALL GETFN (I7)     ! read optional label (following depth parameter)
      LERROR=.FALSE.
      LMISS=.FALSE.
 5    IF (NLS.EQ.NLSMAX) THEN   ! line-sink arrays are at capacity, abort
      WRITE (ILUER,1000) NLSMAX
      RETURN
      ENDIF
      IF (ICODE.EQ.5) THEN  ! (1)
C --------------------------------- string of head specified line sinks (stream command)
      CALL INLINE       ! read ahead
      CTEMP=CVAR(1)
      CTEMP=CTEMP+1.0 ! dummy statement to avoid compiler warning
      IF ((LERROR.OR.LMISS).AND.ICOD.GT.0) THEN ! (2) no coordinates, but may be a legal command
        LERROR=.FALSE.
        LMISS=.FALSE.
        CALL MATCH(AWORD,1,JUMP,LBAD)    ! start looking for command words
        IF (LBAD) THEN  ! (3)
          LBAD=.FALSE.
        ELSE            ! (3)
          IF (JUMP.LT.8.OR.JUMP.GT.10) THEN  ! (4) only "quit", "width", "resistance" and "depth" are legal here.
            WRITE (ILUER,4050) CDUM1,AFILE
            IF (LUCON) READ (ILUIN,1500) ALINE(80)
            ICOD=-5   ! interpret illegal command as quit and use last vertex to end string
          ENDIF ! (4)
        ENDIF   ! (3)
      ENDIF ! (2)
      IF (LFIRST) THEN  ! (5) first vertex, start a new string
        ISTART=NLS+1
      ELSE              ! (5)
        IF (ABS(CLSZS(NLS)-CDUM1).EQ.0.0) THEN  ! (6) zero length line-sink, abort
          WRITE (ILUER,4000) CLSZS(NLS),ALSLAB(NLS),CDUM1
          IOFF=NLS-ISTART+1
          NLS=NLS-IOFF
          NLSH=NLSH-IOFF
          LSEND=.FALSE.
          LSINLETSTRING=.FALSE.
          LSOUTLETSTRING=.FALSE.
          LSLAKESTRING=.FALSE.
          LFIRST=.TRUE.
          ICOD=0  ! signal to LSIN not to accept any further points
          RETURN
        ENDIF     ! (6)
        IF (ABS(CLSZS(NLS)-CDUM1).LT.RTOL) THEN ! (7) very small line-sink, issue a warning
           WRITE (ILUER,4010) CLSZS(NLS),ALSLAB(NLS),CDUM1,AFILE
        ENDIF ! (7)
        CLSZE(NLS)=CDUM1 ! give first point of next line sink
C                        to last point of previous line sink
        CALL PLWIND (CDUM1)  ! upate max. window parameters
      ENDIF  ! (5)
      IF (ICOD.LT.0) THEN  ! (8) end of string, set pointers and link list
        CALL MATCH (APAR,1,JUMP,LBAD) ! check for "end", "overland", "inlet", "outlet" and "lake" commands
        IF (.NOT.LBAD) THEN  ! (9)
          IF (JUMP.EQ.1) LSEND=.TRUE.  ! end stream
          IF (JUMP.EQ.2) THEN  ! (10)        ! overland flow
            ROF=RVAR(2)
            IF (LERROR.OR.LMISS) THEN    ! (11) overland command lacks parameter, set to zero
              WRITE (ILUER,9001) ALINE2
              LERROR=.FALSE.
              LMISS=.FALSE.
              ROF=0.0
            ENDIF  ! (11)
          ENDIF   ! (10)
          IF (JUMP.EQ.3) THEN          ! end stream is inlet for lake
            LSEND=.TRUE.
            LSINLETSTRING=.TRUE.
          END IF                       ! stream is outlet for lake
          IF (JUMP.EQ.4) THEN
            CALL GETFN(2) ! read table filename into AFILE
            ATEMP=AFILE
            IF (LERROR) WRITE (ILUER,9050)
            LSOUTLETSTRING=.TRUE.
          END IF
         IF (JUMP.EQ.5) THEN     ! stream forms lake
         CALL GETFN(2) ! read table filename into AFILE
         ATEMP=AFILE
         IF (LERROR) WRITE(ILUER,9060)
         LSLAKESTRING=.TRUE.
         END IF
         IF (JUMP.EQ.6) THEN   ! read lake levels
         RDUM=RVAR(2)
         IF (LERROR) THEN
         WRITE (ILUER,9060)
         RH0=0.0D0
         ELSE
         RH0=RDUM       ! lake bottom elevation
         ENDIF
         RDUM=RVAR(3)
         IF (LERROR) THEN
         WRITE (ILUER,9060)
         RH2=0.0D0
         ELSE
         RH2=RDUM       ! second head estimate
         ENDIF
         END IF
        ENDIF ! (9)
        NLSTRING=NLSTRING+1
        IF (LSLAKESTRING) THEN
          RLSEVAP(NLSTRING)=RHW
        ELSE
          RLSOVHW(NLSTRING)=RHW
        END IF
        RLSOFST(NLSTRING)=ROF
        IF (LSLAKESTRING) THEN  ! is lake feature, add:
          RLSH0(NLSTRING)=RH0   ! lake bottom elevation
          RLSH1(NLSTRING)=RDUM1 ! first lake stage estimate
          RLSH2(NLSTRING)=RH2   ! second lake stage estimate
        ENDIF
        RLEN=0.0
        DO 8 I=ISTART,NLS                ! calculate length of string
        RLEN=RLEN+ABS(CLSZE(I)-CLSZS(I))
   8    CONTINUE
        DO 9 I=ISTART,NLS             ! calculate overland inflow per line sink
        RLSOFSIG(I)=RLSOFST(NLSTRING)/RLEN
   9    CONTINUE
        IF (ISTART.EQ.NLS) THEN ! string contains only one line sink
          KLSTRING(NLSTRING)=NLS
          KLSTREND(NLSTRING)=NLS
          KLSUP(NLS)=NLS
          KLSDN(NLS)=NLS
          IF (LSEND) KLSDN(NLS)=0
          IF (LSINLETSTRING) LSINLET(NLS)=.TRUE.
          IF (LSOUTLETSTRING) LSOUTLET(NLS)=.TRUE.
          IF (LSLAKESTRING) LSLAKE(NLS)=.TRUE.
          IF (LSOUTLETSTRING) ALSTBLFILENAME(NLSTRING)=ATEMP
          IF (LSLAKESTRING) ALSTBLFILENAME(NLSTRING)=ATEMP
          GOTO 25
        ENDIF   
        IF (RLSH(ISTART).GE.RLSH(NLS)) THEN ! line sinks oriented downstream
          LEQUALHD=RLSH(ISTART).EQ.RLSH(NLS)
          IF (LEQUALHD.AND..NOT.LSLAKESTRING)
     &                        WRITE (ILUER,5000) ALSLAB(ISTART)
          KLSTRING(NLSTRING)=ISTART
          KLSTREND(NLSTRING)=NLS
          DO 10 I=ISTART+1,NLS-1
          KLSDN(I)=I+1
          KLSUP(I)=I-1
  10      CONTINUE
          KLSDN(ISTART)=ISTART+1
          KLSUP(NLS)=NLS-1
          IF (LSEND.OR.LEQUALHD) KLSDN(NLS)=0
          IF (LSINLETSTRING) LSINLET(NLS)=.TRUE.
          IF (LSOUTLETSTRING) LSOUTLET(ISTART)=.TRUE.
          IF (LSOUTLETSTRING) ALSTBLFILENAME(NLSTRING)=ATEMP
          IF (LSLAKESTRING) ALSTBLFILENAME(NLSTRING)=ATEMP
          GOTO 25
        ENDIF
        IF (RLSH(ISTART).LT.RLSH(NLS)) THEN  ! line sinks oriented upstream
          KLSTRING(NLSTRING)=NLS
          KLSTREND(NLSTRING)=ISTART
          DO 20 I=ISTART+1,NLS-1
          KLSUP(I)=I+1
          KLSDN(I)=I-1
  20      CONTINUE
          KLSUP(ISTART)=ISTART+1
          KLSDN(NLS)=NLS-1
          IF (LSEND) KLSDN(ISTART)=0
          IF (LSINLETSTRING) LSINLET(ISTART)=.TRUE.
          IF (LSOUTLETSTRING) LSOUTLET(NLS)=.TRUE.
          IF (LSOUTLETSTRING) ALSTBLFILENAME(NLSTRING)=ATEMP
C         should not occur for a lake feature
        ENDIF
 25     LSEND=.FALSE.           ! reset logicals for next string
c        if (lsoutletstring.or.lslakestring)
c     &  write (iluer,1002) NLSTRING,ALSTBLFILENAME(NLSTRING)
c 1002   format (' RELSPT2: NLSTRING,ALSTBLFILENAME ',i5,2x,a16)
        LSINLETSTRING=.FALSE.
        LSOUTLETSTRING=.FALSE.
        LSLAKESTRING=.FALSE.
        LFIRST=.TRUE.
      ELSE                    ! (8) next line sink, store data
        NLS=NLS+1
        NLSH=NLSH+1
        KLSPTH(NLSH)=NLS
        CLSZS(NLS)=CDUM1
        CALL PLWIND (CDUM1)  ! upate max. window parameters        
        RLSH(NLS)=RDUM1
        RLSWID(NLS)=RDUM2
        ilsbound(nls)=ilsbound0
c        IF (RDUM3.LT.0.01.AND.RDUM3.GT.0.0) THEN  ! surpress low resistance warning
c        WRITE (ILUER,6000) AFILE,RDUM3
c        ENDIF
        IF (RDUM3.GT.10000.0) THEN
        WRITE (ILUER,7000) AFILE,RDUM3
        RDUM3=10000.0
        ENDIF
        RLSRES(NLS)=RDUM3
        RLSDEP(NLS)=RDUM4
        ALSLAB(NLS)=AFILE
        LSGIV(NLS)=.FALSE.
        IF (LSLAKESTRING) LSLAKE(NLS)=.TRUE.
        KLSUP(NLS)=NLS
        KLSDN(NLS)=NLS
        LFIRST=.FALSE.
      ENDIF    ! (8)
        RETURN ! to LSIN for next input line
      ENDIF    ! (1)
c
C     -------------------------- end of ICODE=5 logic (line-sink strings)
c
      IF (ABS(CDUM1-CDUM2).EQ.0.0) THEN     ! zero length line-sink, abort
        WRITE (ILUER,3010) CDUM1,CDUM2,AFILE
        RETURN
      ENDIF
      IF (ABS(CDUM1-CDUM2).LT.RTOL) THEN    ! very small line-sink, issue a warning
        WRITE (ILUER,3000) CDUM1,CDUM2,AFILE
      ENDIF      
      NLS=NLS+1
      CLSZS(NLS)=CDUM1
      CALL PLWIND (CDUM1)           ! update max. window
      CLSZE(NLS)=CDUM2
      CALL PLWIND (CDUM2)           ! update max. window
C      
      IF (ICODE.EQ.1) THEN
C ----------------------------------strength specified line sinks      
      NLSIG=NLSIG+1
      KLSPTS(NLSIG)=NLS
      RLSIG(NLS)=RDUM1
      RLSWID(NLS)=RDUM2
      ilsbound(nls)=ilsbound0
      RLSRES(NLS)=0.0
      ALSLAB(NLS)=AFILE
      KLSUP(NLS)=NLS
      KLSDN(NLS)=NLS
      RETURN
      ENDIF
      IF (ICODE.EQ.2.OR.ICODE.EQ.3.OR.ICODE.EQ.4) THEN
C ---------------------------------head specified line-sink or drain or gallery
      NLSH=NLSH+1
      KLSPTH(NLSH)=NLS
      RLSH(NLS)=RDUM1
      RLSWID(NLS)=RDUM2
      ilsbound(nls)=ilsbound0
c      IF (RDUM3.LT.0.01D0.AND.RDUM3.GT.0.0D0) THEN ! small resistance
c      WRITE (ILUER,6000) AFILE,RDUM3
c      ENDIF
      IF (RDUM3.GT.10000.0) THEN
      WRITE (ILUER,7000) AFILE,RDUM3
      RDUM3=10000.0
      ENDIF
      RLSRES(NLS)=RDUM3
      RLSDEP(NLS)=RDUM4
      ALSLAB(NLS)=AFILE
      LSGIV(NLS)=.FALSE.
      KLSUP(NLS)=NLS
      KLSDN(NLS)=NLS
      if (icode.eq.2) then  ! drain
        lsdrain(nls)=.true.
        lsgiv(nls)=.true.     ! do not include in first iteration
        rlsig(nls)=0.0D0
      end if
      if (icode.eq.3) then   ! gallery
        lsgallery(nls)=.true.
        rlsq(nls)=rqgallery
        rlshmin(nls)=rlsh(nls)
      end if
      IF (rlsdep(nls).EQ.0.0D0.AND.rlsres(nls).GT.0.0D0) THEN
        if (.not.lsgallery(nls).and..not.lsdrain(nls)) then
          WRITE (ILUER,8000) AFILE  ! warn for zero depth on line-sinks (not drains or galleries)
        end if
      ENDIF

      RETURN
      ENDIF
      WRITE (ILUER,9002) ICODE
      RETURN
C      
 1000 FORMAT (' ***ERROR in line sink module:',
     &' too many line sinks (max.=',I3,')',/)
 1500 FORMAT (A1)     
 2000 FORMAT (' ***ERROR in line sink module:',
     &' too many streams (max.=',I3,')',/)
 3000 FORMAT (' ***WARNING in line sink module:',
     &' line sink may be too short.',/,
     &' line sink points: ',2F11.1,2x,2F11.1,/,
     &' line sink label ',A16,/)
 3010 FORMAT (' ***ERROR in line sink module: zero length line sink,',/,
     &' line sink points: ',2F11.1,2X,2F11.1,/,
     &' line sink label ',A16,', line sink ignored!',/)
 4000 FORMAT (' ***ERRROR in line sink module:',
     &' zero length line sink in string.',/,
     &' point ',2F11.1,' with label ',A16,/,' too close to point ',
     &2F11.1,/,' entire string will be ignored.',/,
     &' Remaining points in input file will generate "illegal command"',
     &' messages.',/)
 4050 FORMAT (' ***ERROR in line sink module:',/,
     &' String should end with a coordinate pair without a head.',/,
     &' Point ',2G11.4,' with label ',A16,/,' used as end coordinates.',
     &' Associated head has been ignored.',/,
     &' Press any key to continue.')
 4010 FORMAT (' ***WARNING in line sink module:',/,
     &' point ',2F11.1,' with label ',A16,' may be too close to',/,
     &' point ',2F11.1,' with label ',A16,/)
 5000 FORMAT (' ***WARNING in line sink module:',
     &' no gradient in stream section;',/,
     &' flow assumed into the direction in which the elements were',
     &' entered.',/,
     &' First stream element of string has label: ',A16,/,
     &' Stream is declared and "end stream"; will not be connected.',/)
 6000 FORMAT (' ***WARNING: line sink with label ',A16,' has a ',
     &'resistance ',G11.4,/,
     &' A resistance below 0.01 may lead to an unstable ',
     &'solution procedure!')
 7000 FORMAT (' ***WARNING: line sink with label ',A16,' has a ',
     &'resistance ',G11.4,/,
     &' This line sink will play no significant role in the ',
     &'groundwater flow solution!',/,' Resistance parameter is set to ',
     &'10000 to conform to format of the DATA command.',/)
 8000 FORMAT (' ***WARNING: depth is zero, but resistance is not zero.',
     &/,' Line sink ',A16,' may cause numerical instability!')
 9001 FORMAT(' ***ILLEGAL or MISSING PARAMETER(S) in line sink module'
     &,/,' ',80A1,/)
 9002 FORMAT(' ***ILLEGAL ICODE in RELSPT: ',I10,/)
 9050 FORMAT(' ***ERROR in RELSPT: could not read stage table filename')
 9060 FORMAT (' ***ERROR in RELSPT: illegal or missing parameter')
      END