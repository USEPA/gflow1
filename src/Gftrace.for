C     Last change:  HMH   7 Aug 2013   10:22 am
c     This file contains the following routines and functions
c
c     BLOCK DATA COM3DDAT   initializes COM3D and TRACOM common blocks
c     SUBROUTINE STREAM     driver of streamline tracing (input of commands)
c     function lwellinfo    searches for wells and passes well data (called in STREAM)
c     SUBROUTINE TRACE      streamline tracing routine (executing traces)
c     SUBROUTINE PREDCOR    PREDICTOR CORRECTOR method for integrating velocities in 3-D
c     SUBROUTINE RUNKUT     RUNGE-KUTTA method for integrating velocities in 3-D
c     SUBROUTINE GETSTEP    get the current stepsize
c          ENTRY SETSTEP    set a new stepsize
c          ENTRY ZVAL       returns value of R3DZ (vertical elevation of streamline)
c
C
c
c
c
c
c
c ---------------------------------------------------------------------------------
c
      BLOCK DATA COM3DDAT
c
c ---------------------------------------------------------------------------------
c
      IMPLICIT NONE
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      DATA RETARDATION,RHALFLIFE /1.0,0.0/
      DATA LTRACEOUT,LGRAPHICS,LBOUNDARY /.FALSE.,.TRUE.,.FALSE./
      END
c
c ---------------------------------------------------------------------------------
c      
      SUBROUTINE STREAM (RA,IRA,RSCR)
c
c ---------------------------------------------------------------------------------
c
C
C     Routine drives the tracing of three-dimensional streamlines.
C
      IMPLICIT NONE
      INTEGER(4) IRA,NWORD,NPAR1,NPAR2,NPAR3,NPTMAX,ICOUNT,JUMP,
     &           I,ILINES,J,NWSTR,NST,IVAR
      INTEGER :: ivalues(8)
      LOGICAL LSTEPSELECT,LBAD,LCONT0,LCURS0,
     &        LDARK0,LRET,LSURFER,LREPORT,LCONT,LAY,LDARK,LDBOUTSIDE,
     &        LAY0,LWELLINFO
      REAL(8) RA,RSCR,RCLEVEL,RDTIME,RTIC,RDS,RPI,RDIV,RDUMX,RDUMY,
     &        R1,R2,R3,R4,RTOP,RBAS,RTOPT,RBAST,RDUM,RHEAD,RXW,RYW,RZ,
     &        RQWELL,RWELL,RSWITCH,RQABS,RDTET,RTET,RTET1,RDUM1,RDUM2,
     &        RQI,RXST,REST,SQROOT,RVAR,RFTOP,RFBASE,RHGHT,RF3DSP,
     &        RFHEAD,RFHGHT
      DIMENSION ra(ira,*),rscr(*)
      COMPLEX(8) CI,CZ,CZWELL,CDEL,cvar
      PARAMETER (NWORD=111,NPAR1=8,NPAR2=11,NPAR3=11,NPTMAX=50000)   ! increased from 5,000 on 8/7/2013
      CHARACTER(1) AWORD(NWORD),APAR1(NPAR1),APAR2(NPAR2),APAR3(NPAR3)
      CHARACTER(1) ADUM
      CHARACTER(8) ADATE
      CHARACTER(11) ATIME
      CHARACTER(5) azone
      CHARACTER(16) ATFIL,APATH,AWELL
      INCLUDE 'COM3D.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'DWMN.INC'
      DIMENSION RXST(4,NPTMAX),RQI(3)
      DATA LAY,LCONT,LREPORT,LSURFER /.TRUE.,.FALSE.,.TRUE.,.FALSE./
      DATA LDARK /.FALSE./
      DATA NST,REST,RDS,RTIC,RDTIME,RCLEVEL /0,1.0E30,1.0,0.0,0.0,0.0/
      DATA RPI /3.1415926/
      DATA CI /(0.0,1.0)/
      DATA AWORD /'?',' ',
     &            'P','O','I','N',' ',
     &            'T','I','M','E',' ',
     &            'S','T','E','P',' ',
     &            'G','O',' ',
     &            'Q','U','I','T',' ',
     &            'D','I','R','E',' ',
     &            'S','H','O','W',' ',
     &            'L','A','Y','O',' ',
     &            'C','O','N','T',' ',
     &            'Z','M','A','R',' ',
     &            'C','U','R','S',' ',
     &            'T','M','A','R',' ',
     &            'R','E','P','O',' ',
     &            'G','R','I','D',' ',    ! command added for MASTER project
     &            'P','I','C','T',' ',
     &            'F','I','L','E',' ',
     &            'S','U','R','F',' ',
     &            'S','T','R','E',' ',
     &            'W','E','L','L',' ',
     &            'T','R','A','N',' ',
     &            'C','M','A','R',' ',
     &            'B','O','U','N',' ',
     &            ATERM/
      DATA APAR1 /'O','N',' ',
     &            'O','F','F',' ',
     &            ATERM/
      DATA APAR2 /'F','O','R','W',' ',
     &            'B','A','C','K',' ',
     &            ATERM/
      DATA APAR3 /'D','A','R','K',' ',
     &            'V','I','S','I',' ',
     &            ATERM/
      ICOUNT=0          ! for MASTER project MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      L3DEND=.FALSE.
      LSTEPSELECT=.FALSE.
      LBOUNDARY=.FALSE.
      IF (.NOT.LGRCOM) LCONT=.FALSE.
      RDIV=100                        ! stepsize is 1/100 th of window
      IF (.NOT.L3DPL) THEN
      RDUMX=(RX2-RX1)/RDIV
      RDUMY=(RY2-RY1)/RDIV
      ELSE
      RDUMX=(R3DX2-R3DX1)/RDIV
      RDUMY=(R3DY2-R3DY1)/RDIV
      IF (L3DVER) THEN
      RDUMX=SQROOT(RDUMX*RDUMX+RDUMY*RDUMY)
      RDUMY=(R3DZ2-R3DZ1)/RDIV
      ENDIF
      ENDIF
      RDS=MIN(RDUMX,RDUMY)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
  10  IF (LERROR.OR.LMISS) WRITE (ILUER,2000) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      if (lucon) then
      WRITE (ILUME,1100) NST,NPTMAX
      WRITE (ILUME,1103)
      WRITE (ILUME,1105)
      WRITE (ILUME,1110) REST
      WRITE (ILUME,1120) RDS
      WRITE (ILUME,1125) RTIC,RDTIME
      WRITE (ILUME,1126) RCLEVEL
      WRITE (ILUME,1127) RETARDATION,RHALFLIFE
      IF (L3DREV) WRITE (ILUME,1130)
      IF (.NOT.L3DREV) WRITE (ILUME,1132)
      IF (.NOT.LDARK) WRITE (ILUME,1135)
      IF (LDARK) WRITE (ILUME,1137)
      IF (.NOT.LTRACEOUT) WRITE (ILUME,1138)
      IF (LTRACEOUT) WRITE (ILUME,1139) APATH
      IF (.NOT.LSURFER) WRITE (ILUME,1142)
      IF (LSURFER) WRITE (ILUME,1143) ATFIL
      IF (LREPORT) WRITE (ILUME,1145)
      IF (.NOT.LREPORT) WRITE (ILUME,1146)
      IF (LAY) WRITE (ILUME,1150)
      IF (.NOT.LAY) WRITE (ILUME,1160)
      IF (LCONT) WRITE (ILUME,1162)
      IF (.NOT.LCONT) WRITE (ILUME,1163)
      IF (LCURS) WRITE (ILUME,1165)
      IF (.NOT.LCURS) WRITE (ILUME,1167)
      IF (LGRAPHICS) WRITE (ILUME,1168)
      IF (.NOT.LGRAPHICS) WRITE (ILUME,1169)
      WRITE (ILUME,1170)
      endif
      CALL INLINE
      CALL MATCH (AWORD,1,JUMP,LBAD)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LBAD) GOTO 15
      GOTO (10,13),JUMP
  13  WRITE (ILUER,3000) ALINE2
      LERROR=.FALSE.
      GOTO 10
  15  GOTO (100,200,300,400,500,600,700,800,900,950,960,970,980,985,990,
     &992,994,995,997,996,998,999,993),JUMP  ! label 990 added for MASTER project
C
C     help
C  
 100  AFILE='STRMHLP.HLP    '
C      CALL HELP                       ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     point
C      
  200 IF (NST.EQ.NPTMAX) THEN
        WRITE (ILUER,4990) NPTMAX
        GOTO 10
      ENDIF
      if (lucon) WRITE (ILUME,5000) NST+1
      CALL INLINE
      IF (ALINE(1).EQ.'Q'.OR.ALINE(1).EQ.'q') THEN  ! quit entering points
        IF (LCURS.AND.NST.GT.0.and.lucon) THEN
C          CALL TONE         ! NOT AVAILABLE IN BATCH MODE
          WRITE (ILUME,1900)
        ENDIF
        GOTO 10
      ENDIF
      IF (ALINE(1).EQ.'C'.OR.ALINE(1).EQ.'c') THEN  ! clear the points buffer
        NST=0
        if (lucon) WRITE (ILUME,1901)
        GOTO 10
      ENDIF
      R1=RVAR(1)     ! X-coordinate of streamline starting point
      R2=RVAR(2)     ! Y-coordinate of streamline starting point
      R3=RVAR(3)     ! Z-coordinate of streamline starting point
      IF (LERROR) THEN
        WRITE (ILUER,2000)
        LERROR=.FALSE.
        LMISS=.FALSE.
        GOTO 200
      ENDIF
      IF (LDBOUTSIDE(R1,R2)) THEN ! reject when point is outside the model domain.
      WRITE (ILUER,9100) R1,R2
      GOTO 200
      END IF
      R4=RVAR(4)    ! direction: +1.0=forward,  0.0=use default direction,  -1.0=backward
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      R4=0.0        ! if R4 not provided set to default direction
      ENDIF
      CZ=CMPLX(R1,R2)
      RTOP=RFTOP(CZ) ! changed on 2/9/99 from saturated height to top of confining layer
      RBAS=RFBASE(CZ)
      rhght=rfhght(cz)
      RTOPT=RTOP+0.001*RHGHT
      RBAST=RBAS-0.001*RHGHT
c      IF (R3.LT.RBAST.OR.R3.GT.RTOPT) THEN
      IF (R3.LT.RBAST) THEN  ! changed on 11/21/99; only test for base, drop point down later
        WRITE (ILUER,2005) R3,RBAS,RTOP
        GOTO 200
      ENDIF
      NST=NST+1
      RXST(1,NST)=R1
      RXST(2,NST)=R2
      RXST(3,NST)=R3
      RXST(4,NST)=R4
      GOTO 200
C
C     time
C      
 300  R1=RVAR(2)
      IF (LERROR) GOTO 10
      REST=R1
      GOTO 10
C
C     step
C      
 400  R1=RVAR(2)
      IF (LERROR) GOTO 10
      RDS=R1
      LSTEPSELECT=.TRUE.
      GOTO 10
C
C     go
C      
 500  CONTINUE
C      IF (LCONT) THEN                        ! NOT AVAILABLE IN BATCH MODE
C        CALL PGRPAR (LESCAP)
C      ELSE
C        IF (LGRAPHICS) CALL PLAPAR (LESCAP)
C      ENDIF
C      IF (LESCAP) GOTO 10
      IF (.NOT.L3DPL) THEN
      RDUMX=(RX2-RX1)/RDIV
      RDUMY=(RY2-RY1)/RDIV
      ELSE
      RDUMX=(R3DX2-R3DX1)/RDIV
      RDUMY=(R3DY2-R3DY1)/RDIV
      IF (L3DVER) THEN
      RDUMX=SQROOT(RDUMX*RDUMX+RDUMY*RDUMY)
      RDUMY=(R3DZ2-R3DZ1)/RDIV
      ENDIF
      ENDIF
      RDUM=MIN(RDUMX,RDUMY)
      IF (LSTEPSELECT) THEN
        IF (RDS.GT.2*RDUM.OR.RDS.LT.0.5*RDUM) THEN
C          IF (LGRAPHICS) THEN                      ! NOT AVAILABLE IN BATCH MODE
C            WRITE (ILUME,9070) RDS,RDUM
C          ENDIF
        ENDIF
      ELSE
        RDS=RDUM
      ENDIF
      RDS0=RDS
C      IF (LGRAPHICS) THEN                         ! NOT AVAILABLE IN BATCH MODE
C       IF (LSCRN.AND..NOT.LCONT) CALL PAGE
C       IF (.NOT.LCONT) CALL PLOTON
C       IF (LCONT) CALL CNSTRM (RA,IRA,RSCR)
C       IF (LAY) CALL LAYOUT
C       IF (LSCRN.AND..NOT.LSINGL) CALL PLOTOF
CC
CC   colors: 1=white, 2=red,  3=green, 4=blue, 5=yellow, 6=cyan
CC           7=magenta, 8=dark gray, 9=light gray, 10=light red,
CC           11=light green, 12=light blue, 13=brown, 14=light cyan,
CC           15=light magenta
CC   lines:  1=solid,  2=dash, 3=dots,  4=dash-dot, 5=medium dashed,
CC           6=dash with two dots, 7=short dash
CC
C       IF (LCOLOR) THEN
C        CALL DRWCLL (2)
C        CALL DRWTPL (1)
C       ELSE
C        CALL DRWTPL (2)
C       ENDIF
C       IF (LHPLOT.AND.LCOLOR) THEN
C        WRITE (ILUME,7000)
C        CALL GFNEWPEN
C       ENDIF
C      ENDIF
      do i=1,nst
      cz=CMPLX(rxst(1,i),rxst(2,i))
      rhead=rfhead(cz) ! head at starting point
      rtop=rftop(cz)   ! aquifer top at starting point
      if (rhead.lt.rtop) then ! unconfined flow
        if (rxst(3,i).gt.rhead) rxst(3,i)=rhead ! drop starting point to water table
      end if
      end do
!     now enter trace
      CALL TRACE (RXST,NPTMAX,NST,RDS,REST,RTIC,RDTIME,LREPORT,LSURFER,
     &LDARK,RCLEVEL)
C      IF (LGRAPHICS) THEN                 ! NOT AVAILABLE IN BATCH MODE
C       CALL GRAPHEDIT (RA,IRA,RSCR)
C       IF (LSINGL) CALL PLOTOF
C       IF (.NOT.LSCRN) THEN
C        CALL CLEARSCREEN
C        WRITE (ILUME,9050)
C        GOTO 500
C       ENDIF
C      ENDIF
      GOTO 10
C
C     return
C      
 600  IF (LSURFER) CLOSE (3)
      LSURFER=.FALSE.
      RETURN
C
C     direction
C 
 700  CALL MATCH (APAR2,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      L3DREV=JUMP.EQ.2
      GOTO 10
C
C     show
C      
 800  IF (NST.EQ.0) THEN
        if (lucon) WRITE (ILUME,8005)
        GOTO 10
      ENDIF
      if (lucon) WRITE (ILUME,8000)
      IF (LUOUTFILE) THEN
        call date_and_time (adate,atime,azone,ivalues)
        WRITE (ILUOUT,8007) ADATE,ATIME,ATITLE
      ENDIF
      ILINES=0
      DO 850 I=1,NST
      ILINES=ILINES+1
      if (lucon) WRITE (ILUME,8010) I,(RXST(J,I),J=1,3)
      IF (LUOUTFILE.AND.RXST(4,I).NE.-5.0) 
     &WRITE (ILUOUT,8011) (RXST(J,I),J=1,3)
      IF (LUOUTFILE.AND.RXST(4,I).EQ.-5.0) 
     &WRITE (ILUOUT,8012) (RXST(J,I),J=1,3)
      IF (ILINES.EQ.20) THEN
        ILINES=0
        if (lucon) WRITE (ILUME,1700)
        READ (ILUIN,1600) ADUM
        IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') GOTO 10
        if (lucon) WRITE (ILUME,8000)
      ENDIF
 850  CONTINUE
      IF (LUOUTFILE) WRITE (ILUOUT,8015)
      if (lucon) WRITE (ILUME,1800)
      READ (ILUIN,1600) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') GOTO 10
C      CALL PAGE         ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     layout
C      
 900  CALL MATCH (APAR1,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      LAY=JUMP.EQ.1
      GOTO 10
C
C     contour
C
 950  CALL MATCH (APAR1,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      IF (JUMP.EQ.2) LCONT=.FALSE.
      IF (LGRCOM) THEN
      LCONT=JUMP.EQ.1
      ELSE
      WRITE (ILUER,9000)
      LCONT=.FALSE.
      ENDIF
      GOTO 10
C
C     set depth tic-mark interval
C      
C
 960  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RTIC=RDUM
      GOTO 10
C
C     cursor
C      
 970  CALL MATCH (APAR1,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      IF (JUMP.EQ.1) LCURS=.TRUE.
      IF (JUMP.EQ.2) LCURS=.FALSE.
      IF ((LCURS.AND.LPRINT).OR.(LCURS.AND.LPLOT)) THEN
      WRITE (ILUER,9060)
      LCURS=.FALSE.
      ENDIF
      IF (LCURS.AND.NST.GT.0.and.lucon) THEN
C      CALL TONE                             ! NOT AVAILABLE IN BATCH MODE
      WRITE (ILUME,1900)
      ENDIF
      GOTO 10
C
C     tmark
C      
 980  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RDTIME=RDUM
      GOTO 10      
C
C     report (of residence times)
C
  985 CALL MATCH (APAR1,2,JUMP,LBAD)    
      LMISS=LBAD
      IF (LMISS) GOTO 10
      LREPORT=JUMP.EQ.1
      GOTO 10
C
C     GRID command, added for MASTER project MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C      
 990  CONTINUE                  ! REMOVED, FOR RESEARCH PURPOSES ONLY
C      NNX=IVAR(2)
C      NNY=IVAR(3)
C      IF (LERROR) GOTO 10
C      NSWITCH=IVAR(4) ! >0 end at 0,0   |1| write CZstart and T   |2| also plot CZstart
C
C     |nswitch|=1 trace streamlines starting inside the watershed contour only and when it ends
C                 inside the contour write the starting coordinates and the residence time
C     |nswitch|=2 trace streamlines starting inside the window (all starting points in the grid)
C                 and when it ends inside the watershed contour write as under |nswitch|=1
C                 and plot a dot at the location of the starting point. Streamlines are not plotted
C                 on the screen in this mode.
C
C      IF (LERROR) THEN
C      NSWITCH=0
C      LERROR=.FALSE.
C      ENDIF
C      LMAST=.TRUE.
C      IF (IABS(NSWITCH).EQ.2) LDARK=.TRUE. ! do not draw pathlines
C      IF (LCONT) THEN
C      CALL PGRPAR (LESCAP)
C      ELSE
C      CALL PLAPAR (LESCAP)
C      ENDIF
C      IF (LESCAP) GOTO 10
C      IF (LSCRN.AND..NOT.LCONT) CALL PAGE       ! NOT AVAILABLE IN BATCH MODE
C      IF (.NOT.LCONT) CALL PLOTON               ! NOT AVAILABLE IN BATCH MODE
C      IF (LCONT) CALL CNSTRM (RA,IRA,RSCR)      ! NOT AVAILABLE IN BATCH MODE
C      IF (LAY) CALL LAYOUT                      ! NOT AVAILABLE IN BATCH MODE
C      IF (LSCRN.AND..NOT.LSINGL) CALL PLOTOF    ! NOT AVAILABLE IN BATCH MODE
C      IF (LCOLOR) THEN                          ! NOT AVAILABLE IN BATCH MODE
C      CALL DRWCLL (2)                           ! NOT AVAILABLE IN BATCH MODE
C      CALL DRWTPL (1)                           ! NOT AVAILABLE IN BATCH MODE
C      ELSE                                      ! NOT AVAILABLE IN BATCH MODE
C      CALL DRWTPL (2)                           ! NOT AVAILABLE IN BATCH MODE
C      ENDIF                                     ! NOT AVAILABLE IN BATCH MODE
C      IF (LHPLOT.AND.LCOLOR) THEN               ! NOT AVAILABLE IN BATCH MODE
C      WRITE (ILUME,7000)                        ! NOT AVAILABLE IN BATCH MODE
C      CALL GFNEWPEN                             ! NOT AVAILABLE IN BATCH MODE
C      ENDIF                                     ! NOT AVAILABLE IN BATCH MODE
C      RRDX=(RX2-RX1)/NNX
C      RRDY=(RY2-RY1)/NNY
C      RRX0=RX1+0.5*RRDX
C      RRY0=RY1+0.5*RRDY
C      DO 993 I=1,NNX
C      RRX=RRX0+(I-1)*RRDX
C      DO 991 J=1,NNY
C      RRY=RRY0+(J-1)*RRDY
C      CZ=CMPLX(RRX,RRY)
C      IF (LWSIN(CZ).OR.IABS(NSWITCH).EQ.2) THEN  ! trace pathline
C      NST=1
C      RXST(1,1)=RRX
C      RXST(2,1)=RRY
C      RXST(3,1)=RFBASE(CZ)+RFHGHT(CZ)
C      CALL TRACE (RXST,NPTMAX,NST,RDS,REST,RTIC,RDTIME,LREPORT,LSURFER,
C     &LDARK,RCLEVEL)
C      ENDIF
C 991  CONTINUE
C 993  CONTINUE
C      CALL GRAPHEDIT (RA,IRA,RSCR)          ! NOT AVAILABLE IN BATCH MODE
C      IF (LSINGL) CALL PLOTOF               ! NOT AVAILABLE IN BATCH MODE
C      IF (ICOUNT.GT.0) WRITE (ILUOUT,1002) RESTIM  ! empty buffer
C 1002 FORMAT (G14.7,',',G14.7,',',G14.7,',',G14.7,',',G14.7)
      GOTO 10
C     ------------ end of the logic for the MASTER project MMMMMMMMMMMMMMMMMMMM
c
c
c
 992  CALL MATCH (APAR1,2,JUMP,LBAD)    
      LMISS=LBAD
      IF (LMISS) GOTO 10
      IF (LGRAPHICS) THEN
      LAY0=LAY
      LCONT0=LCONT
      LCURS0=LCURS
      LDARK0=LDARK
      ENDIF
      LGRAPHICS=JUMP.EQ.1
      IF (.not.LGRAPHICS) THEN
      LAY=.FALSE.
      LCONT=.FALSE.
      LCURS=.FALSE.
      LDARK=.true.
      ELSE
      LAY=LAY0
      LCONT=LCONT0
      LCURS=LCURS0
      LDARK=LDARK0
      ENDIF
      GOTO 10
C
C     Boundary (x1) (y1) (x2) (y2)
C
  993 continue
      ctr1=cvar(2)
      ctr2=cvar(4)
      IF (LERROR.OR.LMISS) GOTO 10
      LBOUNDARY=.TRUE.  ! pathlines will stop at (infinitely long) line through ctr1 and ctr2
      ctr2=ctr2-ctr1 ! ctr2 points from first to second point
      GOTO 10
C
C     FILE  FILENAME.PTH  write streamline traces for the GUI
C
  994 CONTINUE
c      CALL MATCH (APAR1,2,JUMP,LBAD) ! on/off not available in batch mode
c      IF (.NOT.LBAD) THEN
c      IF (JUMP.EQ.2) THEN
c      IF (LTRACEOUT) CLOSE (10)
c      LTRACEOUT=.FALSE.
c      GOTO 10
c      ENDIF
c      IF (JUMP.EQ.1) THEN
c      LERROR=.TRUE.
c      GOTO 10
c      ENDIF
c      ENDIF
      CALL GETFN (2)
      IF (LERROR.OR.LMISS) GOTO 10
      IF (LTRACEOUT) CLOSE (10)
      CALL FILECH ('.PTH')
      CALL CRAFIL (10,4,LRET,'.PTH')
      IF (LRET) GOTO 10
      APATH=AFILE
      LTRACEOUT=.TRUE.
      GOTO 10
c
C     surfer option: write boundary line file for 3-D streamlines
C      
  995 CALL MATCH (APAR1,2,JUMP,LBAD)
      IF (.NOT.LBAD) THEN
      IF (JUMP.EQ.2) THEN
      IF (LSURFER) CLOSE (3)
      LSURFER=.FALSE.
      GOTO 10
      ENDIF
      IF (JUMP.EQ.1) THEN
      LERROR=.TRUE.
      GOTO 10
      ENDIF
      ENDIF
      CALL GETFN (2)
      IF (LERROR.OR.LMISS) GOTO 10
      IF (LSURFER) CLOSE (3)
      CALL FILECH ('.BLN')
      CALL CRAFIL (3,3,LRET,'.BLN')
      IF (LRET) GOTO 10
      ATFIL=AFILE
      ISFREC=1
      ISFCNT=0
      LSURFER=.TRUE.
      GOTO 10
C
C     well x y z n
C
 996  RXW=RVAR(2)
      RYW=RVAR(3)
      CZ=CMPLX(RXW,RYW)
      RZ=RVAR(4)
      r3dz=rz
      NWSTR=IVAR(5)
      IF (LERROR.OR.LMISS) GOTO 10
      if (lwellinfo(cz,czwell,rwell,rqwell,awell)) then
       IF (ABS(RQWELL).LT.1.0E-10) THEN  ! Zero pumping rate, abort
         WRITE (ILUER,9091) AWELL,RQWELL
         GOTO 10
       ENDIF
       RHGHT=RFHGHT(CZ)
       RTOP=RFBASE(CZ)+RHGHT
       RBAS=RFBASE(CZ)
       RTOPT=RTOP+0.0001*RHGHT
       RBAST=RBAS-0.0001*RHGHT
       IF (RZ.LT.RBAST.OR.RZ.GT.RTOPT) THEN     ! starting elevation outside aquifer, abort
         WRITE (ILUER,2005) RZ,RBAS,RTOP,AWELL
         GOTO 10
      ENDIF
      RSWITCH=1.0
      IF (RQWELL.GT.0.0) RSWITCH=-1.0 ! trace back in time for pumping well
      CALL DISCH (CZWELL,RQI)  ! find ambient flow
      RQI(3)=0.0
      RQABS=RF3DSP(RQI,RQI)
      IF (RQABS.LT.1.0E-10) THEN ! no ambient flow, radial flow near the well
        RDTET=2.0*RPI/NWSTR
        RTET=0.0
        CDEL=2.0*CMPLX(RWELL,0.0)
      ELSE		         ! orient first streamline opposite ambient flow
        RTET1=50.0*0.18/NWSTR+0.02
        RDTET=(6.2831853-RTET1)/(NWSTR-1)
        CDEL=CMPLX(RQI(1),RQI(2))
        RTET=AIMAG(LOG(CDEL))
        RTET=RTET+0.5*RTET1
        CDEL=2.0*RWELL*EXP(CI*RTET)
      ENDIF
      RTET=RTET-RDTET
      DO I=1,NWSTR  ! add pathline starting points to the buffer
       RTET=RTET+RDTET
       CZ=CZWELL+2.0*RWELL*EXP(CI*RTET)
       NST=NST+1
       IF (NST.GT.NPTMAX) THEN
         WRITE (ILUER,9092) NPTMAX,ALINE2
         GOTO 10
       ENDIF
       RXST(1,NST)=REAL(CZ)
       RXST(2,NST)=AIMAG(CZ)
       RXST(3,NST)=RZ
       RXST(4,NST)=RSWITCH
       END do
      else
        WRITE (ILUER,9090) CZ ! now well was found, issue a warning
      endif
      GOTO 10      
C
C     streamlines visible/dark
C      
 997  CALL MATCH (APAR3,2,JUMP,LBAD)
      IF (LBAD) THEN
      LBAD=.FALSE.
      LERROR=.TRUE.
      GOTO 10
      ENDIF
      LDARK=JUMP.EQ.1
      GOTO 10
C
C     transport
C
 998  continue
C$IF .NOT.LTRANSPORT
      GOTO 13
C$ENDIF
      RDUM1=RVAR(2)            
      RDUM2=RVAR(3)
      IF (LERROR) GOTO 10
      IF (RDUM1.LE.0.0) THEN
      WRITE (ILUER,9001)
      GOTO 10
      ENDIF
      RETARDATION=RDUM1
      IF (RDUM2.LE.0.0) THEN
      WRITE (ILUER,9002)
      GOTO 10
      ENDIF
      RHALFLIFE=RDUM2
      GOTO 10
C
C     cmark
C
 999  CONTINUE
C$IF .NOT.LTRANSPORT
      GOTO 13
C$ENDIF
      RDUM1=RVAR(2)
      IF (LERROR) GOTO 10
      IF (RDUM1.LT.0.0.OR.RDUM1.GT.1.0) THEN
      WRITE (ILUER,9003)
      GOTO 10
      ENDIF
      RCLEVEL=RDUM1
      GOTO 10
C      
 1100 FORMAT ('        --------- TRACE module ---------',/,
     &        ' <F1> = Help',/,
     &        ' POINTS         (',I3,' specified, maximum ',I4,')')
 1103 FORMAT (' WELL  xwell  ywell  zstreamlines  nstreamlines')
 1105 FORMAT (' SHOW           (to display contents of POINTS buffer.)')
 1110 FORMAT (' TIME  ',G14.7,' (max. residence time for a streamline)')
 1120 FORMAT (' STEP  ',G14.7,' (step size along a streamline)')
 1125 FORMAT (' ZMARK ',G14.7,' (streamline depth increment)',/
     &        ' TMARK ',G14.7,' (groundw. res. time marker increment)')
 1126 FORMAT (' CMARK ',G14.7,' (conc. marker incr. between 0 and 1)')
 1127 FORMAT (' TRANSPORT ',G14.7,2X,G14.7,' (retardation) (halflife)')
 1130 FORMAT (' DIRECTION BACKWARD     (forward)')
 1132 FORMAT (' DIRECTION FORWARD      (backward)')
 1135 FORMAT (' STREAMLINES VISIBLE    (dark)')
 1137 FORMAT (' STREAMLINES DARK       (visible)')
 1138 FORMAT (' FILE   OFF             (filename)')
 1139 FORMAT (' FILE   ',A16,' (off)')
 1142 FORMAT (' SURFER OFF             (filename)')
 1143 FORMAT (' SURFER ',A16,' (off)')
 1145 FORMAT (' REPORT ON              (off)')
 1146 FORMAT (' REPORT OFF             (on)')
 1150 FORMAT (' LAYOUT ON              (off)')
 1160 FORMAT (' LAYOUT OFF             (on)')
 1162 FORMAT (' CONTOUR ON             (off)')
 1163 FORMAT (' CONTOUR OFF            (on)')
 1165 FORMAT (' CURSOR  ON             (off)')
 1167 FORMAT (' CURSOR  OFF            (on)')
 1168 FORMAT (' PICTURE ON             (off)')
 1169 FORMAT (' PICTURE OFF            (on)')
 1170 FORMAT (' <F2> or GO',/,' <Esc> or QUIT',/,' >')
 1500 FORMAT (80A1)
 1600 FORMAT (A1)
 1700 FORMAT ('         ------ MORE -------')
 1800 FORMAT (' Press <Enter> to return to menu.')
 1900 FORMAT (' ***** WARNING: points buffer will be cleared! *****',/,
     &' (CURSOR OFF will bring up cursor after tracing points',
     &' in buffer.)',/)
 1901 FORMAT (' Points buffer has been cleared.')
 1950 FORMAT (' ****WARNING: time is set to ',G11.4,/,' Do not perform '
     &'streamline tracing in the presence of transient wells!',/)
 2000 FORMAT (' ***ILLEGAL OR MISSING PARAMETER(S) in trace module:',/,
     &' ',80A1)
 2005 FORMAT (' ***ERROR: z-value (',G11.4,') outside aquifer;',/,
     &' base=',G11.4,' top=',G11.4,' Point not stored.',/)
 2006 FORMAT (' ***ERROR: z-value (',G11.4,') outside aquifer;',/,
     &' base=',G11.4,' top=',G11.4,' No pathlines traced from well:',
     &   A16,/)
 3000 FORMAT (' ***ILLEGAL COMMAND in trace module:',/,' ',80A1)
 4990 FORMAT (' ***BUFFER FULL: currently ',I4,' points stored.')
 5000 FORMAT (' Enter X,Y,Z and RDIRECTION of point ',I4,
     &' type QUIT to stop.'/
     &' >')
 7000 FORMAT (' Mount pen for STREAMLINES, then press LOCAL') 
 8000 FORMAT (' Point #        X               Y               Z',/)
 8005 FORMAT (' No starting points in buffer!')     
 8007 FORMAT (' *Date: ',A8,2X, 'Time: ',A11,/,
     &        ' *Title: ',A16,/,' message nul',/,' error error.log',/,
     &' yes',/,' echo con',/,' quit',/,' trace',/,' points')
 8010 FORMAT (I5,4X,3(F14.4,2X))
 8011 FORMAT (3(F14.4,2X)) 
 8012 FORMAT (3(F14.4,2X),' aborted by Ctrl-Break')  
 8015 FORMAT (' quit',/,'report off',/,' quit',/,' switch',/,
     &' error con',/,' message con',/,' echo off',/,' input con')
 9000 FORMAT (' *** ERROR: no grid or incompatible grid!',/,
     &' You may continue without a contour plot, else return to GRID',
     &' module.')
 9001 FORMAT (' ***ERROR: RETARDATION must be larger than zero.')
 9002 FORMAT (' ***ERROR: HALFLIFE must be larger than zero.')
 9003 FORMAT (' ***ERROR: increment must be between 0 and 1')
 9050 FORMAT (' You may now reset the device to the screen.'//)     
 9060 FORMAT (' ***ERROR in trace module:',/,' cursor cannot ',
     &'be on when device is set to printer or plotter!',/)
 9070 FORMAT(' NOTE: selected stepsize = ',g14.7,' recommended: '
     &,g14.7,/,' Press any key to continue.')
 9071 FORMAT(' WARNING: selected stepsize = ',g14.7,' recommended: '
     &,g14.7)          
 9080 FORMAT (' ***ERROR in trace module:',/,
     &        ' trace file already open, close or continue tracing.')
 9090 FORMAT (' ***ERROR in trace module:',/,
     &        ' No well found at ',E14.7,1X,E14.7)
 9091 FORMAT (' ***ERROR in trace module:',/,
     &        ' Well with label: ',A16,' has near zero discharge: ',
     &        E14.7,' No streamlines are being traced.')
 9092 FORMAT (' ***ERROR in trace module',/,
     &        ' Point buffer exceeded (max:',I5,') while executing:'
     &        ,/,A80)
 9100 FORMAT (' ***ERROR in trace module: point ',2(e11.4),
     &        ' is outside model domain.',/,' Point is ignored.')
      END
c
c ---------------------------------------------------------------------------------
c
      LOGICAL function lwellinfo (cz,czwell,rwell,rqwell,awell)
c
c ---------------------------------------------------------------------------------
c
c
c     Function is true when a 2D or 3D well is found at CZ.
c     If true well information is returned. The function is called by STREAM.
c
      implicit none
      LOGICAL lwlinfo,lw3info
      REAL(8) rwell,rqwell
      COMPLEX(8) cz,czwell
      CHARACTER(16) awell
c
      lwellinfo=.true.
      if (lwlinfo(cz,czwell,rwell,rqwell,awell)) return  ! search 2D wells
      if (lw3info(cz,czwell,rwell,rqwell,awell)) return  ! search 3D wells
      lwellinfo=.false.  ! no well at CZ
      return
      end
c
c ---------------------------------------------------------------------------------
c
      SUBROUTINE TRACE (RXST,NPTMAX,NST,RDS,REST,RTIC,RDTIME,LREPORT,
     & LSURFER,LDARK,RCLEVEL)
c
c ---------------------------------------------------------------------------------
c
C
C     Routine traces streamlines in plan view, when L3DVER=.FALSE.
C     The depth of the streamline is indicated with tick marks.
C     Routine traces streamlines in a vertical plane, when L3DVER=.TRUE.
C     NOTE: Streamlines should lie in the vertical plane!
C
C     RXST    starting points of streamlines
C     NPTMAX  maximum number of points in RXST array
C     NST     number of streamlines
C     RDS     specified stepsize (distance) along each streamline
C     REST    maximum residence time for each streamline
C     RTIC    tick mark interval for streamline depth (horizontal plots only).
C     RDTIME  tick mark interval for travel times
C     LREPORT .TRUE. for writing points and residence times along the 
C             streamlines. When .FALSE. only start and end points and times
C             are reported.
C     LSURFER Echoes coordinates of streamlines to a surfer boundary line file
C             opened in the STREAM subroutine.
C     LDARK   .TRUE. no streamlines are plotted, only tick marks for residence
C             times. Feature is useful when plotting isochrones of travel times
C             during wellhead protection studies.
C     RCLEVEL tick mark interval for concentrations along streamline
C
C
C
C     if the flag LTRACOUT = .TRUE. a pathline file is written for the GUI with the following format:
c
c
c    ---------- start documentation of *.pth file ------------ SEE ALSO "DOC" subfolder
c
c
c   Format of *.pth file generated by GFLOW:
c
c START  x1  y2  z1  t1  c1    iflag   iElementType
c        x1  y2  z1  t1  c1
c        x2  y2  z2  t2  c2
c        x3  y3  z3  t3  c3
c        ......
c        xn  yn  zn  tn  cn
c END    xn  xn  zn  tn  cn    iflag   iElementType  aElementLabel
c
c
c  where x1,y1,z1 is starting point of trace, t1 starting time and c1 starting concentration.
c  x1 .....xn lines are data for all successive points of the trace.
c
c
c
c     iFlag=-1 Error, iFlag not set.
c     iFlag=0  Requested time of travel is exceeded.
c     iFlag=1  Particle left window.
c     iFlag=2  Particle stopped at a boundary (see iElementType).
c     iFlag=3  Particle entered a region where the aquifer is dry.
c     iFlag=4  Particle is trapped in a stagnation point.
c     iFlag=5  Particle crossed aquifer top (entered layer above?).
c     iFlag=6  Particle crossed aquifer bottom (entered layer below?).
c     iFlag=7  Particle stopped by CTR-Break key.
c     iFlag=8  Particle stopped because of infinite velocity (singularity).
c
c     iElementType=-1 No element was reached.
c     iElementType=1  Particle terminated at a well.
c     iElementType=2  Particle terminated at a linesink.
c     iElementType=3  Particle terminated at a sink disc.
c     iElementType=4  Particle terminated at a partially penetrating well.
c     iElementType=5  Particle terminated at a Theis well.
c     iElementType=6  Particle terminated at a 3D sink disc.
c     iElementType=7  particle terminated at stream boundary (Simple WHPA in WhAEM)
c
c     aElementLabel="_NO_ELEMENT_" when iFlag is not equal to 2
c     aElementLabel="_NO_LABEL_"   when iFlag=2, but element has no label
c     aElementLabel= name of element when iFlag=2
c
C       -----------  end documentation if *.pth file -------------
C
C
      IMPLICIT NONE
      INTEGER(4) NPTMAX,NST,ISTATUSBAR,ISTART,IEND,ISTREAM,I,
     &           J,NSWITCH,ICOUNT,ISTAG,NWSTR,idum
      LOGICAL LREPORT,LSURFER,LDARK,L3DREV0,LCURTEMP,LWGEN,
     &        LTIME,LCTIME,LZTOP,LSHORT,LMAST,LTWL
      REAL(8) RXST,RDS,REST,RTIC,RDTIME,RCLEVEL,RTWTIME0,
     &        RXC,RYC,RZ,R1,R2,RFLAGD,RFLAGU,RBETA,RZTIC,
     &        RZTICU,RZTICL,RT,RCT,RTIME,RDUM,RCONCPLOT,RCTIME,
     &        RX,RY,RTOL,RTOP,RBASE,RDTWALL,RDT,RZOLD,RDIST,RZDIS,
     &        RTEST,RCDT,RTWT,RZ0,RESTIM,RXST0,RQI,RPI,RTWTIM,
     &        RFTOP,RFBASE
      COMPLEX(8) CZ,CZOLD,CZ0,CZCROSS,CI
      INCLUDE 'MATCH.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'DWMN.INC'
      INCLUDE 'LUSYS.INC'
      COMMON /MASTER/ RESTIM(5),ICOUNT,NSWITCH,LMAST  ! for MASTER project
      INTEGER(4) IDIRECT
      CHARACTER(1) ACHAR,ADIRECT,AESC
      CHARACTER(4) AVISIBLE
      CHARACTER(11) AZEL
      CHARACTER(16) AWELL
      DIMENSION RXST(4,*),RQI(3),RXST0(3)
      DATA ISTAG /0/
      DATA RPI,CI,NWSTR /3.141592654,(0.0,1.0),10/
      EXTERNAL VELOC
C
C
      AESC=CHAR(27)
      L3DREV0=L3DREV
      LCURTEMP=.FALSE.
      ISTATUSBAR=1
      LWGEN=.FALSE.
      IF (LCURS) THEN         ! cursor set or points from buffer?
        ISTART=0
        IEND=0
        NST=0
      ELSE
        ISTART=1
        IEND=NST
      ENDIF
      ISTREAM=0
      RDS0=RDS
      rtwtime0=rtwtim()
      LTIME=RDTIME.GT.0.0
      LCTIME=RCLEVEL.GT.0.0
      IF (L3DPL) THEN               ! calculate starting position cursor
      RXC=(R3DX1+R3DX2)*0.5
      RYC=(R3DY1+R3DY2)*0.5
      R3DZ=(R3DZ1+R3DZ2)*0.5
      RZ=R3DZ
      ELSE
      RXC=(RX1+RX2)*0.5
      RYC=(RY1+RY2)*0.5
      ENDIF
      LZTOP=.TRUE.
      IF (.NOT.LCURS.AND.NST.EQ.0) THEN  ! no cursor selected, but buffer empty
        ISTART=0
        IEND=0
        LCURTEMP=.TRUE.
        LCURS=.TRUE.
      ENDIF      
C
C ------------------------- START KEY OPTIONS --------------------  ! NOT AVAILABLE IN BATCH MODE
C
C                             IF (LGRAPHICS) THEN
C      CALL PLSETCUR (2,1,6)
C      WRITE (ILUME,5003)
C                              ENDIF
C 1010 CONTINUE
C                              IF (LGRAPHICS) THEN
C   1  IF (LDARK) AVISIBLE='dark'    ! set and print status bar
C      IF (.NOT.LDARK) AVISIBLE='vis.'
C      IF (L3DREV) ADIRECT='<'
C      IF (.NOT.L3DREV) ADIRECT='>'
C      IF (LZTOP) THEN
C      AZEL='aquifer top'
C      ELSE
C      WRITE (AZEL,1500) RZ
C      ENDIF
C      IF (LCURS) THEN
C      CALL PLSETCUR (1,1,1)
C      AESC=CHAR(27)
C      IATR=32
C      CALL PLSETCUR (1,1,1)
C      WRITE (ILUME,1600) AESC,IATR
C      IF (ISTATUSBAR.EQ.1)
C     & WRITE (ILUME,5000) AZEL,ADIRECT,AVISIBLE,RDTIME,RTIC
C      IF (ISTATUSBAR.EQ.2)
C     & WRITE (ILUME,4999) RETARDATION,RHALFLIFE,RCLEVEL,REST
C      IATR=37
C      WRITE (ILUME,1600) AESC,IATR
C 111  IF (LWGEN) GOTO 3
C        LDRAW=.TRUE.
C        IF (L3DPL) R3DZ=RZ
C        CALL CURSOR (ACHAR,RXC,RYC,0,0,2,LDRAW) ! get cursor position and key
C        CZ=CMPLX(RXC,RYC)
C        IF (L3DPL) RZ=R3DZ
C        IF (ACHAR.EQ.AESC) THEN           ! <Esc> leave trace routine
C          IF (LCURTEMP) LCURS=.FALSE.
C          CALL PLSETCUR (1,1,8)
C          WRITE (ILUME,4500)
C          GOTO 500
C        ENDIF
C        IF (ACHAR.EQ.'H'.OR.ACHAR.EQ.'h') THEN ! <H> set streamline elevation
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUME,7000)
C          CALL CURGET (RXC,RYC)
C          ISTATUSBAR=1
C          RDUM=RVAR(1)
C          IF (LERROR) THEN
C            LERROR=.FALSE.
C            LMISS=.FALSE.
C            LZTOP=.TRUE.
C            WRITE (ILUME,7500)
C            GOTO 1
C          ENDIF
C          RZ=RDUM
C          RBASE=RFBASE(CZ)
C          RZTOP=RFBASE(CZ)+RFHGHT(CZ)
C          IF (RZ.LT.RBASE.OR.RZ.GT.RZTOP) THEN
C            CALL PLSETCUR (5,1,4)
C            WRITE (ILUER,9000) RZ,RBASE,RZTOP
C            GOTO 111
C          ENDIF
C          LZTOP=.FALSE.
C          WRITE (ILUME,8000) RZ
C          GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'S'.OR.ACHAR.EQ.'s') GOTO 2 ! <S> start streamline
C        IF (ACHAR.EQ.'<') THEN                   ! < or > direction
C        ISTATUSBAR=1
C        CALL PLSETCUR (5,1,4)
C        L3DREV=.TRUE.
C        L3DREV0=L3DREV
C        GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'>') THEN
C        ISTATUSBAR=1
C        CALL PLSETCUR (5,1,4)
C        L3DREV=.FALSE.
C        L3DREV0=L3DREV
C        GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'D'.OR.ACHAR.EQ.'d') THEN     ! <D> dark streamline
C        ISTATUSBAR=1
C        CALL PLSETCUR (5,1,4)
C        LDARK=.TRUE.
C        GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'V'.OR.ACHAR.EQ.'v') THEN      ! <V> visible streamline
C        ISTATUSBAR=1
C        CALL PLSETCUR (5,1,4)
C        LDARK=.FALSE.
C        GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'/') THEN
C        CALL PLSETCUR (5,1,4)
C        IF (ISTATUSBAR.EQ.1) THEN
C        ISTATUSBAR=2
C        ELSE
C        ISTATUSBAR=1
C        ENDIF
C        GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'T'.OR.ACHAR.EQ.'t') THEN       ! <T> set time markers
C          ISTATUSBAR=1
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUME,7002)
C          CALL CURGET (RXC,RYC)
C          RDUM=RVAR(1)
C          IF (LERROR) THEN
C            LERROR=.FALSE.
C            LMISS=.FALSE.
C            RDTIME=0.0
C            LTIME=.FALSE.
C            WRITE (ILUME,7502)
C            GOTO 1
C          ENDIF
C          RDTIME=RDUM
C          IF (RDTIME.LE.0.0) THEN
C            RDTIME=0.0
C            LTIME=.FALSE.
C            WRITE (ILUME,7502)
C          ELSE
C            LTIME=.TRUE.
C          ENDIF
C          GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'Z'.OR.ACHAR.EQ.'z') THEN     ! <Z> set elevation markers
C          ISTATUSBAR=1
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUME,7004)
C          CALL CURGET (RXC,RYC)
C          RDUM=RVAR(1)
C          IF (LERROR) THEN
C            LERROR=.FALSE.
C            LMISS=.FALSE.
C            RTIC=0.0
C            WRITE (ILUME,7504)
C            GOTO 1
C          ENDIF
C          RTIC=RDUM
C          GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'C'.OR.ACHAR.EQ.'c') THEN ! <C> set concentration marker
C          ISTATUSBAR=2
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUME,7005)
C          CALL CURGET (RXC,RYC)
C          RDUM=RVAR(1)
C          IF (LERROR) THEN
C            LERROR=.FALSE.
C            LMISS=.FALSE.
C            RCLEVEL=0.0
C            LCTIME=.FALSE.
C            WRITE (ILUME,7505)
C            GOTO 1
C          ENDIF
C          IF (RDUM.LT.0.0.OR.RDUM.GT.1.0) THEN
C            WRITE (ILUER,7509)
C            GOTO 1
C          ENDIF
C          RCLEVEL=RDUM
C          IF (RCLEVEL.GT.0.0) THEN
C          LCTIME=.TRUE.
C          ELSE
C          LCTIME=.FALSE.
C          ENDIF
C          GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'E'.OR.ACHAR.EQ.'e') THEN ! <E> set endtime
C          ISTATUSBAR=2
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUME,7006)
C          CALL CURGET (RXC,RYC)
C          RDUM=RVAR(1)
C          IF (LERROR) THEN
C            LERROR=.FALSE.
C            LMISS=.FALSE.
C            REST=1.0E30
C            WRITE (ILUME,7506)
C            GOTO 1
C          ENDIF
C          REST=RDUM
C          GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'R'.OR.ACHAR.EQ.'r') THEN ! <R> retardation/halflife
C          ISTATUSBAR=2
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUME,7007) RETARDATION,RHALFLIFE
C          CALL CURGET (RXC,RYC)
C          RDUM1=RVAR(1)
C          RDUM2=RVAR(2)
C          IF (LERROR) THEN
C           LERROR=.FALSE.
C           LMISS=.FALSE.
C           WRITE (ILUER,6000)
C           GOTO 1
C          ENDIF
C          IF (RDUM1.LE.0.0) THEN
C            WRITE (ILUER,7507)
C            GOTO 1
C          ENDIF
C          IF (RHALFLIFE.LT.0.0) THEN
C            WRITE (ILUER,7508)
C            GOTO 1
C          ENDIF
C          RETARDATION=RDUM1
C          RHALFLIFE=RDUM2
C          GOTO 1
C        ENDIF
C        IF (ACHAR.EQ.'W'.OR.ACHAR.EQ.'w') THEN ! <W> generate NWSTR streamlines
CC                                              ! starting from nearest well
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUME,7010) NWSTR
C          CALL CURGET (RXC,RYC)
C          IDUM=IVAR(1)
C          IF (.NOT.LERROR.AND..NOT.LMISS) THEN
C            NWSTR=IDUM
C          ELSE
C            LERROR=.FALSE.
C            LMISS=.FALSE.
C          ENDIF
C        IF (.not.lwellinfo(cz,czwell,rwell,rqwell,awell)) THEN
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUER,6010)
C          GOTO 111
C        ENDIF
C        IF (RQWELL.LT.0.0.AND.L3DREV) THEN
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUER,6020) AWELL
C          GOTO 111
C        ENDIF
C        IF (RQWELL.GT.0.0.AND..NOT.L3DREV) THEN
C          CALL PLSETCUR (5,1,4)
C          WRITE (ILUER,6030) AWELL
C          GOTO 111
C        ENDIF
C        CALL DISCH (CZWELL,RQI)
C        RQI(3)=0.0
C        RQABS=RF3DSP(RQI,RQI)
C        IF (RQABS.LT.1.0E-10) THEN
C        RDTET=2.0*RPI/NWSTR
C        CDEL=2.0*CMPLX(RWELL,0.0)
C        ELSE
C        RTET1=50.0*0.18/NWSTR+0.02
C        RDTET=(6.2831853-RTET1)/(NWSTR-1)
C        CDEL=CMPLX(RQI(1),RQI(2))
C        RTET=AIMAG(LOG(CDEL))
C        RTET=RTET+0.5*RTET1
C        CDEL=2.0*RWELL*CEXP(CI*RTET)
C        ENDIF
C        LWGEN=.TRUE.
C        IWSTR=1
C        CZ=CZWELL+CDEL
C        RXC=REAL(CZ)
C        RYC=AIMAG(CZ)
C        GOTO 2
C      ENDIF
C  3   IF (LWGEN) THEN                     ! all streamlines for well completed
C        IF (IWSTR.EQ.NWSTR) THEN
C          LWGEN=.FALSE.
C          GOTO 111
C        ENDIF
C        IWSTR=IWSTR+1                     ! select next streamline for well
C        RTET=RTET+RDTET
C        CZ=CZWELL+2.0*RWELL*CEXP(CI*RTET)
C        RXC=REAL(CZ)
C        RYC=AIMAG(CZ)
C        GOTO 2
C      ENDIF
C        CALL PLSETCUR (5,1,4)
C        WRITE (ILUER,6000)                ! illegal key stroke
C        GOTO 111
C
C       ----------------------------- END OF KEY OPTIONS ---------------
C        
C   2    IF (NST.LT.NPTMAX) THEN
C          NST=NST+1                       ! add point to buffer
C          ISTART=ISTART+1
C          IEND=IEND+1
C        ELSE                        ! Buffer full, shift points downward and
C          DO 222 I=1,NPTMAX-1       ! add new point to last storage location
C          RXST(1,I)=RXST(1,I+1)     ! in the buffer. Hence, last NPTMAX points
C          RXST(2,I)=RXST(2,I+1)     ! will be retained in the buffer at all
C          RXST(3,I)=RXST(3,I+1)     ! times.
C          RXST(4,I)=RXST(4,I+1)
C          NST=NPTMAX
C          ISTART=NST
C          IEND=NST
C 222      CONTINUE
C        ENDIF
C        RXST(1,ISTART)=RXC          ! set starting point from cursor position
C        RXST(2,ISTART)=RYC
C        RXST(4,ISTART)=0.0
C        IF (LZTOP.AND..NOT.L3DPL) THEN
C          RXST(3,ISTART)=RFBASE(CZ)+RFHGHT(CZ)
C        ELSE
C          RXST(3,ISTART)=RZ
C        ENDIF
C      ENDIF
C      ENDIF        ! ----------- end of LGRAPHICS condition -------------------------
      R1=RX2-RX1                    ! calculate z-marker length and sign
      R2=RY2-RY1
      RFLAGD=MAX(R1,R2)/200.
      RFLAGU=-RFLAGD
      RBETA=RHALFLIFE/0.693147181
C
C     ----------------------------- LOOP OVER ALL STREAMLINES --------
C
      iticks_trace1=0
      iticks_trace2=0
      iticks_near1=0
      iticks_near2=0
      DO 10 I=ISTART,IEND
c      if (i.gt.1) call timer (iticks_trace2)
c      iticks_trace=iticks_trace2-iticks_trace1
c      if (iticks_trace.gt.0) then ! this will not report the last trace!
c       idum=iticks_trace-iticks_near
c       write (ilume,1001) inear_call,iticks_near,idum,iticks_trace
c 1001  format ('number of "near" subroutine calls=',i7,/,
c     &         'time for "near" subroutine calls=',i7,'E-2 seconds.',/,
c     &         'time for rest of tracing logic=  ',i7,'E-2 seconds.',/,
c     &         'total trace time=                ',i7,'E-2 seconds.')
c       iticks_trace1=0
c       iticks_trace2=0
c       iticks_near1=0
c       iticks_near2=0
c       iticks_trace=0
c       iticks_near=0
c       inear_call=0
c      end if
c      call timer (iticks_trace1)
      if (i.eq.0) GOTO 10   ! avoid trace when a GO command is given without any points specified. (6/1/2000)
      LCLOSE=.FALSE.
      ISTAG=0
c			initialize variables for BATCH mode
      iFlag=-1
      iElementType=-1
      aElementLabel="_NO_ELEMENT_"
c
      ISTREAM=ISTREAM+1
      RSTEP=RDS
      RZTIC=RXST(3,I)
      RZTICU=RZTIC+RTIC
      RZTICL=RZTIC-RTIC
      RT=0.0            ! groundwater residence time
      RCT=0.0           ! contaminant residence time
      RTIME=RDTIME      ! GW res. time at which tick mark is plotted
      if (ltwl()) then
      rdum=rtwtim()
      call time(rtwtime0)  ! reset time for transient wells
      end if
      RCONCPLOT=1.0-RCLEVEL
      IF (RCONCPLOT.GT.0.0) THEN
      RCTIME=-RBETA*LOG(RCONCPLOT) ! Cont. res. time at which tick mark plotted
      ELSE
      LCTIME=.FALSE.
      ENDIF
      L3DEND=.FALSE.
      CZ=CMPLX(RXST(1,I),RXST(2,I))       ! select starting point from array
      if (lboundary) then  ! set up for SBNEAR
       ctr0=ctr1-cz ! vector from starting point of trace to point 1 of boundary
       ctr0=ctr0/ctr2 ! as above, but divided by vector from point 1 to 2 on boundary
      end if
      R3DZ=RXST(3,I)
      RXST0(1)=RXST(1,I)
      RXST0(2)=RXST(2,I)
      RXST0(3)=RXST(3,I)
      IF (INT(RXST(4,I)+0.1).EQ.0) L3DREV=L3DREV0
      IF (INT(RXST(4,I)+0.1).EQ.1) L3DREV=.FALSE.
      IF (INT(RXST(4,I)-0.1).EQ.-1) L3DREV=.TRUE.
      LSHORT=.FALSE.
C      IF (LGRAPHICS) CALL PLSETCUR (5,1,4)           ! NOT AVAILABLE IN BATCH MODE
      IF (LSINGL) THEN
      WRITE (ILUME,1000) ISTREAM          ! print streamline start
      ELSE
      WRITE (ILUME,1000) I
      ENDIF
C      IF (LGRAPHICS) CALL PLSETCUR (6,1,2)           ! NOT AVAILABLE IN BATCH MODE
      WRITE (ILUME,2000) (RXST(J,I),J=1,3),RT
      RCONCTR=1.0
c
      IF (LTRACEOUT) THEN                     ! write "start record"   BATCH mode
      ILAYER=1
      IDIRECT=1
      IF (L3DREV) IDIRECT=-1
      WRITE (10,2100) (RXST0(J),J=1,3),
     &                    RT,RCONCTR,ILAYER,IDIRECT
      ENDIF
c
C      IF (LGRAPHICS) CALL PLSETCUR (7,1,1)           ! NOT AVAILABLE IN BATCH MODE
      IF (LUOUTFILE) THEN    ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      WRITE (ILUOUT,1000) I                 ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      WRITE (ILUOUT,2000) (RXST(J,I),J=1,3),RT ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      ENDIF                                 ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C      IF (.NOT.LDARK) CALL DRAW (CZ,0)       ! set pen at start of pathline   ! NOT AVAILABLE IN BATCH MODE
C      IF (LSURFER) THEN
C      LSURF=.TRUE.
C      IF (.NOT.LDARK) CALL DRAW (CZ,0)
C      LSURF=.FALSE.
C      ENDIF
c									BATCH mode
      IF (LTRACEOUT)WRITE (10,2200)           ! write starting "vertex record"
     &              (RXST(J,I),J=1,3),RT,RCONCTR
c
  5   CONTINUE                     ! start new pathline increment
      RX=REAL(CZ)                  ! check for starting point outside window
      RY=AIMAG(CZ)
      L3DEND=.TRUE.
      RTOL=0.1*RDS0
      IF (L3DPL) THEN
      iFlag=1
      IF (RX.LT.R3DX1-RTOL) GOTO 8
      IF (RX.GT.R3DX2+RTOL) GOTO 8
      IF (RY.LT.R3DY1-RTOL) GOTO 8
      IF (RY.GT.R3DY2+RTOL) GOTO 8
      iFlag=6
      IF (R3DZ.LT.R3DZ1) GOTO 8
      iFlag=5
      IF (R3DZ.GT.R3DZ2) GOTO 8
      iFlag=-1
      ELSE
      iFlag=1
      IF (RX.LT.RX1) GOTO 8
      IF (RX.GT.RX2) GOTO 8
      IF (RY.LT.RY1) GOTO 8
      IF (RY.GT.RY2) GOTO 8
      RTOP=RFTOP(CZ)
      RBASE=RFBASE(CZ)
      iFlag=5
      IF (R3DZ.GT.RTOP+RTOL) GOTO 8
      iFlag=6
      IF (R3DZ.LT.RBASE-RTOL) GOTO 8
      iFlag=-1
      ENDIF
      L3DEND=.FALSE.
      CZ0=CZ
      RZ0=R3DZ
C      IF (.NOT.LHPLOT.AND..NOT.LDARK) CALL DRAW (CZ,0)  ! NOT AVAILABLE IN BATCH MODE
  6   CALL PREDCOR (RDT,CZ,VELOC,RDS,RDTWALL,CZCROSS) ! generate new point (SET NEW STEP)
      IF (ISTAG.EQ.0) THEN    ! check for being at a stagnation point
        CZOLD=CZ
        RZOLD=R3DZ
        ENDIF
        ISTAG=ISTAG+1
        IF (ISTAG.EQ.25) THEN   ! compare current point with 25 points back
          ISTAG=0
          RDIST=ABS(CZ-CZOLD)
          RZDIS=ABS(R3DZ-RZOLD)
          RTEST=3.0*RSTEP
          IF (RDIST.LT.RTEST.AND.RZDIS.LT.RTEST) THEN  ! no progress,
             L3DEND=.TRUE.                       !  hung in stagnation point
             iFlag=4		
             ISTAG=0
             RT=1.0E20
C             IF (LGRAPHICS) CALL PLSETCUR (7,1,1)                        ! NOT AVAILABLE IN BATCH MODE
             if (lucon) WRITE (ILUME,4000) CZ,R3DZ,RT
           ENDIF
        ENDIF
C
      RT=RT+RDT                    ! increment groundwater residence time
      RCDT=RDT*RETARDATION
      RCT=RCT+RCDT                 ! increment contaminant residence time 
      IF (L3DEND) GOTO 8
      IF (LREPORT.and.lucon) WRITE (ILUME,4000) CZ,R3DZ,RT
      IF (RT.GT.REST) THEN         ! past endtime: abort the trace
        L3DEND=.TRUE.
        iFlag=0
        IF (RDTWALL.GT.0.0) THEN ! crossed a slurry wall with finite resistance, STOP trace
          RT=REST    ! select endtime as ending time (particle ended in the wall)
          CZ=CZCROSS ! set endpoint at wall, leave R3DZ as calculated (includes jump across wall)
          RCT=RT*RETARDATION
        ELSE         ! did not end in a slurry wall, finish trace with one last step.
          RT=RT-RDT
          RCT=RT*RETARDATION
          RDT=REST-RT
          IF (RDT.LE.0.0) GOTO 8
          CZ=CZ0
          R3DZ=RZ0
          CALL RUNKUT (RDT,CZ,VELOC,0.0)      ! make one last step to complete
        END IF
      ENDIF                               ! streamline up to endtime
C      IF (LTIME.AND.RT.GT.RTIME) THEN     ! draw T-tic marks if appropriate  ! NOT AVAILABLE IN BATCH MODE
C  7   CALL TMARK (CZ0,RZ0,CZ,RDT,RT,RTIME,LDARK,RXST0,LSHORT,.FALSE.)
C      RTIME=RTIME+RDTIME
C      IF (RT.GT.RTIME) GOTO 7
C      ENDIF
C      IF (LCTIME.AND.RCT.GT.RCTIME) THEN
C  71  IF (RCONCPLOT.LE.0.0) GOTO 72
C      CALL TMARK(CZ0,RZ0,CZ,RCDT,RCT,RCTIME,.TRUE.,RXST0,.TRUE.,.TRUE.)
C      RCONCPLOT=RCONCPLOT-RCLEVEL
C      IF (RCONCPLOT.LE.0.0) GOTO 72
C      RCTIME=-RBETA*LOG(RCONCPLOT)
C      IF (RCT.GT.RCTIME) GOTO 71
C      ENDIF
C  72  CONTINUE
C      IF (.NOT.L3DVER.AND.RTIC.GT.0.0.AND..NOT.LDARK) THEN ! draw Z-tic marks
C      LDOWN=R3DZ.LT.RZTIC
C      LUP=R3DZ.GE.RZTIC
C      IF (LDOWN.AND.R3DZ.LE.RZTICL) THEN
C      CALL ZMARK (CZ0,CZ,RZ0,RZTICL,RTIC,RFLAGD)      ! set downward tic mark
C      RZTIC=RZTICL+RTIC
C      RZTICU=RZTIC+RTIC
C      ENDIF
C      IF (LUP.AND.R3DZ.GE.RZTICU) THEN
C      CALL ZMARK (CZ0,CZ,RZ0,RZTICU,RTIC,RFLAGU)      ! set upward tic mark
C      RZTIC=RZTICU-RTIC
C      RZTICL=RZTIC-RTIC
C      ENDIF
C      ENDIF
c
      IF (LTRACEOUT)  ! write "vertex record"                 BATCH mode
     &        WRITE (10,2200) CZ,R3DZ,RT,RCONCTR
c                                                                ! NOT AVAILABLE IN BATCH MODE
C      IF (.NOT.LDARK) CALL DRAW (CZ,1)                          ! draw next pathline section
C      IF (LSURFER) THEN
C      LSURF=.TRUE.
C      IF (.NOT.LDARK) CALL DRAW (CZ,1)
C      LSURF=.FALSE.
C      ENDIF
      RSTEP=RDS0
c                  rtwt      is the absolute time for the transient wells
c                  rtwtime0  is the starting time for the trace (0, unless set with the TIME command)
c                  RT        is the residence time of the particle
      if (L3DREV) then
        rtwt=rtwtime0-RT   ! tracing back in time
      else
        rtwt=rtwtime0+RT   ! tracing froward in time
      endif
      rtwt=MAX(rtwt,0.0)   ! prevent negative time for transient wells
      call time (rtwt)                    ! update time for transient wells
      IF (.NOT.L3DEND) GOTO 5             ! continue with next pathline increment
C
C     ---------------------- END OF STREAMLINE -----------------------
C      
  8   CONTINUE
c
      IF (LTRACEOUT) THEN  !						BATCH mode
      WRITE (10,2200) CZ,R3DZ,RT,RCONCTR  ! write ending "vertex record"  
      WRITE (10,2300) CZ,R3DZ,         ! write the "end record" 
     &     RT,RCONCTR,iFlag,iElementType,aElementLabel
      ENDIF
c                                                   ! NOT AVAILABLE IN BATCH MODE
C      IF (.NOT.LDARK) CALL DRAW (CZ,1)
C      IF (LGRAPHICS) CALL GFPLOT (-2,0,0)
C      IF (LGRAPHICS) CALL PLSETCUR (7,1,1)
c      IF (LSINGL) THEN
      WRITE (ILUME,3000) CZ,R3DZ,RT
c      write (iluer,3001) RT           ! special request Steve
c      ENDIF
      IF (LUOUTFILE) WRITE (ILUOUT,3000) CZ,R3DZ,RT
C      IF (LUOUTFILE.AND.LMAST.AND.LWSIN(CZ)) THEN  ! MASTER project MMMMMMMM
C       IF (NSWITCH.GT.0) THEN         ! test to see if streamline ends at well
C        IF (ABS(CZ).GT.10.0) GOTO 10  ! well at (0,0)
C       ENDIF                          ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C       IF (IABS(NSWITCH).GT.0) THEN ! write residence starting point and time M
C         WRITE (ILUOUT,1004) (RXST(JJ,1),JJ=1,3),RT
C 1004    FORMAT (4(G14.7,2X))
C        IF (IABS(NSWITCH).EQ.2) THEN
C         LDARK=.FALSE.
C         CZ0=CMPLX(RXST(1,1),RXST(2,1))
C        IF (LGRAPHICS) CALL DRAW (CZ0,0)           ! NOT AVAILABLE IN BATCH MODE
C        IF (LGRAPHICS) CALL DRAW (CZ0,1)           ! NOT AVAILABLE IN BATCH MODE
C         LDARK=.TRUE.
C        ENDIF
C       ENDIF
C      ENDIF                          ! end code for MASTER project MMMMMMMMMM
      IF (LREPORT.AND..NOT.LCURS) THEN
        if (lucon) WRITE (ILUME,5002)
      ENDIF
   10 CONTINUE
C
C     ----------------------- END OF STREAMLINE TRACING LOOP ------------
C 
C                                   IF (LGRAPHICS)  THEN            ! NOT AVAILABLE IN BATCH MODE
C      IF (LCURS) THEN
C      GOTO 1010
C      ELSE
C      IF (LMAST) then
C        call time (rtwtime0)
C        RETURN
C      endif
C  15   IF (LSCRN) THEN        ! now activate temporary cursor
C        ISTART=NST            ! for adding more streamlines
C        IEND=NST
C        LCURTEMP=.TRUE.
C        LCURS=.TRUE.
C        LWGEN=.FALSE.
C        CALL PLSETCUR (2,1,7)
C        WRITE (ILUME,5001)
C        GOTO 1010
C       ELSE
C        call time (rtwtime0)
C        RETURN
C       ENDIF
C      ENDIF
C                                   ENDIF
C      
  500 IF (LSURFER) THEN
      LSURF=.TRUE.
      CZ=(0.0,0.0)
C      IF (.NOT.LDARK) CALL DRAW (CZ,-2)            ! NOT AVAILABLE IN BATCH MODE
      LSURF=.FALSE.
      ENDIF
      call time (rtwtime0)
      RETURN
C      
 1000 FORMAT ('+Streamline',I3)
 1500 FORMAT (E14.7)
 1600 FORMAT ('+',A1,'[',I2,'m')       
 2000 FORMAT ('+start xyz ',3(E14.7,1X),' start time ',E14.7)
 2100 FORMAT ('START',1X,5(E13.6,1X),2I10)
 2200 FORMAT (6X,5(E13.6,1X))
 2300 FORMAT ('END',3X,5(E13.6,1X),2I10,1X,A16)
 3000 FORMAT ('+end   xyz ',3(E14.7,1X),' end   time ',E14.7)
 3001 format (d14.7)
 3050 FORMAT ('+End of streamline.')
 4000 FORMAT ('+xyz ',3(E14.7,1X),'   residence time ',E14.7)
 4500 FORMAT ('+                                                     ',
     &'                          ')
C
 4999 FORMAT ('+Retard.=',G11.4,' Halfl.=',G11.4,' Cmark=',G11.4,
     &' Endtime=',G11.4)      
 5000 FORMAT ('+H=',A11,'  Direct.=',A1,'  Str.line=',A4,
     &'  Tmark=',G11.4,' Zmark=',G11.4)
 5001 FORMAT('+KEYS: S=1 streamline, W=# streamlines from well,',
     &' H=elevation start pnt,',/,
     &' T=Tmark, Z=Zmark, </>=backward/forward, D/V=dark/visible,',/,
     &' Ctrl-Break=cancel tracing, Esc =exit, edit graph.')
 5002 FORMAT (' Press any key to continue.')
 5003 FORMAT('+KEYS: S=1 streamline, W=# streamlines from well,',
     &' H=elevation start pnt,',/,
     &' T=Tmark, Z=Zmark, C=Cmark, E=Endtime, </>=backw./forw.,',
     &' D/V=dark/vis.,',/,
     &' R=retard./halfl.,  / =status, Ctrl-Break =end trace,',
     &' Esc =edit graph.')
 6000 FORMAT ('+***ILLEGAL KEY STROKE: try again.')
 6010 FORMAT ('+***ERROR: no (steady state) wells available.')
 6020 FORMAT ('+***ERROR: well ',A16,' is an injection well.',/,
     &' Cannot trace backward.')
 6030 FORMAT ('+***ERROR: well ',A16,' is a pumping well.',/,
     &' Cannot trace forward.')
 7500 FORMAT ('+Starting elevation set at water table or aquifer top,',
     &' until new H-key stroke.    ')
 7502 FORMAT('+No markers for travel time.                            ')
 7504 FORMAT('+No tick marks for streamline depth.                    ')
 7505 format('+No tick marks for concentrations.                      ')
 7506 FORMAT('+Endtime is set to "infinity".                          ')
 7507 FORMAT('+***ERROR: Retardation must be larger than 0.0   ')
 7508 FORMAT('+***ERROR: Halflife must be larger or equal to 0.0      ')
 7509 FORMAT('+***ERROR: Increment must be between 0.0 and 1.0        ')
 7000 FORMAT(' Enter z-value (press <Enter> for water table or aquifer',
     &' top) >')
 7002 FORMAT (' Enter time interval for travel time markers. >')
 7004 FORMAT(' Enter depth increment for streamline depth tick marks.>')
 7005 FORMAT(' Enter concentr. incr. (0-1) for concentr. tick marks. >')
 7006 FORMAT (' Enter endtime (Press <Enter> for "infinity".) >')
 7007 FORMAT (' Retard.=',G11.4,' halfl.=',G11.4,' Enter new values >')
 7010 FORMAT (' Number of streamlines is ',I3,'. Enter new number, or ',
     &'press <Enter> >')     
 8000 FORMAT ('+Starting elevation set at ',G14.7,' until new H-key ',
     &'stroke.      ') 
 9000 FORMAT ('+***ERROR: RZ (',E14.7,') outside aquifer:',/
     &' base is at ',E14.7,' and top is at ',E14.7)
      END
c
c ---------------------------------------------------------------------------------
c
      SUBROUTINE PREDCOR (RDT,CZ,VELOC,RDS,RDTWALL,CZCROSS)
c
c ---------------------------------------------------------------------------------
c
C
C     PREDICTOR CORRECTOR method for integrating velocities in 3-D.
C
C Input:
C     CZ      current location (x,y) along the streamline (R3DZ contains
C             current z- location)
C     VELOC   external velocity function
C     RSTEP   in common /grid/ is the step size [L],
C             it may be adjusted by a specific discharge function.
C     RDS     originally specified step size
C
C Output:
C     RDTWALL residence time inside the wall, only larger than zero when crossing a
C             slurry wall with a finite resistance.
C     CZCROSS intersection of trace with slurry wall.
C     RDT     time step for current step RSTEP [L] along pathline
C     CZ      new location (x,y) along the pathline (R3DZ contains new
C             z- location).
C     RSTEP   may have been adjusted near stagnation point or near inhomogeneity
C     L3DEND  may be set .TRUE. in case a stagnation point is encountered.
C
C     Currently no measures to avoid integration across inhomogeneity 
C     boundaries!
C
C
      IMPLICIT NONE
      LOGICAL LREDO,LEULER,lfinterface
      REAL(8) RXI1,RXI2,RVI1,RVI2,RV1,RV2,RZ0,RDT,RDS,K0,B0,H0,
     &     RTEST,RHGHT,RFHGHT,RFPER,RFBASE,RBASE,RZNEW,RZEULER,
     &     RF3DSP,RDELTA,RDTWALL,SQROOT,K1,B1,RZMIN,RZMAX,RFPERM,
     &     rfinterface,rbloc,z0,z1,deltaz
      COMPLEX(8) CZ,CZNEW,CZEULER,CZCROSS,cz1
      INCLUDE 'GRID.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION RXI1(3),RXI2(3),RVI1(3),RVI2(3)
      RSTEP=RDS
      RDTWALL=0.0
      CZCROSS=(1.0E+21,1.0E+21) ! flag for CZCROSS not set in DBNEAR
      RXI1(1)=REAL(CZ)
      RXI1(2)=AIMAG(CZ)
      RXI1(3)=R3DZ
      RZ0=R3DZ
      RDT=0.0
      k0=rfperm(cz)
      b0=rfbase(cz)
      h0=rfhght(cz)
      if (lfinterface()) z0=rfinterface(cz)
      CALL VELOC (CZ,RVI1)           ! calculate first velocity
      RV1=RF3DSP(RVI1,RVI1)
      IF (RV1.LT.1.0E-12*K0) THEN       ! zero velocity, abort
      WRITE (ILUER,1000)
      L3DEND=.TRUE.
      iFlag=4
      RETURN
      ENDIF
      IF (RV1.GT.1.E15) THEN         ! infinite velocity, abort
      WRITE (ILUER,2000)
      L3DEND=.TRUE.
      iFlag=8
      RETURN
      ENDIF
      RV1=SQROOT(RV1)
      IF (LCLOSE) THEN
        RV1=RF3DSP(RVI1,RDIRI)
        RVI1(1)=RV1*RDIRI(1)
        RVI1(2)=RV1*RDIRI(2)
      ENDIF
  5   lredo=.false.
      RDELTA=RSTEP
      IF (L3DREV) RDELTA=-RSTEP
      RXI2(1)=RXI1(1)+RVI1(1)*RDELTA/RV1  ! first estimate of new point
      RXI2(2)=RXI1(2)+RVI1(2)*RDELTA/RV1
      RXI2(3)=RXI1(3)+RVI1(3)*RDELTA/RV1
      CZNEW=CMPLX(RXI2(1),RXI2(2))
      R3DZ=RXI2(3)
      rznew=r3dz
      CZEULER=CZNEW
      RZEULER=RZNEW
      k1=rfperm(cznew)
      b1=rfbase(cznew)
      if (k0.ne.k1.or.b0.ne.b1) then    ! crossing a transmissivity inhomogeneity use Euler
        rvi2(1)=rvi1(1)
        rvi2(2)=rvi1(2)
        rvi2(3)=rvi1(3)
        goto 15
      endif
      CALL VELOC (CZNEW,RVI2)                 ! calculate second velocity
        RV2=RF3DSP(RVI2,RVI2)
        IF (RV2.LT.1.0E-12*K0) THEN       ! zero velocity, abort
        WRITE (ILUER,1000)
        L3DEND=.TRUE.
        iFlag=4
        RETURN
        ENDIF
        IF (RV2.GT.1.E15) THEN         ! infinite velocity, abort
          WRITE (ILUER,2000)
          L3DEND=.TRUE.
          iFlag=8
          RETURN
        ENDIF
      rtest=rf3dsp(rvi1,rvi2)
      if (rtest.lt.0.0) then !  opposite velocities, apply euler' method
        rvi2(1)=rvi1(1)
        rvi2(2)=rvi1(2)
        rvi2(3)=rvi1(3)
      else                    ! apply predictor corrector method        
        RV2=SQROOT(RV2)
        IF (LCLOSE) THEN
          RV2=RF3DSP(RDIRI,RVI2)
          RVI2(1)=RV2*RDIRI(1)
          RVI2(2)=RV2*RDIRI(2)
        ENDIF
        RVI2(1)=0.5*(RVI1(1)+RVI2(1))  ! calculate final velocity
        RVI2(2)=0.5*(RVI1(2)+RVI2(2))
        RVI2(3)=0.5*(RVI1(3)+RVI2(3))
      endif
  15  RV2=RF3DSP(RVI2,RVI2)
      RV2=SQROOT(RV2)
      RDELTA=RSTEP
      IF (L3DREV) RDELTA=-RSTEP
      RXI2(1)=RXI1(1)+RVI2(1)*RDELTA/RV2  ! calculate new point
      RXI2(2)=RXI1(2)+RVI2(2)*RDELTA/RV2
      RXI2(3)=RXI1(3)+RVI2(3)*RDELTA/RV2
      CZNEW=CMPLX(RXI2(1),RXI2(2))
      RZNEW=RXI2(3)
C     ---------------------- adjust new point or end streamline
C                            when near an analytic element
c      inear_call=inear_call+1
c      call timer(iticks_near1)
c      GOTO 666 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BLOCKED FOR TIMING!!!!!!
      CALL WLNEAR (CZ,CZNEW,RZ0,RZNEW,LREDO) ! wells
      CALL LSNEAR (CZ,CZNEW,RZ0,RZNEW,CZEULER,RZEULER,LEULER,LREDO) ! line sinks
      CALL DBNEAR (CZ,CZNEW,RZ0,RZNEW,RDTWALL,CZCROSS,LREDO) ! slurry walls
      CALL W3NEAR (CZ,CZNEW,RZ0,RZNEW,LREDO) ! ppwells
      CALL TWNEAR (CZ,CZNEW,RZ0,RZNEW,LREDO) ! Theis wells
      CALL DINEAR (CZ,CZNEW,RZ0,RZNEW,LREDO) ! sinkdiscs (3D)
      CALL PDNEAR (CZ,CZNEW,RZ0,RZNEW,LREDO) ! sinkdiscs (2D)
      CALL SBNEAR (CZ,CZNEW,RZ0,RZNEW,LREDO) ! stream boundary for Simple WHPA in WhAEM  (routine at end of this file)
c 666  continue
c      call timer(iticks_near2)                          ! temporary timer logic
c      iticks_near=iticks_near+iticks_near2-iticks_near1 ! temporary timer logic
      IF (LREDO) GOTO 5
      RHGHT=RFHGHT(CZNEW)  ! don't let streamline get out of aquifer
      RBASE=RFBASE(CZNEW)
      if (lfinterface()) THEN ! interface flow
        z1=rfinterface(cznew)
        rbloc=MAX(rbase,z1)
        RZMAX=Rbloc+RHGHT
        RZMIN=Rbloc-0.001*rhght
      else
        RZMAX=RBASE+RHGHT
        RZMIN=RBASE-0.001*rhght
      endif
      if (rbase.ne.b0) THEN ! crossed jump in aquifer base, adjust pathline elevation
      if (lfinterface()) then
        b0=MAX(b0,z0)
        b1=MAX(rbase,z1)
      end if
      deltaz=rvi2(3)*rdelta/rv2
      rznew=(rxi1(3)-b0)*rhght/h0+b1+deltaz  ! this logic is slightly inaccurate since b0 and b1
c                                                          are not calculated exactly at the boundary. For instance,
c           if the interface is above both aquifer bases still a small jump in aquifer base is seen.
c      write (iluer,1001) rxi1(3),b0,b1,rhght,h0,deltaz,rznew
c 1001 format (' predcor1: rxi1(3),b0,b1,rhght,h0,deltaz,rznew '
c     &          ,/,7(d14.7))
      end if
      IF (RZNEW.GT.RZMAX) THEN
      L3DEND=.TRUE.
      iFlag=5
      ENDIF
      IF (RZNEW.LT.RZMIN) THEN
      L3DEND=.TRUE.
      iFlag=6
      ENDIF
C
      IF (LEULER) THEN
c      rv1=  ! for 50:1 anisotropy ratio
c     &SQRT(50.0d0*(rvi1(1)*rvi1(1)+rvi1(2)*rvi1(2))+rvi1(3)*rvi1(3))
        RDT=RSTEP/RV1+RDTWALL
      ELSE
c      rv2=   ! for 50:1 anisotropy ratio
c     &SQRT(50.0d0*(rvi2(1)*rvi2(1)+rvi2(2)*rvi2(2))+rvi2(3)*rvi2(3))
        RDT=RSTEP/RV2+RDTWALL
      ENDIF
      CZ=CZNEW
      R3DZ=RZNEW
      RETURN
 1000 FORMAT ('***ZERO VELOCITY in pathline tracing.                  ')
 2000 FORMAT ('***INFINITE VELOCITY in pathline tracing.              ')
      END
c
c ---------------------------------------------------------------------------------
c
      SUBROUTINE RUNKUT (RDT,CZ,VELOC,RDELTA)
c
c ---------------------------------------------------------------------------------
c
C
C     RUNGE-KUTTA method of integrating velocities in 3-D.
C     If RDELTA is larger than zero, the time step will be
C     adjusted to approximately match the path length step size
C     RDELTA.
C
C     RDT    time step used in runga kutta calculation
C     CZ     current location (x,y) along the streamline (R3DZ contains
C            current z- location)
C     VELOC  external velocity function
C     RDELTA switch to 1) RDELTA>0 adjust RDT to stepsize RSTEP
C                      2) RDELTA=0 use specified timestep RDT
C     Note: RSTEP may be adjusted by a specific discharge function, however,
C     values below 0.01 times the specified value (RDELTA) will be overruled.
C
C     Note: step size means distance along the streamline [L].
C
C     AVOIDING INTEGRATION ACROSS INHOMOGENEITIES (DISCONTINUOUS VELOCITY)
C     Function checks if permeability of porosity changes when calculating
C     the velocity at the four different points (4th order method).
C     If so, the point in question, which is apparently inside an inhomogeneity
C     is returned as the end point of the step. The associated time step
C     is calculated based on Euler's Method for this step.
C
C     OVERESTIMATION OF TIME STEP
C     If the actual step size appears much larger (more than twice) than the
C     intended step size, the integration is repeated with half the residence
C     time. This correction is only made once per subroutine call!
C
C
      IMPLICIT NONE
      INTEGER(4) I,J
      LOGICAL LFIRST
      REAL(8) RDT,RDELTA,RFAC,RSUM,RK,RVI,RDUM,RFPERM,RDELT,
     &        R2,RF3DSP,R,SQROOT,RZ0,RPERM0,RPERM,RHDUM,RVDUM
      COMPLEX(8) CZ,CZ0
      INCLUDE 'GRID.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION RVI(3),RK(3),RSUM(3),RFAC(4)
      SAVE
      DATA RFAC /1.0,2.0,2.0,1.0/
      LFIRST=.TRUE.           ! used to invoke step size test at end of routine
      IF (RDELTA.GT.0.0) THEN
        CALL VELOC (CZ,RVI)           ! estimate RDT from Euler's Method
        RDUM=RFPERM(CZ)
        IF (L3DEND) RETURN
        RSTEP=MAX(RSTEP,0.01*RDELTA)  ! do not tolerate a step less than
        RDELT=RSTEP                     ! 0.01 of the specified step
        R2=RF3DSP (RVI,RVI)
        R=SQROOT(R2)
        IF (R.LT.1.E-10) THEN
C          CALL PLOTOF          ! NOT AVAILABLE IN BATCH MODE
          WRITE (ILUER,1000)
          L3DEND=.TRUE.
          RETURN
        ENDIF
        IF (R.GT.1.E20) THEN
C          CALL PLOTOF         ! NOT AVAILABLE IN BATCH MODE
          WRITE (ILUER,2000)
          L3DEND=.TRUE.
          RETURN
        ENDIF
        RDT=RDELT/R
      ENDIF
      CZ0=CZ
      RZ0=R3DZ
  8   RSUM(1)=0.0
      RSUM(2)=0.0
      RSUM(3)=0.0
      RPERM0=RFPERM(CZ)
      DO 20 I=1,4
      CALL VELOC (CZ,RVI)
      RPERM=RFPERM(CZ)
      IF (RDELTA.GT.0.0.AND.RPERM.NE.RPERM0) THEN ! crossed inhomogeneity
        RDELT=ABS(CZ0-CZ)        ! use current point as endpoint of this step
        RDELT=SQROOT(RDELT*RDELT+(R3DZ-RZ0)*(R3DZ-RZ0))
        RDT=RDELT/R               ! calculate timestep based on Euler's method
      RETURN
      ENDIF
      DO 10 J=1,3
      RK(J)=RDT*RVI(J)
      IF (L3DREV) RK(J)=-RK(J)
      RSUM(J)=RSUM(J)+RFAC(I)*RK(J)
  10  CONTINUE
      CZ=CZ0+0.5*CMPLX(RK(1),RK(2))
      R3DZ=RZ0+0.5*RK(3)
  20  CONTINUE
      CZ=CZ0+CMPLX(RSUM(1),RSUM(2))/6.0
      R3DZ=RZ0+RSUM(3)/6.0
      IF (.NOT.LFIRST) RETURN
      LFIRST=.FALSE.
      RHDUM=ABS(CZ0-CZ)                  ! test for RDS too large
      RVDUM=RZ0-R3DZ
      RDUM=SQROOT(RHDUM*RHDUM+RVDUM*RVDUM)
      IF (RDUM.GT.2.0*RDELT) THEN         ! redo with reduced time step
      RDT=0.5*RDT
      CZ=CZ0
      R3DZ=RZ0
      GOTO 8
      ENDIF
      RETURN
 1000 FORMAT ('***ZERO VELOCITY IN <RUNKUT>')
 2000 FORMAT ('***INFINITE VELOCITY IN <RUNKUT>')
      END
c
c ---------------------------------------------------------------------------------
c
      SUBROUTINE GETSTEP (RSTP0,RSTP,LEND,LREV)
c
c ---------------------------------------------------------------------------------
c
C
C     Routine is called by specific discharge functions to terminate
C     streamline plotting when at singularity.
C
C     RSTP   stepsize along streamline given to spec. disch. function
C     LEND   used to pass L3DEND to specific discharge function
C     LREV   used to pass L3DREV to specific discharge function
C
      IMPLICIT NONE
      LOGICAL LEND,LREV
      REAL(8) RSTP0,RSTP,RZ
      INCLUDE 'GRID.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'LUSYS.INC'
      RSTP0=RDS0
      RSTP=RSTEP
      LEND=L3DEND
      LREV=L3DREV
      RETURN
C     -------------
      ENTRY SETSTEP (RSTP,LEND)
C     -------------
C
C     RSTP  used to pass new stepsize to RSTEP (changed when LEND set TRUE)
C     LEND  used to update L3DEND by spec. disch. function
C
      RSTEP=RSTP
      L3DEND=LEND
      RETURN
C     -------------
      ENTRY ZVAL(RZ)
C     -------------      
C
C     Returns value of R3DZ            
C
      RZ=R3DZ
      RETURN
      END
C
C --------------------------------------------------------------------
C
      SUBROUTINE SBNEAR (CZ,CZNEW,RZ0,RZNEW,LREDO)
c
c --------------------------------------------------------------------
C
C     Routine terminates the particle trace if CZ and CZNEW on different sides of
C     infinitely long line
c
c     Procedure: (1) calculate vector from first point on line to CZ and
c                    calculate dot product with vector from first point to second point
c                (2) calculate vector from first point on line to CZNEW and
c                    calculate dot product with vector from first point to second point
c                (3) check to see if dot products differ in sign (means that particle crossed line)
c                    note: if one of the dot products is zero, still stop, since particle is on the line
c
C
      IMPLICIT NONE
      LOGICAL LREDO,L3DEND,L3DREV
      REAL(8) RZ0,RZNEW,RSD0,RSTEP,rtest
      COMPLEX(8) CZ,CZNEW
      INCLUDE 'TRACOM.INC'
      include 'lusys.inc'
C
      IF (LBOUNDARY) THEN
        call getstep (rsd0,rstep,l3dend,l3drev)
        IF (L3DEND) RETURN  ! done, pathline ended elsewhere
        ctrnew=ctr1-cznew
        ctrnew=ctrnew/ctr2
        rtest=AIMAG(ctr0)*AIMAG(ctrnew)
c        write (ilume,1001) ctr1,ctr2,ctr0,ctrnew,rtest
c 1001 format ('sbnear1: ctr1,ctr2,ctr0,ctrnew,rtest ',
c     &         4(d11.4,1x),1x,d11.4)
        IF (rtest.LE.0.0d0) THEN ! crossing boundary stop pathline trace
          L3DEND=.TRUE.
          CALL SETSTEP (RSTEP,L3DEND)
          iFlag=2
          iElementType=7
          aElementLabel='_STREAM_BOUNDARY_'
          write (ilume,1002) l3dend
 1002 format ('sbnear2: l3dend ',2x,l3)
        END IF
      ENDIF
      return
      END SUBROUTINE

