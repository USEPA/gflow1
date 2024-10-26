C     Last change:  HMH   6 Jan 2001   11:20 am
c       This file contains the following routines and functions:
c
C	BLOCK DATA TWDAT    initializes transient well common blocks
C	SUBROUTINE TWIN     handles input of transient wells as defined in the .dat file
C
c -------------------
c no   twmat.for   file, transient wells are superimposed on a steady state solution
c -------------------
c
c -------------------
c twfun.for contains:
c -------------------
c
C     REAL FUNCTION RFTWPT    contribution to the discharge potential due to all transient wells
C     SUBROUTINE TWQI         returns the discharge vector components due to all transient wells
c
c ---------------------
c no  twnflow.for   file, transient wells should not be used close to any boundary conditions
c ---------------------
c
c
c ------------------
c twio.for contains:
c ------------------
c
C     SUBROUTINE twextract        writes data for transient wells to .xtr file
C     SUBROUTINE TWDATIO          writes all transient well data to a ".dat" file
C     SUBROUTINE TWIO             reads or writes contents of TWCOM.INC from or to .sol file
C     SUBROUTINE TWOUT            writes contents of TWCOM formatted to ILUOUT for debugging purposes
c
c ---------------------
c no    twcheck.for    file, transient wells are not part of the (steady state) solution procedure
c ---------------------
c
c
c --------------------
c twnear.for contains:
c --------------------
c
C     SUBROUTINE TWNEAR  check for streamline trace near well and slows down trace or ends trace
C
c --------------------
c twserv.for contains:
c --------------------
c
C     SUBROUTINE TWTIME      input routine for current time
C     SUBROUTINE TIME        called by TWTIME
C     REAL FUNCTION RTWTIM   reports the current time setting
C     LOGICAL FUNCTION LTWL  TRUE if transient wells are present and time is larger than zero
c
c ---------------------
c twextra.for contains:
c ---------------------
c
C     SUBROUTINE TWPLOT     plots transient wells on analytic element layout
C
c
C
C -------------------------------------------------------------------------
C
      BLOCK DATA TWDAT
C
C -------------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'twcom.inc'
      INCLUDE 'lusys.inc'
      DATA RPI4,NTW /12.56637061,0/
      DATA RTWH0,RTWQ /NTWMX*0.0,NTWMX*0.0 /
      data RTWST / NTWMX*0.0 /
      DATA RTWR0,RTWSTO /NTWMX*0.2,NTWMX*0.2/
      DATA CTWZ /NTWMX*(0.0,0.0)/
      DATA RTWRD0,RTWCT,RTWST0 /0.2,0.0,0.2/
      END
C
C -------------------------------------------------------------------------
C
      SUBROUTINE TWIN
C
C -------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) JUMP
      REAL(8) RDUM,RVAR,RDUM1,RDUM2,RDUM3,RDUM4,RDUM5
      LOGICAL LENTRY,LBAD
      COMPLEX(8) CDUM,CVAR
      INCLUDE 'gvcom.inc'      
      INCLUDE 'twcom.inc'      
      INCLUDE 'match.inc'
      INCLUDE 'lusys.inc'
      CHARACTER(1) AWORD(23)
      DATA AWORD /'?',' ',
     &            'S','T','O','R',' ',
     &            'R','A','D','I',' ',
     &            'D','I','S','C',' ',
     &            'Q','U','I','T',' ',
     &            ATERM/
      LENTRY=.FALSE.
      LERROR=.FALSE.
      LMISS=.FALSE.
C      CALL CLEARSCREEN    ! NOT AVAILABLE IN BATCH MODE
  10  IF (LENTRY) GOTO 400
      IF (LERROR) WRITE (ILUER,2000) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      if (lucon) WRITE (ILUME,1000) NTWMX,RTWST0,RTWRD0
      CALL INLINE
  11  CALL MATCH (AWORD,1,JUMP,LBAD)
C      CALL CLEARSCREEN    ! NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LBAD) GOTO 15
      GOTO (10,13), JUMP
  13  WRITE (ILUER,3000) ALINE2
      LERROR=.FALSE.
      GOTO 10
  15  GOTO (100,200,300,400,500),JUMP
C
C     help
C
 100  AFILE='TWHLP.HLP       '
C      CALL HELP                  ! NOT AVAILABLE IN BATCH MODE
      LENTRY=.FALSE.
      GOTO 10
C
C     storage coefficient (default)
C
 200  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RTWST0=RDUM
      GOTO 10
C
C     radius (default)
C
 300  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      IF (RDUM.LT.1.0E-10) THEN
        WRITE (ILUER,3500) RDUM,RTWRD0
        GOTO 10
      ENDIF
      RTWRD0=RDUM
      GOTO 10
C
C     discharge specified transient wells
C
 400  if (lucon) WRITE (ILUME,4000)
      CALL INLINE
      CDUM=CVAR(1)
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 11
      ENDIF      
      RDUM1=RVAR(3)
      RDUM2=RVAR(4)
      RDUM3=RVAR(5)
      IF (LERROR) GOTO 10
      RDUM4=RVAR(6)
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      RDUM4=RTWST0
      RDUM5=RTWRD0
      CALL GETFN(6)
      GOTO 405
      ENDIF
      RDUM5=RVAR(7)
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      RDUM5=RTWRD0
      CALL GETFN(7)
      ELSE
      CALL GETFN(8)
      ENDIF
C     ----------------check if still room for more wells and store data
 405  IF (NTW.EQ.NTWMX) THEN
      WRITE (ILUER,5000) NTWMX
      GOTO 10
      ENDIF
      NTW=NTW+1
      CTWZ(NTW)=CDUM
      RTWQ(NTW)=RDUM1
      RTWST(NTW)=RDUM2
      RTWH0(NTW)=RDUM3
      RTWSTO(NTW)=RDUM4
      IF (RDUM5.LT.1.0E-10) THEN
        WRITE (ILUER,3500) RDUM5,RTWRD0
        RTWR0(NTW)=RTWRD0
      ELSE
        RTWR0(NTW)=RDUM5        
      ENDIF      
      ATWLAB(NTW)=AFILE
      LENTRY=.TRUE.
      CDUM=CTWZ(NTW)+RTWR0(NTW)
      CALL PLWIND (CDUM)
      CDUM=CTWZ(NTW)-RTWR0(NTW)
      CALL PLWIND (CDUM)
      CDUM=CTWZ(NTW)+CMPLX(0.0,RTWR0(NTW))
      CALL PLWIND (CDUM)
      CDUM=CTWZ(NTW)-CMPLX(0.0,RTWR0(NTW))
      CALL PLWIND (CDUM)
      GOTO 400
C
C     return
C
 500  LERROR=.FALSE.
      LMISS=.FALSE.
      RETURN
C
 1000 FORMAT ('          --------- TRANSIENT WELL module ----------'/
     &        ' Maximum number of transient wells:',I4,/
     &        ' Available commands:'/
     &        ' <F1> = Help'/
     &        ' STORAGE  ',G11.4,/
     &        ' RADIUS   ',G11.4,/
     &        ' DISCHARGE '/
     &        ' <Esc> or QUIT'/' >')
 2000 FORMAT (' ***ILLEGAL or MISSING PARAMETERS in twell module:',/,
     &        ' ',80A1)
 3000 FORMAT (' ***ILLEGAL COMMAND in twell module:',/,' ',80A1)
 3500 FORMAT (' ***ERROR in well module: radius of ',G11.4,
     &' is too small.',/,' Current default value used: ',G14.7,/)
 4000 FORMAT (' x  y  discharge  start-time initial-head [storage] ',
     &        ' [radius] [label]',/)
 5000 FORMAT (' ***ERROR in TWELL module: too many wells. (max.=',I3,
     &        ')')
      END
      
