c       This file contains the following routines and functions:
c
C	BLOCK DATA WLDAT    initializes well common blocks
C	SUBROUTINE WLIN     handles input of wells as defined in the .dat file
C
c -------------------
c wlmat.for contains:
c -------------------
c
C    WLCZC   generates control point vector for wells
C    WLMAT   generates matrix coefficients for wells
C    WLKNO   generates known vector for wells
C    WLSUB   substitutes solution vector into strengths parameters for wells  
c
c -------------------
c wlfun.for contains:
c -------------------
c
C     SUBROUTINE WELQI         returns the discharge vector due to all (2D) wells
C     COMPLEX FUNCTION CFWLQC  returns the "W-function" due to well i
C     COMPLEX FUNCTION CFWEOM  returns the complex potential due to all (2D) wells   
C     COMPLEX FUNCTION CFWLOMC returns the complex potential for well i with strength 1.  
c
c ---------------------
c wlnflow.for contains:
c ---------------------
c
C     REAL FUNCTION RFNFWL    returns the flow across a line due to all wells
C     REAL FUNCTION RFNFWLCO  coefficient function for the flow across a line due to well i.
c
c ------------------
c wlio.for contains:
c ------------------
c
C     SUBROUTINE wlextracthead    writes data for head specified wells to .xtr file
C          ENTRY wlextractdisch   writes data for discharge specified wells to .xtr file
C     SUBROUTINE WLDATIO          rites all well data to a ".dat" file
C     SUBROUTINE WLIO             reads or writes contents of WLCOM.INC from or to .sol file
C     SUBROUTINE WLOUT            writes contents of WLCOM formatted to ILUOUT for debugging purposes
c
c ---------------------
c wlcheck.for contains:
c ---------------------
c
C     SUBROUTINE WLERROR   maximum error in the boundary conditions specified at wells
c
c --------------------
c wlnear.for contains:
c --------------------
c
C     SUBROUTINE WLNEAR  check for streamline trace near well and slows down trace or ends trace
c
c --------------------
c wlserv.for contains:
c --------------------
c
C     Logical FUNCTION lwlinfo    routine returns information for the 2D well
c
c ---------------------
c wlextra.for contains:
c ---------------------
c
C     SUBROUTINE WEPLOT        plots all (2D) wells
C     SUBROUTINE WLCUR         facilitates graphical display or modification of well attributes
C     SUBROUTINE WLCHEK        check of boundary conditions
C
c
C-------------------------------------------------------------------------------------------------------
C
C     Last change:  HMH   2 Jan 2001    5:33 pm
      BLOCK DATA WLDAT
      IMPLICIT NONE
      INCLUDE 'wlcom.inc'
      DATA CWLZ /NWLMX*(0.0,0.0)/
      DATA RWLQ,RWLH,RWLR /NWLMX*0.0,NWLMX*0.0,NWLMX*0.0/
      DATA NWL,RPI2/0,6.2831853/
      DATA LWLH /NWLMX*.FALSE./
      DATA RWLD0 /NWLMX*0.0,NWLMX*0.0,NWLMX*0.0/
      DATA AWLAB /NWLMX*'                '/
      DATA RWLRD0 /0.1/
      END
C
C --------------------------------------------------------------------
C
      SUBROUTINE WLIN (LSOL,RA,IRA,RSCR)
      IMPLICIT NONE
      INTEGER(4) IRA,JUMP,ICODE
      LOGICAL LSOL,LENTRY,LBAD
      REAL(8) RA,RSCR,RDUM1,RVAR,RDUM2,RDUM
      COMPLEX(8) CDUM,CVAR
      CHARACTER(1) AWORD(28)
      INCLUDE 'wlcom.inc'
      INCLUDE 'lusys.inc'
      INCLUDE 'match.inc'
      DIMENSION RA(IRA,*),RSCR(*)
      DATA AWORD /'?',' ',
     .            'H','E','A','D',' ',
     .            'D','I','S','C',' ',
     .            'R','A','D','I',' ',
     .            'C','U','R','S',' ',
     .            'Q','U','I','T',' ',
     .            ATERM/
      LENTRY=.FALSE.
      LERROR=.FALSE.
      LMISS=.FALSE.
C      CALL CLEARSCREEN                ! NOT AVAILABLE IN BATCH MODE
  10  IF (LENTRY) THEN
      JUMP=ICODE
      GOTO 15
      ENDIF
      IF (LERROR) WRITE (ILUER,2000) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      if (lucon) then
      WRITE (ILUME,1000) NWLMX
      WRITE (ILUME,1010) RWLRD0
      WRITE (ILUME,1020)
      endif
      CALL INLINE
  11  CALL MATCH (AWORD,1,JUMP,LBAD)
C      CALL CLEARSCREEN               ! NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LBAD) GOTO 15
      GOTO (10,13),JUMP
  13  WRITE (ILUER,3000) ALINE2
      LERROR=.FALSE.
      GOTO 10
  15  GOTO (100,200,300,400,500,600),JUMP
C  
 100  AFILE='                '
      AFILE(1:9)='WLHLP.HLP'
C      CALL HELP                     ! NOT AVAILABLE IN BATCH MODE
      LENTRY=.FALSE.
      GOTO 10
C
C     head specified wells
C
 200  if (lucon) WRITE (ILUME,4000)
      CALL INLINE
      CDUM=CVAR(1)
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 11
      ENDIF
      RDUM1=RVAR(3)
      IF (LERROR) GOTO 10
      LENTRY=.TRUE.
      ICODE=2
      RDUM2=RVAR(4)
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      RDUM2=RWLRD0
      CALL GETFN (4)
      ELSE
      CALL GETFN (5)
      ENDIF
C -----------------------------check if still space for more wells      
      IF (NWL.EQ.NWLMX) THEN
      WRITE (ILUER,1800) NWLMX
      RETURN
      ENDIF
C ------------------------------store data      
      NWL=NWL+1
      CWLZ(NWL)=CDUM
      LWLH(NWL)=.TRUE.
      RWLH(NWL)=RDUM1
      IF (RDUM2.LT.1.0E-10) THEN
        WRITE (ILUER,3500) RDUM2,RWLRD0
        RWLR(NWL)=RWLRD0
      ELSE
        RWLR(NWL)=RDUM2
      ENDIF
      AWLAB(NWL)=AFILE
      RWLQ(NWL)=0.0
      LSOL=.FALSE.
      CDUM=CWLZ(NWL)+RWLR(NWL)
      CALL PLWIND (CDUM)
      CDUM=CWLZ(NWL)-RWLR(NWL)
      CALL PLWIND (CDUM)
      CDUM=CWLZ(NWL)+CMPLX(0.0,RWLR(NWL))
      CALL PLWIND (CDUM)
      CDUM=CWLZ(NWL)-CMPLX(0.0,RWLR(NWL))
      CALL PLWIND (CDUM)
      GOTO 200
C
C     discharge specified wells
C
 300  if (lucon) WRITE (ILUME,6000)
      CALL INLINE
      CDUM=CVAR(1)
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 11
      ENDIF
      RDUM1=RVAR(3)
      IF (LERROR) GOTO 10
      LENTRY=.TRUE.
      ICODE=3
      RDUM2=RVAR(4)
      IF (LERROR) THEN
      LERROR=.FALSE.
      RDUM2=RWLRD0
      CALL GETFN (4)
      ELSE
      CALL GETFN (5)
      ENDIF      
C -----------------------------check if still space for more wells      
      IF (NWL.EQ.NWLMX) THEN
      WRITE (ILUER,1800) NWLMX
      RETURN
      ENDIF
C ------------------------------store data      
      NWL=NWL+1
      CWLZ(NWL)=CDUM
      LWLH(NWL)=.FALSE.
      RWLQ(NWL)=RDUM1
      IF (RDUM2.LT.1.0E-10) THEN
        WRITE (ILUER,3500) RDUM2,RWLRD0
        RWLR(NWL)=RWLRD0
      ELSE
        RWLR(NWL)=RDUM2
      ENDIF
      AWLAB(NWL)=AFILE
      RWLH(NWL)=0.0
      LSOL=.FALSE.
      CDUM=CWLZ(NWL)+RWLR(NWL)
      CALL PLWIND (CDUM)
      CDUM=CWLZ(NWL)-RWLR(NWL)
      CALL PLWIND (CDUM)
      CDUM=CWLZ(NWL)+CMPLX(0.0,RWLR(NWL))
      CALL PLWIND (CDUM)
      CDUM=CWLZ(NWL)-CMPLX(0.0,RWLR(NWL))
      CALL PLWIND (CDUM)
      GOTO 300
C
C     Default radius
C
 400  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      IF (RDUM.LT.1.0E-10) THEN
        WRITE (ILUER,3500) RDUM,RWLRD0
        GOTO 10
      ENDIF
      RWLRD0=RDUM
      GOTO 10
C
C     Cursor
C
 500  CONTINUE
      IF (.NOT.LSOL) THEN
C      CALL TONE                          ! NOT AVAILABLE IN BATCH MODE
      WRITE (ILUER,1900)
      LERROR=.TRUE.
      ENDIF
      LENTRY=.FALSE.
C      CALL WLCUR (RA,IRA,RSCR,LSOL)      ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     Return
C      
 600  LERROR=.FALSE.
      LMISS=.FALSE.
      RETURN
C
 1000 FORMAT ('                       -------- WELL module --------'/
     &        ' Maximum number of wells:',I4,/
     &        ' Available commands:'/' <F1> = Help'/' HEAD'/
     &        ' DISCHARGE'/' CURSOR')
 1010 FORMAT (' RADIUS   ',G11.4)     
 1020 FORMAT (' <Esc> or QUIT'/' >')
 1500 FORMAT (80A1)
 1800 FORMAT (' ***ERROR in well module: too many wells. (max=',I3,')')
 1900 FORMAT (' ***WARNING: no valid solution, data likely in ERROR !'/)
 2000 FORMAT (' ***ILLEGAL or MISSING PARAMETER(S) in well module:',/,
     &' ',80A1)
 3000 FORMAT (' ***ILLEGAL COMMAND in well module:',/,' ',80A1)
 3500 FORMAT (' ***ERROR in well module: radius of ',G11.4,
     &' is too small.',/,' Current default value used: ',G14.7,/)
 4000 FORMAT (' x  y   head    [radius]  [label]',/)
 6000 FORMAT (' x  y   discharge   [radius]  [label]',/)
C
      END

