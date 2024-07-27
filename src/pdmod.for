C     Last change:  HMH  16 Jan 2001    8:26 pm
c     This file contains the following routines and functions
c
c     BLOCK DATA PDDAT    initialize common blocks
c     SUBROUTINE PDIN     input routine for 2D sink discs
c     SUBROUTINE PDENTER  input service routine called by PDIN
c
c -------------------
c PDMAT.FOR contains:
c -------------------
c
c     SUBROUTINE PDCZC       generates control points
c     SUBROUTINE PDMAT       generates matrix coefficients
c     SUBROUTINE PDKNO       generates known vector
c     SUBROUTINE PDSUB       substitutes solution
c
c ---------------------
c PDCHECK.FOR contains:
c ---------------------
c
c     SUBROUTINE PDERROR     determines max. error at control points
c
c
c ------------------
c PDIO.FOR contains:
c ------------------
c
c      SUBROUTINE PDIO           reads and writes contents of PDCOM to .SOL file
c      SUBROUTINE PDDATIO        writes (2D) sink disc data to .DAT file
c      subroutine pdextracthead  writes head specified disc data to .XTR file
c          entry pdextractdisch  writes discharge specified disc data to .XTR file
c      SUBROUTINE PDOUT          writes common blocks (formatted) for debugging purposes
c
c -------------------
c PDFUN.FOR contains:
c -------------------
c
c     COMPLEX FUNCTION CFPDOM     returns complex potential when outside disc
c     COMPLEX FUNCTION CFPDCO     coefficient function for complex potential when outside disc
c     REAL FUNCTION RFPDPT        real potential for all discs
c     REAL FUNCTION RFPDCO        coefficient function for discharge potential when inside disc
c     SUBROUTINE PDQI             discharge vector for all discs
c
c ---------------------
c PDNFLOW.FOR contains:
c ---------------------
c
c     REAL FUNCTION RFNFPD      returns flow across a straight line element
c     REAL FUNCTION RFNFPDCO    normal flow coefficient function
c     REAL FUNCTION RFDPPDI     normal flow coefficient function when inside disc
c     REAL FUNCTION RFDPPDO     normal flow coefficient function when outside disc
c
c
c --------------------
c PDSERV.FOR contains:
c --------------------
c
c     REAL FUNCTION RFPDS       returns exfiltration rate due to all (2D) sink discs
c     REAL FUNCTION RFPDSBOTTOM returns bottom exfiltration due to discharge specified discs
c     REAL FUNCTION RFPDH       returns exfiltration rate due to all head specified discs
c
c
c --------------------
c PDNEAR.FOR contains:
c --------------------
c
c     SUBROUTINE PDNEAR      aborts streamline tracing when at sink disc
c
c ---------------------
c PDEXTRA.FOR contains:
c ---------------------
c
c     SUBROUTINE PDCHEK      check of boundary conditions
c     REAL FUNCTION RFBCPD   branch cut correction for sink discs (not used)
c     SUBROUTINE PDCUR       for interactive graphics data retrieval and entry
c     SUBROUTINE PDPLOT      adds (2D) sink discs to the layout
c
c
c ----------------------------------------------------------------------
c
      BLOCK DATA PDDAT
c
c ----------------------------------------------------------------------
c
      IMPLICIT NONE
      INCLUDE 'PDCOM.INC'
      DATA NPD,NPDHD,NPDRC /3*0/
      DATA LPDHD,LPDRC,LPDRCH,LPDPER /2*.TRUE.,2*.FALSE./
      DATA RPI2,RPDC0,RPDD0 /6.283185308,2*0.0/
      DATA IPDHD,IPDRC /NPDMX*0,NPDMX*0/
      DATA LPDGIV /NPDMX*.FALSE./
      DATA APDLAB /NPDMX*'                '/
      DATA CPDZ /NPDMX*(0.0,0.0)/
      DATA RPDS,RPDSB,RPDH,RPDR /NPDMX*0.,NPDMX*0.,NPDMX*0.,NPDMX*0./
      DATA RPDRES,RPDEP /NPDMX*0.0,NPDMX*0.0/
      END
c
c ----------------------------------------------------------------------
c
      SUBROUTINE PDIN (LSOL,RA,IRA,RSCR)
c
c ----------------------------------------------------------------------
c
C
C   Discsink function
C
      IMPLICIT NONE
      INTEGER(4) IRA,ICODE,JUMP
      LOGICAL LSOL,LBAD
      REAL(8) RA,RSCR,RDUM,RVAR
      CHARACTER(1) AWORD(43),APAR(25)
      INCLUDE 'PDCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
      DIMENSION RA(IRA,*),RSCR(*)
      DATA AWORD /'D','I','S','C',' ',
     &            'H','E','A','D',' ',
     &            '?',' ',
     &            'R','E','S','I',' ',
     &            'D','E','P','T',' ',
     &            'P','L','O','T',' ',
     &            'H','I','G','H',' ',
     &            'Q','U','I','T',' ',
     &            'C','U','R','S',' ',
     &            ATERM/
      DATA APAR  /'A','L','L',' ',
     &            'H','E','A','D',' ',
     &            'R','E','C','H',' ',
     &            'P','E','R','C',' ',
     &            'N','O','N','E',' ',
     &            ATERM/ 
      ICODE=0
      LERROR=.FALSE.
      LMISS=.FALSE.      
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
 10   IF (LERROR.OR.LMISS) WRITE (ILUER,1000) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      IF (ICODE.NE.0) THEN
      JUMP=ICODE
      GOTO 15
      ENDIF
      if (lucon)  then
      WRITE (ILUME,2000) NPDMX,RPDC0,RPDD0
      IF (LPDHD.AND..NOT.LPDRC) WRITE (ILUME,2010)
      IF (LPDRC.AND..NOT.LPDHD) WRITE (ILUME,2020)
      IF (LPDHD.AND.LPDRC)      WRITE (ILUME,2030)
      IF (LPDRCH.AND.LPDPER)      WRITE (ILUME,2040)
      IF (.NOT.LPDRCH.AND.LPDPER) WRITE (ILUME,2050)
      IF (.NOT.LPDRCH.AND..NOT.LPDPER) WRITE (ILUME,2060)
      WRITE (ILUME,2070)
      endif
      CALL INLINE
  11  CALL MATCH (AWORD,1,JUMP,LBAD)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE      IF (.NOT.LBAD) GOTO 15
      GOTO (10,13),JUMP
  13  WRITE (ILUER,1010) ALINE2
      LERROR=.FALSE.
      GOTO 10
  15  GOTO (100,200,300,400,500,600,700,800,900),JUMP
C
C     Enter discharge specified sinkdiscs (lakes, swamps, or 
C                                      recharge due to precipitation)
C
 100  if (lucon) WRITE (ILUME,3000)
      CALL INLINE
      RDUM=RVAR(1)
      ICODE=1
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 11
      ENDIF
      CALL PDENTER (ICODE)
      IF (LERROR) GOTO 10
      LSOL=.FALSE.
      GOTO 100  
C
C     Enter head specified sinkdiscs (lakes or swamps)
C
 200  if (lucon) WRITE (ILUME,4000)
      CALL INLINE
      RDUM=RVAR(1)
      ICODE=2
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 11
      ENDIF
      CALL PDENTER (ICODE)
      IF (LERROR) GOTO 10
      LSOL=.FALSE.
      GOTO 200
C
C     Help
C
 300  AFILE='PDHLP.HLP       '
C      CALL HELP                    ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     Enter default resistance of layer underneath sinkdisc
C
 400  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RPDC0=RDUM
      GOTO 10
C
C     Enter default depth from water table to bottom of resistance layer
C
 500  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RPDD0=RDUM
      GOTO 10
C
C     Set plotting of sinkdiscs as part of layout
C 
 600  CALL MATCH (APAR,2,JUMP,LBAD)
      LERROR=LBAD
      ICODE=0
      IF (LERROR) GOTO 10
      GOTO (610,620,630),JUMP
 610  LPDHD=.TRUE.
      LPDRC=.TRUE.
      GOTO 10
 620  LPDHD=.TRUE.
      LPDRC=.FALSE.
      GOTO 10
 630  LPDHD=.FALSE.
      LPDRC=.TRUE.
      GOTO 10
C
C     Highlight recharging or percolating sinkdiscs that are head specified
C
 700  CALL MATCH (APAR,2,JUMP,LBAD)
      ICODE=0
      LERROR=LBAD
      IF (LERROR) GOTO 10
      LERROR=JUMP.LT.3
      IF (LERROR) GOTO 10
      IF (JUMP.EQ.3) THEN  ! highlight recharging sinkdiscs
      LPDRCH=.TRUE.
      LPDPER=.TRUE.
      ENDIF
      IF (JUMP.EQ.4) THEN  ! highlight percolating sinkdiscs
      LPDRCH=.FALSE.
      LPDPER=.TRUE.
      ENDIF
      IF (JUMP.EQ.5) THEN  ! highlight none
      LPDRCH=.FALSE.
      LPDPER=.FALSE.
      ENDIF
      GOTO 10
C
C     Return
C
 800  RETURN
C
C     Cursor
C 
  900 IF (.NOT.LSOL) THEN
C      CALL TONE                           ! NOT AVAILABLE IN BATCH MODE
      WRITE (ILUER,1020)
      LERROR=.TRUE.
      ENDIF
C      CALL PDCUR (RA,IRA,RSCR,LSOL)       ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
 1000 FORMAT (' ***ILLEGAL or MISSING PARAMETERS in sinkdisc module:',/,
     &        ' ',80A1) 
 1010 FORMAT (' ***ILLEGAL COMMAND in sinkdisc module:',/,' ',80A1)
 1020 FORMAT (' ***WARNING: no valid solution, data likely in ERROR !'/)
 1500 FORMAT (80A1)
 2000 FORMAT ('            -------- SINKDISC module --------'/
     &' Maximum number of sinkdiscs:',I4,/,
     &' Available commands:',/,
     &' <F1> = Help ',/,
     &' HEAD ',/,
     &' DISCHARGE ',/,
     &' RESISTANCE ',G11.4,/,
     &' DEPTH      ',G11.4,/,
     &' CURSOR')
 2010 FORMAT (' PLOT head      (all / discharge)')
 2020 FORMAT (' PLOT discharge (all / head)')     
 2030 FORMAT (' PLOT all       (head / discharge)')
 2040 FORMAT (' HIGHLIGHT RECHARGING (percolating/none)')
 2050 FORMAT (' HIGHLIGHT PERCOLATING (recharging/none)')
 2060 FORMAT (' HIGHLIGHT NONE  (recharging/percolating)')
 2070 FORMAT (' <Esc> or QUIT',/,' >')
 3000 FORMAT (' x0 y0  xc yc  extraction rate top [extr. rate bottom]',
     &' [label]',/,' >') 
 4000 FORMAT (' x0 y0  xc yc   head  [resistance]  [depth]  [label]',/,
     &        ' >')
      END 
c
c ----------------------------------------------------------------------
c
      SUBROUTINE PDENTER (ICODE)
c
c ----------------------------------------------------------------------
c
      IMPLICIT NONE
      INTEGER(4) ICODE
      REAL(8) RDUM1,RDUM2,RVAR,RDUM3,RDUM4,RTOL,RFGRTOL
      COMPLEX(8) CDUM1,CVAR,CDUM2,CZ
      INCLUDE 'PDCOM.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
C -------------------------- check for space in arrays
      IF (NPD.EQ.NPDMX) THEN
      WRITE (ILUER,1000) NPDMX
      LERROR=.TRUE.
      RETURN
      ENDIF
C -------------------------retrieve data from ALINE
      CDUM1=CVAR(1)
      CDUM2=CVAR(3)
      RDUM1=ABS(CDUM1-CDUM2)
      RDUM2=RVAR(5)
      IF (LERROR) RETURN
      RDUM3=RVAR(6)
      IF (LERROR) THEN
      IF (ICODE.EQ.1) RDUM3=0.0     ! default aquifer bottom exfiltration rate
      IF (ICODE.EQ.2) RDUM3=RPDC0   ! default resistance of sinkdisc bottom
      RDUM4=RPDD0
      CALL GETFN (6)
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 10
      ENDIF
      RDUM4=RVAR(7)
      IF (LERROR) THEN
      RDUM4=RPDD0
      CALL GETFN (7)
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 10
      ENDIF
      CALL GETFN (8)
C -----------------------substitute data in arrays
  10  RTOL=RFGRTOL()
      IF (RDUM1.LT.RTOL) THEN       ! radius may be too small
        IF (RDUM1.EQ.0.0) THEN
          WRITE (ILUER,3000) AFILE
          RETURN
        ELSE
          WRITE (ILUER,3100) AFILE,RDUM1
        ENDIF
      ENDIF
      NPD=NPD+1
      CPDZ(NPD)=CDUM1
      RPDR(NPD)=RDUM1
C ----------------------------   update max. window
      CZ=CDUM1+RDUM1
      CALL PLWIND (CZ)
      CZ=CDUM1-RDUM1
      CALL PLWIND (CZ)
      CZ=CDUM1+CMPLX(0.0,RDUM1)
      CALL PLWIND (CZ)
      CZ=CDUM1-CMPLX(0.0,RDUM1)
      CALL PLWIND (CZ)
      IF (ICODE.EQ.1) THEN
C ---------------------------discharge specified sinkdiscs
      NPDRC=NPDRC+1
      IPDRC(NPDRC)=NPD
      RPDS(NPD)=RDUM2+RDUM3
      RPDSB(NPD)=RDUM3
      RPDH(NPD)=0.0
      RPDRES(NPD)=0.0
      RPDEP(NPD)=0.0
      APDLAB(NPD)=AFILE
      RETURN
      ENDIF
      IF (ICODE.EQ.2) THEN
C ----------------------------head specified sinkdiscs
      NPDHD=NPDHD+1
      IPDHD(NPDHD)=NPD
      RPDS(NPD)=0.0
      RPDSB(NPD)=0.0
      RPDH(NPD)=RDUM2
      IF (RDUM3.LT.1.0.AND.RDUM3.GT.0.0) THEN
      WRITE (ILUER,4000) AFILE,RDUM3
      RDUM3=0.0
      ENDIF
      IF (RDUM3.GT.10000.0) THEN
      WRITE (ILUER,5000) AFILE,RDUM3
      RDUM3=10000.0
      ENDIF
      RPDRES(NPD)=RDUM3
      RPDEP(NPD)=RDUM4
      APDLAB(NPD)=AFILE
      RETURN
      ENDIF      
      WRITE (ILUER,2000) ICODE
      LERROR=.TRUE.
      RETURN
 1000 FORMAT (' ***TOO MANY SINKDISCS (max.=',I3,')')
 2000 FORMAT (' ***ERROR in PDENTER: ICODE out of range (ICODE=',I3,')')
 3000 FORMAT (' ***ERROR: sinkdisc ',A16,' has a zero radius!',/,
     &' Sinkdisc has been ignored.',/)
 3100 FORMAT (' ***WARNING: sinkdisc ',A16,' has a radius ',G11.4,/,
     &' which may be too small.',/)     
 4000 FORMAT (' ***WARNING: sinkdisc with label ',A16,' has a ',
     &'resistance ',G11.4,/,' A resistance below 1.0 may lead to an ',
     &'unstable solution procedure!',/,
     &' Resistance parameter is set to 0.0 , ',
     &'meaning no resistance is assumed.',/)
 5000 FORMAT (' ***WARNING: sinkdisc with label ',A16,' has a ',
     &'resistance ',G11.4,/,
     &' This sinkdisc will play no significant role in the ',
     &'groundwater flow solution!',/,' Resistance parameter is set to ',
     &'10000 to conform to format of the DATA command.',/)
      END
