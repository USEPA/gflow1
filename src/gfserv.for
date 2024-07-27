C     Last change:  HMH  18 Dec 2006    5:12 pm
c     This file contains the following routines and functions
c
c     SUBROUTINE PLAPAR      set up parameters for plotting a layout
c          ENTRY PGRPAR      for plotting a grid (contour plot).
C          ENTRY GGRPAR      for generating a grid.
C          ENTRY WVPAR       for setting a window or viewport.
c     SUBROUTINE PLWIND      updates window parameters (to zoom to extent)
c     REAL FUNCTION RFGRTOL  returns tolerance factor based on current windows parameters
c     SUBROUTINE CIRCLINE    intersection points of a horizontal circle with a vertical plane (Called in DIPLOT and W3PLOT)
c     subroutine GetMainData returns parameters in common main (main.inc)
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE PLAPAR (LESCAP)
c
c -------------------------------------------------------------------------------
c
C     All entries have the argument LESCAP to indicate a return
C     without further action.
C
C     Command words and parameters:
C   #     HELP       [help]
C   1     VIEWPORT   (x1, y1, x2, y2)
C   2     WINDOW     (x1, y1, x2, y2) or (x1,y1,z1,x2,y2,z2) or
C                    all  save  delete  select cursor
C   3     AQUIFER    (ilayer)
C   4     HORIZONTALPOINTS (nx)
C   5     PLOT       (heads / potentials / streamlines / flow net / discharge)
C   6     DEVICE     (screen / postscript / HPGL / hiplot / slides)
C   7     CHART      (large / small)
C   8     SPEED      (1-3 inch/sec)
C   9     STEPSIZE   (0.001",0.005" and 0.1 mm)
C   10    LINE       (0,1,....,9,:)
C   11    MONITOR    (single / dual)
C   12    DISPLAY    (monochrome / color)
C   13    PALETTE    (1,2,3,4)
C   14    MINUSGRID  (filename)
C   15    TMARK      (1,2,3,4,5,6) {.,+,*,square,cross,diamond)
C   16    SAVE       (filename)
C   17    LOAD       (filename)
C   18    SURFER     (filename)
C   19    GENERATE
C   20    GO
C   21    DOTMAP     (on / off)
C         RETURN
C
C     SPECIAL PARAMETERS SET IN COMMON /GRID/
C
C     LRESP       output array to indicate which commands have
C                 been used.
C     ILAYER      input/output parameter indicating the current
C                 aquifer to which commands apply.
C
      IMPLICIT NONE
      INTEGER(4) NWORD,NPAR1,NPAR2,NPAR3,NPAR4,
     &           NPAR5,NPAR6,NPAR7,NPAR8,NPAR9,
     &           NCMPAR,IVAR,ISIZE,I,JUMP,IRGBL,
     &           IXX1,IYY1,IXX2,IYY2,IW,J,IERR,IL
      LOGICAL LESCAP,LCMD,LCMPL,LCMPG,LCMGG,
     &        LCMWV,LDSCR,LDPRI,LDPLO,LDHIP,LNOPLT,
     &        LPLAY,LPGRID,LGGRID,LWIVI,LGO,LBAD,L1,L2,LRET,LRESS
      REAL(8) RVAR,RA,R1,R2,R3,R4,R5,R6,RDUM,RTWTIM
      DIMENSION RA(ISIZE,*)
      COMPLEX(8) CZ1,CZ2
      CHARACTER(1) AWORD,APAR1,APAR2,APAR3,APAR4,
     &             APAR5,APAR6,APAR7,APAR8,APAR9,
     &             BLCHART
      PARAMETER (NWORD=111,NPAR1=77,NPAR2=26,NPAR3=11,NPAR4=23)
      PARAMETER (NPAR5=11,NPAR6=11,NPAR7=30,NPAR8=25,NPAR9=8)
      DIMENSION APAR1(NPAR1),APAR2(NPAR2),APAR3(NPAR3),APAR4(NPAR4)
      DIMENSION APAR5(NPAR5),APAR6(NPAR6),APAR7(NPAR7),APAR8(NPAR8)
      DIMENSION APAR9(NPAR9)
      DIMENSION AWORD(NWORD)
      INCLUDE 'LUSYS.INC'
      INCLUDE 'GRID.INC' 
      INCLUDE 'COM3D.INC'
      INCLUDE 'DWMN.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'PLOT.INC'
      PARAMETER (NCMPAR=21)
      DIMENSION LCMD(NCMPAR),LCMPL(NCMPAR),LCMPG(NCMPAR),LCMGG(NCMPAR),
     .LCMWV(NCMPAR),LDSCR(NCMPAR),LDPRI(NCMPAR),LDPLO(NCMPAR),
     .LDHIP(NCMPAR)
      DATA LCMPL /.FALSE.,.TRUE.,3*.FALSE.,.FALSE.,4*.FALSE.,2*.TRUE.,
     &            .FALSE.,6*.FALSE.,.TRUE.,.FALSE./
      DATA LCMPG /5*.FALSE.,.FALSE.,4*.FALSE.,2*.TRUE.,.FALSE.,
     &6*.FALSE.,.TRUE.,.FALSE./
      DATA LCMGG /.FALSE.,.TRUE.,.FALSE.,2*.TRUE.,8*.FALSE.,.TRUE.,
     &            .FALSE.,3*.TRUE.,.FALSE.,2*.TRUE./
      DATA LCMWV /.FALSE.,.TRUE.,19*.FALSE./
      DATA LDSCR /.FALSE.,.TRUE.,3*.FALSE.,5*.FALSE.,2*.TRUE.,.FALSE.,
     &            6*.FALSE.,2*.TRUE./
      DATA LDPRI /.FALSE.,.TRUE.,4*.FALSE.,5*.FALSE.,.TRUE.,.FALSE.,
     &            6*.FALSE.,.TRUE.,.FALSE./
      DATA LDPLO /.FALSE.,.TRUE.,4*.FALSE.,5*.FALSE.,.TRUE.,.FALSE.,
     &            6*.FALSE.,.TRUE.,.FALSE./
      DATA LDHIP /.FALSE.,.TRUE.,4*.FALSE.,4*.TRUE.,.FALSE.,.TRUE.,
     &            7*.FALSE.,.TRUE.,.FALSE./
      DATA LNOPLT /.TRUE./
      DIMENSION IRGBL(9,4)
      DATA IRGBL / 600, 600,   0, 600,   0,   0,   0, 600,   0,
     .            1000,1000, 400,1000, 400, 400, 400,1000, 400,
     .             600, 600, 600, 600,   0, 600,   0, 600, 600,
     .            1000,1000,1000,1000, 400,1000, 400,1000,1000/
      DATA AWORD /'?',' ',
     &            'V','I','E','W',' ',
     &            'W','I','N','D',' ',
     &            'A','Q','U','I',' ',
     &            'H','O','R','I',' ',
     &            'P','L','O','T',' ',
     &            'D','E','V','I',' ',
     &            'C','H','A','R',' ',
     &            'S','P','E','E',' ',
     &            'S','T','E','P',' ',
     &            'L','I','N','E',' ',
     &            'M','O','N','I',' ',
     &            'D','I','S','P',' ',
     &            'P','A','L','E',' ',
     &            'M','I','N','U',' ',
     &            'T','M','A','R',' ',
     &            'S','A','V','E',' ',
     &            'L','O','A','D',' ',
     &            'S','U','R','F',' ',
     &            'G','E','N','E',' ',
     &            'Q','U','I','T',' ',
     &            'G','O',' ',
     &            'D','O','T','M',' ',
     &            ATERM/
      DATA APAR1 /'H','E','A','D',' ',
     &            'P','O','T','E',' ',
     &            'S','T','R','E',' ',
     &            'D','I','S','C',' ',
     &            'F','L','O','W',' ',
     &            'L','E','A','K','A','N','C','E',' ',
     &            'L','E','A','K','A','G','E',' ',
     &            'B','A','S','E',' ',
     &            'T','O','P',' ',
     &            'H','E','I','G',' ',
     &            'C','O','N','D',' ',
     &            'P','O','R','O',' ',
     &            'R','E','C','H',' ',
     &            'I','N','T','E',' ',
     &            ATERM/
      DATA APAR2 /'S','C','R','E',' ',
     &            'P','R','I','N',' ',
     &            'P','L','O','T',' ', 
     &            'H','I','P','L',' ',
     &            'S','L','I','D',' ',
     &            ATERM/
      DATA APAR3 /'L','A','R','G',' ',
     &            'S','M','A','L',' ',
     &            ATERM/
      DATA APAR4 /'0',' ','1',' ','2',' ','3',' ','4',' ','5',' ',
     &            '6',' ','7',' ','8',' ','9',' ',':',' ',ATERM/
      DATA APAR5 /'S','I','N','G',' ',
     &            'D','U','A','L',' ',
     &            ATERM/
      DATA APAR6 /'M','O','N','O',' ',
     &            'C','O','L','O',' ',
     &            ATERM/
      DATA APAR7 /'D','O','T',' ',
     &            'P','L','U','S',' ',
     &            'A','S','T','E',' ',
     &            'S','Q','U','A',' ',
     &            'C','R','O','S',' ',
     &            'D','I','A','M',' ',
     &            ATERM/
      DATA APAR8 /'A','L','L',' ',
     &            'C','U','R','S',' ',
     &            'S','A','V','E',' ',
     &            'D','E','L','E',' ',
     &            'S','E','L','E',' ',
     &            ATERM/
      DATA APAR9 /'O','N',' ',
     &            'O','F','F',' ',
     &            ATERM/
C
C--------------------------
C           PLAPAR
C--------------------------
C
C                          "entry" used for layout generation
C
      LPLAY=.TRUE.
      LPGRID=.FALSE.
      LGGRID=.FALSE.
      LWIVI=.FALSE.
      GOTO 8
C -------------------------
      ENTRY PGRPAR (LESCAP)
C -------------------------
C                            entry used for contour plot generation
      LPLAY=.FALSE.
      LPGRID=.TRUE.
      IF (LNOPLT) THEN 
       WRITE (ILUER,999) 
 999   FORMAT (' Generate grid before plotting!') 
       RETURN
      ENDIF
      LGGRID=.FALSE.
      LWIVI=.FALSE.
      GOTO 8
C ----------------------------------
      ENTRY GGRPAR (LESCAP,RA,ISIZE)
C ----------------------------------
C                                    entry used for grid generation
      LPLAY=.FALSE.
      LPGRID=.FALSE.
      LGGRID=.TRUE.
      LWIVI=.FALSE.
      GOTO 8
C ------------------------
      ENTRY WVPAR (LESCAP)
C ------------------------
C                             entry used to change window &/or viewport
C --------------------------
      LPLAY=.FALSE.
      LPGRID=.FALSE.
      LGGRID=.FALSE.
      LWIVI=.TRUE.
      GOTO 8
C
   8  DO 9 I=1,15
      LRESP(I)=.FALSE.
   9  CONTINUE
      ADISPG='DISPLAY2' ! default display to screen
      lploton=.false.
      lsingl=.false.
C      IF (.NOT.LERROR) CALL CLEARSCREEN    ! NOT AVAILABLE IN BATCH MODE
      LESCAP=.FALSE.
      LERROR=.FALSE.
      LMISS=.FALSE.
  10  IF (LERROR.OR.LMISS) THEN
      IF (LPLAY) WRITE (ILUER,2000) ALINE2
      IF (LPGRID) WRITE (ILUER,2001) ALINE2
      IF (LGGRID) WRITE (ILUER,2002) ALINE2
      IF (LWIVI) WRITE (ILUER,2003) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      ENDIF
      DO 19 I=1,NCMPAR
      IF (LPLAY) LCMD(I)=LCMPL(I)
      IF (LPGRID) LCMD(I)=LCMPG(I)
      IF (LGGRID) LCMD(I)=LCMGG(I)
      IF (LWIVI) LCMD(I)=LCMWV(I)
  19  CONTINUE
      LGO=.FALSE.
      if (lucon) then
      IF (LGGRID) THEN
      WRITE (ILUME,1001)
      ELSE
      WRITE (ILUME,1000)
      ENDIF
C         NOTE; this menu has not been updated since GFLOW1 runs only in
C         batch mode. For instance, interface contouring is not added to the menu.
      IF (LCMD(1)) WRITE (ILUME,1010) IX1,IY1,IX2,IY2,NDOTX,NDOTY
      IF (LCMD(2).AND..NOT.L3DPL)WRITE (ILUME,1020) RX1,RY1,RX2,RY2
      IF (LCMD(2).AND.L3DPL)WRITE(ILUME,1021)R3DX1,R3DY1,R3DZ1,
     &                                      R3DX2,R3DY2,R3DZ2
      IF (LCMD(3)) WRITE (ILUME,1030) ILAYER
      IF (LCMD(4)) WRITE (ILUME,1040) NX
      IF (LCMD(5).AND.LEAKAGE) WRITE (ILUME,1045)
      IF (LCMD(5).AND.LEAKANCE) WRITE (ILUME,1046)
      IF (LCMD(5).AND.LHEAD) WRITE (ILUME,1051)
      IF (LCMD(5).AND.LPHI) WRITE (ILUME,1052)
      IF (LCMD(5).AND.LPSI) WRITE (ILUME,1053)
      IF (LCMD(5).AND.LDISCH) WRITE (ILUME,1054)
      IF (LCMD(5).AND.LFLOWNET) WRITE (ILUME,1055)
      IF (LCMD(21).AND.LDOTMAP) WRITE (ILUME,1057)
      IF (LCMD(21).AND..NOT.LDOTMAP) WRITE (ILUME,1058)
      IF (LCMD(6).AND.LSCRN.AND..NOT.LSLIDES) WRITE (ILUME,1061)
      IF (LCMD(6).AND.LPRINT) WRITE (ILUME,1062)
      IF (LCMD(6).AND.LPLOT)  WRITE (ILUME,1063)
      IF (LCMD(6).AND.LHPLOT) WRITE (ILUME,1064)
      IF (LCMD(6).AND.LSLIDES) WRITE (ILUME,1065)      
      IF (LCMD(7).AND.BLCHRT.EQ.'F') WRITE (ILUME,1071)
      IF (LCMD(7).AND.BLCHRT.EQ.'H') WRITE (ILUME,1072)
      IF (LCMD(8)) WRITE (ILUME,1080) BLSPD
      IF (LCMD(9)) WRITE (ILUME,1090) BLRES
      IF (LCMD(10)) WRITE (ILUME,1100) BLINE
      IF (LCMD(12).AND..NOT.LCOLOR) WRITE (ILUME,1130)
      IF (LCMD(12).AND.LCOLOR) WRITE (ILUME,1140)
      IF (LCMD(13).AND.IPALET.EQ.1) WRITE (ILUME,1150) IPALET
      IF (LCMD(13).AND.IPALET.EQ.2) WRITE (ILUME,1152) IPALET
      IF (LCMD(13).AND.IPALET.EQ.3) WRITE (ILUME,1154) IPALET
      IF (LCMD(13).AND.IPALET.EQ.4) WRITE (ILUME,1156) IPALET
      IF (LCMD(14)) WRITE (ILUME,1160)
      IF (LCMD(15).AND.IPMTP.EQ.1) WRITE (ILUME,1170)
      IF (LCMD(15).AND.IPMTP.EQ.2) WRITE (ILUME,1175)
      IF (LCMD(15).AND.IPMTP.EQ.3) WRITE (ILUME,1180)
      IF (LCMD(15).AND.IPMTP.EQ.4) WRITE (ILUME,1185)
      IF (LCMD(15).AND.IPMTP.EQ.5) WRITE (ILUME,1190)
      IF (LCMD(15).AND.IPMTP.EQ.6) WRITE (ILUME,1191)
      IF (LCMD(16)) WRITE (ILUME,1192)
      IF (LCMD(17)) WRITE (ILUME,1193)
      IF (LCMD(18)) WRITE (ILUME,1194)
      IF (LCMD(19)) WRITE (ILUME,1195)
      WRITE (ILUME,1197)
      endif
C      CALL SCALE                           ! NOT AVAILABLE IN BATCH MODE
      if (lucon) WRITE (ILUME,2010)
      CALL INLINE
      CALL MATCH (AWORD,1,JUMP,LBAD)
c      CALL CLEARSCREEN                     !  NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LBAD) GOTO 25
      GOTO (10,23),JUMP
  23  IF (LPLAY) WRITE (ILUER,3000)  ALINE2
      IF (LPGRID) WRITE (ILUER,3001) ALINE2
      IF (LGGRID) WRITE (ILUER,3002) ALINE2
      IF (LWIVI) WRITE (ILUER,3003) ALINE2
      LERROR=.FALSE.
      GOTO 10
  25  GOTO (100,195,200,260,300,
     &      350,400,450,500,550,
     &      600,650,700,750,800,
     &      850,860,870,880,900,
     &      925,950,975),JUMP
C
C     HELP
C
 100  IF (LPLAY)  AFILE='PLCMPL.HLP      '
      IF (LPGRID) AFILE='PLCMPG.HLP      '
      IF (LGGRID) AFILE='PLCMGG.HLP      '
      IF (LWIVI)  AFILE='PLCMWV.HLP      '
      CALL HELP
      GOTO 10
C
C     VIEWPORT
C
 195  IF (.NOT.LCMD(1)) GOTO 23
      IXX1=IVAR(2)
      IYY1=IVAR(3)
      IXX2=IVAR(4)
      IYY2=IVAR(5)
      IF (LERROR) GOTO 10
      IX1=IXX1
      IY1=IYY1
      IX2=IXX2
      IY2=IYY2
      LRESP(1)=.TRUE.
      GOTO 10
C
C     WINDOW
C
 200  IF (.NOT.LCMD(2)) GOTO 23
      CALL MATCH (APAR8,2,JUMP,LBAD)
      IF (.NOT.LBAD) THEN
      GOTO (210,220,230,240,250),JUMP
C -------------------window all      
 210  L3DPL=.FALSE.
      L3DHOR=.FALSE.
      L3DVER=.FALSE.
      IF (RXMN.LT.RXMX.AND.RYMN.LT.RYMX) THEN
      RX1=RXMN
      RY1=RYMN
      RX2=RXMX
      RY2=RYMX
      LRESP(2)=.TRUE.
      GOTO 256
      ELSE
      WRITE (ILUER,2015)
      GOTO 10
      ENDIF
C -------------------window cursor      
  220 IF (.NOT.LSCRN) THEN
       WRITE (ILUER,5000)
       GOTO 10
      ENDIF
C      RXX1=0.5*(RX1+RX2)
C      RYY1=0.5*(RY1+RY2)
C      CALL PLOTON            ! COMMENTED OUT STATEMENTS NOT AVAILABLE IN BATCH MODE
C      CALL LAYOUT
C      IF (.NOT.LSINGL) CALL PLOTOF
C      WRITE (ILUME,4000)
C      CALL CURSOR (ACHAR,RXX1,RYY1,15,2,2,.TRUE.)      ! cross hair for X1,Y1
C      IF (ACHAR.EQ.'Q'.OR.ACHAR.EQ.'q') THEN
C        CALL PLOTOF
C        GOTO 10
C      ENDIF
C      RXX2=RXX1
C      RYY2=RYY1
C      WRITE (ILUME,4050)
C      CALL CURSOR (ACHAR,RXX2,RYY2,15,2,2,.TRUE.)      ! cross hair for X2,Y2
C      IF (ACHAR.EQ.'Q'.OR.ACHAR.EQ.'q') THEN
C        CALL PLOTOF
C        GOTO 10
C      ENDIF
C      CALL CURSOR (ACHAR,RXX2,RYY2,0,2,2,.TRUE.)      ! cancel cross hairs
C      CALL DRWTPL (1)                                 ! plot window
C      CALL DRWCLL (2)
C      CZ=CMPLX(RXX1,RYY1)
C      CALL DRAW (CZ,0)
C      CZ=CMPLX(RXX2,RYY1)
C      CALL DRAW(CZ,1)
C      CZ=CMPLX(RXX2,RYY2)
C      CALL DRAW(CZ,1)
C      CZ=CMPLX(RXX1,RYY2)
C      CALL DRAW(CZ,1)
C      CZ=CMPLX(RXX1,RYY1)
C      CALL DRAW(CZ,1)
C      CALL GFPLOT (0,0,0)
C      RX1=RXX1                                        ! set window
C      RY1=RYY1
C      RX2=RXX2
C      RY2=RYY2
C      IF (LSINGL) THEN
C      CALL PLSETCUR (1,1,3)
C      CALL PLOTOF
C      ENDIF
      GOTO 256
C -----------------window save      
  230 IF (NGRWIN.EQ.NGRWMX) THEN
      WRITE (ILUER,2030) NGRWMX
      GOTO 10
      ENDIF
      NGRWIN=NGRWIN+1
      CZWIN(1,NGRWIN)=CMPLX(RX1,RY1)
      CZWIN(2,NGRWIN)=CMPLX(RX2,RY1)
      CZWIN(3,NGRWIN)=CMPLX(RX2,RY2)
      CZWIN(4,NGRWIN)=CMPLX(RX1,RY2)
      WRITE (ILUME,4070) NGRWIN,NGRWMX
      GOTO 10
C ----------------window delete  
  240 CZ1=CMPLX(RX1,RY1)
      CZ2=CMPLX(RX2,RY2)
      DO 243 IW=1,NGRWIN
      IF (CZ1.EQ.CZWIN(1,IW).AND.CZ2.EQ.CZWIN(3,IW)) THEN
      NGRWIN=NGRWIN-1
      WRITE (ILUME,4080) NGRWIN
      IF (IW.LT.NGRWIN) GOTO 10
      DO 242 I=IW,NGRWIN
      DO 241 J=1,4
      CZWIN(J,IW)=CZWIN(J,IW+1)
  241 CONTINUE
  242 CONTINUE
      GOTO 10
      ENDIF
  243 CONTINUE
      WRITE (ILUER,2040) CZ1,CZ2
      GOTO 10      
C ----------------window select  
  250 IF (.NOT.LSCRN) THEN
      WRITE (ILUER,5000) 
      GOTO 10
      ENDIF
      IF (NGRWIN.EQ.0) THEN
      WRITE (ILUER,2051)
      GOTO 10
      ENDIF
C      CALL PLOTON
C      CALL LAYOUT
C      CALL DRWCLL (2)
C      DO 252 IW=1,NGRWIN
C      CALL DRAW(CZWIN(4,IW),0)
C      DO 251 J=1,4
C      CALL DRAW (CZWIN(J,IW),1)
C 251  CONTINUE
C 252  CONTINUE
C      CALL GFPLOT (0,0,0)
C      IF (.NOT.LSINGL) CALL PLOTOF
C      WRITE (ILUME,2052)
C      CALL CURSOR (ACHAR,RXX1,RYY1,0,0,2,.TRUE.)
C      IF (ACHAR.EQ.'Q'.OR.ACHAR.EQ.'q') GOTO 10
C      RDIS=1.0E21
C      CZ=CMPLX(RXX1,RYY1)
C      DO 254 I=1,NGRWIN
C      DO 253 J=1,4
C      R=ABS(CZ-CZWIN(J,I))
C      IF (R.LT.RDIS) THEN
C        RDIS=R
C        IW=I
C      ENDIF
C 253  CONTINUE
C 254  CONTINUE
C      RX1=REAL(CZWIN(1,IW))
C      RY1=AIMAG(CZWIN(1,IW))
C      RX2=REAL(CZWIN(3,IW))
C      RY2=AIMAG(CZWIN(3,IW))
C      IF (LSINGL) THEN
C       CALL PLSETCUR (1,1,3)
C       CALL PLOTOF
C      ENDIF
      GOTO 256
      ENDIF
C ----------------window x1,y2,(z1)  x2,y2,(z2)      
      LBAD=.FALSE.
      R1=RVAR(2)
      R2=RVAR(3)
      R3=RVAR(4)
      R4=RVAR(5)
      IF (LERROR) GOTO 10
      R5=RVAR(6)
      IF (LERROR) THEN
      L3DPL=.FALSE.
      L3DHOR=.FALSE.
      L3DVER=.FALSE.
      LERROR=.FALSE.
      LMISS=.FALSE.
      RX1=R1
      RY1=R2
      RX2=R3
      RY2=R4
      LRESP(2)=.TRUE.
      GOTO 256
      ENDIF
      R6=RVAR(7)
      IF (LERROR) GOTO 10
      L3DPL=.TRUE.
      R3DX1=R1
      R3DY1=R2
      R3DZ1=R3
      R3DX2=R4
      R3DY2=R5
      R3DZ2=R6
      RX1=R1
      RY1=R2
      RX2=R4
      RY2=R5
      L3DHOR=R3DZ1.EQ.R3DZ2
      L3DVER=R3DZ1.NE.R3DZ2
      LRESP(2)=.TRUE.
 256  L1=RX1.EQ.RWX1.AND.RX2.EQ.RWX2
      L2=RY1.EQ.RWY1.AND.RY2.EQ.RWY2
      LGRCOM=L1.AND.L2
      GOTO 10      
C
C     AQUIFER
C
 260  IF (.NOT.LCMD(3)) GOTO 23
      IXX1=IVAR(2)
      IF (LERROR) GOTO 10
      ILAYER=IXX1
      LRESP(3)=.TRUE.
      GOTO 10
C
C     HORIZONTAL POINTS (of grid)
C
 300  IF (.NOT.LCMD(4)) GOTO 23
      IXX1=IVAR(2)
      IF (LERROR) GOTO 10
      IF (IXX1.GT.NXMAX) THEN
        WRITE (ILUER,5070) NXMAX
        GOTO 10
      ENDIF
      NX=IXX1
      LRESP(4)=.TRUE.
      GOTO 10
C
C     PLOT
C
 350  IF (.NOT.LCMD(5)) GOTO 23
      CALL MATCH (APAR1,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      LHEAD=JUMP.EQ.1
      LPHI=JUMP.EQ.2
      LPSI=JUMP.EQ.3
      LDISCH=JUMP.EQ.4
      LFLOWNET=JUMP.EQ.5
      LEAKANCE=JUMP.EQ.6
      LEAKAGE=JUMP.EQ.7
      LBASE=JUMP.EQ.8
      LTOP=JUMP.EQ.9
      LHEIGHT=JUMP.EQ.10
      LCONDUCTIVITY=JUMP.EQ.11
      LPOROSITY=JUMP.EQ.12
      LRECHARGE=JUMP.EQ.13
      LINTELV=JUMP.EQ.14
      LRESP(5)=.TRUE.
      IF (LPSI.OR.LFLOWNET) THEN
        RDUM=RTWTIM()
        IF (RDUM.GT.0.0) THEN
C          CALL TONE                   !  NOT AVAILABLE IN BATCH MODE
          WRITE (ILUME,5050) RDUM
        ENDIF
      ENDIF
      GOTO 10
C
C     DEVICE
C
 400  IF (.NOT.LCMD(6)) GOTO 23
      CALL MATCH (APAR2,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      LSCRN=.FALSE.
      LPRINT=.FALSE.
      LPLOT=.FALSE.
      LHPLOT=.FALSE.
      LSLIDES=.FALSE.
      NDOTX=32767
      NDOTY=32767
      GOTO (405,410,415,420,425),JUMP
C ----------------------------------screen  (GSS)    
 405  ADISPG='DISPLAY2'      
      LSCRN=.TRUE.
      NDOTY=23900
      IF (.NOT.LSLIDES) THEN
      IY2=MIN(IY2,NDOTY)
      ENDIF    
      IF (LCOLOR) ADISPG='DISPLAY2'
      DO 408 I=7,15
      LCMPL(I)=LDSCR(I)
      LCMPG(I)=LDSCR(I)
 408  CONTINUE      
      IF (.NOT.LCOLOR) THEN
      LCMPL(13)=.FALSE.
      LCMPG(13)=.FALSE.
      ENDIF
      GOTO 430
C ----------------------------------printer   (GSS)
 410  ADISPG='DISPLAY3'
      LPRINT=.TRUE.
      DO 411 I=7,15
      LCMPL(I)=LDPRI(I)
      LCMPG(I)=LDPRI(I)
 411  CONTINUE
      GOTO 430
C -----------------------------------plotter   (GSS)
 415  ADISPG='PLOTTER '
      LPLOT=.TRUE.
      DO 416 I=7,15
      LCMPL(I)=LDPLO(I)
      LCMPG(I)=LDPLO(I)
 416  CONTINUE
      GOTO 430
C ------------------------------------hiplot      
 420  LHPLOT=.TRUE.
      DO 421 I=7,15
      LCMPL(I)=LDHIP(I)
      LCMPG(I)=LDHIP(I)
 421  CONTINUE
C ------------------------------------- slides
 425  LSLIDES=.TRUE.
      GOTO 405
C --------------------------------------      
 430  CONTINUE
C      CALL DOTS                          ! NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LSLIDES) THEN
      IX2=MIN(IX2,NDOTX)
      IY2=MIN(IY2,NDOTY)
      ENDIF
      LRESP(6)=.TRUE.
      IF (LGO) GOTO 950
      GOTO 10   
C
C     CHART
C
 450  IF (.NOT.LCMD(7)) GOTO 23
      CALL MATCH(APAR3,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      IF (JUMP.EQ.1) BLCHRT='F'
      IF (JUMP.EQ.2) BLCHRT='H'
      LRESP(7)=.TRUE.
      GOTO 10
C
C     SPEED
C
 500  IF (.NOT.LCMD(8)) GOTO 23
      IXX1=IVAR(2)
      IF (LERROR) GOTO 10
      IF (IXX1.EQ.1) BLSPD='1'
      IF (IXX1.EQ.2) BLSPD='2'
      IF (IXX1.EQ.3) BLSPD='3'
      IF (IXX1.GE.1.AND.IXX1.LE.3) GOTO 501
      WRITE (ILUER,2000)
      GOTO 10
 501  LRESP(8)=.TRUE.
      GOTO 10
C
C     STEPSIZE
C
 550  IF (.NOT.LCMD(9)) GOTO 23
      CALL MATCH (APAR4,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      IF (JUMP.GT.1.AND.JUMP.LT.5) THEN
      IF (JUMP.EQ.2) BLRES='1'
      IF (JUMP.EQ.3) BLRES='2'
      IF (JUMP.EQ.4) BLRES='M'
C      CALL DOTS                 ! NOT AVAILABLE IN BATCH MODE
      IX1=0
      IY1=0
      IX2=NDOTX
      IY2=NDOTY
      LRESP(9)=.TRUE.
      GOTO 10
      ELSE
      WRITE (ILUER,2000)
      GOTO 10
      ENDIF
C
C     LINE
C
 600  IF (.NOT.LCMD(10)) GOTO 23
      CALL MATCH (APAR4,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      IF (JUMP.EQ.1) BLINE='0'
      IF (JUMP.EQ.2) BLINE='1'
      IF (JUMP.EQ.3) BLINE='2'
      IF (JUMP.EQ.4) BLINE='3'
      IF (JUMP.EQ.5) BLINE='4'
      IF (JUMP.EQ.6) BLINE='5'
      IF (JUMP.EQ.7) BLINE='6'
      IF (JUMP.EQ.8) BLINE='7'
      IF (JUMP.EQ.9) BLINE='8'
      IF (JUMP.EQ.10)BLINE='9'
      IF (JUMP.EQ.11)BLINE=':'
      LRESP(10)=.TRUE.
      GOTO 10
C
C     MONITOR (single / dual)
C
 650  IF (.NOT.LCMD(11)) GOTO 23
      CALL MATCH (APAR5,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      LSINGL=JUMP.EQ.1
      LRESP(11)=.TRUE.
      IF (LSINGL) NDOTY=23900
      IF (.NOT.LSINGL) NDOTY=32767
      GOTO 10
C
C     DISPLAY (monochrome / color)
C      
 700  IF (.NOT.LCMD(12)) GOTO 23
      CALL MATCH (APAR6,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      IF (JUMP.EQ.1) THEN     ! graphics mono
      LCOLOR=.FALSE.
      IF (LSCRN) ADISPG='DISPLAY2'
      IF (LPLAY) LCMPL(13)=.FALSE.
      IF (LPGRID) LCMPG(13)=.FALSE.
      ELSE
      LCOLOR=.TRUE.           ! graphics color
      IF (LSCRN) ADISPG='DISPLAY2'
      ENDIF
      LRESP(13)=.TRUE.
      GOTO 10
C
C     PALETTE
C
 750  IF (.NOT.LCMD(13)) GOTO 23
      IPALET=IVAR(2)
      IF (LERROR) GOTO 10
      DO 751 I=1,9
 751  IRGBIN(I)=IRGBL(I,IPALET)
      LPALET=.TRUE.
      LRESP(13)=.TRUE.
      GOTO 10
C
C     MINUSGRID (subtract a grid from the current grid)
C 
 800  IF (.NOT.LCMD(14)) GOTO 23
      IF (LSUBTRACT.OR.LFLOWNET.OR.L3DVER.OR..NOT.LGRCOM) THEN
      IF (LFLOWNET) WRITE (ILUER,5061)
      IF (L3DVER) WRITE (ILUER,5062)
      IF (LSUBTRACT) WRITE (ILUER,5063)
      IF (.NOT.LGRCOM) WRITE (ILUER,5064)
      GOTO 10
      ENDIF
      CALL OPFILU (2,1,LRET,'.GRI')
      IF (LRET) GOTO 10
      CALL BUFIN (73,2,IERR)
      CALL GRTEMPIO (73,2,LRET)
      IF (LRET) THEN
        WRITE (ILUER,5060)
        CALL BUFEX (73,2,IERR)
        CLOSE (2)        
        GOTO 10
      ENDIF
      CALL RATEMPIO (73,2)
      CALL BUFEX (73,2,IERR)
      CLOSE (2)
      CALL DIFFGRID (RA,ISIZE)      
      LRESP(14)=.TRUE.
      LSUBTRACT=.TRUE.
      LGRID=.TRUE.
      WRITE (ILUME,5065) AFILE
      GOTO 10
C
C     TMARKER
C
 850  IF (.NOT.LCMD(15)) GOTO 23
      CALL MATCH (APAR7,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      IPMTP=JUMP
      LRESP(15)=.TRUE.
      GOTO 10
C
C     SAVE a grid 
C      
 860  IF (.NOT.LCMD(16)) GOTO 23
      IF (LFLOWNET) THEN
        WRITE (ILUER,6000)
        GOTO 10
      ENDIF
      CALL CRFILU (2,1,LRET,'.GRI')
      IF (LRET) GOTO 10
      CALL BUFIN (41,2,IERR)
      CALL GRIO (41,2)
c      write (iluer,1002) isize
c 1002 format (' GFSERV2: save grid--isize=',i5)
c      call gridout(1) ! debugging, remove
      CALL RAIO (41,2,RA,ISIZE)
      CALL BUFEX (41,2,IERR)
      CLOSE (2)
      GOTO 10
C
C     LOAD a grid
C
 870  IF (.NOT.LCMD(17)) GOTO 23
      CALL OPFILU (2,1,LRET,'.GRI')
      IF (LRET) GOTO 10
      NLEVELS=0
      NTEXTS=0
      LFLOWNET=.FALSE.
      CALL BUFIN (73,2,IERR)
      CALL GRIO (73,2)
      CALL RAIO (73,2,RA,ISIZE)
      CALL BUFEX (73,2,IERR)
      CLOSE (2)
      LSUBTRACT=.FALSE.
      IF (ILAYER.NE.IGRIDCODE) THEN ! verify grid is compatible with solution
        WRITE (ILUER,5055)
        LGRCOM=.FALSE.
        GOTO 10
      ENDIF
      LGRID=.TRUE.
      GOTO 10
C
C     SURFER
C
 880  IF (.NOT.LCMD(18)) GOTO 23
      CALL GFSURF (RA,ISIZE)
      IF (LERROR.OR.LMISS) GOTO 10
      GOTO 10
C
C     GO
C
 900  IF (.NOT.LCMD(19)) GOTO 23
      LESCAP=.FALSE.
      GOTO 951
C
C     RETURN
C      
 925  LESCAP=.TRUE.
      GOTO 951
C
C     GO
C
 950  LGO=.TRUE.
      IF (LCURS.AND..NOT.LSCRN) THEN
        WRITE (ILUER,3500)
        CALL INLINE
        IF (ALINE(1).EQ.'Q'.OR.ALINE(1).EQ.'q') GOTO 925
        LPRINT=.FALSE.
        LPLOT=.FALSE.
        LHPLOT=.FALSE.
        GOTO 405
      ENDIF
      IF (.NOT.LCMD(20)) GOTO 23
 951  CONTINUE
C      CALL SCALE                           ! NOT AVAILABLE IN BATCH MODE
      IF (LERROR) GOTO 10
      DO 952 I=0,7
      J=11+I
      IWKING(J)=ICHAR(ADISPG(I+1:I+1))
 952  CONTINUE
      IF (LGGRID) THEN
        LNOPLT=.FALSE.
        IF (LESCAP) THEN
          LRESS=.FALSE.
          DO 953 IL=2,5
          IF (LRESP(IL)) LRESS=.TRUE.
 953      CONTINUE
          IF (LRESS) LNOPLT=.TRUE.
           ELSE
            LGRID=.TRUE.
            LSUBTRACT=.FALSE.
        ENDIF
      ENDIF               
      RETURN
C
C     dotmap on/off
C      
C
 975  CALL MATCH(APAR9,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      IF (JUMP.EQ.1) LDOTMAP=.TRUE.
      IF (JUMP.EQ.2) LDOTMAP=.FALSE.
      GOTO 10
C      
 1000 FORMAT ('      --------- GRAPHICS module --------',/,
     &        ' <F1> = Help')
 1001 FORMAT (/,'      ----------- GRID module ----------',/,
     &        ' <F1> = Help')     
 1010 FORMAT (' VIEWPORT ',4(I5,2X)/
     &        ' max. coordinates: ',I5,',',I5)
 1020 FORMAT (' WINDOW ',4(F14.3,2X),
     &'                                  (all,cursor,select,save,delete)
     &')
 1021 FORMAT (' WINDOW ',2(F11.1,1X),F11.3,1X,2(F11.1,1X),F11.3)
 1030 FORMAT (' AQUIFER ',I1)
 1040 FORMAT (' HORIZONTALPOINTS ',I4)
 1045 format (' PLOT LEAKAGE     (potentials, streamlines, discharge,',
     &' flownet,',/,
     &        '                   heads, leakance)')
 1046 format (' PLOT LEAKANCE    (potentials, streamlines, discharge,',
     &' flownet,',/,
     &        '                   heads, leakage)')     
 1051 FORMAT (' PLOT HEADS       (potentials, streamlines, discharge,',
     &' flownet,',/,
     &        '                   leakage, leakance)')
 1052 FORMAT (' PLOT POTENTIALS  (heads, streamlines, discharge,',
     &' flownet,',/,
     &        '                   leakage, leakance)')
 1053 FORMAT (' PLOT STREAMLINES (heads, potentials, discharge,',
     &' flownet,',/,
     &        '                   leakage, leakance)')
 1054 FORMAT (' PLOT DISCHARGE   (heads, potentials, streamlines,',
     &' flownet,',/,
     &        '                   leakage, leakance)')
 1055 FORMAT (' PLOT FLOWNET     (heads, potential, streamlines,',
     &' discharge,',/,
     &        '                   leakage, leakance)')     
 1057 FORMAT (' DOTMAP ON        (off)')
 1058 FORMAT (' DOTMAP OFF       (on)')     
 1061 FORMAT (' DEVICE SCREEN    (printer, plotter, slides)')
 1062 FORMAT (' DEVICE PRINTER   (screen, plotter, slides)',/,
     &' (You will be prompted for POSTSCRIPT filename at the start of',
     &' plotting.)')
 1063 FORMAT (' DEVICE PLOTTER   (screen, printer, slides)',/,
     &' (You will be prompted for HPGL filename at the start of',
     &' plotting.)') 
 1064 FORMAT (' DEVICE HIPLOT    (DMP40 driver)')
 1065 FORMAT (' DEVICE SLIDES    (screen, printer, plotter)',/,
     &' (Shoot slides from graphics on the screen.)')
 1071 FORMAT (' CHART LARGE')
 1072 FORMAT (' CHART SMALL')
 1080 FORMAT (' SPEED ',A1,'  inches/second')
 1090 FORMAT (' STEPSIZE ',A1,' (1,2,3...0.001",0.005",0.1mm)')
 1100 FORMAT (' LINE ',A1)
 1110 FORMAT (' MONITOR SINGLE           (dual)')
 1120 FORMAT (' MONITOR DUAL             (single)')
 1130 FORMAT (' DISPLAY MONOCHROME       (color)')
 1140 FORMAT (' DISPLAY COLOR            (monochrome)')
 1150 FORMAT (' PALETTE ',I2,' (1-4) YELLOW')
 1152 FORMAT (' PALETTE ',I2,' (1-4) RED')
 1154 FORMAT (' PALETTE ',I2,' (1-4) GREEN')
 1156 FORMAT (' PALETTE ',I2,' (1-4) BLUE')
 1160 FORMAT (' MINUSGRID                (filename)')
 1170 FORMAT (' TMARKER  DOT             (plus,asterisk,square,cross,',
     &                                                   'diamond)')
 1175 FORMAT (' TMARKER  PLUS            (dot,asterisk,square,cross,',
     &                                                   'diamond)')
 1180 FORMAT (' TMARKER  ASTERISK        (dot,plus,square,cross,',
     &                                                   'diamond)')
 1185 FORMAT (' TMARKER  SQUARE          (dot,plus,asterisk,cross,',
     &                                                   'diamond)')
 1190 FORMAT (' TMARKER  CROSS           (dot,plus,asterisk,square,',
     &                                                   'diamond)')
 1191 FORMAT (' TMARKER  DIAMOND         (dot,plus,asterisk,square,',
     &                                                   'cross)')
 1192 FORMAT (' SAVE                     (filename)')
 1193 FORMAT (' LOAD                     (filename)')
 1194 FORMAT (' SURFER                   (filename)')
 1195 FORMAT (' <Esc> or QUIT',/,' <F2> or GENERATE')
 1197 FORMAT (' <F2> or GO',/,' <Esc> or QUIT')
 2000 FORMAT (' ***ILLEGAL OR MISSING PARAMETER(S) in layout:',/,
     &        ' ',80A1)
 2001 FORMAT (' ***ILLEGAL OR MISSING PARAMETER(S) in contour graphics:'
     &       ,/,' ',80A1)
 2002 FORMAT (' ***ILLEGAL OR MISSING PARAMETER(S) in grid module:',/,
     &        ' ',80A1)
 2003 FORMAT (' ***ILLEGAL OR MISSING PARAMETER(S) in window module:',/,
     &        ' ',80A1)
 2010 FORMAT (' >')
 2015 FORMAT (' ***ERROR: no features to adjust window to, ',
     &'command ignored!',/)
 2030 FORMAT (' ***ERROR: no place in buffer for another window (max.='
     &,I2,').'/
     &'  Delete one or more windows before saving this window.',/)
 2040 FORMAT (' ***ERROR: window: ',4F11.1,/,
     &' not found in buffer, no window has been deleted.',/)
 2051 FORMAT (' ***ERROR: no window available in buffer.',/)
 2052 FORMAT (' Position the cursor on a corner of the window to be',
     &' selected and press <Enter>.')          
 2500 FORMAT (80A1)
 3000 FORMAT (' ***ILLEGAL COMMAND in layout:',/,
     &        ' ',80A1)
 3001 FORMAT (' ***ILLEGAL COMMAND in contour graphics:',/,
     &        ' ',80A1)
 3002 FORMAT (' ***ILLEGAL COMMAND in grid module:',/,
     &        ' ',80A1)
 3003 FORMAT (' ***ILLEGAL COMMAND in window module:',/,
     &        ' ',80A1)
 3500 FORMAT (' ***DEVICE must be SCREEN in CURSOR mode!',/,
     &' Press <Enter> to select screen and continue, type QUIT to', 
     &' abort.')
 4000 FORMAT ('+Position cursor at LOWER LEFT corner of window; ',
     &'press <Enter>')
 4050 FORMAT ('+Position cursor at UPPER RIGHT corner of window; ',
     &'press <Enter>')
 4070 FORMAT (' Current window has been saved: ',
     &I2,' windows in buffer.',/,' Buffer capacity is ',I2,/)
 4080 FORMAT (' Current window has been deleted from the buffer: ',
     &I3,' windows left.'/)
 5000 FORMAT (' ***ERROR: device must be set to screen to use cursor.')     
 5050 FORMAT (' ****WARNING: time is set to ',G11.4,/,
     &' Do not grid streamlines when transient wells are present!',/)
 5055 FORMAT (' ***ERROR: grid is not compatible with solution!',/)
 5060 FORMAT (' ***ERROR: grids incompatable, cannot subtract.',/)     
 5061 FORMAT (' ***ERROR: cannot subtract flownets.',/)
 5062 FORMAT (' ***ERROR: substraction of vertical sections is not',
     &        ' supported.',/)
 5063 FORMAT (' ***ERROR: you have already subtracted a grid!',/,
     &' First LOAD or GENERATE a new grid before subtracting a grid.',/)     
 5064 FORMAT (' ***ERROR: no grid to subtract from!',/,
     &' First LOAD or GENERATE a new grid before subtracting a grid.',/)      
 5065 FORMAT (' Grid in file ',A16,' subtracted from current grid.')
 5070 FORMAT (' ***ERROR in grid module: max. hor. points is ',I4,/)
 6000 FORMAT (' ***ERROR in grid module: flownets cannot be saved.',/)
      END
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE PLWIND (CZ)
c
c -------------------------------------------------------------------------------
c
C
C     Routine compares CZ to existing minimum and maximum window
C     parameters and updates these if CZ is outside the maximum window
c     (to zoom to extent)
C
      IMPLICIT NONE
      REAL(8) RXIN,RYIN
      COMPLEX(8) CZ
      INCLUDE 'GRID.INC'
      INCLUDE 'LUSYS.INC'
C
      RXIN=REAL(CZ)
      RYIN=AIMAG(CZ)
      RXMN=MIN(RXMN,RXIN)
      RXMX=MAX(RXMX,RXIN)
      RYMN=MIN(RYMN,RYIN)
      RYMX=MAX(RYMX,RYIN)
      RETURN
      END
c
c -------------------------------------------------------------------------------
c
      REAL(8) FUNCTION RFGRTOL()
c
c -------------------------------------------------------------------------------
c
C
C     Function returns a tolerance factor which is 0.001 (1/1000)
C     of the minimum window coordinates.
C     It is important, therfore, that the window parameters are set
C     early on (preferably first) when entering data.
C
      IMPLICIT NONE
      REAL(8) RDX,RDY,RDEL,AMIN1
      INCLUDE 'GRID.INC'
      INCLUDE 'LUSYS.INC'
C
      RFGRTOL=0.0
      RDX=RX2-RX1
      RDY=RY2-RY1
      IF (RDX.LT.0.OR.RDY.LT.0) THEN
        WRITE (ILUER,1000)
        RETURN
      ENDIF
      RDEL=MIN(RDX,RDY)
      RFGRTOL=RDEL/1000.0
      IF (RFGRTOL.EQ.0.0) WRITE (ILUER,2000)
      RETURN
 1000 FORMAT (' ***ERROR in RFGRTOL: illegal window parameters.')
 2000 FORMAT (' ***WARNING: RFGRTOL (spatial tolerance) equal to zero.')
      END
c
c -------------------------------------------------------------------------------
c
      SUBROUTINE CIRCLINE (RX0,RY0,RAD0,CZ01,CZ02,CZ1,CZ2,LFAIL)
c
c -------------------------------------------------------------------------------
c
C
C     Routine calculates the intersection points of a horizontal
C     circle with a vertical plane (the vertical window)
C     Called in DIPLOT and W3PLOT
C
      IMPLICIT NONE
      LOGICAL LFAIL
      REAL(8) RX0,RY0,RAD0,RML,RY0L,RAL,RBL,RCL,ROOTL,SQROOT,
     &        RXX1,RYY1,RXX2,RYY2
      COMPLEX(8) CZ01,CZ02,CZ1,CZ2,CMPLX
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      LFAIL=.FALSE.
      IF (.NOT.L3DVER) THEN
        WRITE (ILUER,1000)
        LFAIL=.TRUE.
        RETURN
      ENDIF
      RML=(R3DY2-R3DY1)/(R3DX2-R3DX1)
      RY0L=R3DY2-RML*R3DX2
      RAL=1.0+RML*RML
      RBL=RML*(RY0L-RY0)-2.0*RX0
      RCL=RX0*RX0+(RY0L-RY0)*(RY0L-RY0)-RAD0*RAD0
      ROOTL=RBL*RBL-4.0*RAL*RCL
      IF (ROOTL.LE.0.0) THEN
        LFAIL=.TRUE.
        RETURN
      ENDIF
      RXX1=(-RBL-SQROOT(ROOTL))/(2.0*RAL)
      RXX2=(-RBL+SQROOT(ROOTL))/(2.0*RAL)          
      RYY1=RML*RXX1+RY0L
      RYY2=RML*RXX2+RY0L
      CZ1=CMPLX(RXX1,RYY1)
      CZ2=CMPLX(RXX2,RYY2)
C                       generate "shadow" points of circle on plane
      RX0=-RBL/(2.0*RAL)
      RY0=RML*RX0+RY0L
      RBL=RML*(RY0L-RY0)-2.0*RX0
      RCL=RX0*RX0+(RY0L-RY0)*(RY0L-RY0)-RAD0*RAD0
      ROOTL=RBL*RBL-4.0*RAL*RCL
      RXX1=(-RBL-SQROOT(ROOTL))/(2.0*RAL)
      RXX2=(-RBL+SQROOT(ROOTL))/(2.0*RAL)          
      RYY1=RML*RXX1+RY0L
      RYY2=RML*RXX2+RY0L
      CZ01=CMPLX(RXX1,RYY1)
      CZ02=CMPLX(RXX2,RYY2)      
      RETURN
 1000 FORMAT (' ***ERROR in CIRCLINE: no vertical window set.')            
      END
c
c ------------------------------------------------------------------------------
c
      subroutine GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
c
c     Routine returns parameters from common main
c
      implicit none
      LOGICAL lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut
      CHARACTER(8) aBasenameOut
      CHARACTER(16) aDateTimeOut
      integer nsolOut
      include 'main.inc'
c
      lsolOut=lsol
      loadsolOut=loadsol
      linalreadyOut=linalready
      lErrorReportOut=lErrorReport
      lDirectfromDiskOut=lDirectFromDisk
      aBasenameOut=aBasename
      aDateTimeOut=aDateTime
      nsolOut=nsol
      return
      end subroutine



