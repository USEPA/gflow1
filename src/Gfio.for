C     Last change:  HMH  17 Nov 2008    8:40 am
c     This file contains the following routines or functions
c
c     BLOCK DATA SETLUSYS   set logical units
c     SUBROUTINE LUOUT      writes common block LUSYS (formatted) for debugging purposes
c     SUBROUTINE CHKOUT     writes common block CHCK (formatted) for debugging purposes
c     SUBROUTINE CM3DOUT    writes common block COM3D (formatted) for debugging purposes
c     SUBROUTINE INITIO     reads and writes some initialization data to and from disc
c     SUBROUTINE DATIO      driver for creating a GFLOW input script file (.dat)
c     SUBROUTINE STDATIO    writes header to a .dat file (called in DATIO))
c     SUBROUTINE ENDATIO    writes trailer to a .dat file (called in DATIO))
c     SUBROUTINE WINDATIO   writes windows parameters to a .DAT file
c     subroutine extract    driver to generate the *.xtr file (extract file)
c     SUBROUTINE GFIO       reads or writes contents of MAIN.INC to .sol file
c     SUBROUTINE CM3DIO     reads or writes contents of COM3D.INC common blocks to .SOL file
c     SUBROUTINE MNOUT      writes the contents of the common block MAIN (formatted) for debugging purposes
c     SUBROUTINE DEBUG      driver to write common blocks (formatted) for debugging purposes
C     SUBROUTINE MATOUT     writes the matrix and associated arrays (formatted) for debugging purposes
c     SUBROUTINE HALT       prints a diagnostic in GFLOW.OPS and stops program execution
c     SUBROUTINE BIO        driver for generating the binary .sol file (solution file)
c     SUBROUTINE CRFILU     create and/or open a file of a specified type (in argument) with prompt for name
c          ENTRY CRAFIL     create and/or open a file of a specified type (in argument) without prompt for name
c     SUBROUTINE OPFILU     open an existing file with prompting for name
c          ENTRY OPAFIL     open an existing file without prompting for name
c     SUBROUTINE FILECH     check for proper file extension
c     subroutine trywrite   try write operation to ensure that file is write enabled
c     SUBROUTINE GETFNUL    get filename out of ALINE maintaining upper and lower characters
c          ENTRY GETFN      get filename out of ALINE
c     SUBROUTINE GFSURF     write .BLN and .GRD file
c     SUBROUTINE HELP       for viewing help files when using GFLOW interactively (calls VIEW)
c     SUBROUTINE VIEW       views GFLOW help files (called by HELP)
c     SUBROUTINE MATCH      match input string with command words
c     SUBROUTINE VAR        driver for IVAR and RVAR routines
c     FUNCTION IVAR         read an integer out of the input line
c     FUNCTION RVAR         read a real out of the input line
c     SUBROUTINE TIDY       remove blanks and separators from input line and setup pointer array
c     SUBROUTINE TIDYUL     like TIDY, but now maintaining upper and lower case
c     LOGICAL FUNCTION LSEP .TRUE. if ACHARR is a separator. Converts all uppercases to lower cases
c              ENTRY LSEPUL like LSEP, but maintains upper and lower case
c     FUNCTION AFUPC        changes lower case to upper case
c     SUBROUTINE INLINE     reads ALINE and sets pointers (calls TIDY)
c     COMPLEX FUNCTION CVAR return a complex number from ALINE
c     SUBROUTINE SWITCH     driver for SWTCHS
c     SUBROUTINE SWTCHS     reassign logical units for input and output functions
c     SUBROUTINE LUSWAP     swap lugical units, called by SWTCHS
c     SUBROUTINE WriteMatrix  write coefficient matrix to disk (no coefficient corrections)
c     SUBROUTINE LoadMatrix   load matrix from disk, still to be corrected
c     SUBROUTINE WriteDecompMatrix  write decomposed matrix to disk (full system)
c     SUBROUTINE LoadDecompMatrix   load decomposed matrix from disk (reduce with Sherman Morrison)
c
c
c
C
C     NOTE: grid related IO routines are in GFGRID.FOR
c           trace related IO routines are in GFTRACE.FOR
C
c
C
C --------------------------------------------------------------------------
C    
      BLOCK DATA SETLUSYS
C
C --------------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'LUSYS.INC'
      DATA ILUIN,ILUOUT,ILUME,ILUER,ILUTMP,ILUPL,ILUECH,ILUDUM
     .    /5,    8,     13,   7,    4,     1,    9,     2*0/
C
C     NOTE: logical unit 2 is reserved for general file io
c           logical unit 6 is reserved for default output, do not assign
c           logical unit 10 is reserved for coefficient matrix io
c           logical unit 11 is reserved for decomposed matrix io
c
C     
      DATA LECHO,LUCON,LUOUTFILE /.FALSE.,.TRUE.,.FALSE./
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE LUOUT (ICALL)
C
C --------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) ICALL
      INCLUDE  'LUSYS.INC'
      WRITE (ILUOUT,500) ICALL
      WRITE (ILUOUT,1000) LECHO,LUCON,LUOUTFILE
      WRITE (ILUOUT,2000) ILUIN,ILUOUT,ILUME,ILUER
      WRITE (ILUOUT,3000) ILUTMP,ILUPL,ILUECH,ILUDUM
      RETURN
C
 500  FORMAT (' LUOUT: ICALL=',I3)
 1000 FORMAT (' LUOUT: LECHO,LUCON,LUOUTFILE ',3L3)
 2000 FORMAT (' LUOUT: ILUIN,ILUOUT,ILUME,ILUER ',4I3)
 3000 FORMAT (' LUOUT: ILUTMP,ILUPL,ILUEC,ILUDUM(2) ',5I3)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE CHKOUT (ICALL)
C
C --------------------------------------------------------------------------
C
C
C     Routine outputs data in common block /CHCK/
C
      IMPLICIT NONE
      INTEGER(4) I,ICALL,ICOUNT
      CHARACTER(1) ADUM
      INCLUDE 'CHKCOM.INC'
      INCLUDE 'LUSYS.INC'
C
      WRITE (ILUOUT,1000) NPIEZ,NPZMX
      IF (NPIEZ.GT.0) THEN
      WRITE (ILUOUT,1500)
      ICOUNT=0
      DO 10 I=1,NPIEZ
      ICOUNT=ICOUNT+1
      IF (ICOUNT.EQ.20) THEN
      WRITE (ILUOUT,3000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,3000) ICALL
      READ (ILUIN,4000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN      
      ICOUNT=0
      ENDIF
      WRITE (ILUOUT,2000) CZPIEZ(I),RHPIEZ(I),APZLAB(I)
 10   CONTINUE
      ENDIF
      WRITE (ILUOUT,3000) ICALL
      IF (LUOUTFILE) WRITE (ILUME,3000) ICALL
      READ (ILUIN,4000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      RETURN
 1000 FORMAT (' CHKOUT: NPIEZ,NPZMX ',2I5)
 1500 FORMAT (/'     X         Y             PHI           LABEL',/)
 2000 FORMAT (2G11.4,2X,G14.7,2X,A16)
 3000 FORMAT (' CHKOUT: ICALL = ',I3,
     &' Press <Enter> to continue or type R to abort.')
 4000 FORMAT (A1)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE CM3DOUT (ICALL)
C
C --------------------------------------------------------------------------
C
C
C     Routine writes the contents of common /COM3D/
C
      IMPLICIT NONE
      INTEGER(4) ICALL,I
      CHARACTER(1) ADUM
      INCLUDE 'COM3D.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'LUSYS.INC'
C
      WRITE (ILUOUT,1000) R3DH,R3DK,R3DPOR
      WRITE (ILUOUT,1001) R3DSTO,R3DZ
      WRITE (ILUOUT,1010) (R3DCPZ(I),I=1,NY)
      WRITE (ILUOUT,1011) RETARDATION,RHALFLIFE
      WRITE (ILUOUT,2000) ICALL      
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      WRITE (ILUOUT,1020) L3DREV,L3DEND,L3DPL,L3DHOR,L3DVER
      WRITE (ILUOUT,1030) R3DX1,R3DY1,R3DZ1
      WRITE (ILUOUT,1040) R3DX2,R3DY2,R3DZ2
      WRITE (ILUOUT,2000) ICALL      
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      RETURN
C
 1000 FORMAT (' CM3DOUT: R3DH,R3DK,R3DPOR',3(2X,G11.4))
 1001 FORMAT (' CM3DOUT: R3DPOR,R3DZ',2(2X,G11.4))
 1010 FORMAT (20(' CM3DOUT: R3DCPZ(I)',5(2X,G11.4)/))
 1011 FORMAT (' CM3DOUT: RETARDATION,RHALFLIFE ',2G14.7)
 1020 FORMAT (' CM3DOUT: L3DREV,L3DEND,L3DPL,L3DHOR,L3DVER',5(2X,L2))
 1030 FORMAT (' CM3DOUT: R3DX1,R3DY1,R3DZ1',3(2X,G11.4))
 1040 FORMAT (' CM3DOUT: R3DX2,R3DY2,R3DZ2',3(2X,G11.4))
 2000 FORMAT (' CM3DOUT: ICALL=',I3,' press <Enter> to continue.')
 3000 FORMAT (A1)
      END         
C
C --------------------------------------------------------------------------
C
      SUBROUTINE INITIO (ICODE,ILU)
C
C --------------------------------------------------------------------------
C
C
C     Routine reads and writes some initialization data to and from disc
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR,NWORD
      INCLUDE 'LUSYS.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'DWMN.INC'
C
      IF (ICODE.EQ.41) WRITE (ILUME,1000)
      IF (ICODE.EQ.73) WRITE (ILUME,2000)
      CALL BUFIN (ICODE,ILU,IERR)
      CALL BUFIO4 (NGRWIN,1,ILU,ICODE,IERR)
      IF (NGRWIN.GT.0) THEN
      NWORD=4*NGRWIN
      CALL BUFIOC (CZWIN,NWORD,ILU,ICODE,IERR)
      ENDIF
      CALL BUFIOL (LCOLOR,2,ILU,ICODE,IERR)
      CALL BUFIO4 (IX1,4,ILU,ICODE,IERR)
      CALL BUFEX (ICODE,ILU,IERR)
      IF (IERR.NE.0) WRITE (ILUER,2000) IERR
C     Check for legal viewport parameters and correct when necessary
      IF (IX2.GT.NDOTX) IX2=NDOTX
      IF (IY2.GT.NDOTY) IY2=NDOTY      
      RETURN
 1000 FORMAT (' Writing initialization file.')      
 2000 FORMAT (' Reading initialization file.')
 3000 FORMAT (' ***ERROR: in INITIO, IERR=',I4)      
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE DATIO
C
C --------------------------------------------------------------------------
C
C
C     Routine generates an input file ------.dat based on the current
C     data in the program
C
      IMPLICIT NONE
      INTEGER(4) ILU
      LOGICAL LRET
      INCLUDE 'GRID.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
C
      ILU=2
C
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      WRITE (ILUME,1000)
      CALL CRFILU (ILU,2,LRET,'.DAT')
      IF (LRET) THEN
       WRITE (ILUME,4000)
       RETURN
      ENDIF
      WRITE (ILU,2500)
      CALL STDATIO (ILU)
      CALL WINDATIO(ILU)
      CALL GVDATIO (ILU)
      CALL PDDATIO (ILU)
      CALL DBDATIO (ILU)
      CALL WLDATIO (ILU)
      CALL LSDATIO (ILU)
      CALL TWDATIO (ILU)      
      CALL W3DATIO (ILU)
      CALL DIDATIO (ILU)
      CALL ENDATIO (ILU)
      CLOSE (ILU)
      WRITE (ILUME,3000) AFILE
      RETURN
 1000 FORMAT ('     ----- Writing Input Data File -----')
 2000 FORMAT (' *Date: ',A8,2X, 'Time: ',A11,/,
     &        ' *Title: ',A16)
 2500 FORMAT (' * Written by GFLOW1 version 1.1 ')     
 3000 FORMAT (' Input data file ',A16,' has been written.')
 4000 FORMAT (' No input data file has been written.')
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE STDATIO (ILU)
C
C --------------------------------------------------------------------------
C
C
C     Routine writes header to a ".dat" file.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      INTEGER(4) ILU,I,ICHAR
      INCLUDE 'GRID.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
C
      WRITE (ILU,1000)
      WRITE (ILU,2000) RMODORIGX,RMODORIGY,IMODORCODE
      DO 10 I=1,16
      IF (ICHAR(ATITLE(I:I)).GT.32) GOTO 15
  10  CONTINUE
      RETURN      
  15  WRITE (ILU,3000) ATITLE
      RETURN
C     
 1000 FORMAT (' error error.log',/,
     &        ' yes',/,
     &        ' message nul',/,
     &        ' echo con',/,
     &        ' quit')
 2000 FORMAT (' modelorigin ',G21.14,1X,G21.14,1X,I2)
 3000 FORMAT (' title ',A16)      
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE ENDATIO (ILU)
C
C --------------------------------------------------------------------------
C
C
C     Routine writes trailer to a ".dat" file.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      INTEGER(4) ILU
      INCLUDE 'LUSYS.INC'
C
      WRITE (ILU,1000)
      RETURN
C     
 1000 FORMAT (' switch',/,
     &        ' error con',/,
     &        ' message con',/,
     &        ' echo off',/,
     &        ' input con')
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE WINDATIO (ILU)
C
C --------------------------------------------------------------------------
C
C
C     Routine writes window to a ".dat" file.
C     Current 3D window is written (if L3DPL).
C     Maximum window is written looking at window all and current window.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      INTEGER(4) ILU
      INCLUDE 'GRID.INC'
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
C
      IF (RXMN.GE.RXMX.OR.RYMN.GE.RYMX) RETURN  ! zero area in window
      IF (L3DPL) THEN         ! first check for valid window
      IF (R3DX1.GT.R3DX2.OR.R3DY1.GT.R3DY2.OR.R3DZ1.GT.R3DZ2) RETURN
      IF (R3DX1.EQ.R3DX2.AND.R3DY1.EQ.R3DY2.AND.R3DZ1.EQ.R3DZ2) RETURN
      WRITE (ILU,1000) R3DX1,R3DY1,R3DZ1,R3DX2,R3DY2,R3DZ2
      ELSE  
      WRITE (ILU,2000) RX1,RY1,RX2,RY2 ! write current window coordinates
      ENDIF
      RETURN
C     
 1000 FORMAT (' layout',/,
     &        ' window ',2(2F11.1,F11.3),/,
     &        ' quit')
 2000 FORMAT (' layout',/,
     &        ' window ',4G14.7,/,
     &        ' quit')
      END
C
C --------------------------------------------------------------------------
C
      subroutine extract (lsol)
C
C --------------------------------------------------------------------------
C
c
c     Routine allows all aquifer parameters plus head and velocity
c     at x,y,z to be written to the screen or a specified file: FILENAME.XTR
c     Routine also allows relevant data from the well, line sink and inhomogeneity modules
c     to be extracted.
c
      implicit none
      INTEGER(4) ilu,jump,ipar,ipoint
      logical lfile,lret,lbad,lsol,lfirst,lflux,ldboutside
      REAL(8) rfpor,rfperm,rfhght,rfdbgm,rfpds,rfpdsbottom,rb,
     &     rp,rk,rn,rnb,rh,rdum,rvi(3),rfbase,rfhead,rvar,
     &     rnormalflow,rfnormalflow,rnormalflownum,rfnumnf,
     &     rx,ry,rfinterface,rflsnearflow,rflsnearlake,
     &     rflkleakage,rflkresistance,rflklowerhead,rhl,rc,
     &     rflkrecharge,rflksubcellleakage
      COMPLEX(8) cz,cdum,cvar,cz1,cz2
      character(1) aword(73),apar(11)
      character(16) afilename,alabel

      include 'com3d.inc'
      include 'lusys.inc'
      include 'match.inc'
      data aword /'?',' ',
     &            'F','I','L','E',' ',
     &            'Q','U','I','T',' ',
     &            'L','I','N','E',' ',
     &            'S','I','N','K',' ',
     &            'S','D','3','D',' ',
     &            'W','E','L','L',' ',
     &            'P','P','W','E',' ',
     &            'T','H','E','I',' ',
     &            'I','N','H','O',' ',
     &            'L','E','A','K',' ',
     &            'F','L','U','X',' ',
     &            'I','N','T','E',' ',
     &            'G','A','G','E',' ',
     &            'L','A','K','E',' ',
     &            ATERM/
      data apar  /'H','E','A','D',' ',
     &            'D','I','S','C',' ',
     &           ATERM/
c
      lerror=.false.
      lmiss=.false.
      lfile=.false.
      lfirst=.true.
      lflux=.false.
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
   10 if (lerror.or.lmiss) write (iluer,1000) aline2
      lerror=.false.
      lmiss=.false.
      if (.not.lsol) write (iluer,1500)
      if (lfile) then
        if (lucon) write (ilume,2001) afilename
      else
        if (lucon) write (ilume,2000)
      endif
      call inline
   11 call match (aword,1,jump,lbad)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      if (lbad.and..not.lflux) goto 100
      if (lbad.and.lflux) GOTO 810
      lfirst=.true.
   12 goto (200,300,400,500,550,580,600,650,675,700,710,
     &      800,850,900,910), jump
c
c     read coordinates and label
c      
 100  if (lflux) GOTO 400
      lbad=.false.
      cdum=cvar(1)
      if (lerror) goto 10
      cz=cdum
      rdum=rvar(3)
      if (lerror) then
        rdum=rfbase(cz)+rfhght(cz)
        call getfn(3)
      else
        call getfn(4)
      endif
      lerror=.false.
      lmiss=.false.
      r3dz=rdum
c                 write extract data
      rx=REAL(cz)
      ry=AIMAG(cz)
      if (ldboutside(rx,ry)) then
       rh=-9999.0
      else
       rh=rfhead(cz)
       r3dz=MIN(r3dz,rh)   ! drop point down to water table.
      end if
      rn=rfpds(cz)+rfdbgm(cz)+rflkrecharge(cz)
      rnb=rflkleakage(cz)+rflksubcellleakage(cz)    ! replaced rfpdbottom of module PD (ponds)
      rk=rfperm(cz)
      rp=rfpor(cz)
      rb=rfbase(cz)
      rhl=rflklowerhead(cz)
      rc=rflkresistance(cz)
      call veloc(cz,rvi)
      if (lfile) then
       if (lfirst) then
         write (ilu,5000)
         lfirst=.false.
       endif
       write (ilu,3000) cz,r3dz,rp,rk,rb,rn,rnb,rh,rhl,rc,rvi,afile              ! file
       if (lucon)
     &  write (ilume,4000) cz,r3dz,rp,rk,rb,rn,rnb,rh,rhl,rc,rvi,afile  ! screen
       else
       if (lucon)
     &  write (ilume,4000) cz,r3dz,rp,rk,rb,rn,rnb,rh,rhl,rc,rvi,afile             ! screen
      endif
      goto 10
c
c     help
c
 200  if (lflux) GOTO 400
      AFILE='                '
      afile(1:11)='extract.hlp'
C      call help                     ! NOT AVAILABLE IN BATCH MODE
      goto 10
c
c     file
c
 300  if (lflux) GOTO 400
      call getfn(2)
      if (lerror.or.lmiss) goto 10
      ilu=11
      call crafil (ilu,-5,lret,'.XTR')  ! record length = 256
      if (lret) goto 10
      lfile=.true.
      afilename=afile
      if (.not.lsol) write (ilu,1501)
      goto 10
c
c     quit
c
 400  continue
      if (lflux) then ! quit to terminate coordinate pairs
      lflux=.false.
        if (ipoint.gt.1) then
          write (ilu,8000)
          write (ilu,8010) alabel,rnormalflow                            ! ****** reactivate to supress debuging
c          write (ilu,8011) alabel,rnormalflow,rnormalflownum              ! ****** block to supress debugging
        else
          write (iluer,1600) ipoint
        END if
        if (jump.eq.3) then
        GOTO 10  ! was quit command, read a new command.
        else
        GOTO 12  ! was not quit command, return to proper command.
        end if
      end if
      write (ilu,9000)
      close (ilu)
      lfile=.false.
      return
c
C     Extract line sink data      <linesink head  /  linesink discharge>
C
 500  if (lflux) GOTO 400
      if (.not.lfile) then
        write (iluer,1502)
        goto 10
      endif
      call match (apar,2,ipar,lbad)
      lerror=lbad
      if (lerror) then
        lbad=.false.
        goto 10
      endif
      if (ipar.eq.1) call lsextracthead (ilu)
      if (ipar.eq.2) call lsextractdisch (ilu)
      if (lucon) write (ilume,6000) afilename
      goto 10
c
C     Extract two dimensional sink disc data  (2D)    <si  head / si discharge>
C
 550  if (lflux) GOTO 400
      if (.not.lfile) then
        write (iluer,1502)
        goto 10
      endif
      call match (apar,2,ipar,lbad)
      lerror=lbad
      if (lerror) then
        lbad=.false.
        goto 10
      endif
      if (ipar.eq.1) call pdextracthead (ilu)
      if (ipar.eq.2) call pdextractdisch (ilu)
      if (lucon) write (ilume,6050) afilename
      goto 10
c
C     Extract three dimensional sink disc data  (3D)       <sd3d  head / sd3d discharge>
C
 580  if (lflux) GOTO 400
      if (.not.lfile) then
        write (iluer,1502)
        goto 10
      endif
      call match (apar,2,ipar,lbad)
      lerror=lbad
      if (lerror) then
        lbad=.false.
        goto 10
      endif
      if (ipar.eq.1) call diextracthead (ilu)
      if (ipar.eq.2) call diextractdisch (ilu)
      if (lucon) write (ilume,6080) afilename
      goto 10
c
c     Extract well data                          <well head / well discharge>
c
 600  if (lflux) GOTO 400
      if (.not.lfile) then
        write (iluer,1502)
        goto 10
      endif
      call match (apar,2,ipar,lbad)
      lerror=lbad
      if (lerror) then
        lbad=.false.
        goto 10
      endif
      if (ipar.eq.1) call wlextracthead (ilu)
      if (ipar.eq.2) call wlextractdisch (ilu)
      if (lucon) write (ilume,7000) afilename
      goto 10
c
c     Extract partially penetrating well data            ppwell head /
c
 650  if (lflux) GOTO 400
      if (.not.lfile) then
        write (iluer,1502)
        goto 10
      endif
      call match (apar,2,ipar,lbad)
      lerror=lbad
      if (lerror) then
        lbad=.false.
        goto 10
      endif
      if (ipar.eq.1) call w3extracthead (ilu)
      if (ipar.eq.2) call w3extractdisch (ilu)
      if (lucon) write (ilume,7000) afilename
      goto 10
c
c     Extract Theis well data
c
 675  if (lflux) GOTO 400
      if (.not.lfile) then
        write (iluer,1502)
        goto 10
      endif
      call twextract (ilu)
      if (lucon) write (ilume,7000) afilename
      goto 10
c
c     Extract inhomogeneity data
c
 700  if (lflux) GOTO 400
      if (.not.lfile) then
        write (iluer,1502)
        goto 10
      endif
      call inhomextract(ilu)
      GOTO 10
c
c     Extract leakage data
c
 710  if (lflux) GOTO 400
      if (.not.lfile) then
        write (iluer,1502)
        goto 10
      endif
      call lkextract(ilu)
      GOTO 10
c
c     Extract flux across polyline (flux command followed by coordinate pairs)
c
 800  if (lflux) GOTO 400
      if (.not.lfile) then
        write (iluer,1502)
        goto 10
      endif
      lflux=.true.
      rnormalflow=0.0
      rnormalflownum=0.0
      ipoint=0
      alabel='NO_LABEL_GIVEN   '
      call getfn(2)
      if (lerror.or.lmiss) goto 10
      alabel=afile
      GOTO 10
 810  lbad=.false.
      cdum=cvar(1)
      if (lerror) goto 10
      ipoint=ipoint+1
      if (ipoint.gt.1) then
      cz1=cz2
      cz2=cdum
      rnormalflow=rnormalflow+rfnormalflow(cz1,cz2)
      rnormalflownum=rnormalflownum+rfnumnf(cz1,cz2)        ! ****** block to surpress debugging
      else
      cz2=cdum
      end if
      GOTO 10
c
c     write interface elevation
c
 850  cz=cvar(2)
      if (lerror) GOTO 10
      rdum=rfinterface (cz)
      write (ilu,8500)
      write (ilu,8501) cz,rdum
      if (lucon) write (ilume,8500)
      if (lucon) write (ilume,8501) cz,rdum
      GOTO 10
C
C     write gauging station flow
C
 900  cdum=cvar(2)
      call getfn(4)
      if (lerror) GOTO 10
      rdum=rflsnearflow (cdum)
      write (ilu,8600) cdum,rdum,afile
      GOTO 10
C
C     write lake stage
C
 910  cdum=cvar(2)
      call getfn(4)
      if (lerror) GOTO 10
      rdum=rflsnearlake (cdum)
      write (ilu,8700) cdum,rdum,afile
      GOTO 10
c
 1000 format (' ***ILLEGAL or MISSING PARAMETERS.'/,' ',80a1)
 1500 format (' ***WARNING: no solution, data may be in error!',/)
 1501 format ('* warning: no solution, data may be in error!')
 1502 format (' ***ERROR: no .XTR file open, use FILE command first!')
 1600 format (' ***ERROR:',i3,'point in FLUX command:',
     &        ' no flow calculated!')
 2000 format ('                  ---------- EXTRACT module ---------'/
     &        ' Available commands:',/,
     &        ' <F1> = Help',/,' FILE  (filename)',/,
     &        ' LINESINK (head/discharge)',/,
     &        ' WELL     (head/discharge)',/,
     &        ' INHOMOGENEITY            ',/,
     &        ' FLUX  (label)            ',/,
     &        ' (x)  (y)  [z]  [label]',/,' <Esc> or QUIT',/,' >')
 2001 format ('                  ---------- EXTRACT module ---------'/
     &        ' Available commands:',/,
     &        ' <F1> = Help',/,' FILE  (',a,')',/,
     &        ' LINESINK (head/discharge)',/,
     &        ' WELL     (head/discharge)',/,     
     &        ' (x)  (y)  [z]  [label]',/,' <Esc> or QUIT',/,' >')
 3000 format ('@',10(e14.7,','),4(e21.14,','),a)
 4000 format (' x,y,z ',3(e14.7,2x),/,
     &        ' porosity     = ',e14.7,' hydraulic conductivity = '
     &        ,e14.7,/,
     &        ' net recharge = ',e14.7,' leakage (bottom)       = '
     &        ,e14.7,/,
     &        ' head         = ',e14.7,/,
     &        ' lower head   = ',e14.7,/,
     &        ' resistance   = ',e14.7,/,
     &        ' Vx,Vy,Vz     = ',3(e14.7,2x),/,
     &        ' label        = ',a,//)
 5000 format ('*      x              y              z       ',
     &         '   porosity    hydr. conduct.   base elevation',
     &         ' net recharge  leakage (bottom)      head     ',
     &         ' lower head     resistance             ',
     &         '      Vx                    Vy                    Vz  ',
     &         '            label')
 6000 format (' Line sink data written to file: ',/,a,/)
 6050 format (' Sink Disc (2D) data to file: ',/,a,/)
 6080 format (' Sink Disc (3D) data to file: ',/,a,/)
 7000 format (' Well data written to file: ',/,a,/)
 7050 format(' Partially Penetrating Well data written to file: ',/,a)
 8000 format ('! flux_inspection_line_label   normal_flow      ',
     &        '  numerical_nf')
 8010 format (3x,a16,13x,',',e14.7)
 8011 format (3x,a16,13x,',',e14.7,e14.7)
 8500 format (' interface ')
 8501 format ('@',3(d14.7,','))
 8600 format ('@ gage, ',3(d14.7,',',2x),a16)
 8700 format ('@ lake_stage, ',3(d14.7,',',2x),a16)
 9000 format ('*  end of extract file')
c 8010 format (3x,a16,13x,',',e14.7,',',e14.7)
      end
C
C --------------------------------------------------------------------------
C
      SUBROUTINE GFIO (ICODE,ILU,RVERSION,ierr)
C
C --------------------------------------------------------------------------
C
C
C     Routine reads or writes contents of MAIN.INC to external file
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR
      REAL(8) RVERSION      
      INCLUDE 'MAIN.INC'
      INCLUDE 'LUSYS.INC'
C
      if (ierr.ne.0) return
      CALL BUFIOR (RVERSION,1,ILU,ICODE,IERR)
      CALL BUFIOL (LDISPL,2,ILU,ICODE,IERR)
      CALL BUFIO4 (NSOL,1,ILU,ICODE,IERR)
      CALL BUFIOR (RERMAX,1,ILU,ICODE,IERR)
      IF (RVERSION.EQ.1.0) RETURN
      IF (RVERSION.EQ.2.0) RETURN
      IF (RVERSION.EQ.3.0) RETURN
      IF (RVERSION.EQ.4.0) RETURN
      IF (RVERSION.EQ.5.0) RETURN
      IF (RVERSION.EQ.6.0) RETURN
      IF (RVERSION.EQ.7.0) RETURN
      CALL BUFIO4 (NOUTERLOOP,1,ILU,ICODE,IERR)
      IF (RVERSION.EQ.8.0) RETURN
      IF (RVERSION.EQ.9.0) RETURN
      CALL BUFIOA (aDateTime,1,ILU,ICODE,IERR)
C
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE CM3DIO (ICODE,ILU,RVERSION,ierr)
C
C --------------------------------------------------------------------------
C
C
C     Routine reads or writes contents of COM3D.INC common blocks to
C     an external file.
C     ICODE=41 write
C     ICODE=73 read
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR
      REAL(8) RVERSION
      INCLUDE 'COM3D.INC'
C
      if (ierr.ne.0) return
      CALL BUFIOR (R3DH,300,ILU,ICODE,IERR)
      CALL BUFIO4 (IMODORCODE,1,ILU,ICODE,IERR)
      CALL BUFIOR (RMODORIGX,2,ILU,ICODE,IERR)
      CALL BUFIOL (L3DREV,2,ILU,ICODE,IERR)
      CALL BUFIOR (R3DX1,6,ILU,ICODE,IERR)
      CALL BUFIOL (L3DPL,3,ILU,ICODE,IERR)
      IF (RVERSION.EQ.1.0) RETURN
      IF (RVERSION.EQ.2.0) RETURN
C
      RETURN
      END
C
C --------------------------------------------------------------------------
C     
      SUBROUTINE MNOUT (ICALL)
C
C --------------------------------------------------------------------------
C
C
C     Routine writes the contents of the common block /MAIN/
C
      IMPLICIT NONE
      INTEGER(4) ICALL
      CHARACTER(1) ADUM
      INCLUDE 'MAIN.INC'
      INCLUDE 'LUSYS.INC'
C
      WRITE (ILUOUT,1000) RPI2,RO2PI,ROPI
      WRITE (ILUOUT,2000) CPI,C2PI
      WRITE (ILUOUT,3000) CI,CIM
      WRITE (ILUOUT,4000) LSOL,NSOL,RERMAX,LDISPL
      WRITE (ILUME,5000) ICALL
      READ (ILUIN,1500) ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      RETURN
 1000 FORMAT (' MN: RPI2,RO2PI,ROPI ',3G14.7)
 1500 FORMAT (A1)
 2000 FORMAT (' MN: CPI,C2PI ',4G14.7)
 3000 FORMAT (' MN: CI,CIM   ',4G14.7)
 4000 FORMAT (' MN: LSOL,NSOL,RERMAX,LDISPL ',L2,2X,I3,2X,G11.4,2X,L2)
 5000 FORMAT (' MN: ICALL=',I3,' press <Enter> to continue.')
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE DEBUG      
C
C --------------------------------------------------------------------------
C
C
C     Subroutine writes all common variables 
C
      IMPLICIT NONE
      INTEGER(4) JUMP
      LOGICAL LBAD
      CHARACTER(1) APAR(81)
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
      DATA APAR /         
     .             'Q','U','I','T',' ', 
     .             'G','R','I','D',' ',
     .             'P','P','W','E',' ',
     .             'C','O','M','3',' ',
     .             'G','I','V','E',' ',
     .             'S','I','N','K',' ',
     .             'L','I','N','E',' ',
     .             'P','L','O','T',' ',
     .             'S','D','3','D',' ',
     .             'T','W','E','L',' ',
     .             'W','E','L','L',' ',
     .             'I','N','H','O',' ',
     .             'M','A','I','N',' ',
     .             'L','U','S','Y',' ',
     .             'C','H','E','C',' ',
     .             'L','E','A','K',' ',
     .             ATERM/   
C
      LERROR=.FALSE. 
      LMISS=.FALSE.
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
  10  IF (LERROR.OR.LMISS) WRITE (ILUER,400) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      WRITE (ILUME,100)
      CALL INLINE
      CALL MATCH(APAR,1,JUMP,LBAD)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LBAD) GOTO 15
      GOTO (10,13), JUMP
  13  WRITE (ILUME,300) ALINE2
      GOTO 10       
  15  GOTO (900,994,995,996,997,998,999,1000,1002,1004,1006,1008,
     &      1009,1010,1020,1030),JUMP
 900  RETURN
 994  CALL GRIDOUT (1)
      GOTO 10
 995  CALL W3OUT (1)
      GOTO 10
 996  CALL CM3DOUT (1)
      GOTO 10
 997  CALL GVOUT (1)
      GOTO 10  
 998  CALL PDOUT (1)
      GOTO 10
 999  CALL LSNKOUT (1)
      GOTO 10
 1000 CONTINUE
C      CALL PLOUT (1)     ! NOT AVAILABLE IN BATCH MODE
      GOTO 10 
 1002 CALL DIOUT (1)
      GOTO 10
 1004 CALL TWOUT (1)
      GOTO 10
 1006 CONTINUE
      CALL WLOUT (1)
      GOTO 10    
 1008 CALL DBOUT (1)
      GOTO 10
 1009 CALL MNOUT (1)
      GOTO 10      
 1010 CALL LUOUT (1)
      GOTO 10      
 1020 CALL CHKOUT (1)
      GOTO 10
 1030 CONTINUE
C      CALL LKOUT (1)          ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C            
 100  FORMAT (' Debug module commands:'/
     &' GRID'/' PPWEL'/' COM3D'/' GIVEN'/' SINKDISC'/' LINESINK'/' PLOT'
     &/' SD3D'/' TWELL'/' WELL'/' INHOMOGENEITY'/' MAIN'/' LUSYS'/
     &' CHECK'/' LEAKAGE'/' <Esc> or QUIT'/' >')
 300  FORMAT (' ***ILLEGAL COMMAND in debug module:',/,
     &        ' ',80A1)
 400  FORMAT (' ***ILLEGAL or MISSING PARAMETER(S) in debug module:',/,
     &        ' ',80A1)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE MATOUT (ICALL,RA,RB,CZI,CALPH,DRFAC,M,N,ITYPE)
C
C --------------------------------------------------------------------------
C
C
C     Subroutine prints matrix solution arrays
C
      IMPLICIT NONE
      INTEGER(4) ICALL,M,N,ITYPE,I,J
      REAL(8) RA,RB,DRFAC
      COMPLEX(8) CZI,CALPH
      CHARACTER(1) ADUM
      INCLUDE 'LUSYS.INC'
      DIMENSION RA(M,N),CZI(*),CALPH(*),RB(*),ITYPE(*),DRFAC(4,*)
      if (n.gt.10) then
        write (ilume,1001)
 1001 format (' Skipped MATOUT, matrix too large (larger than M X 10).')
        return
      end if
      WRITE (ilume,1100) ICALL
      IF (LUOUTFILE) WRITE (ILUME,1100) ICALL
      WRITE (ilume,1000) M,N
      IF (M.NE.0.AND.N.NE.0) THEN
      WRITE (ilume, 1500) (CZI(I),I=1,M)
      WRITE (ilume,1100) ICALL
      IF (LUOUTFILE) WRITE (ILUME,1100) ICALL
      WRITE (ilume,1700) (CALPH(I),I=1,M)
      WRITE (ilume,1100) ICALL
      WRITE (ilume,1800) ((DRFAC(J,I),J=1,4),I=1,M)
      WRITE (ilume,1100) ICALL
      IF (LUOUTFILE) WRITE (ILUME,1100) ICALL
      ENDIF
      IF (M.NE.0.AND.N.NE.0) THEN
c      WRITE (ilume,1200) (ITYPE(I),I=1,M)
c      DO 20 I=1,M
c      WRITE (ilume,2000)
c      DO 10 J=1,N
c      WRITE (ilume,2100) I,J  ! to verify matrix structure in output
c  10  CONTINUE
c  20  CONTINUE
c      write (ilume,2001) ((ra(i,j),j=1,m),i=1,n)
c 2001 format (7(D14.7,1x))
      write (ilume,2000)
      select case (n)
          case (1)
          do i=1,m
            write (ilume,2001) ra(i,1:n)
 2001      format (d14.7)
          end do
          case (2)
          do i=1,m
            write (ilume,2002) ra(i,1:n)
 2002      format (2(d14.7))
          end do
          case (3)
          do i=1,m
            write (ilume,2003) ra(i,1:n)
 2003      format (3(d14.7))
          end do
          case (4)
          do i=1,m
            write (ilume,2004) ra(i,1:n)
 2004      format (4(d14.7))
          end do
          case (5)
          do i=1,m
            write (ilume,2005) ra(i,1:n)
 2005      format (5(d14.7))
          end do
          case (6)
          do i=1,m
            write (ilume,2006) ra(i,1:n)
 2006      format (6(d14.7))
          end do
          case (7)
          do i=1,m
            write (ilume,2007) ra(i,1:n)
 2007      format (7(d14.7))
          end do
          case (8)
          do i=1,m
            write (ilume,2008) ra(i,1:n)
 2008      format (8(d14.7))
          end do
          case (9)
          do i=1,m
            write (ilume,2009) ra(i,1:n)
 2009      format (9(d14.7))
          end do
          case (10)
          do i=1,m
            write (ilume,2010) ra(i,1:n)
 2010      format (10(d14.7))
          end do
          case default
          write (ilume,2011) n
 2011 format (' number of columns out of range: ',i10)
      end select
c      DO 40 I=1,M
c      WRITE (ilume,2000)
c      DO 30 J=1,N
c      WRITE (ilume,2200) RA(I,J)  ! print matrix coefficients
c  30  CONTINUE
c  40  CONTINUE
      WRITE (ilume,1100) ICALL
      IF (LUOUTFILE) WRITE (ILUME,1100) ICALL
      WRITE (ilume,3000) (RB(J),J=1,M)
      WRITE (ilume,1100) ICALL
      IF (LUOUTFILE) WRITE (ILUME,1100) ICALL
      ENDIF
      WRITE (ilume,6000) ICALL     
 1000 FORMAT ( ' M,N=',2I3/)
 1100 FORMAT ('   MATOUT: ICALL=',I3)
 1200 FORMAT (10('  ITYPE(I)=',10(I3)/))
 1500 FORMAT (10('     CZI(I)=',(2G11.4)/))
 1700 FORMAT (10('   CALPH(I)=',(2G11.4)/))
 1800 FORMAT (10(' DRFAC(1-4,I)=',(4G11.4)/))
 2000 FORMAT (' Matrix RA:',/)
 2100 FORMAT ('&RA(',I3,',',I3,') ')
 2200 FORMAT ('&',E11.4,1X)
 3000 FORMAT (10('      RB(J)=',(G14.7)/))
 4000 FORMAT (A1)
 6000 FORMAT (' LEAVING ICALL',I1/)
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE HALT (NLINE)
C
C --------------------------------------------------------------------------
C
C
C     Routine prints a diagnostic in GFLOW.OPS and stops program execution.
C
      IMPLICIT NONE
      INTEGER(4) I,NLINE,itic,itic0
      CHARACTER*1 adum
      INCLUDE 'LUSYS.INC'
      WRITE (*,1000)
      WRITE (*,2000) (AMESS(I),I=1,NLINE)
      write (*,4000)
      !call timer  (itic0)
      do
      !call timer  (itic)
      if (itic-itic0.gt.500) GOTO 10
      end do
   10 OPEN (UNIT=11,FILE='GFLOW.OPS',STATUS='REPLACE')
      WRITE (11,1000)
      WRITE (11,3000) (AMESS(I),I=1,NLINE)
      CLOSE (11)
      STOP
 1000 FORMAT (//,' FATAL ERROR IN GFLOW1.EXE:',/)
 2000 FORMAT (' ',A60)
 3000 FORMAT (A60)
 4000 format (/,' WAIT for Box to close...')
      END SUBROUTINE
C
C --------------------------------------------------------------------------
C
      SUBROUTINE BIO (ILU,AIO)
C
C --------------------------------------------------------------------------
C
C
C     Routine performs a binary read or write operation to
C     an external file.
C      
      IMPLICIT NONE
      INTEGER(4) ILU,ICODE,IERR
      LOGICAL LREAD,LWRITE,LRET
      REAL(8) RCURRENTVERSION,RVERSION
      CHARACTER(5) AIO
      INCLUDE 'LUSYS.INC'
      include 'match.inc'
C           version 1.0  changed into 2.0  on 11/28/95
C           version 2.0  changed into 3.0  on 10/15/99
C           version 3.0  changed into 4.0  in 2000
C           version 4.0  changed into 5.0  on 06/19/01
C           version 5.0  changed into 6.0  on 12/29/01
C           version 6.0  changed into 7.0  on 3/14/02
C           VERSION 7.0  changed into 8.0  on 5/21/02
C           VERSION 8.0  changed into 9.0  on 7/22/02
c           version 9.0  changed into 10.0 on 9/1/02
C           version 10.0 changed into 11.0 on 7/15/03
c           version 11.0 changed into 12.0 on 2/8/05
c           version 12.0 changed into 13.0 on 7/21/05
c           version 13.0 changed into 14.0 on 3/14/06
c           version 14.0 changed into 15.0 on 4/26/06
c           version 15.0 changed into 16.0 on 7/18/06
c           version 16.0 changed into 17.0 on 11/7/06
      DATA RCURRENTVERSION /17.0/
C
      rversion=0.0
      LREAD=AIO.EQ.'READ '
      LWRITE=AIO.EQ.'WRITE'
      IF (.NOT.LREAD.AND..NOT.LWRITE) THEN
       WRITE (ILUER,1000) AIO
       RETURN
      ENDIF
      ICODE=0
      IF (LWRITE) THEN
        CALL CRFILU (ILU,1,LRET,'.SOL')
        ICODE=41
      ELSE
        CALL OPFILU (ILU,1,LRET,'.SOL')
        ICODE=73
      ENDIF
      IF (LRET) RETURN
      ierr=0
      CALL BUFIN (ICODE,ILU,IERR)
C ---------------------------begin program spec. statements
      IF (LWRITE) RVERSION=RCURRENTVERSION
      CALL GFIO (ICODE,ILU,RVERSION,ierr)
      if (rversion.ne.rcurrentversion) then
        WRITE (ILUER,2000) afile
      AMESS(1)='Solution file (*.sol) is incompatible with this solver.'
      AMESS(2)='Delete the *.sol file and run the solver again.'
        CALL HALT(2)  ! stop program execution for batch version
      end if
      CALL WLIO (ICODE,ILU,RVERSION,ierr)
      CALL LSIO (ICODE,ILU,RVERSION,ierr)
      CALL PDIO (ICODE,ILU,RVERSION,ierr)
      CALL DBIO (ICODE,ILU,RVERSION,ierr)
      CALL W3IO (ICODE,ILU,RVERSION,ierr)
      CALL GVIO (ICODE,ILU,RVERSION,ierr)
      CALL TWIO (ICODE,ILU,RVERSION,ierr)
      CALL DIIO (ICODE,ILU,RVERSION,ierr)
      CALL LKIO (ICODE,ILU,RVERSION,ierr)
      CALL GRIDIO (ICODE,ILU,RVERSION,ierr)
      CALL CM3DIO (ICODE,ILU,RVERSION,ierr)
C ---------------------------end program spec. statements
      CALL BUFEX (ICODE,ILU,IERR)
      CLOSE (ILU)
      if (ierr.ne.0) then
        write (iluer,3000)
      end if
      RETURN
  900 FORMAT (A1)      
 1000 FORMAT (' ***ERROR in BIO: illegal argument AIO=',A5)
 2000 format (' ***ERROR in BIO: ',a16,'.sol is incompatible ',/
     &        'with the current solver. It is an older solution file.')
 3000 format (' ***ERROR encountered in reading or writing *.sol file!')
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE CRFILU(ILU,ITYPE,LRET,AEXT)
C
C --------------------------------------------------------------------------
C
C
C     |ITYPE| = 1 unformatted file
C     |ITYPE| = 2 formatted file
C     |ITYPE| = 3 formatted random access, recl=31
C     |ITYPE| = 4 formatted file of record length 132
C     |ITYPE| = 5 formatted file of record length 256
C     ITYPE < 0 do not report IO
C
      IMPLICIT NONE
      INTEGER(4) ILU,ITYPE,ITYP,IERR
      LOGICAL LRET,LREP,LFILEOPEN
      CHARACTER(4) AEXT
      CHARACTER (LEN=7)  AWRITE
      INCLUDE 'DWMN.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
C
      LRET=.TRUE.
      CALL GETFN (2)
  5   IF (LERROR) RETURN
      IF (LMISS) THEN
      LMISS=.FALSE.
 10   IF (LMISS.OR.LERROR) THEN
      IF (LPLOTON.AND.LSINGL) THEN
      WRITE (ILUME,9009)
      WRITE (ILUER,9101)
      ELSE
      WRITE(ILUER,9001)
      ENDIF
      ENDIF
      LERROR=.FALSE.
      LMISS=.FALSE.
      IF (LPLOTON.AND.LSINGL) THEN
      WRITE (ILUME,9009)
      WRITE (ILUME,9006)
      ELSE
      WRITE (ILUME,9005)
      ENDIF
      CALL INLINE
      CALL GETFN(1)
      IF (LMISS) THEN
      LMISS=.FALSE.
      LERROR=.FALSE.
      LRET=.TRUE.
      RETURN
      ENDIF
      IF(LERROR) GOTO 10
      ENDIF
C     ----------------------------------
      ENTRY CRAFIL (ILU,ITYPE,LRET,AEXT)
C     ----------------------------------  
      LERROR=.FALSE.
      LMISS=.FALSE.
      LRET=.FALSE.
C
C     Check for proper extension
C
      CALL FILECH (AEXT)
C
      LREP=ITYPE.GT.0
      ITYP=IABS(ITYPE)
c
c     testing write enable assuming batch operation, no devices will be opened.
c
      call trywrite (ilu,afile)
      IF (ITYP.EQ.1) THEN
      OPEN(UNIT=ILU,FILE=AFILE,STATUS='NEW',IOSTAT=IERR,
     .     FORM='UNFORMATTED',ERR=5000)
      IF (LREP) WRITE (ILUME,9002) AFILE,ILU
      ENDIF
      IF (ITYP.EQ.2) THEN
      IF(INDEX(AFILE,'PRN').NE.0.OR.INDEX(AFILE,'LPT1').NE.0.OR.
     .   INDEX(AFILE,'CON').NE.0.OR.INDEX(AFILE,'NUL').NE.0) THEN
        OPEN(UNIT=ILU,FILE=AFILE)
      ELSE
        OPEN(UNIT=ILU,ERR=5000,FILE=AFILE,STATUS='NEW',IOSTAT=IERR)
      ENDIF      
      ENDIF
      IF (ITYP.EQ.3) THEN
      OPEN(UNIT=ILU,ERR=5000,ACCESS='DIRECT',FILE=AFILE,
     .FORM='FORMATTED',STATUS='NEW',RECL=46,IOSTAT=IERR)
      ENDIF
      IF (ITYP.EQ.4) THEN
        OPEN(UNIT=ILU,ERR=5000,FILE=AFILE,STATUS='NEW',RECL=132,
     &       IOSTAT=IERR)
      ENDIF
      IF (ITYP.EQ.5) THEN
        OPEN(UNIT=ILU,ERR=5000,FILE=AFILE,STATUS='NEW',RECL=256,
     &       IOSTAT=IERR)
      ENDIF
      GOTO 8000
 5000 IF (LPLOTON.AND.LSINGL) THEN
      WRITE (ILUME,9009)
      WRITE(ILUME,9110) AFILE
      ELSE
      WRITE (ILUME,9010) AFILE
      ENDIF
 5001 IF (LPLOTON.AND.LSINGL) THEN
C      CALL CURGET (RXC,RYC)         ! NOT AVAILABLE IN BATCH MODE
      ELSE
      CALL INLINE
      ENDIF
      IF (ALINE(1).EQ.'Q'.OR.ALINE(1).EQ.'q') THEN
      LRET=.TRUE.
      LERROR=.FALSE.
      LMISS=.FALSE.
      RETURN
      ENDIF
      IF (ALINE(1).EQ.'Y'.OR.ALINE(1).EQ.'y') THEN
      IF (ITYP.EQ.1) THEN
        OPEN(UNIT=ILU,FILE=AFILE,STATUS='UNKNOWN',
     .    FORM='UNFORMATTED',ERR=6000,IOSTAT=IERR)
      ENDIF
      IF (ITYP.EQ.2) THEN
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='UNKNOWN',ERR=6000,IOSTAT=IERR)
      ENDIF
      IF (ITYP.EQ.3) THEN
      OPEN(UNIT=ILU,ERR=6000,ACCESS='DIRECT',FILE=AFILE,
     .FORM='FORMATTED',STATUS='UNKNOWN',RECL=46,IOSTAT=IERR)
      CLOSE (UNIT=ILU,STATUS='DELETE',ERR=7000,IOSTAT=IERR)
      OPEN(UNIT=ILU,ERR=6000,ACCESS='DIRECT',FILE=AFILE,
     .FORM='FORMATTED',STATUS='NEW',RECL=46,IOSTAT=IERR)
      ENDIF
      IF (ITYP.EQ.4) THEN
      OPEN(UNIT=ILU,ERR=5000,FILE=AFILE,STATUS='UNKNOWN',RECL=132,
     &       IOSTAT=IERR)
      ENDIF 
      IF (ITYP.EQ.5) THEN
      OPEN(UNIT=ILU,ERR=5000,FILE=AFILE,STATUS='UNKNOWN',RECL=256,
     &       IOSTAT=IERR)
      ENDIF      
      IF (LREP) WRITE (ILUME,9003) AFILE,ILU      
      GOTO 8000
      ENDIF
      IF (ALINE(1).EQ.'N'.OR.ALINE(1).EQ.'n') THEN
      LERROR=.FALSE.
      LMISS=.TRUE.
      GOTO 5
      ENDIF
      IF (LPLOTON.AND.LSINGL) THEN
      WRITE (ILUME,9009)
      WRITE (ILUER,9122)
      ELSE
      WRITE (ILUER,9022)
      WRITE (ILUER,9010) AFILE
      ENDIF
      GOTO 5001
 6000 WRITE (ILUER,9011) MOD(IERR,256),ITYPE
      AMESS(1)='Failure to allocate or access file.'
      AMESS(2)='Make sure the directory or file are write enabled.'
      CALL HALT(2)  ! stop program execution for batch version
c ----------------------- end of program execution -----------------------
 7000 WRITE (ILUER,9012) MOD(IERR,256),ITYPE
      LERROR=.TRUE.
      RETURN
C
 8000 INQUIRE (UNIT=ILU,OPENED=LFILEOPEN,READWRITE=AWRITE)            ! ensure file is connected
      IF (AWRITE.EQ.'NO') THEN
      AMESS(1)='No write access to file.'
      AMESS(2)='Make sure the directory and file are write enabled.'
      CALL HALT(2)  ! stop program execution for batch version
c ----------------------- end of program execution -----------------------
      END IF
      IF (LFILEOPEN) THEN
        LRET=.FALSE.
        LERROR=.FALSE.     
      ELSE
      AMESS(1)='Failure to access file.'
      AMESS(2)='Make sure the file is write enabled.'
      CALL HALT(2)  ! stop program execution for batch version
c---------------------- end of program execution --------------------------
      ENDIF
      RETURN
C      
 9002 FORMAT (' File ',A16,' has been opened and assigned to LU:',I3,/)
 9102 FORMAT ('+File ',A16,' has been opened and assigned to LU:',I3) 
 9003 FORMAT(' File ',A16,' has been reopened and assigned to LU:',I3,/)
 9103 FORMAT('+File ',A16,' has been reopened and assigned to LU:',I3) 
 9001 FORMAT(' ***ILLEGAL OR MISSING PARAMETERS IN CRFILE'/)
 9101 FORMAT('+***ILLEGAL OR MISSING PARAMETERS IN CRFILE') 
 9005 FORMAT(' Enter name of file to be created. (CR to abort IO) >')
 9006 FORMAT('+Enter name of file to be created. (CR to abort IO) >') 
 9009 FORMAT ('+                                                     ',
     &'                                        ')
 9010 FORMAT(' File ',A16,' already exists; REPLACE? (Yes/No/Quit) >')     
 9110 FORMAT('+File ',A16,' already exists; REPLACE? (Yes/No/Quit) >')
 9011 FORMAT(' *** ERROR ',I5,' while allocating file (ITYPE=',I2,').'/)
 9111 FORMAT('+*** ERROR ',I5,' while allocating file (ITYPE=',I2,').') 
 9012 FORMAT(' *** ERROR ',I5,' while closing file (ITYPE=',I2,').'/)
 9112 FORMAT('+*** ERROR ',I5,' while closing file (ITYPE=',I2,').') 
 9022 FORMAT(' *** INCORRECT RESPONSE:')
 9122 FORMAT('+*** INCORRECT RESPONSE: REPLACE? (Yes/No/Quit) >') 
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE OPFILU(ILU,ITYPE,LRET,AEXT)
C
C --------------------------------------------------------------------------
C
c     open an existing file with prompting for a name
c
C
C     |ITYPE| = 1 unformatted file
C     |ITYPE| = 2 formatted file
C     |ITYPE| = 3 formatted, direct access file (recl=31)
C     ITYPE < 0 do not report IO
C      
      IMPLICIT NONE
      INTEGER(4)ILU,ITYPE,ITYP,IERR
      LOGICAL LRET,LREP,LFILEOPEN
      CHARACTER(4) AEXT
      INCLUDE 'DWMN.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      LRET=.TRUE.
      CALL GETFN (2)
  5   IF (LERROR) THEN    
      LERROR=.FALSE.
      LMISS=.FALSE.
 10   IF(LMISS.OR.LERROR) THEN
      IF (LPLOTON.AND.LSINGL) THEN
      WRITE (ILUME,9009)
      WRITE (ILUER,9101)
      ELSE
      WRITE(ILUER,9001)
      ENDIF
      ENDIF
      LERROR=.FALSE.
      LMISS=.FALSE.
      IF (LPLOTON.AND.LSINGL) THEN
      WRITE (ILUME,9009)
      WRITE (ILUME,9105)
C      CALL CURGET (RXC,RYC)         ! NOT AVAILABLE IN BATCH MODE
      ELSE
      WRITE(ILUME,9005)
      CALL INLINE
      ENDIF
      CALL GETFN(1)
      IF (LMISS) THEN
      LMISS=.FALSE.
      LERROR=.FALSE.
      LRET=.TRUE.
      RETURN
      ENDIF
      IF(LERROR) GOTO 10
      ENDIF
C
C     Check for proper file extension
C
      CALL FILECH (AEXT)
C      
C     -----------------------------
      ENTRY OPAFIL (ILU,ITYPE,LRET)
C     -----------------------------
c
c     open an existing file without prompting for a name
c
C
C     |ITYPE| = 1 unformatted file
C     |ITYPE| = 2 formatted file
C     |ITYPE| = 3 formatted, direct access file (recl=31)
C     ITYPE < 0 do not report IO
C      
      LRET=.FALSE.
      LERROR=.FALSE.
      LMISS=.FALSE.
      LREP=ITYPE.GT.0
      ITYP=IABS(ITYPE) 
      IF (ITYP.EQ.1) THEN
      OPEN(UNIT=ILU,FILE=AFILE,STATUS='OLD',FORM='UNFORMATTED',
     .     ERR=5000,IOSTAT=IERR)
      ENDIF
      IF (ITYP.EQ.2) THEN
      OPEN (UNIT=ILU,FILE=AFILE,STATUS='OLD',ERR=5000,IOSTAT=IERR)
      ENDIF
      IF (ITYP.EQ.3) THEN
      OPEN(UNIT=ILU,ERR=5000,ACCESS='DIRECT',FILE=AFILE,
     .STATUS='OLD',RECL=31,IOSTAT=IERR)
      ENDIF
      IF (LREP) WRITE (ILUME,9007) AFILE,ILU
      INQUIRE (UNIT=ILU,OPENED=LFILEOPEN)          ! ensure file is connected
      IF (LFILEOPEN) THEN
        LRET=.FALSE.
        LERROR=.FALSE.     
      ELSE
        LRET=.TRUE.
        LERROR=.TRUE.
      ENDIF
      RETURN
 5000 IF (LPLOTON.AND.LSINGL) THEN
      WRITE (ILUME,9009)
      WRITE (ILUER,9110) AFILE,MOD(IERR,256)
      ELSE
      WRITE(ILUER,9010) AFILE,MOD(IERR,256)
      ENDIF
      LERROR=.TRUE.
      GOTO 5
 9000 FORMAT(80A1)
 9001 FORMAT(' ***ILLEGAL OR MISSING PARAMETERS IN OPFILU'/)
 9101 FORMAT('+***ILLEGAL OR MISSING PARAMETERS IN OPFILU') 
 9005 FORMAT(' Enter name of file to be opened. (CR to abort IO) >')
 9105 FORMAT('+Enter name of file to be opened. (CR to abort IO) >') 
 9007 FORMAT(' File ',A16,' assigned to logical unit ',I2,/)
 9107 FORMAT('+File ',A16,' assigned to logical unit ',I2) 
 9010 FORMAT(' ***File ',A16,' does not exist. (IOSTAT=',I5,')')
 9110 FORMAT('+***File ',A16,' does not exist. (IOSTAT=',I5,')') 
 9009 FORMAT ('+                                                     ',
     &'                                        ') 
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE FILECH (AEXT)
C
C --------------------------------------------------------------------------
C
C
C     Routine ensures that the filename extension is AEXT
C
      IMPLICIT NONE
      INTEGER(4) I
      CHARACTER(4) AEXT
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'DWMN.INC'
C
      DO 10 I=2,13
      IF (AFILE(I:I).EQ.' ') THEN
      AFILE(I:I+3)=AEXT
      RETURN
      ENDIF
      IF (AFILE(I:I).EQ.'.') THEN
      IF (AFILE(I:I+3).EQ.AEXT) RETURN
      IF (LPLOTON.AND.LSINGL) THEN
      WRITE (ILUME,9009)
      WRITE (ILUER,1001) AFILE(I:I+3),AEXT
      ELSE
      WRITE (ILUER,1000) AFILE(I:I+3),AEXT
      ENDIF
      AFILE(I:I+3)=AEXT
      RETURN
      ENDIF
  10  CONTINUE
      AFILE(13:16)=AEXT
      IF (LPLOTON.AND.LSINGL) THEN
      WRITE (ILUME,9009)
      WRITE (ILUER,2001) AFILE
      ELSE
      WRITE (ILUER,2000) AFILE
      ENDIF
      RETURN
 1000 FORMAT (' ***WARNING: extension "',A4,'" repaced by "',A4,'"!')
 1001 FORMAT ('+***WARNING: extension "',A4,'" repaced by "',A4,'"!') 
 2000 FORMAT (' ***WARNING: filename too long, replaced by ',A16)
 2001 FORMAT ('+***WARNING: filename too long, replaced by ',A16) 
 9009 FORMAT ('+                                                     ',
     &'                                        ')
      END
C
C --------------------------------------------------------------------------
C
      subroutine trywrite (ilu,afil)
C
C --------------------------------------------------------------------------
C
c
c     routine checks if write operation to the file is possible
c
      implicit none
      INTEGER(4) ilu
      CHARACTER*16 afil
      include 'lusys.inc'
c
      open (UNIT=ilu,FILE=afil,ERR=100,STATUS='replace')
      write (ilu,FMT="('*')",ERR=200)
      close (ilu)
      return
 100  AMESS(1)='Failure to replace file:'
      AMESS(2)=afil
      AMESS(3)='Make sure directory and file are write enabled'
      CALL HALT(3)  ! stop program execution for batch version
      stop
 200  AMESS(1)='Failure to write to file:'
      AMESS(2)=afil
      AMESS(3)='Make sure directory and file are write enabled'
      CALL HALT(3)  ! stop program execution for batch version
c
      end subroutine
C
C --------------------------------------------------------------------------
C
      SUBROUTINE GETFNUL (NPAR)
C
C --------------------------------------------------------------------------
C
C
C  Parameter NPAR in input line is interpreted as a file-name
C  and placed, left-justified, in AFILE.
C  NCHAR gives the number of characters found in the name.
C
      IMPLICIT NONE
      INTEGER(4) NPAR,I,I1,I2
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'DWMN.INC'
      DO 10 I=1,ilinelength         ! restore ALINE prior to TIDY call
      ALINE(I)=ALINE2(I)
  10  CONTINUE
      CALL TIDYUL          ! new TIDYUL call maintaining upper and lower case
C     ------------------
      ENTRY GETFN (NPAR)
C     ------------------            
      AFILE = '                '
      NCHAR = 0
      IF (NPAR.LE.0) THEN
            LMISS=.TRUE.
            RETURN
      ENDIF            
      I1=ILPNT(NPAR)
      NCHAR=MAX0(0,ILPNT(NPAR+1)-I1)
      IF (NCHAR.EQ.0)  THEN
            LMISS=.TRUE.
            RETURN
      ENDIF
      IF (NCHAR.GT.16) THEN
            IF (LOPEN) THEN
C            IF (LOPEN) CALL PLSETCUR (7,1,1)   ! NOT AVAILABLE IN BATCH MODE
            ELSE
            LERROR=.TRUE.
            ENDIF
            WRITE (ILUER,2000) NCHAR
      ENDIF
      NCHAR=MIN0(NCHAR,16)
      I2=I1+NCHAR-1
      WRITE (AFILE,1000) (ALINE(I),I=I1,I2)
      RETURN
 1000 FORMAT (16A1)
 2000 FORMAT (' ***WARNING: string has ',I3,
     &        ' characters, truncated to maximum length of 16.')
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE GFSURF (RA,ISIZE)
C
C --------------------------------------------------------------------------
C
C
C     Write .BLN and .GRD file
C
      IMPLICIT NONE
      INTEGER(4) ISIZE,NROW,NXROW,NXREST,I,J,K
      LOGICAL LRET
      REAL(8) RA
      CHARACTER(16) AFILT,AFILESURFER1,AFILESURFER2
      DIMENSION RA(ISIZE,*)
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'TRACOM.INC'
C
      CALL GETFN (2)
      IF (LERROR.OR.LMISS) RETURN
      AFILT=AFILE
      CALL FILECH ('.GRD')
      CLOSE (3)
      CALL CRAFIL (3,2,LRET,'.GRD')
      AFILESURFER1=AFILE
      IF (LRET) RETURN
      WRITE (3,1000)
      WRITE (3,1010) NX,NY
      WRITE (3,1020) RX1,RX2
      WRITE (3,1020) RY1,RY2
      CALL MINMAX (RA,ISIZE)
      WRITE (3,1020) RMIN,RMAX
      NROW=NX/5
      NXROW=NROW*5
      NXREST=NX-NXROW        
      DO 20 J=1,NY
      IF (NXROW.GT.0) THEN
      DO 10 I=1,NXROW,5
      WRITE (3,1030) (RA(I+K,J),K=0,4)
  10  CONTINUE
      ENDIF
      IF (NXREST.EQ.1) WRITE (3,1032) RA(NXROW+1,J)
      IF (NXREST.EQ.2) WRITE (3,1033) (RA(NXROW+K,J),K=0,1)
      IF (NXREST.EQ.3) WRITE (3,1034) (RA(NXROW+K,J),K=0,2)
      IF (NXREST.EQ.4) WRITE (3,1035) (RA(NXROW+K,J),K=0,3)            
  20  CONTINUE
      CLOSE (3)
      if (lucon) WRITE (ILUME,2000) AFILESURFER1
      IF (.NOT.LGRAPHICS) RETURN
C      AFILE=AFILT                         ! NOT AVAILABLE IN BATCH MODE
C      CALL FILECH ('.BLN')
C      CALL CRAFIL (3,3,LRET,'.BLN')
C      AFILESURFER2=AFILE
C      IF (LRET) RETURN
C      LSURF=.TRUE.
C      ISFREC=1
C      ISFCNT=0
C      call ploton
C      CALL LAYOUT
C      call plotof
C      CZ=(0.0,0.0)
C      CALL DRAW (CZ,-2)
C      LSURF=.FALSE.
C      CLOSE (3)
C      if (lucon) WRITE (ILUME,3000) AFILESURFER2
      RETURN
 1000 FORMAT ('DSAA')
 1010 FORMAT (I5,2X,I5)
 1020 FORMAT (G14.7,2X,G14.7)
 1030 FORMAT (5(G14.7,1X))
 1032 FORMAT (G14.7)
 1033 FORMAT (2(G14.7,1X))
 1034 FORMAT (3(G14.7,1X))
 1035 FORMAT (4(G14.7,1X))  
 2000 FORMAT (' Grid has been written to SURFER grid file: ',A16)
 3000 Format (' Layout has been written to SURFER boundary line file: ',
     &A16,/)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE HELP
C
C --------------------------------------------------------------------------
C
C
C     Routine views help files
C
      IMPLICIT NONE
      INCLUDE 'MAIN.INC'
      CALL SYSTEM (AHELPDRIVE)
      CALL SYSTEM (AHELPPATH)
      CALL VIEW
      CALL SYSTEM (ACURRENTDRIVE)
      CALL SYSTEM (ACURRENTPATH)
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE VIEW
C
C --------------------------------------------------------------------------
C
C
C     Routine reads and scrolls or pages through ASCII files
C
      IMPLICIT NONE
      INTEGER(4) ILU,NLINETOTAL,ISTARTLINE,IENDLINE,ILINE,IKEY,IXKEY,
     &  NLINEMAX,IMXLN,I,J,IATR,IBUFST,IBUFEN,ILINE1,ILINE2,NLINEMX2
      LOGICAL LUPLINE,LDOWNLINE,LUPAGE,LDOWNPAGE,LRETURN
      PARAMETER (NLINEMAX=100,IMXLN=24)        ! 100 test value
      CHARACTER(1) AESC
      CHARACTER(80) ALINEBUFFER(NLINEMAX)
      CHARACTER(80) APAGEBUFFER
      DIMENSION APAGEBUFFER(IMXLN)
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
C
      DO 1 I=1,NLINEMAX                   ! empty buffers
      ALINEBUFFER(I)='                                                 
     &                              '
  1   CONTINUE
      DO 2 I=1,IMXLN
      APAGEBUFFER(I)=ALINEBUFFER(1)
  2   CONTINUE      
C  
      AESC=CHAR(27)
      ILU=2
      NLINEMX2=NLINEMAX/2
      OPEN (ILU,ERR=100,FILE=AFILE,STATUS='OLD',POSITION='REWIND')     
      NLINETOTAL=0
  5   READ (ILU,END=6,FMT=1000) ALINEBUFFER(1)
      NLINETOTAL=NLINETOTAL+1           ! determine length of file
      GOTO 5
  6   IF (NLINETOTAL.EQ.0) THEN         ! empty file, return
        WRITE (ILUER,3000)
        RETURN
      ENDIF
      IBUFST=1
      IBUFEN=MIN(NLINETOTAL,NLINEMAX)
      ISTARTLINE=1
      IENDLINE=MIN(IMXLN,NLINETOTAL)
      REWIND ILU
  7   CONTINUE
      I=0
      DO 10 J=IBUFST,IBUFEN               ! read lines into linebuffer
      I=I+1
      READ (ILU,END=200,FMT=1000) ALINEBUFFER(I)
  10  CONTINUE
  17  CONTINUE      
      IF (ISTARTLINE.LT.IBUFST) THEN      ! hit top of linebuffer, refill
        IBUFST=IBUFST-NLINEMX2
        IBUFST=MAX(IBUFST,1)
        IBUFEN=IBUFST+NLINEMAX-1
        REWIND ILU
        IF (IBUFST.EQ.1) GOTO 7
        DO 18 I=1,IBUFST-1
        READ (ILU,END=200,FMT=1000) ALINEBUFFER(1)
  18    CONTINUE
        GOTO 7
      ENDIF
      IF (IENDLINE.GT.IBUFEN) THEN        ! hit bottom of linebuffer, refill
        IBUFEN=IBUFEN+NLINEMX2
        IBUFEN=MIN(IBUFEN,NLINETOTAL)
        IBUFST=IBUFEN-NLINEMAX+1
        REWIND ILU
        IF (IBUFST.EQ.1) GOTO 7
        DO 20 I=1,IBUFST-1
        READ (ILU,END=200,FMT=1000) ALINEBUFFER(1)
  20    CONTINUE
        GOTO 7
      ENDIF          
      I=0
      ILINE1=ISTARTLINE-IBUFST+1
      ILINE2=IENDLINE-IBUFST+1
      DO 21 ILINE=ILINE1,ILINE2           ! read linebuffer into page buffer
      I=I+1
      APAGEBUFFER(I)=ALINEBUFFER(ILINE)
  21  CONTINUE   
C            IF (LOPEN) CALL PLSETCUR (7,1,1)   ! NOT AVAILABLE IN BATCH MODE
      WRITE (ILUME,2000) APAGEBUFFER      ! write page buffer to screen
      IATR=32
      WRITE (ILUME,1600) AESC,IATR        ! set text color green
      WRITE (ILUME,4000)                  ! write menu bar at screen bottom
      IATR=37
      WRITE (ILUME,1600) AESC,IATR        ! set text color back to white
  22  continue
      ikey=0 ! to avoid compiler warning
      LUPLINE=IKEY.EQ.1072
      LDOWNLINE=IKEY.EQ.1080
      LUPAGE=IKEY.EQ.1073
      LDOWNPAGE=IKEY.EQ.1081
      LRETURN=IKEY.EQ.27
      IF (LUPLINE) THEN
      IF (ISTARTLINE.EQ.1) GOTO 22        ! at top of file, get new key
      ISTARTLINE=ISTARTLINE-1             ! move up one line
      IENDLINE=IENDLINE-1
      GOTO 17
      ENDIF
      IF (LDOWNLINE) THEN
      IF (IENDLINE.EQ.NLINETOTAL) GOTO 22 ! at bottom of file, get new key
      ISTARTLINE=ISTARTLINE+1             ! move down one line
      IENDLINE=IENDLINE+1
      GOTO 17
      ENDIF
      IF (LUPAGE) THEN
      IF (ISTARTLINE.LE.IMXLN) THEN       ! move to top of file
        ISTARTLINE=1
        IENDLINE=MIN(IMXLN,NLINETOTAL)
      ELSE
        ISTARTLINE=ISTARTLINE-IMXLN         ! move up one page
        IENDLINE=IENDLINE-IMXLN
      ENDIF
      GOTO 17
      ENDIF
      IF (LDOWNPAGE) THEN
      IF (IENDLINE.GT.NLINETOTAL-IMXLN) THEN ! move to bottom of file
        IENDLINE=NLINETOTAL
        ISTARTLINE=MAX(NLINETOTAL-IMXLN+1,1)
      ELSE
        ISTARTLINE=ISTARTLINE+IMXLN            ! move down one page
        IENDLINE=IENDLINE+IMXLN
      ENDIF
      GOTO 17
      ENDIF
      IF (LRETURN) THEN       ! clear screen, close file, and return
      CLOSE (ILU)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      RETURN
      ENDIF
      GOTO 22
C      
  100 CONTINUE
c      call inline
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      WRITE (ILUER,5000) AFILE
      RETURN
  200 CONTINUE
c      call inline
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      WRITE (ILUER,6000)
      RETURN
C
 1000 FORMAT (A80)
 1600 FORMAT ('+',A1,'[',I2,'m')        
 1999 FORMAT (' ')
 2000 FORMAT (23('+',A80,/),'+',A80)
 3000 FORMAT (' ***ERROR in view routine, file is empty.',/)
 4000 FORMAT ('+ Press UP/DOWN ARROW keys or PAGE keys to scan text,',
     &        ' press <Esc> to exit')
 5000 FORMAT (' ***ERROR in view routine, cannot open file:',A16,/)     
 6000 FORMAT (' ***ERROR in view routine, try to read beyond end of',
     &        ' file.',/)
      END      
C
C --------------------------------------------------------------------------
C

      SUBROUTINE MATCH (NTABB,NPAR,JUMP,LBAD)   ! *** BRUTE FORCE FIX OF SUBSCRIPT OUT OF BOUNDS ERROR
C
C --------------------------------------------------------------------------
C
C
C TO MATCH INPUT STRING (IN ALINE) TO KEYWORD IN TABLE
C
C INPUT:  NTAB  TABLE OF KEYWORDS IN A1 FORMAT, SEPARATED
C -----         BY BLANKS AND TERMINATED WITH $.
C               ALL KEYWORDS MUST BE OF AT LEAST TWO
C               CHARACTERS, UNLESS A SINGLE-LETTER
C               KEYWORD IS REALLY INTENDED.
C         NPAR  PARAMETER NUMBER IN INPUT LINE.
C
C OUTPUT: LBAD   .TRUE. IF UNSUCCESSFUL
C ------  JUMP  SEQUENCE NUMBER OF KEYWORD IF SUCCESSFUL.
C               = 1 FOR MISSING PARAMETER (LBAD=.TRUE.)
C               = 2 IF WORD NOT FOUND (LBAD=.TRUE.)
C
      IMPLICIT NONE
      INTEGER(4) NPAR,JUMP,I,IL1,IL2,NT,IL
      LOGICAL LBAD,LONG
      CHARACTER(1) IBLK,NTI,NTABB(*),NTAB(1000)
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
      DATA IBLK /' '/
C
C      WRITE(ILUME,*)'NTAB = ',NTAB
C      WRITE(ILUME,*)'NTDIM = ',NTDIM
C      WRITE(ILUME,*)'NPAR = ',NPAR
      I=0
   1  I=I+1
      NTAB(I)=NTABB(I)
      IF (NTAB(I).NE.ATERM) GOTO 1    ! store key words in temporary array to avoid array overflow
      NTAB(I+1)=ATERM
      LONG=.FALSE.
      IL1=ILPNT(NPAR)
      IL2=ILPNT(NPAR+1)-1
      IF(IL2.GT.0) GOTO 20
C---MISSING PARAMETER---
      JUMP=1
      GOTO 150
C-----------------------
   20 I=0
      NT=0
   40 NT=NT+1
      IF(NTAB(I+2).EQ.IBLK.AND.IL2.GT.IL1) GOTO 80
      DO 50 IL=IL1,IL2
      I=I+1
      NTI=NTAB(I)
      IF(NTI.NE.ALINE(IL)) GOTO 70
   50 CONTINUE
      IF (.NOT.LONG) GOTO 60
      IF (NTAB(I+1).NE.IBLK) GOTO 80
   60 JUMP=NT
      LBAD=.FALSE.
      RETURN
   65 IF (NTI.EQ.IBLK) GOTO 40
      IL2=IL2-1
      LONG=.TRUE.
      IF (IL2.LT.IL1) GOTO 100
      GOTO 20
   70 IF(NTI.EQ.IBLK.OR.NTI.EQ.ATERM) GOTO 65
   80 I=I+1
      NTI=NTAB(I)
      IF(NTI.EQ.IBLK) GOTO 40
      IF(NTI.NE.ATERM) GOTO 80
      GOTO 65
C---NOT FOUND---
  100 JUMP=2
  150 LBAD=.TRUE.
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE VAR(NPAR)
C
C --------------------------------------------------------------------------
C
C
C Common Routine for IVAR and RVAR
C
      IMPLICIT NONE
      INTEGER(4) NPAR,K20,NP,I,ILL,NUM,IL,N1
      CHARACTER(1) IBUF(20),IBL
      CHARACTER(20) BUF
      COMMON /CVARR/ BUF
      INCLUDE 'MATCH.INC'
      DATA IBL,K20 /' ',20/
      LERROR=.FALSE.
      NP=NPAR
      NPBAD=NPAR-1
      DO 10 I=1,K20
   10 IBUF(I)=IBL
      ILL=ILPNT(NP)
      NUM=ILPNT(NP+1)-ILL
      NUM=MIN0(NUM,K20)
      IF(NUM.GT.0) GOTO 15
      LMISS=.TRUE.
      LERROR=.TRUE.
C     NERR=7
      WRITE (BUF,100) IBUF
      RETURN
   15 DO 20 IL=1,NUM
      N1=K20-NUM+IL
      IBUF(N1)=ALINE(ILL)
   20 ILL=ILL+1
      WRITE (BUF,100) IBUF
      LMISS=.FALSE.
  100 FORMAT(20A1)
      RETURN
      END
C
C --------------------------------------------------------------------------

      INTEGER(4) FUNCTION IVAR(NPAR)
C
C --------------------------------------------------------------------------
C
C
C Return Integer Value of Parameter NPAR
C LMISS is set .TRUE. if missing.
C LERROR is set .TRUE. if format bad and missing.
C
      IMPLICIT NONE
      INTEGER(4) NPAR,IV
      CHARACTER BUF*20
      COMMON /CVARR/ BUF
      INCLUDE 'MATCH.INC'
      IF(LERROR) GOTO 100
      CALL VAR(NPAR)
      READ (BUF,200,ERR=100) IV
      IVAR=IV
      RETURN
  100 LERROR=.TRUE.
      IVAR=0
      RETURN
  200 FORMAT(I20)
      END
C
C --------------------------------------------------------------------------
C
      REAL(8) FUNCTION RVAR(NPAR)
C
C --------------------------------------------------------------------------
C
C
C RETURN REAL VALUE OF PARAMETER NPAR
C
      IMPLICIT NONE
      INTEGER(4) NPAR
      REAL(8) RV
      CHARACTER BUF*20
      COMMON /CVARR/ BUF
      INCLUDE 'MATCH.INC'
      IF(LERROR) GOTO 100
      CALL VAR(NPAR)
      READ (BUF,200,ERR=100) RV
      RVAR=RV
      RETURN
  100 LERROR=.TRUE.
      RVAR=0.0
      RETURN
  200 FORMAT(F20.0)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE TIDY
C
C --------------------------------------------------------------------------
C
C
C ELIMINATE BLANKS, ETC. FROM INPUT LINE AND
C MAKE INDEX TO LOCATION OF PARAMETERS
C
      IMPLICIT NONE
      INTEGER(4) IL1,IL2,NPAR,I
      LOGICAL LSEP
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      SAVE
      LERROR=.FALSE.
      IL1=1
      IL2=1
      NPAR=0
      DO 5 I=1,40
    5 ILPNT(I)=0
      GOTO 20
C---NOW WITHIN A STRING---
   10 IF(LSEP(ALINE(IL2))) GOTO 30
   15 ALINE(IL1)=ALINE(IL2)
      IL1=IL1+1
      IL2=IL2+1
      IF(IL2.LE.ilinelength) GOTO 10
      GOTO 50
C---NOW IN A GAP---
   20 IF(.NOT.LSEP(ALINE(IL2))) GOTO 40
   30 IL2=IL2+1
      IF(IL2.LE.ilinelength) GOTO 20
      GOTO 50
C---START OF A STRING---
   40 NPAR=NPAR+1
      ILPNT(NPAR)=IL1
      GOTO 15
   50 ILPNT(NPAR+1)=IL1
      RETURN
 1000 FORMAT ('+',132A1)
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE TIDYUL
C
C --------------------------------------------------------------------------
C
C
C ELIMINATE BLANKS, ETC. FROM INPUT LINE AND
C MAKE INDEX TO LOCATION OF PARAMETERS
C
      IMPLICIT NONE
      INTEGER(4) IL1,IL2,NPAR,I
      LOGICAL LSEPUL
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      SAVE
      LERROR=.FALSE.
      IL1=1
      IL2=1
      NPAR=0
      DO 5 I=1,40
    5 ILPNT(I)=0
      GOTO 20
C---NOW WITHIN A STRING---
   10 IF(LSEPUL(ALINE(IL2))) GOTO 30
   15 ALINE(IL1)=ALINE(IL2)
      IL1=IL1+1
      IL2=IL2+1
      IF(IL2.LE.ilinelength) GOTO 10
      GOTO 50
C---NOW IN A GAP---
   20 IF(.NOT.LSEPUL(ALINE(IL2))) GOTO 40
   30 IL2=IL2+1
      IF(IL2.LE.ilinelength) GOTO 20
      GOTO 50
C---START OF A STRING---
   40 NPAR=NPAR+1
      ILPNT(NPAR)=IL1
      GOTO 15
   50 ILPNT(NPAR+1)=IL1
      RETURN
 1000 FORMAT ('+',132A1)
      END      
C
C --------------------------------------------------------------------------
C
      LOGICAL FUNCTION LSEP(ACHARR)
C
C --------------------------------------------------------------------------
C
C
C Returns .TRUE. if ACHARR is a separator.
C Converts all uppercases to lower cases
C
      IMPLICIT NONE
      INTEGER(4) I
      LOGICAL LSEPUL
      CHARACTER(1) ACHARR,CSEP(2),AFUPC
      SAVE
      DATA CSEP /' ',','/
C
      ACHARR=AFUPC(ACHARR)      ! switch to upper case
C     --------------------
      ENTRY LSEPUL (ACHARR)    ! maintain upper and lower case
C     --------------------
      DO 10 I=1,2
      IF(ACHARR.EQ.CSEP(I)) THEN
      LSEP=.TRUE.
      LSEPUL=.TRUE.
      RETURN
      ENDIF
   10 CONTINUE
      LSEP=.FALSE.
      LSEPUL=.FALSE.
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      CHARACTER(1) FUNCTION AFUPC(ACHARR)
C
C --------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) I
      CHARACTER(1) ACHARR,AFLWC
      CHARACTER(26) AUPPER,ALOWER
      SAVE
      DATA AUPPER /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA ALOWER /'abcdefghijklmnopqrstuvwxyz'/
      AFUPC=ACHARR
      DO 100 I=1,26
      IF (ACHARR.EQ.ALOWER(I:I)) AFUPC=AUPPER(I:I)
  100 CONTINUE
      RETURN
      ENTRY AFLWC (ACHARR)
      AFUPC=ACHARR
      DO 200 I=1,26
      IF (ACHARR.EQ.AUPPER(I:I)) AFLWC=ALOWER(I:I)
  200 CONTINUE
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE INLINE
C
C --------------------------------------------------------------------------
C
C
C     Routine reads ALINE and sets pointers 
C
      IMPLICIT NONE
      INTEGER(4) I,IKEY,ILENGTH,IERR
      CHARACTER(256) AMESSAGE
      CHARACTER(86) ASTRING
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MAIN.INC'
C
   5  IF (LUCON) THEN                    ! <<< reading from the keyboard or file >>>
        DO 6 I=1,ilinelength
        ALINE(I)=' '
   6    CONTINUE
        I=0
   7    ikey=0
        IF (IKEY.EQ.1059) THEN                    ! help
          ALINE(1)='?'
          ALINE(2)=' '
          GOTO 20
        ENDIF
        IF (IKEY.EQ.1060) THEN                    ! go or generate
          ALINE(1)='G'
          ALINE(2)='O'
        GOTO 20
        ENDIF
        IF (IKEY.EQ.27.OR.IKEY.EQ.1068) THEN     ! quit
          ALINE(1)='Q'
          ALINE(2)=' '
        GOTO 20
        ENDIF
        IF (IKEY.EQ.8) THEN                      ! backspace
          IF (I.EQ.0) GOTO 5                     ! line buffer was empty
          ALINE(I)=' '
          WRITE (ILUME,2000) CHAR(IKEY)          ! send backspace to screen
          WRITE (ILUME,2000) ALINE(I)            ! blank out last character
          WRITE (ILUME,2000) CHAR(IKEY)          ! backspace to overwrite blank
          I=I-1                                  ! reduce line buffer by 1
          GOTO 7                                 ! look for next key
         ENDIF
         IF (IKEY.EQ.13) GOTO 20
           I=I+1
           IF (I.GT.ilinelength) GOTO 20
           ALINE(I)=CHAR(IKEY)
           WRITE (ILUME,2000) ALINE(I)          ! echo key to the screen
           GOTO 7
      ELSE
        READ (ILUIN,1000,END=45) ALINE          ! <<<< reading from a file >>>>
      ENDIF
  20  IF (LECHO.AND..NOT.LUCON) WRITE (ILUECH,1500) ALINE
      IF (LECHO.AND.LUCON) WRITE (ILUECH,1000) ALINE
      DO 30 I=1,ilinelength
      ALINE2(I)=ALINE(I)
  30  CONTINUE      
  35  CALL TIDY
      IF (ALINE(1).EQ.'*') GOTO 5
      if (aline(1).eq.'$') then
        PRINT*,' Press any key to continue.'
        goto 5
      endif
      RETURN
C     --------- hit end of file without proper trailer statements      
  45  WRITE (ILUER,3000)
      PRINT 3000
      PRINT*,' Press any key to continue.'
C      CALL TONE              ! NOT AVAILABLE IN BATCH MODE
      CLOSE(ILUIN)
      DO 50 ILENGTH=ilinelength,1,-1
      IF (AHELPPATH(ILENGTH:ILENGTH).NE.' ') GOTO 52
  50  CONTINUE
  52  ILENGTH=ILENGTH-3
      ASTRING(1:ILENGTH)=AHELPPATH(4:ILENGTH+3)
      ASTRING(ILENGTH+1:ILENGTH+10)='\trail.dat'
      OPEN (ILUIN,FILE=ASTRING,STATUS='OLD',IOSTAT=IERR,ERR=55)
      ALINE(1)='*'
      RETURN
C  55  CALL IOSTAT_MSG (IERR,AMESSAGE)
   55 PRINT*, IERR
      PRINT*,' Press any key to continue.'
      CLOSE (ILUIN)
      OPEN (ILUIN,FILE='CON')
      LUCON=.TRUE.
      ALINE(1)='*'
      RETURN
C      
 1000 FORMAT (132a1)
 1500 format ('+',132a1)
 2000 FORMAT ('&',A1)
 3000 FORMAT (' *** ERROR: Premature end of file encountered,',/,
     &' Data may be incomplete!! Attempt to use file "trail.dat".')
 4000 FORMAT (A152)
      END
C
C --------------------------------------------------------------------------
C
      COMPLEX(8) FUNCTION CVAR(IPAR)
C
C --------------------------------------------------------------------------
C
c     return a complex variable from ALINE
c
      IMPLICIT NONE
      INTEGER(4) IPAR
      REAL(8) RX1,RX2,RVAR
      INCLUDE 'MATCH.INC'
      CVAR=(.0,.0)
      RX1=RVAR(IPAR)
      RX2=RVAR(IPAR+1)
      IF(LERROR) RETURN
      CVAR=CMPLX(RX1,RX2)
      RETURN
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE SWITCH
C
C --------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) ILUI,ILUM,ILUE,ILUO,ILUP,ILUC
      INCLUDE 'LUSYS.INC'
      SAVE      
      ILUI=ILUIN
      ILUM=ILUME
      ILUO=ILUOUT
      ILUE=ILUER
      ILUP=ILUPL
      ILUC=ILUECH
      CALL SWTCHS(ILUI,ILUM,ILUE,ILUO,ILUP,ILUC)
      RETURN
      END
C
C --------------------------------------------------------------------------
C      
      SUBROUTINE SWTCHS(ILUI,ILUM,ILUE,ILUO,ILUP,ILUC)
C
C --------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) ILUI,ILUM,ILUE,ILUO,ILUP,ILUC,
     &           NWORD,INREC,JUMP,IDOM,IVAR,
     &           ILUNEW,ILUOLD,IERR,IDUM
      LOGICAL LERREPORT,LBAD,LLU,LTEMP
      CHARACTER(1) AWORD,APAR
      CHARACTER(16) ADEF(6),AFERROR
      PARAMETER (NWORD=47)
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MAIN.INC'
      INCLUDE 'TRACOM.INC'
      DIMENSION AWORD(NWORD),APAR(8)
      SAVE
      DATA ADEF /4*'CON             ','COM1            ',
     .             'NUL             '/
      DATA AWORD /'I','N','P','U',' ',
     .            'O','U','T','P',' ',
     .            'M','E','S','S',' ',
     .            'E','R','R','O',' ',
     .            'P','L','O','T',' ',
     .            'E','C','H','O',' ',
     .            '?',' ',
     .            'Q','U','I','T',' ',
     .            'Y','E','S',' ',
     .            'P','I','C','T',' ',
     .             ATERM/
      DATA APAR  /'O','N',' ',
     .            'O','F','F',' ',
     .             ATERM/
      LERREPORT=INDEX(ADEF(4),'CON').EQ.0
      IF (LERREPORT) THEN
        INQUIRE (ILUE,NEXTREC=INREC)
        AFERROR=ADEF(4)
      ENDIF
      LERROR=.FALSE.
      LMISS=.FALSE.
      IF (ILPNT(3).NE.0) GOTO 100       ! switch filename = input filename
  10  IF(LERROR) WRITE (ILUE,9001) ALINE2
      IF (LMISS) WRITE (ILUE,9002) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      if (lucon) then
      WRITE(ILUM,9500)
      WRITE (ILUM,9510) ADEF(1)
      WRITE (ILUM,9520) ADEF(2)
      WRITE (ILUM,9530) ADEF(3)
      WRITE (ILUM,9540) ADEF(4)
      IF (LECHO) WRITE (ILUM,9555) ADEF(6)
      IF (.NOT.LECHO) WRITE (ILUM,9556)
      IF (LGRAPHICS) WRITE (ILUM,9557)
      IF (.NOT.LGRAPHICS) WRITE (ILUM,9558)
      WRITE (ILUM,9600)
      endif
  11  CALL INLINE
      CALL MATCH(AWORD,1,JUMP,LBAD)
      IF (.NOT.LBAD) GOTO 15
      GOTO (10,13),JUMP
  13  WRITE (ILUE,9003) ALINE2
      LERROR=.FALSE.
      GOTO 10
C
C     Check for new logical unit number
C      
  15  LLU=ILPNT(4).NE.0
      IF(LLU) IDUM=IVAR(3)
      IF (IDUM.EQ.1.OR.IDUM.EQ.2) THEN
        IDUM=0
        WRITE (ILUER,9013)
        GOTO 10
      ENDIF
      IF(LMISS.OR.LERROR) GOTO 10
      ILUNEW=IDUM
      GOTO (100,200,300,400,500,550,600,700,800,900), JUMP
C
C     Change logical unit numbers if LLU=.TRUE.
C      
 100  IF (LUCON.AND.LINALREADY) THEN
      if (lucon) WRITE (ILUME,9015)
      ENDIF
      LINALREADY=.TRUE.
      CALL LUSWAP(LLU,ILUOLD,ILUNEW,ILUI)
      GOTO 5000
 200  CALL LUSWAP(LLU,ILUOLD,ILUNEW,ILUO)
      GOTO 6000
 300  CALL LUSWAP(LLU,ILUOLD,ILUNEW,ILUM)
      GOTO 6000
 400  CALL LUSWAP(LLU,ILUOLD,ILUNEW,ILUE)
      GOTO 6000
 500  CALL LUSWAP(LLU,ILUOLD,ILUNEW,ILUP)      
      GOTO 6000
 550  CALL LUSWAP(LLU,ILUOLD,ILUNEW,ILUC)
      GOTO 6000
C
C     Help
C      
 600  AFILE='SWTCHHLP.HLP   '
      GOTO 10
C
C     Return
C
 700  ILUIN=ILUI
      ILUOUT=ILUO
      ILUME=ILUM
      ILUER=ILUE
      ILUPL=ILUP
      ILUECH=ILUC
      IF (LERREPORT.AND.INREC.GT.1.and.lucon) THEN
        WRITE (ILUME,9990)
        CALL INLINE
        AFILE=AFERROR
      ENDIF
      return
C
C     dummy YES to avoid errors in reading instructions from a file
C     whereby a "yes" is added in anticipation of the need to delete
C     a file
C     
 800  GOTO 10
C
C
 900  CALL MATCH (APAR,2,JUMP,LBAD)    
      LMISS=LBAD
      IF (LMISS) GOTO 10
      LGRAPHICS=JUMP.EQ.1
      GOTO 10
C
C
C     Reassign input logical unit
C      
 5000 CALL GETFN (2)
      LUCON=AFILE(1:4).EQ.'CON '.OR.AFILE(1:4).EQ.'con ' 
      CALL FILECH ('.DAT')
      IF (LMISS) GOTO 10
C      LUCON=INDEX(AFILE,'CON').NE.0
      CLOSE(ILUOLD,ERR=6400,STATUS='KEEP')
      OPEN(ILUNEW,ERR=6500,FILE=AFILE,STATUS='OLD')
      GOTO 10
C
C     Reassign output logical units
C      
 6000 CALL GETFN (2)
      IF (LMISS) GOTO 10
C ----------------------------check for echo command to set LECHO      
      IF (JUMP.EQ.6) THEN
      IF (INDEX(AFILE,'OFF').NE.0) THEN
      CALL FLUSH (ILUOLD)
      CLOSE(ILUOLD,ERR=6400,STATUS='KEEP')      
      LECHO=.FALSE.
      GOTO 10
      ELSE
      LECHO=.TRUE.
      ENDIF
      ENDIF
C ----------------------------check for OUTPUT command to set LUOUTFILE
      IF (JUMP.EQ.2) THEN
      IF (AFILE(1:4).EQ.'CON '.OR.AFILE(1:4).EQ.'con ') THEN
      LUOUTFILE=.FALSE.
      ELSE
      LUOUTFILE=.TRUE.
      ENDIF
      ENDIF
C ---------------------------------------------------------------      
      CALL FLUSH (ILUOLD)
      CLOSE(ILUOLD,ERR=6400,STATUS='KEEP')
      IF (AFILE(1:4).EQ.'PRN '.OR.AFILE(1:4).EQ.'prn '.OR.
     &    AFILE(1:5).EQ.'LPT1 '.OR.AFILE(1:5).EQ.'lpt1 '.OR.
     &    AFILE(1:4).EQ.'CON '.OR.AFILE(1:4).EQ.'con '.OR.
     &    AFILE(1:4).EQ.'NUL '.OR.AFILE(1:4).EQ.'nul ') THEN
        OPEN(ILUNEW,FILE=AFILE)
      ELSE
        call trywrite(ilunew,afile)
        OPEN(ILUNEW,ERR=6600,FILE=AFILE,STATUS='NEW')  ! opening output file or device
      ENDIF
      ADEF(JUMP)=AFILE
      GOTO 10
C
C     Errors in closing or opening logical units.
C      
 6300 WRITE (ILUE,9009) AFILE,MOD(IERR,256)
      LECHO=.FALSE.
      AMESS(1)='Failure to create or access file'
      AMESS(2)='Make sure directory or file is write enabled'
      CALL HALT(2)  ! stop program execution for batch version
c ------------------------- end of program execution ------------------------
      GOTO 10
 6400 WRITE(ILUE,9010)
      AMESS(1)='Failure to access file'
      AMESS(2)='Make sure file is write enabled'
      CALL HALT(2)  ! stop program execution for batch version
c ------------------------- end of program execution ------------------------
      GOTO 10
 6500 if (lucon) WRITE(ILUM,9011)
      LUCON=.TRUE.
      GOTO 10      
 6600 if (lucon) WRITE(ILUM,9012)
      LTEMP=LECHO
      LECHO=.FALSE.
      CALL INLINE
      LECHO=LTEMP
      IF(ALINE(1).EQ.'Y'.OR.ALINE(1).EQ.'y') THEN
        CLOSE(ILUNEW,ERR=6400,STATUS='DELETE')
        OPEN(ILUNEW,FILE=AFILE,ERR=6300,IOSTAT=IERR)
      ELSE
        GOTO 10
      ENDIF
      ADEF(JUMP)=AFILE
      GOTO 10
C      
 9001 FORMAT(' ***ILLEGAL PARAMETER(S) in switch module:',/,80A1)
 9002 FORMAT(' ***MISSING PARAMETER(S) in switch module:',/,80A1)
 9003 FORMAT(' ***ILLEGAL COMMAND in switch module:',/,80A1)
 9009 FORMAT(' ***ERROR in opening file ',A16,' error code=',I3)
 9010 FORMAT(
     .' ***ERROR IN CLOSING FILE; correct problem and type return',/)
 9011 FORMAT(
     .' ***FILE DOES NOT EXIST; try again.')
 9012 FORMAT(
     .' ***FILE ALREADY EXISTS; do you want to delete it? Y/N',/,
     .' >')
 9013 FORMAT (' ***ERROR: logical unit 1 or 2 is reserved!',/)
 9015 FORMAT (' WARNING: input data will be ADDED to current data!',/,
     &' Press <Esc> to cancel, or any other key to continue.')
 9500 FORMAT (' ---------- SWITCH module ----------',/,
     &' Current assignments in SWITCH:'/
     &' IO function (filename, con, prn, lpt1, nul) ')
 9510 FORMAT (' INPUT       (',A16,'             ) ') 
 9520 FORMAT (' OUTPUT      (',A16,'             ) ') 
 9530 FORMAT (' MESSAGE     (',A16,'             ) ') 
 9540 FORMAT (' ERROR       (',A16,'             ) ') 
 9550 FORMAT (' PLOT        (',A16,'             ) ')
 9555 FORMAT (' ECHO        (',A16,'             ) ')
 9556 FORMAT (' ECHO OFF  (file/device name) [logical unit]')
 9557 FORMAT (' PICTURE ON    (off)')
 9558 FORMAT (' PICTURE OFF   (on)')
 9600 FORMAT (/' <F1>=Help',/,' <Esc> or QUIT '/' >') 
 9990 FORMAT (//' ********* ERRORS DETECTED DURING INPUT!',//,
     &        ' Press <Enter> for error report.')
      END
C
C --------------------------------------------------------------------------
C
      SUBROUTINE LUSWAP(LLU,ILUOLD,ILUNEW,ILU)
C
C --------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) ILUOLD,ILUNEW,ILU
      LOGICAL LLU
      SAVE
      ILUOLD=ILU
      IF(LLU) THEN
        ILU=ILUNEW
      ELSE
        ILUNEW=ILUOLD
      ENDIF
      RETURN
      END
c
c ---------------------------------------------------------------------------
c
      subroutine WriteMatrix (dra,isize,n,lErrorReport,ltimer)
c
c ---------------------------------------------------------------------------
c
c     routine writes the coefficient matrix to disk
c
c     dra          coefficient matrix of linear set of equations for unknown strength parameters
c     isize        number of rows
c     n            number of columns
c     lErrorReport .true. for reporting activity
c     ltimer       .true. fro reporting execution time
c
c    uses logical unit 10
c
      implicit none
      LOGICAL ltimer,lErrorReport,lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut
      INTEGER(4) n,istatus,i,j,isize,iticks,iticks1,iticks2,nsolOut
      REAL(8) dra
      CHARACTER(8) aBasenameOut
      CHARACTER(16) aDateTimeOut
      CHARACTER(256) amessage
      DIMENSION dra(isize,n)
      include 'lusys.inc'
c
      if (lErrorReport) write (ilume,1001) n
 1001 format (' writing coefficient matrix to disk, size=',i4)
      !call timer (iticks1)
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      REWIND (UNIT=10) ! make sure we are at start of file
      write (UNIT=10,IOSTAT=istatus,ERR=200) aDateTimeOut,n
      do i=1,isize
      write (UNIT=10,IOSTAT=istatus,ERR=200)dra(i,1:n) ! changed to row by row write on 4-19-06
      end do
      !call timer  (iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,2001) iticks
 2001 format(' Written coefficient matrix. Execution time=        ',i10,
     &' E-2 seconds.')

      return
c  200 call iostat_msg(istatus,amessage)
 200  write (*,2000) istatus ! error in accessing file
      write (ilume,2000) istatus
 2000 format (' Error in WRITE routine in WriteMatrix: IOSTAT=',/,a80,/,
     &        ' program execution aborted.')
      AMESS(1)='Error when writing to matrix file *.mtr'
      AMESS(2)='See IO error in Message.log file. Execution halted.'
      CALL HALT(2) ! stop program execution for batch version
c
      end subroutine
c
c ---------------------------------------------------------------------------
c
      subroutine LoadMatrix (dra,isize,n,lErrorReport,ltimer)
c
c ---------------------------------------------------------------------------
c
c     routine reads the coefficient matrix from disk
c
c     dra          coefficient matrix of linear set of equations for unknown strength parameters
c     isize        number of rows
c     n            number of columns
c     lErrorReport .true. for reporting activity
c     ltimer       .true. fro reporting execution time
c
c     using logical unit 10
c
      implicit none
      LOGICAL ltimer,lErrorReport
      INTEGER(4) n,istatus,i,j,isize,iticks,iticks1,iticks2,nEquation
      REAL(8) dra
      CHARACTER(256) amessage
      CHARACTER(16) aDateTimeStamp
      DIMENSION dra(isize,n)
      include 'lusys.inc'
c
      if (lErrorReport) write (ilume,1001) n
 1001 format (' loading coefficient matrix from disk, size=',i4)
      !call timer (iticks1)
      REWIND (UNIT=10)
      read (UNIT=10,IOSTAT=istatus,ERR=100) aDateTimeStamp,nEquation  ! not used here for speed
      do i=1,isize
      read (UNIT=10,IOSTAT=istatus,ERR=100)dra(i,1:n) ! changed to row by row read on 4-19-06
      end do
      !call timer  (iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,2001) iticks
 2001 format(' Loaded coefficient matrix. Execution time=         ',i10,
     &' E-2 seconds.')

      return
c
c  100 call IOSTAT_MSG (istatus,amessage)
 100  write (*,1000) istatus ! error in reading file
      write (ilume,1000) istatus
 1000 format (' Error in LoadMatrix:',/,a80,/,
     &        ' program execution aborted.')
      AMESS(1)='Error when reading matrix file *.mtr'
      AMESS(2)='See IO error in Message.log file. Execution halted.'
      CALL HALT(2) ! stop program execution for batch version
c
      end subroutine
c
c ---------------------------------------------------------------------------
c
      subroutine WriteDecompMatrix(dra,ipiv,n,lErrorReport,ltimer)
c
c ---------------------------------------------------------------------------
c
c     routine writes the decomposed matrix and pivot vector to disk
c
c     dra          coefficient matrix of linear set of equations for unknown strength parameters
c     ipiv         pivot vector
c     n            number of columns
c     lErrorReport .true. for reporting activity
c     ltimer       .true. fro reporting execution time
c
c     Using logical unit 11
c
c
      implicit none
      LOGICAL ltimer,lErrorReport,lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut
      INTEGER(4) ipiv,n,istatus,i,j,iticks,iticks1,iticks2,nsolOut
      REAL(8) dra
      CHARACTER(8) aBasenameOut
      CHARACTER(16) aDateTimeOut
      CHARACTER(256) amessage
      DIMENSION dra(n,n),ipiv(n)
      include 'lusys.inc'
c
      if (lErrorReport) write (ilume,1001) n
 1001 format (' writing decomposed matrix to disk, size=',i4)
      !call timer (iticks1)
      call GetMainData (lsolOut,loadsolOut,linalreadyOut,
     &           lErrorReportOut,lDirectfromDiskOut,
     &           aBasenameOut,aDateTimeOut,nsolOut)
      REWIND (UNIT=11)
      write (UNIT=11,IOSTAT=istatus,ERR=200) aDateTimeOut,n
      write (UNIT=11,IOSTAT=istatus,ERR=200) ((dra(i,j),i=1,n),j=1,n)
      write (UNIT=11,IOSTAT=istatus,ERR=200) (ipiv(i),i=1,n)
      !call timer  (iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,2001) iticks
 2001 format(' Written decomposed matrix. Execution time=         ',i10,
     &' E-2 seconds.')
      return
c  200 call iostat_msg(istatus,amessage)
 200  write (*,2000) amessage ! error in accessing file
      write (ilume,2000) amessage
 2000 format (' Error in WriteDecompMatrix: IOSTAT=',/,a80,/,
     &        ' program execution aborted.')
      AMESS(1)='Error writing to decomposed matrix file *.dec'
      AMESS(2)='See IO error in Message.log file. Execution halted.'
      CALL HALT(2) ! stop program execution for batch version
c
      end subroutine
c
c ---------------------------------------------------------------------------
c
      subroutine LoadDecompMatrix (dra,ipiv,n,lErrorReport,ltimer)
c
c ---------------------------------------------------------------------------
c
c     routine reads the decomposed matrix and pivot vector from disk
c
c     dra          coefficient matrix of linear set of equations for unknown strength parameters
c     ipiv         pivot vector
c     isize        number of rows
c     n            number of columns
c     lErrorReport .true. for reporting activity
c     ltimer       .true. fro reporting execution time
c
c     using logical unit 11
c
c
      implicit none
      LOGICAL ltimer,lErrorReport
      INTEGER(4) ipiv,n,istatus,i,j,iticks,iticks1,iticks2,
     &           nEquation
      REAL(8) dra
      CHARACTER(256) amessage
      CHARACTER(16) aDateTimeStamp
      DIMENSION dra(n,n),ipiv(n)
      include 'lusys.inc'
c
      if (lErrorReport) write (ilume,1001) n
 1001 format (' loading decomposed matrix from disk, size=',i4)
      !call timer (iticks1)
      REWIND (UNIT=11)
      read (UNIT=11,IOSTAT=istatus,ERR=100) aDateTimeStamp,nEquation  ! not used here for speed
      read (UNIT=11,IOSTAT=istatus,ERR=100) ((dra(i,j),i=1,n),j=1,n)
      read (UNIT=11,IOSTAT=istatus,ERR=100) (ipiv(i),i=1,n)
      !call timer  (iticks2)
      iticks=iticks2-iticks1
      if (ltimer) write (ilume,2001) iticks
 2001 format(' Loaded decomposed matrix. Execution time=          ',i10,
     &' E-2 seconds.')

      return
c
c  100 call IOSTAT_MSG (istatus,amessage)
  100 write (*,1000) istatus ! error in reading file
      write (ilume,1000) istatus
 1000 format (' Error in LoadDecompMatrix:',/,a80,/,
     &        ' program execution aborted.')
      AMESS(1)='Error reading decomposed matrix file *.dec'
      AMESS(2)='See IO error in Message.log file. Execution halted.'
      CALL HALT(2) ! stop program execution for batch version
c
      end subroutine
