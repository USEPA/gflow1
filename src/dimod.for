C     Last change:  HMH   6 Jan 2001    4:34 pm
c     This file contains the following routines and functions
c
c     BLOCK DATA DIDAT      initializes common blocks
C     SUBROUTINE DISCIN     input of (3D) sink discs
c
c -------------------
c DIMAT.FOR contains:
c -------------------
c
c     SUBROUTINE DICZC     generates control points
c     SUBROUTINE DIMAT     generates matrix coefficients
c     SUBROUTINE DIKNO     generates known vector
c     SUBROUTINE DISUB     substitutes solution vector
c
c ---------------------
c DICHECK.FOR contains:
c ---------------------
c
c     SUBROUTINE DIERROR      returns max. error at control points
c
c -------------------
c DIFUN.FOR contains:
c -------------------
c
c     REAL FUNCTION RFDIPT     contribution to disch. pot. of all (3D) discs
c     REAL FUNCTION RFDIGP     same as RFDIPT, but for filling grids
c     SUBROUTINE DISCQ2D       generates discharge vector due to all (3D) discs (in 2D zone)
c          ENTRY DISCQ3D       generates discharge vector in 3D zone
c     REAL FUNCTION RFDIPC     coefficient function for disc + images
c     COMPLEX FUNCTION CDIOM   complex potential outside 3D zone
c     COMPLEX FUNCTION CDICOMC coefficient functions for complex potential
c     REAL FUNCTION RFDISC     coefficient function for disch. pot. for single disc
c     REAL FUNCTION RFDIOM     solid angle due to a disc
c     SUBROUTINE DIQIC         coefficient function for discharge vector due to a single disc
c     SUBROUTINE DIQTZ         returns radial and vertical discharge components
c
c
c ------------------
c DIIO.FOR contains:
c ------------------
c
c     SUBROUTINE DIIO         reads or writes common block data to .SOL file
c     SUBROUTINE DIDATIO      writes 3D disc data to .DAT file
c     SUBROUTINE DIOUT        writes common block data (formatted) to disk for debugging purposes
c     subroutine diextracthead  writes head specified 3D disc data to .XTR file
c         entry diextractdisch  writes discharge specified 3D disc data to .XTR file
c
c ---------------------
c DINFLOW.FOR contains:
c ---------------------
c
c     REAL FUNCTION RFNFDI    returns flow across line due to all 3D discs 
c     REAL FUNCTION RFNFDICO  coefficient functions for RFNFDI
c     REAL  FUNCTION RFDPDIO  coefficient function when line is outside the disc
c     REAL  FUNCTION RFDPDII  coefficient function when line is inside the disc 
c
c --------------------
c DINEAR.FOR contains:
c --------------------
c
c    SUBROUTINE DINEAR    aborts streamline tracing when at 3D sink disc
c
c --------------------
c DISERV.FOR contains:
c --------------------
c
c  **** no file DISERV.FOR ****
c
c ---------------------
c DIEXTRA.FOR contains:
c ---------------------
c
c     SUBROUTINE DIPLOT     adds 3D sink discs to the layout
c     SUBROUTINE DICHEK     interactive routine for checking boundary conditions
c     SUBROUTINE DICUR      graphical data querries and data entry (not developed)
c
c
c
c
c ------------------------------------------------------------------
c
      BLOCK DATA DIDAT
c
c ------------------------------------------------------------------
c
      IMPLICIT NONE
      INCLUDE 'dicom.inc'
      DATA RDIZ,RDIR,RDIS,RDIH,RDICPT,RDIA /N_DI_INIT*0.0/
C
C     NDIRH is the number of aquifer thicknesses away from the
C     rim of the disc where a two-dimensional solution is used.
C     NDIMAG is the number of times that the aquifer and its
C     lower or upper image are repeated (imaged) upward or downward.
C
      DATA NDIS,NDIRH,NDIMAG /0,5,25/
      DATA RTAU,RTAU2,RDISR,RDISR2,RTMR,RTPR,RTDN,ROOT /8*0.0/
      DATA RTARA,RK,RK2,REPS /3*0.0,0.001/
      DATA RPI,RPI2,RPI4 /3.141592654,6.283185308,12.56637062/
      END
c
c ------------------------------------------------------------------
c
      SUBROUTINE DISCIN (LSOL,RA,IRA,RSCR)
c
c ------------------------------------------------------------------
c
C
C     input routine for three dimensional sinkdiscs
C      
      IMPLICIT NONE
      INTEGER(4) IRA,JUMP
      LOGICAL LSOL,LBAD,LHEADSPEC
      REAL(8) RA,RSCR,RTOP,RFTOP,RBASE,RFBASE,RFPOR,RFPERM,
     &        RDUM1,RVAR,RDUM2,RDUM3,RDUM4,RDUM5,RDUM6
      COMPLEX(8) CZ
      CHARACTER(1) AWORD(23)
      INCLUDE 'dicom.inc'
      INCLUDE 'com3d.inc'
      INCLUDE 'match.inc'
      INCLUDE 'lusys.inc'
      DIMENSION RA(IRA,*),RSCR(*)
      DATA AWORD /'?',' ',
     .            'D','I','S','C',' ',
     .            'H','E','A','D',' ',
     .            'Q','U','I','T',' ',
     .            'C','U','R','S',' ',
     .            ATERM/
      LERROR=.FALSE.
      LMISS=.FALSE.
      CZ=(0.0,0.0)
      RTOP=RFTOP(CZ)
      RBASE=RFBASE(CZ)
      R3DH=RTOP-RBASE
      R3DPOR=RFPOR(CZ)
      R3DK=RFPERM(CZ)      
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
  10  IF (LERROR.AND.LMISS) WRITE (ILUER,2000) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.  
      if (lucon) WRITE (ILUME,1000) NDIMAX
      CALL INLINE
  11  CALL MATCH (AWORD,1,JUMP,LBAD)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LBAD) GOTO 15
      GOTO (10,13) JUMP
  13  WRITE (ILUER,3000) ALINE2
      LERROR=.FALSE.
      GOTO 10
  15  GOTO (100,200,300,400,600),JUMP
C
C     help
C  
 100  AFILE='DISCHLP.HLP    '
      CALL HELP
      GOTO 10
C
C     discharge specified disc
C      
 200  if (lucon) WRITE (ILUME,5000)
      LHEADSPEC=.FALSE.
      GOTO 500
C
C     head specified disc
C      
 300  if (lucon) WRITE (ILUME,6000)
      LHEADSPEC=.TRUE.
      GOTO 500
C
C     return
C      
 400  RETURN
C
C     enter disc data
C 
 500  IF (NDIS.EQ.NDIMAX) THEN
      WRITE (ILUER,4000) NDIMAX
      GOTO 10
      ENDIF
      if (lucon) WRITE (ILUME,6500)
      CALL INLINE
      RDUM1=RVAR(1)
      IF (LERROR) THEN
      LERROR=.FALSE.
      LMISS=.FALSE.
      GOTO 11
      ENDIF
      RDUM2=RVAR(2)
      RDUM3=RVAR(3)
      RDUM4=RVAR(4)
      IF (RDUM4.LT.1.0E-10) THEN
        WRITE (ILUER,3500) RDUM4
        IF (LHEADSPEC) THEN
          GOTO 300
        ELSE
          GOTO 200
        ENDIF
      ENDIF
      RDUM5=RVAR(5)
      IF (RDUM3.LT.RBASE.OR.RDUM3.GT.RTOP) THEN
        WRITE (ILUER,7000) RDUM3,RBASE,RTOP
        IF (LHEADSPEC) THEN
          GOTO 300
        ELSE
          GOTO 200
        ENDIF
      ENDIF
      IF (LHEADSPEC) THEN
        RDUM6=RVAR(6)
        IF (LERROR) GOTO 10
        CALL GETFNUL (7)
      ELSE
        CALL GETFNUL (6)
      ENDIF
      NDIS=NDIS+1
      IF (LERROR) THEN
        LERROR=.FALSE.
        ADILAB(NDIS)="                "
      ELSE
        ADILAB(NDIS)=AFILE
      ENDIF
      RDIZ(1,NDIS)=RDUM1
      RDIZ(2,NDIS)=RDUM2
      RDIZ(3,NDIS)=RDUM3
      RDIR(NDIS)=RDUM4
      CZ=CMPLX(RDUM1+RDUM4,RDUM2)   ! update maximum window parameters
      CALL PLWIND (CZ)
      CZ=CMPLX(RDUM1-RDUM4,RDUM2)
      CALL PLWIND (CZ)
      CZ=CMPLX(RDUM1,RDUM2+RDUM4)
      CALL PLWIND (CZ)
      CZ=CMPLX(RDUM1,RDUM2-RDUM4)
      CALL PLWIND (CZ)
      IF (.NOT.LHEADSPEC) RDIS(NDIS)=RDUM5
      IF (LHEADSPEC) RDIH(NDIS)=RDUM5
      IF (LHEADSPEC) RDICPT(NDIS)=RDUM6
      LDIH(NDIS)=LHEADSPEC
      RDIA(NDIS)=RPI*RDUM4*RDUM4
      LSOL=.FALSE.
      IF (LHEADSPEC) GOTO 300
      GOTO 200
C
C     cursor
C
 600  IF (.NOT.LSOL) THEN
      WRITE (ILUER,8000)
C      CALL TONE                          ! NOT AVAILABLE IN BATCH MODE
      ENDIF
C      CALL DICUR (RA,IRA,RSCR,LSOL)      ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
 1000 FORMAT (' -------- SD3D module ----------',/
     &        '  (Three-dimensional sinkdiscs) '//,
     &' Maximum number of sinkdiscs:    ',I3,/,
     &' <F1> = Help',/,
     &' HEAD',/,
     &' DISCHARGE',/,
     &' <Esc> or QUIT',/,
     &' >')
 2000 FORMAT (' ***ILLEGAL OR MISSING PARAMETERS in SD3D module:',/,
     &        ' ',80A1)
 3000 FORMAT (' ***ILLEGAL COMMAND in SD3D module:',/,' ',80A1)
 3500 FORMAT (' ***ERROR in SD3D module: radius ',G11.4,
     &' is too small.',/,' Disc is ignored!',/)
 4000 FORMAT (' *** ERROR: Too many discs (max=',I3,
     &'), disc not stored.')
 5000 FORMAT (' Enter X,Y,Z of center   Radius    Exfiltration rate',
     &' [label]',/)
 6000 FORMAT (' Enter X,Y,Z of center   Radius    Head     Distance',
     &' [label]',/,
     &' (distance from the center where head is specified)',/)
 6500 FORMAT (' >')     
 7000 FORMAT (' ***DISC OUTSIDE AQUIFER: elevation of disc=',G11.4,/,
     &        ' aquifer base=',G11.4,' aquifer top=',G11.4,/)
 8000 FORMAT (' ***WARNING: no valid solution, data likely in error.')
C
      END
