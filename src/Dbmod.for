C     Last change:  HMH  20 Jan 2016    3:23 pm
c ------------------------------------------------------------------------------
c
c     This file contains the following routines or functions:
c
c     DBDAT     block data 
c     DBIN      input of inhomogeneities or horizontal barriers
c     DBDOUBLPT removes vertices that occur on top of another vertex
c     DBRECH    routine generates coefficients to add specified exfiltration rate to inhomogeneties.
c     DBORIEN   ensure that line doublets are ordered with the domain to the left
c     DBSORT    sort nested closed domains from outside to inside
c     ldominsidedom  true if one domain inside the other
C
c
c                 ----------------------------------------
c ---------------------
c dbmat.for contains:
c ---------------------
c     DBKFAC generates the jump conditions for the transmissivity
c     DBCZC  generates control points and sets equation type
c     DBMAT  generates matrix equations
c     DBKNO  generates the known vector
c     DBSUB  substitutes the solution vector in the common blocks
c     RFINTLABLEFT  integral for linear strength along line doublet
c     RFINTLABRIGHT integral for linear strength along line doublet
c     RFINTLABMU    integral for quadratic strength along line doublet
c     RFDBNFCONDITION generates normal flow condition for matrix coefficients
c     RFDBCONDUCTANCE calculates the conductance for a slurry wall
c
c                 ----------------------------------------
c ---------------------
c dbfun.for contains:
c ---------------------
c     CDBOM  complex potential due to all line doublets
c     DBQI   discharge vector due to all line doublets
c     CFDBF  coefficient function F for discharge potential
c     CFDBG  coefficient function G for discharge potential
c     CFDBS  coefficient function S for discharge potential
c     CFDBFD coefficient function F for derivative of discharge potential
c     CFDBGD coefficient function G for derivative of discharge potential
c     CFDBSD coefficient function S for derivative of discharge potential
c     RFDBPT driver to add all contributions due to added exfiltration to discharge potential
c     RFDBRP returns contribution to the discharge potential due to added exfiltration for one inhomogeneity
c                 ----------------------------------------
c ----------------------
c dbnflow.for
c ----------------------
C     RFNFDB       returns flow across a line due to all line doublets and line dipoles
c     RFNFDBGAM    coefficient function for flow across line due to added exfiltration rate for inhom.
c     RFNFDBGAMI   coefficient function for flow across line when inside inhomogeneity
c     RFNFDBGAMO   coefficient function for flow across line when outside inhomogeneity
c     DBINTERSECT  calculates intersection of line with line doublet and the flow in the branch cut
c     DFDBAREA     calculates area cut of inhomogeneity by a line
c                ------------------------------------------
c ---------------------
c dbio.for contains:
c ---------------------
c     DBIO          write/reads contents of common blocks to disk
c     INHOMEXTRACT  writes all line doublet data to file *.xtr
c     DBDATIO       writes all inhomogeneity and slurry wall (horizontal barrier) data to *.dat file
c     DBOUT         formatted write of all commomblock data to a specified logical unit (debugging) 
c                 ----------------------------------------
c ---------------------
c dbcheck.for contains:
c ---------------------
c     DBERROR    driver for reporting errors at control points
c     RFDBER     returns the error at control points for a string
c                 ----------------------------------------
c ---------------------
c dbnear.for contains:
c ---------------------
c     DBNEAR    handles pathline trace near a inhomogeneity or slurry wall
c                 ----------------------------------------
c ---------------------
c dbserv.for contains:
c ---------------------
c     DBPREP    driver for assigning data to common block for line doublet functions
c     DBFILL    fills common blocks for line doublets on either side of a single node
c     LINDOM    true if inside a inhomogeneity, passes relevant data for that inhomogeneity
c     LFDBIN    true is inside inhomogeneity string
c     RFDBGM    returns exfiltration rate at point due to all inhomogeneities
c     LDBOUTSIDE    used to determine if point is in the inactive zone for a closed slurry wall (horizontal barrier)
c     DBCHECKBOTTOM stops program if bottom elevation of inhomogeneity is at or above aquifer top 
c                 ----------------------------------------
c ---------------------
c dbextra.for contains:
c ---------------------
c     DBCHECK  interactive routine to display errors at comtrol points
c     DBPLOT   plot all line doublets
c     DBCUR    interactive graphics routine to display and manipulate date on screen
c     ELSERACH finds line doublet nearest to point
c     PTSEARCH finds vertex nearest to point 
c --------------------------------------------------------------------------------------------------------
c
      BLOCK DATA DBDAT
c
c --------------------------------------------------------------------------------------------------------
c
      IMPLICIT NONE
      INCLUDE 'DBCOM.INC'
      DATA NDB,NDBSTR /2*0/
      DATA DPI,D2PI /.3141592653589793D+1,.6283185307179586D+1/
      DATA CDI,DBEPS  /(0.0D+1,1.0D0),1.0D-3/
      DATA LDBPLT,ldbmatrix /.TRUE.,.false./
      DATA DBSTR,DBQSTR,DBSEN /NDBZMX*0.0D0,NDBZMX*0.0D0,NDBZMX*0.0D0/
      DATA CDDBSS,DBQSTS /NDBZMX*(0.0D+00,0.0D+00),NDBZMX*0.0D+00/
      DATA rdberrstrt,rdberrcntr /ndbzmx*0.0,ndbzmx*0.0/ ! initialze errors to zero
      DATA DBAVS  /NDBSMX*0.0D+00/
      DATA RDBGAM /NDBSMX*0.0D+00/
      DATA RDBRAD /NDBSMX*0.0D+00/
      DATA DBAREA /NDBSMX*0.0D+00/
      DATA CDBZ0  /NDBSMX*(0.0D0,0.0D0)/
      DATA RDBOFFSETN,RDBOFFSETC /0.5,0.05/ ! defines location of integration intervals on slurry walls
      DATA rdboffsetinhom /0.05/
      DATA ISTR1OLD,ISTR2OLD /0,0/  ! starting values for arguments in LDOMINSIDEDOM
      END
c
C -------------------------------------------------------------------------------------------------------      
c
      SUBROUTINE DBIN (LSOL,RA,IRA,RSCR)
c
c --------------------------------------------------------------------------------------------------------
C
C     Input routine for line doublet/dipole strings which define areal
C     inhomogeneity regions (difference in the hydraulic conductivity,
C     porosity, and areal recharge).
C     Domains may occur inside of each other, but may not overlap.
C     Domains with different hydraulic conductivity may not have common boundaries.
C
      IMPLICIT NONE
      INTEGER(4) IRA,JUMP,I,INODEX
      LOGICAL LSOL,LBAD,LWRONG,LPOR,LABELG,LCLOSED,lporositydefault
      REAL(8) RA,RSCR,RDUM,RVAR,RFPOR,R4,R3,R2,R1,RFTOP
      COMPLEX(8) CZEX,CZ,CVAR,C1
      CHARACTER(1) AWORD(33),APAR(8),APAR2(28),ALINEDB(80)
      CHARACTER(16) AFILEG
      INCLUDE 'DBCOM.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      DIMENSION RA(IRA,*),RSCR(*)
      DATA AWORD /'?',' ',
     &            'I','N','H','O',' ',
     &            'T','R','A','N',' ',
     &            'P','L','O','T',' ',
     &            'Q','U','I','T',' ',
     &            'C','U','R','S',' ',
     &            'S','L','U','R',' ',
     &            ATERM/
      DATA APAR  /'O','N',' ','O','F','F',' ',ATERM/
      DATA APAR2 /'O','P','E','N',' ',
     &            'C','L','O','S','E','D',' ',
     &            'I','N','S','I','D','E',' ',
     &            'O','U','T','S','I','D','E',' ',
     &            ATERM/
C
      LERROR=.FALSE.
      LMISS=.FALSE.
C      CALL CLEARSCREEN                          ! NOT AVAILABLE IN BATCH MODE
  10  IF (LERROR) WRITE (ILUER,1300) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      if (lucon) WRITE (ILUME,1500) NDBSMX,NDBZMX
      IF (LDBPLT.and.lucon) WRITE (ILUME,1501)
      IF (.NOT.LDBPLT.and.lucon) WRITE (ILUME,1502)
      if (lucon) WRITE (ILUME,1503)
      CALL INLINE
      CALL MATCH (AWORD,1,JUMP,LBAD)
C      CALL CLEARSCREEN                         ! NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LBAD) GOTO 15
      GOTO (10,13),JUMP
  13  WRITE (ILUER,1400) ALINE2
      GOTO 10
  15  GOTO (50,100,200,300,400,500,600),JUMP
C
C     help
C
 50   AFILE='DBHLP.HLP       '
C      CALL HELP                               ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
c
c     enter inhomogeneity domain (only k may jump) ! command is for compatibilty with older GFLOW data files
c
c     inhomogeneity  k  N  [n]  [label]
c
 100  IF (NDBSTR.EQ.NDBSMX) THEN
C ------------no more place for strings (domains) 
      WRITE (ILUER,2000) NDBSMX
      GOTO 10
      ENDIF
      RDUM=RVAR(2)
      IF (LERROR) GOTO 10
C ------------found k value, create new string (domain)      
      LSOL=.FALSE.
      NDBSTR=NDBSTR+1
      RDBK(NDBSTR)=RDUM
      IDOMAINTYPE(NDBSTR)=1
      RDUM=RVAR(3)
      IF (LERROR) GOTO 10
C ------------found added recharge value, store
      RDBGAM(NDBSTR)=RDUM
      RDUM=RVAR(4)            
C ------------store porosity if found, else store later,
C             also check for presence of a label
      LWRONG=(RDUM.LE.0.0.OR.RDUM.GT.1.0).AND..NOT.(LERROR.OR.LMISS)
      IF (LWRONG) WRITE (ILUER,7000)
      IF (LERROR.OR.LMISS.OR.LWRONG) THEN
        LPOR=.FALSE.
        IF (LWRONG) THEN
          CALLGETFN (5)
        ELSE
          CALL GETFN(4)
        ENDIF
        LABELG=.NOT.LMISS
        LERROR=.FALSE.
      ELSE
        RDBP(NDBSTR)=RDUM
        LPOR=.TRUE.
        CALL GETFN(5)
        LABELG=.NOT.LMISS
      ENDIF
      AFILEG=AFILE
      NDBSTI(NDBSTR)=0
      IDBSTA(NDBSTR)=NDB+1
C -------------read new node (coordinate pair)      
 120  if (lucon) WRITE (ILUME,1800) NDBSTR
      CALL INLINE
      DO I=1,80
      ALINEDB(I)=ALINE2(I)
      END DO
      CALL MATCH (AWORD,1,JUMP,LBAD)
      IF (.NOT.LBAD) THEN
C -------------command word found instead of node      
        IF (NDBSTI(NDBSTR).LT.3) THEN
C -------------however, no string completed yet; do not enter string
        WRITE (ILUER,3000) NDBSTI(NDBSTR)
        NDB=IDBSTA(NDBSTR)-1
        NDBSTR=NDBSTR-1
        goto 15
        ENDIF
C -------------remove points that are too close to other points
  131   CALL DBDOUBLPT (NDBSTR,INODEX,CZEX)
        IF (INODEX.GT.0) THEN
          WRITE (ILUER,6000) CZEX
          GOTO 131
        ENDIF
        IF (INODEX.EQ.-1) THEN  ! domain too small: do not enter string
          WRITE (ILUER,6010) NDBSTR,ADBLAB(NDB)
          NDB=NDB-3
          NDBSTR=NDBSTR-1
          goto 15
        ENDIF
C -------------correct orientation of new string (if necessary)      
        CALL DBORIEN
C -------------generate recharge coefficients      
        CALL DBRECH (NDBSTR)      
C -------------sort strings to ensure that when strings are nested,
C              they are ordered from outside to inside.
        CALL DBSORT
        GOTO 15
      ENDIF
      LBAD=.FALSE.
      CZ=CVAR(1)
      IF (LERROR) THEN
      WRITE (ILUER,4000) ALINEDB
      NDB=IDBSTA(NDBSTR)-1
      NDBSTR=NDBSTR-1
      LERROR=.FALSE.
      LMISS=.FALSE.      
      return
      ENDIF
C -------------coordinate pair found, now see if there is still place
      IF (NDB.EQ.NDBZMX) THEN
      WRITE (ILUER,5000) NDBZMX
      NDBSTR=NDBSTR-1
      RETURN
      ENDIF
C --------------get porosity if not yet stored
      IF (.NOT.LPOR) RDUM=RFPOR(CZ)
C --------------add node to string (domain)      
      NDB=NDB+1
      NDBSTI(NDBSTR)=NDBSTI(NDBSTR)+1
      CDBZ(NDB)=CZ
      CALL PLWIND (CZ)                    ! update max. window
      IF (.NOT.LPOR) THEN
C --------------store porosity if not yet done      
      RDBP(NDBSTR)=RDUM
      LPOR=.TRUE.
      ENDIF
C --------------store global base elevation (for old data file with no base jumps)
      CALL GVPAR (R1,R2,R3,RDUM,R4,C1)
      RDBB(NDBSTR)=RDUM
C ---------------store label if present
      CALL GETFN (3)
      IF (LMISS) THEN
        LMISS=.FALSE.
        IF (LABELG) ADBLAB(NDB)=AFILEG
      ELSE
        ADBLAB(NDB)=AFILE
      ENDIF            
      GOTO 120
C
C      enter transmissivity domain (both k and b may jump)
c
c     transmissivity  k  b  h  N  [n]  [label]
C      
 200  IF (NDBSTR.EQ.NDBSMX) THEN
C ------------no more place for strings (domains) 
      WRITE (ILUER,2000) NDBSMX
      GOTO 10
      ENDIF
      RDUM=RVAR(2)
      IF (LERROR) GOTO 10
C ------------found k value, create new string (domain)      
      LSOL=.FALSE.
      NDBSTR=NDBSTR+1
      RDBK(NDBSTR)=RDUM
      IDOMAINTYPE(NDBSTR)=1
      RDUM=RVAR(3)
      IF (LERROR) GOTO 10
C ------------found bottom elevation value, store
      RDBB(NDBSTR)=RDUM
      RDUM=RVAR(4)
      IF (LERROR) GOTO 10
C ------------found estimated average head for the domain, store
      RDBH(NDBSTR)=RDUM
      RDUM=RVAR(5)
      IF (LERROR) GOTO 10
C ------------found added recharge value, store
      RDBGAM(NDBSTR)=RDUM
      RDUM=RVAR(6)
      lporositydefault=rdum.EQ.-9999.0   ! added 5/18/04 because of change in GUI, see also next statement
C ------------store porosity if found, else store later,
C             also check for presence of a label
      LWRONG=(RDUM.LE.0.0.OR.RDUM.GT.1.0).AND..NOT.(LERROR.OR.LMISS)
     &.AND..not.lporositydefault ! don't check if default porosity is used
      IF (LWRONG) WRITE (ILUER,7000)
      IF (LERROR.OR.LMISS.OR.LWRONG.or.lporositydefault) THEN
        LPOR=.FALSE.
        IF (LWRONG) THEN
          CALLGETFN (7)
        ELSE
          CALL GETFN(6)
        ENDIF
        LABELG=.NOT.LMISS
        LERROR=.FALSE.
      ELSE
        RDBP(NDBSTR)=RDUM
        LPOR=.TRUE.
        CALL GETFN(7)
        LABELG=.NOT.LMISS
      ENDIF
      AFILEG=AFILE
      NDBSTI(NDBSTR)=0
      IDBSTA(NDBSTR)=NDB+1
C -------------read new node (coordinate pair)      
 220  if (lucon) WRITE (ILUME,1800) NDBSTR
      CALL INLINE
      DO I=1,80
      ALINEDB(I)=ALINE2(I)
      END DO
      CALL MATCH (AWORD,1,JUMP,LBAD)
      IF (.NOT.LBAD) THEN
C -------------command word found instead of node      
        IF (NDBSTI(NDBSTR).LT.3) THEN
C -------------however, no string completed yet; do not enter string
        WRITE (ILUER,3000) NDBSTI(NDBSTR)
        NDB=IDBSTA(NDBSTR)-1
        NDBSTR=NDBSTR-1
        goto 15
        ENDIF
C -------------remove points that are too close to other points
  231   CALL DBDOUBLPT (NDBSTR,INODEX,CZEX)
        IF (INODEX.GT.0) THEN
          WRITE (ILUER,6000) CZEX
          GOTO 231
        ENDIF
        IF (INODEX.EQ.-1) THEN  ! domain too small: do not enter string
          WRITE (ILUER,6010) NDBSTR,ADBLAB(NDB)
          NDB=NDB-3
          NDBSTR=NDBSTR-1
          goto 15
        ENDIF
C -------------correct orientation of new string (if necessary)
        CALL DBORIEN
C -------------generate recharge coefficients
        CALL DBRECH (NDBSTR)      
C -------------sort strings to ensure that when strings are nested,
C              they are ordered from outside to inside.
        CALL DBSORT
        GOTO 15
      ENDIF
      LBAD=.FALSE.
      CZ=CVAR(1)
      IF (LERROR) THEN
      WRITE (ILUER,4000) ALINEDB
      NDB=IDBSTA(NDBSTR)-1
      NDBSTR=NDBSTR-1
      LERROR=.FALSE.
      LMISS=.FALSE.      
      return
      ENDIF
C -------------coordinate pair found, now see if there is still place
      IF (NDB.EQ.NDBZMX) THEN
      WRITE (ILUER,5000) NDBZMX
      NDBSTR=NDBSTR-1
      RETURN
      ENDIF
C --------------get porosity if not yet stored
      IF (.NOT.LPOR) RDUM=RFPOR(CZ)
C --------------check for valid aquifer base (when k>0), else abort entering domain
      IF (RDBK(NDBSTR).GT.0.0.AND.RDBB(NDBSTR).GE.RFTOP(CZ)) THEN
      WRITE (ILUER,5010) RDBB(NDBSTR),NDBSTR
      NDB=IDBSTA(NDBSTR)-1
      NDBSTR=NDBSTR-1
      RETURN
      END IF
C --------------add node to string (domain)      
      NDB=NDB+1
      NDBSTI(NDBSTR)=NDBSTI(NDBSTR)+1
      CDBZ(NDB)=CZ
      CALL PLWIND (CZ)                    ! update max. window
      IF (.NOT.LPOR) THEN
C --------------store porosity if not yet done      
      RDBP(NDBSTR)=RDUM
      LPOR=.TRUE.
      ENDIF
C ---------------store label if present
      CALL GETFN (3)
      IF (LMISS) THEN
        LMISS=.FALSE.
        IF (LABELG) ADBLAB(NDB)=AFILEG
      ELSE
        ADBLAB(NDB)=AFILE
      ENDIF            
      GOTO 220
C
C     set plotting of inhomogeneities on or off
C
 300  CALL MATCH (APAR,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
      LDBPLT=JUMP.EQ.1
      GOTO 10
C
C     return to main program
C
 400  CONTINUE
c -------------set line-doublet tolerance
      call setdbeps()
C -------------calculate RDBTFACNOD, RDBTFACCTR, RDBBFACNOD and RDBBFACCTR for each node and element center (used in matrix solution)
      CALL DBKFAC
      RETURN 
C
C     cursor
C
 500  IF (.NOT.LSOL) THEN
C      CALL TONE              ! NOT AVAILABLE IN BATCH MODE
      WRITE (ILUER,1450)
      LERROR=.TRUE.
      ENDIF
c      CALL DBCUR (RA,IRA,RSCR,LSOL)              ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C    Enter a slurry wall
C
 600  IF (NDBSTR.EQ.NDBSMX) THEN
C ------------no more place for strings (domains) 
      WRITE (ILUER,2000) NDBSMX
      GOTO 10
      ENDIF
      CALL MATCH (APAR2,2,JUMP,LBAD)
      LERROR=LBAD
      IF (LERROR) GOTO 10
C ----------- found open or closed slurry wall, start a string
      LSOL=.FALSE.
      NDBSTR=NDBSTR+1
      IF (JUMP.EQ.1) THEN
      IDOMAINTYPE(NDBSTR)=2
      LCLOSED=.FALSE.
      ENDIF
      IF (JUMP.EQ.2) THEN
      IDOMAINTYPE(NDBSTR)=3
      RDBT(NDBSTR)=0.0    ! solution on inside and outside (parameter "closed")
      LCLOSED=.TRUE.
      ENDIF
      IF (JUMP.EQ.3) THEN
      IDOMAINTYPE(NDBSTR)=3
      RDBT(NDBSTR)=-9999.0  ! solution on inside only (parameter "inside")
      LCLOSED=.TRUE.
      ENDIF
      IF (JUMP.EQ.4) THEN
      IDOMAINTYPE(NDBSTR)=3
      RDBT(NDBSTR)=+9999.0  ! solution on outside only (parameter "outside")
      LCLOSED=.TRUE.
      ENDIF
      RDBGAM(NDBSTR)=0.0 ! add no recharge to slurry walls
      RDUM=RVAR(3)
      IF (LERROR) GOTO 10
C ------------found k value, store
      IF (ABS(RDBT(NDBSTR)).EQ.9999.0) THEN ! one sided model domain, k must be zero
        IF (RDUM.NE.0.0) THEN
          WRITE (ILUER,8000) NDBSTR,RDUM
          RDUM=0.0
        END IF
      END IF
      RDBK(NDBSTR)=RDUM
      RDUM=RVAR(4)
      IF (LERROR) GOTO 10
C ------------found width, store
      RDBW(NDBSTR)=RDUM
      RDUM=RVAR(5)
      IF (LERROR) GOTO 10
C ------------ found porosity
      LWRONG=(RDUM.LE.0.0.OR.RDUM.GT.1.0)
      IF (LWRONG) THEN       ! out of range, store back ground porosity later
      IF (ABS(RDBT(NDBSTR)).NE.9999.0) WRITE (ILUER,7000) ! if not inside/ouside domian issue warning
      ELSE
      RDBP(NDBSTR)=RDUM      ! store porosity
      ENDIF
      RDUM=RVAR(6)
      IF (LERROR) GOTO 10
C ------------ found bottom elevation of slurry wall, store
      RDBB(NDBSTR)=RDUM
      CALL GETFN(7)
      LABELG=.NOT.LERROR
      IF (LERROR) GOTO 10
C ------------ found label, store
      AFILEG=AFILE
      NDBSTI(NDBSTR)=0
      IDBSTA(NDBSTR)=NDB+1
C -------------read new node (coordinate pair)      
 620  if (lucon) WRITE (ILUME,1800) NDBSTR
 621  CALL INLINE
      DO  I=1,80
       ALINEDB(I)=ALINE2(I)
      END DO
      CALL MATCH (AWORD,1,JUMP,LBAD)
      IF (.NOT.LBAD) THEN
C -------------command word found instead of coordinate pair
       IF((LCLOSED.AND.NDBSTI(NDBSTR).LT.3).OR.NDBSTI(NDBSTR).LT.2) THEN
C -------------however, no string completed yet; do not enter string
        WRITE (ILUER,3000) NDBSTI(NDBSTR)
        NDB=IDBSTA(NDBSTR)-1
        NDBSTR=NDBSTR-1
        GOTO 15
        ENDIF
C -------------remove points that are too close to other points
  631   CALL DBDOUBLPT (NDBSTR,INODEX,CZEX)
        IF (INODEX.GT.0) THEN
          WRITE (ILUER,6000) CZEX
          GOTO 631
        ENDIF
        IF (INODEX.EQ.-1) THEN  ! domain too small: do not enter string
          WRITE (ILUER,6010) NDBSTR,ADBLAB(NDB)
          NDB=NDB-3
          NDBSTR=NDBSTR-1
          GOTO 15
        ENDIF
C --------------- now store background porosity if necessary
        IF (LWRONG) THEN
        RDBP(NDBSTR)=RFPOR(CZ)
        END IF
        IF (.NOT.LCLOSED) NDBSTI(NDBSTR)=NDBSTI(NDBSTR)-1
C -------------correct orientation of new string (if necessary)      
        IF (LCLOSED) CALL DBORIEN
C -------------sort strings to ensure that when strings are nested,
C              they are ordered from outside to inside.
        IF (LCLOSED) CALL DBSORT
        GOTO 15
      ENDIF
      LBAD=.FALSE.
      CZ=CVAR(1)
      IF (LERROR) THEN ! Error in reading coordinates, abort entering of string
      WRITE (ILUER,4000) ALINEDB
      NDB=IDBSTA(NDBSTR)-1
      NDBSTR=NDBSTR-1
      LERROR=.FALSE.
      LMISS=.FALSE.      
      RETURN
      ENDIF
C -------------coordinate pair found, now see if there is still place
      IF (NDB.EQ.NDBZMX) THEN
      WRITE (ILUER,5000) NDBZMX
      NDBSTR=NDBSTR-1
      RETURN
      ENDIF
C --------------add node to string (domain)
      NDB=NDB+1
      NDBSTI(NDBSTR)=NDBSTI(NDBSTR)+1
      CDBZ(NDB)=CZ
      CALL PLWIND (CZ)                    ! update max. window
C ---------------store label if present
      CALL GETFN (3)
      IF (LMISS) THEN
        LMISS=.FALSE.
        IF (LABELG) ADBLAB(NDB)=AFILEG
      ELSE
        ADBLAB(NDB)=AFILE
      ENDIF            
      GOTO 620
C 
 1000 FORMAT (80A1)
 1300 FORMAT (' ***ILLEGAL or MISSING PARAMETERS in inhomogeneity ',
     &'module:',/,' ',80A1)
 1400 FORMAT (' ***ILLEGAL COMMAND in inhomogeneity module:',/,
     &        ' ',80A1)
 1450 FORMAT (' ***WARNING: no valid solution, data likely in ERROR !'/)
 1500 FORMAT ('                -------- INHOMOGENEITY module --------'/
     &' Maximum number of strings:',I4,/,
     &' Maximum number of doublets:',I4,/,
     &' Available commands:',/,
     &' <F1> = Help',/,
     &' INHOMOGENEITY (k) (added recharge if <0) [porosity] [label]',/,
     &' SLURRYWALL (open/closed) (k) (w) (n) (b) (label',/,
     &' CURSOR')
 1501 FORMAT (' PLOT on  (off)')
 1502 FORMAT (' PLOT off (on)')
 1503 FORMAT (' <Esc> or QUIT',/,' >')
 1800 FORMAT (' Enter x, y, and [label] of the next node of string ',
     &          I3,/,' >')
 2000 FORMAT (' ***ERROR in INHOMOGENEITY: too many string (max.=',I3,
     &        ')',/,' this string is ignored!',/)
 3000 FORMAT (' ***ERROR in INHOMOGENEITY: only ',I1,' coordinate ',
     &        'pair(s) entered.',/,' this string is ignored!',/)
 4000 FORMAT (' ***ERROR in INHOMOGENEITY: expected a coordinate pair, '
     &        ,'but got:',/,80A1,/,' this string is ignored!',/)
 5000 FORMAT (' ***ERROR in INHOMOGENEITY: too many elements (max.=',
     &     I5,')',/,' this string is ignored!',/)
 5010 FORMAT (' ***ERROR in INHOMOGENEITY: base elevation ',E14.7,
     & ' of string ',I3,' is at or above aquifer top.',/,
     & ' Domain not entered.')
 6000 FORMAT(' ***WARNING in INHOMOGENEITY: point ',2F11.1,/,
     &' coincides with another point and has been removed.',/)     
 6010 FORMAT (' ***WARNING in INHOMOGENEITY: string ',I3,
     &' with label ',a16,' not entered.',/)
 7000 FORMAT (' ***ERROR in INHOMOGENEITY: porosity out of range;',/,
     &' must be between 0 and 1 (background porostity used).')
 8000 FORMAT (' ***WARNING in SLURRY WALLS: string number ',I3,/,
     &' is declared an INSIDE or OUTSIDE domain and should have a zero',
     &' conductivity.',/,
     &' Found k=',E14.7,', but k has been set to zero.',/)
      END
c
c -------------------------------------------------------------------------------------------------
c

      SUBROUTINE DBDOUBLPT (ISTR,INODEX,CZ)
c
C --------------------------------------------------------------------------------------------------
C
C     Routine checks for double points (points too close) in
C     string ISTR. If point is found, it is removed, while its
C     address is stored in INODEX and its coordinates are stored 
C     in CZ. If no point is found INODEX=0 
C     If domain has only three points left still contains a line doublet
C     that is too small INODEX=-1
C     Note: only one point is removed per call! Recall routine to
C     check for additional double points.
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INODEX,INOD1,INODL,INOD,II1,II2,I
      REAL(8) RTEST,RFGRTOL,RLEN
      COMPLEX(8) CZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      INODEX=0
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      RTEST=RFGRTOL()
      DO 30 INOD=INOD1,INODL
      ii1=inod
      IF (INOD.EQ.INODL.AND.IDOMAINTYPE(ISTR).NE.2) THEN
      ii2=inod1
      ELSE
      ii2=inod+1
      ENDIF
      rlen=abs(cdbz(ii1)-cdbz(ii2))
      if (rlen.le.rtest.and.rlen.gt.0.0) then
        if (lucon) then
          write (iluer,2000) ii1,ii2,adblab(ii1),adblab(ii2)
          call inline
        else
          write (iluer,2001) ii1,ii2,adblab(ii1),adblab(ii2)
        endif
      endif
      IF (rlen.eq.0.0) THEN ! points coincide!!!!
      IF (NDBSTI(ISTR).EQ.3) THEN   ! domain too small, don't remove point
      WRITE (ILUER,1000) ISTR
      inodex=-1
      RETURN
      ENDIF      
      CZ=CDBZ(INOD)
      INODEX=INOD
      write (iluer,1001) inod,adblab(inod)
 1001 format (' DBDOUBLPT1: adblab(',i4,')=',a16)
      IF (INOD.LT.INODL-1) THEN
      INODL=INODL-1
      DO 20 I=INOD+1,INODL  ! leave first point intact; no update of IDBSTA
      CDBZ(I)=CDBZ(I+1)
      ADBLAB(I)=ADBLAB(I+1)
   20 CONTINUE
      ENDIF
      NDBSTI(ISTR)=NDBSTI(ISTR)-1
      NDB=NDB-1
      RETURN
      ENDIF
   30 CONTINUE
      RETURN
 1000 FORMAT ('+***ERROR in DBDOUBLPT: cannot remove double point,',
     &' domain ',I3,' too small.')
 2000 FORMAT ('+***WARNING point ',i3,' may be too close to point '
     &         ,i3,'                           ',/,
     &' labels are: ',a16,2x,a16,/,' Press <Enter> to continue.')
 2001 FORMAT(' ***WARNING in INHOMOGENEITY: point ',i3,
     &' may be too close to point ',i3,/,
     &' labels are: ',a16,2x,a16,/)
         
      END
c
c -------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBRECH (ISTR)
c
c --------------------------------------------------------------------------------------------------
C
C     Routine generates coefficients to add recharge to inhomogeneities
C      
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD1,INODL,NODM1,NOD,INOD
      REAL(8) RV1,RV2,RU1,RU2,RDUM,RFDBRP
      COMPLEX(8) CZZ,CZ1,CZ2,CZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
C
      IF (ISTR.LT.0.OR.ISTR.GT.NDBSTR) THEN
      WRITE (ILUER,1000) ISTR,NDBSTR
      RETURN
      ENDIF      
      DBAREA(ISTR)=0.0
      INOD1=IDBSTA(ISTR)
      NOD=INOD1+NDBSTI(ISTR)-1
      NODM1=NOD-1
      DO 10 INOD=INOD1,NOD          ! calculate center of string and initialize
      CDBZ0(ISTR)=CDBZ0(ISTR)+CDBZ(INOD)
      CDDBSS(INOD)=CMPLX(0.0D+00,0.0D+00)
      DBQSTS(INOD)=0.0D+00
  10  CONTINUE
      CDBZ0(ISTR)=CDBZ0(ISTR)/NDBSTI(ISTR)
      DO 20 INOD=INOD1,NODM1          ! calculate DBAREA, aimag(CDDBSS), RDBRAD
      CZZ=CDBZ(INOD)-CDBZ0(ISTR)
      RU1=REAL(CZZ)
      RU2=AIMAG(CZZ)
      CZZ=CDBZ(INOD+1)-CDBZ0(ISTR)
      RV1=REAL(CZZ)
      RV2=AIMAG(CZZ)
      DBAREA(ISTR)=DBAREA(ISTR)+0.5*(RU1*RV2-RU2*RV1)
      CDDBSS(INOD+1)=CDDBSS(INOD+1)-CDI*RDBGAM(ISTR)*DBAREA(ISTR)  ! imaginary strength parameter at INOD+1
  20  CONTINUE
      RU1=RV1
      RU2=RV2
      CZZ=CDBZ(INOD1)-CDBZ0(ISTR)
      RV1=REAL(CZZ)
      RV2=AIMAG(CZZ)
      DBAREA(ISTR)=DBAREA(ISTR)+0.5*(RU1*RV2-RU2*RV1)
      RDBRAD(ISTR)=DBAREA(ISTR)/DPI
      CZ2=CDBZ(INOD1)         ! calculate real(CDDBSS), DBQSTS
      RDUM=RFDBRP(CZ2,ISTR)   ! added inside potential at vertex INOD1
      DO 30 INOD=INOD1,NODM1
      CZ1=CZ2
      CDDBSS(INOD)=CDDBSS(INOD)-RDUM  ! real (linear) strength parameter for vertex INOD
      CZ2=CDBZ(INOD+1)
      RDUM=RFDBRP(CZ2,ISTR)   ! added inside potential at vertex INOD+1
      CZ=0.5*(CZ1+CZ2)
      DBQSTS(INOD)=-RFDBRP(CZ,ISTR)
     &             -0.5*(REAL(CDDBSS(INOD))-RDUM) ! real (parabolic) sterngth parameter at line-doublet center
  30  CONTINUE
      CDDBSS(NOD)=CDDBSS(NOD)-RDUM
      CZ=0.5*(CDBZ(NOD)+CDBZ(INOD1))
      DBQSTS(NOD)=-RFDBRP(CZ,ISTR)
     &             -0.5*REAL(CDDBSS(NOD)+CDDBSS(INOD1))  
      RETURN
 1000 FORMAT (' *** ERROR in DBRECH: ISTR out of range.',/,
     &' ISTR is ',I3,' Current number of strings is ',I3)      
      END
c
c ---------------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBORIEN
c
c ---------------------------------------------------------------------------------------------------------------
C
C     Check and correct orientation of last doublet string
C     When traversing along the line doublets from z1 to z2 the
C     domain must be on the LEFT-HAND side.
C      
      IMPLICIT NONE
      INTEGER(4) INOD,INOD1,NOD,IEN,I
      REAL(8) RDUM
      COMPLEX(8) CZ,C1,C2,CDUM
      CHARACTER(16) ADBDUM
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      IF (NDB.EQ.0) RETURN
      IF (NDBSTR.EQ.0) THEN
      WRITE (ILUER,1000) NDB
      RETURN
      ENDIF
      IF (IDOMAINTYPE(NDBSTR).EQ.2) RETURN ! open string
      INOD1=IDBSTA(NDBSTR)
      NOD=INOD1+NDBSTI(NDBSTR)-1
      C1=0.5*(CDBZ(INOD1+1)+CDBZ(INOD1))
      C2=0.5*(CDBZ(INOD1+1)-CDBZ(INOD1))
      CZ=C2*(0.0,1.0E-3)+C1
      CALL DBPREP (CZ)
      RDUM=0.0
      DO 10 INOD=INOD1,NOD
      RDUM=RDUM+AIMAG(CDDBLN(INOD))
  10  CONTINUE
C ------------if inside RDUM~2pi, else RDUM~0, test for >1!  
      IF (RDUM.GT.1.0) RETURN
C ------------wrong orientation, correct!      
      IEN=NDBSTI(NDBSTR)/2
      DO 20 I=1,IEN
      CDUM=CDBZ(NOD-I+1)
      ADBDUM=ADBLAB(NOD-I+1)
      CDBZ(NOD-I+1)=CDBZ(INOD1+I-1)
      ADBLAB(NOD-I+1)=ADBLAB(INOD1+I-1)
      CDBZ(INOD1+I-1)=CDUM
      ADBLAB(INOD1+I-1)=ADBDUM
  20  CONTINUE
C ------------reset CZ0 in DBPREP, since data in DBCOM has changed
      CZ=(1.0E21,1.0E21)
      CALL DBPREP(CZ)      
      RETURN
 1000 FORMAT (' ***ERROR in DBORIEN: ',I4,' doublet nodes, but no ',
     &        'domains encountered!')
      END
c
C -----------------------------------------------------------------------------------------------------------
c
      SUBROUTINE DBSORT
c
c -----------------------------------------------------------------------------------------------------------
C
C     If last domain is inside another, it is sequenced to follow
C     the domain in which it is located.
c     This routine is called after the completion of the entry of a new
c     "transmissivity inhomogeneity,"
!     Added 9/3/99: Ignore recharge only domains.
C      
      IMPLICIT NONE
      INTEGER(4) INODL1,ISTREN,ISTR0,ISTR,ISTRL,I1,ISTNEW,NSTNEW,
     &           IDOMNEW,IEN,I,istrlast
      LOGICAL L1,ldominsidedom
      REAL(8) RKNEW,RPNEW,RBNEW,RTNEW,RWNEW,RHNEW,DAVNEW,RDGNEW,
     &        DBANEW,RDRNEW
      COMPLEX(8) CZ,C1,C2,CDZNEW
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      IF (NDB.EQ.0) RETURN
      if (ndbstr.eq.1) return ! first string don't try to sort
      IF (NDBSTR.EQ.0) THEN
      WRITE (ILUER,1000) NDB
      RETURN
      ENDIF
      
!     Test if last string is a "recharge only" domain
      if (rdbk(ndbstr).eq.-9999.0.and.rdbb(ndbstr).eq.-9999.0)
     & return
C     ..........................................................
C     Test if last string occurs inside an existing string
C     ..........................................................
C ---------test for being inside any of the previous strings
      ISTREN=NDBSTR-1
      ISTR0=0
      istrlast=ndbstr
      DO 20 ISTR=1,ISTREN
      IF (IDOMAINTYPE(ISTR).NE.1) GOTO 20 ! not a transmissivity inhomogeneity, skip
c      write (ilume,1003) istren,ndbstr,istrlast,istr
c 1003 format ('dbsort3: istren,ndbstr,istrlast,istr ',4(i5,2x))
      L1=ldominsidedom(istrlast,istr).and..not.! added 9/3/99: ignore "recharge only" domains
     &(rdbk(istr).eq.-9999.0.and.rdbb(istr).eq.-9999.0)
      IF (L1) ISTR0=ISTR
  20  CONTINUE
      IF (ISTR0.EQ.0) GOTO 100
C ---------not inside any other domain, but may be on the outside
C          of existing domains
      IF (ISTR0.EQ.ISTREN) RETURN
C ---------inside last domain, thus already in proper sequence: return
      GOTO 200
C ---------reorder the strings by placing last string after string ISTR0
C     ...........................................................
C     Test for any of the existing strings to occur inside of the
C     last string
C     ...........................................................
  100 ISTRL=NDBSTR
      DO 110 ISTR=1,ISTREN
      ISTR0=ISTR-1
      I1=IDBSTA(ISTR)
      C1=0.5*(CDBZ(I1+1)+CDBZ(I1))
      C2=0.5*(CDBZ(I1+1)-CDBZ(I1))
      CZ=C2*(0.0,1.0E-3)+C1
      CALL DBPREP (CZ)
      L1=LDBINS (ISTRL).and..not.
     &(rdbk(istr).eq.-9999.0.and.rdbb(istr).eq.-9999.0) !  ignore "recharge only" domains
      IF (L1) GOTO 200
C -----------string ISTR occurs inside last string: reorder      
 110  CONTINUE
C -----------no existing strings occur inside last string: return
      do i=1,ndbstr
      inodl1=idbsta(i)
c      write (ilume,1001) i,adblab(inodl1)
c 1001 format (' dbsort1: string ',i3,' has starting label ',a16)
      end do
      RETURN
C     ...................................................................
C     Reorder pointers and arrays to place last string after string ISTR0
C     ...................................................................
  200 ISTNEW=IDBSTA(NDBSTR)
      NSTNEW=NDBSTI(NDBSTR)
      IDOMNEW=IDOMAINTYPE(NDBSTR)
      RKNEW=RDBK(NDBSTR)
      RPNEW=RDBP(NDBSTR)
      RBNEW=RDBB(NDBSTR)
      RTNEW=RDBT(NDBSTR)
      RWNEW=RDBW(NDBSTR)
      RHNEW=RDBH(NDBSTR)
      DAVNEW=DBAVS(NDBSTR)
      RDGNEW=RDBGAM(NDBSTR)
      DBANEW=DBAREA(NDBSTR)
      CDZNEW=CDBZ0(NDBSTR)
      RDRNEW=RDBRAD(NDBSTR)
      IEN=NDBSTR-ISTR0-1
      DO 230 I=1,IEN
      IDBSTA(NDBSTR-I+1)=IDBSTA(NDBSTR-I)
      NDBSTI(NDBSTR-I+1)=NDBSTI(NDBSTR-I)
      IDOMAINTYPE(NDBSTR-I+1)=IDOMAINTYPE(NDBSTR-I)
      RDBK(NDBSTR-I+1)=RDBK(NDBSTR-I)
      RDBP(NDBSTR-I+1)=RDBP(NDBSTR-I)
      RDBB(NDBSTR-I+1)=RDBB(NDBSTR-I)
      RDBT(NDBSTR-I+1)=RDBT(NDBSTR-I)
      RDBW(NDBSTR-I+1)=RDBW(NDBSTR-I)
      RDBH(NDBSTR-I+1)=RDBH(NDBSTR-I)
      DBAVS(NDBSTR-I+1)=DBAVS(NDBSTR-I)
      RDBGAM(NDBSTR-I+1)=RDBGAM(NDBSTR-I)
      DBAREA(NDBSTR-I+1)=DBAREA(NDBSTR-I)
      CDBZ0(NDBSTR-I+1)=CDBZ0(NDBSTR-I)
      RDBRAD(NDBSTR-I+1)=RDBRAD(NDBSTR-I)
 230  CONTINUE
      IDBSTA(ISTR0+1)=ISTNEW
      NDBSTI(ISTR0+1)=NSTNEW
      IDOMAINTYPE(ISTR0+1)=IDOMNEW
      RDBK(ISTR0+1)=RKNEW
      RDBP(ISTR0+1)=RPNEW
      RDBB(ISTR0+1)=RBNEW
      RDBT(ISTR0+1)=RTNEW
      RDBW(ISTR0+1)=RWNEW
      RDBH(ISTR0+1)=RHNEW
      DBAVS(ISTR0+1)=DAVNEW
      RDBGAM(ISTR0+1)=RDGNEW
      DBAREA(ISTR0+1)=DBANEW
      CDBZ0(ISTR0+1)=CDZNEW
      RDBRAD(ISTR0+1)=RDRNEW
C ----------reset CZ0 in DBPREP, since data in DBCOM has changed.
      CZ=(1.0E21,1.0E21)      
      CALL DBPREP(CZ)
c
c     debugging
c
      do i=1,ndbstr
      inodl1=idbsta(i)
c      write (ilume,1002) i,adblab(inodl1)
c 1002 format (' dbsort2: string ',i3,' has starting label ',a16)
      end do
      RETURN
C      
 1000 FORMAT (' ***ERROR in DBSORT: ',I4,' doublet nodes, but no ',
     &        'domains encountered!')
      END
c
c
C -----------------------------------------------------------------------------------------------------------
c
      LOGICAL function ldominsidedom(istr1,istr2)
c
c -----------------------------------------------------------------------------------------------------------
C
C    True if domain istr1 is inside istr2. In testing it ensures that the vertex
c    of domain istr1, used in the logical function ldbins is not shared by domain
c    istr2!
C      
      IMPLICIT NONE
      INTEGER(4) istr1,istr2,i,istart,iend,j,jstart,jend
      COMPLEX(8) CZ
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      save
c
c     Find vertex of domain istr1 that is not shared with a vertex of domain istr2
c
      if (istr1old.eq.istr1.and.istr2old.eq.istr2) return ! we did this already
      istr1old=istr1
      istr2old=istr2
      ldominsidedom=.true.
      if (istr1.eq.istr2) RETURN ! sure, domain is inside itself
      istart=idbsta(istr1)
      iend=istart+ndbsti(istr1)-1
      do 10 i=istart,iend
        jstart=idbsta(istr2)
        jend=jstart+ndbsti(istr2)-1
        cz=cdbz(i)
        do j=jstart,jend
          if (cz.eq.cdbz(j)) GOTO 10 ! coincides, try other CZ
        end do
        GOTO 20 ! cz does not coincide, now see is inside domain istr2
 10   continue
c     all vertices of domain istr1 coincide with vertices of domain istr2
      return
 20   continue  ! check if cz is inside or not
      call dbprep(cz)
      ldominsidedom=ldbins(istr2)
      return
      end
c
c -------------------------------------------------------------------------------------------------
c

      SUBROUTINE setdbeps
c
C --------------------------------------------------------------------------------------------------
C
C     Routine sets a line-doublet tolerance for the physical domain
c     based of the smallest line-doublet length
C
      IMPLICIT NONE
      INTEGER(4) ISTR,INOD1,INODL,INOD,II1,II2
      REAL(8) RLEN
      INCLUDE 'DBCOM.INC'
      INCLUDE 'LUSYS.INC'
      if (ndbstr.eq.0) return
      rlen=1.0d+31
      do 20 istr=1,ndbstr
      INOD1=IDBSTA(ISTR)
      INODL=INOD1+NDBSTI(ISTR)-1
      DO 10 INOD=INOD1,INODL
      ii1=inod
      IF (INOD.EQ.INODL.AND.IDOMAINTYPE(ISTR).NE.2) THEN
      ii2=inod1
      ELSE
      ii2=inod+1
      ENDIF
      rlen=MIN(rlen,abs(cdbz(ii1)-cdbz(ii2)))
   10 CONTINUE
   20 continue
      dbeps=0.001D+0*rlen
      RETURN
      END


