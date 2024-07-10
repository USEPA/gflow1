C     Last change:  HMH  24 Jan 2018    5:32 pm
c     This file contains the following routines and functions
c
c     GFLOW1 main program
c     SUBROUTINE SETUP      input of analytic elements and progDRAm execution instructions as defined in .dat file
c     SUBROUTINE INITS      initialization of constants (called by SETUP)
c
c -------------------  
c GFMAT.FOR contains:
c -------------------
c
c     SUBROUTINE SOLUT           driver for matrix solution process
c     subroutine condition       normalizes the matrix to diagonal elements of 1
c     subroutine adjust          adjusts the known vector based on normalized matrix
c     SUBROUTINE DECOMP          decomposes the matrix for gaussian elimination
c     SUBROUTINE SOLVE           forward elimination and back substitution (use after DECOMP)
c     subroutine Gauss_Seidel    double sweep Gauss-Seidel matrix solver
c     subroutine gfread_convergence_file
c                                reads file with GW solution convergence criteria 
c
c -------------------
c GFFUN.FOR contains:
c -------------------
c
c     COMPLEX FUNCTION COMEGA       all functions that contribute to the complex potential
c     COMPLEX FUNCTION CFBIGZ       returns Z for complex line integDRAls
c     REAL FUNCTION RFPOT           returns the real part of the complex potential: the discharge potential
c     REAL FUNCTION RFPSI           returns the imaginary part of the complex potential: the stream function
c     REAL FUNCTION RFHEAD          returns the head measured with respect to the datum mean sea level (msl)
c     REAL FUNCTION RFDSCH          returns the absolute value of the discharge vector
c     SUBROUTINE DISCH              returns the components if the discharge vector and the vertical specific discharge
c          ENTRY DISCHSPEC          returns specific discharge vector due to 2D functions (called by SDISCH)
c     SUBROUTINE SDISCH             returns the complete specific discharge vector (incl. approx. vertical discharge)
c     SUBROUTINE VELOC              returns the velocity vector (incl. approx. vertical velocity)
c     REAL FUNCTION RFK             complete elliptic integral of the first kind K(k)
c     REAL FUNCTION RFE             complete elliptic integral of the second kind E(k)
c     REAL FUNCTION RFPI            complete elliptic integral of the third kind Pi(k,alpha)
c     REAL FUNCTION ERFC            complementary error function Erfc(x)
c     SUBROUTINE LEGEND             Legendre functions: Pn(0), Pn(x) and P'n(x)
c     REAL FUNCTION RE1             Exponential integral function E1(x)
c     FUNCTION RF3DSP               scalar product function (3D vectors)
c     REAL FUNCTION SQROOT          shell for intrinsic square root function
c     REAL FUNCTION RFSCALAR        returns the scalar product of 2D vectors (in complex numbers)
c     FUNCTION RFSP3D               scalar product function (three dimensional)
c     REAL FUNCTION RFVP2D          two dimensional vector product
c     SUBROUTINE VPRO3D             three dimensional vector product
c
c ---------------------
c GFNFLOW.FOR contains:
c ---------------------
c
c     REAL function rfnormalflow      driver for calculating the flow across a straight line element
c     REAL FUNCTION RFBRANCH          returns jump in PSI across line for a unit withdrawal at a point
c     SUBROUTINE BRANCHCUT            returns jump in PSI across line for linesink or linedipole
c     REAL FUNCTION RFAREATRIANGLE    returns area for triangle
c     REAL FUNCTION RFNUMNF           returns a numerical approximation of the flow across a straight line element
c     LOGICAL FUNCTION LINSECTCIRCLE  true if line intersects a circle
c     LOGICAL FUNCTION LINSECTLINE    true if two lines intersect
c     REAL*8 FUNCTION DFA             returns signed area of triangle 
c
c ------------------
c GFIO.FOR contains:
c ------------------
c
c     BLOCK DATA SETLUSYS   set logical units
c     SUBROUTINE LUOUT      writes common block LUSYS (formatted) for debugging purposes
c     SUBROUTINE CHKOUT     writes common block CHCK (formatted) for debugging purposes
c     SUBROUTINE INITIO     reads and writes some initialization data to and from disc
c     SUBROUTINE DATIO      driver for creating a GFLOW input script file (.dat)
c     SUBROUTINE STDATIO    writes header to a .dat file (called in DATIO))
c     SUBROUTINE ENDATIO    writes trailer to a .dat file (called in DATIO))
c     SUBROUTINE WINDATIO   writes windows parameters to a .DAT file
c     subroutine extract    driver to generate the *.xtr file (extract file)
c     SUBROUTINE GFIO       reads or writes contents of MAIN.INC to .sol file
c     SUBROUTINE MNOUT      writes the contents of the common block MAIN (formatted) for debugging purposes
c     SUBROUTINE DEBUG      driver to write common blocks (formatted) for debugging purposes
C     SUBROUTINE MATOUT     writes the matrix and associated arrays (formatted) for debugging purposes
c     SUBROUTINE HALT       prints a diagnostic in GFLOW.OPS and stops program execution
c     SUBROUTINE BIO        driver for generating the binary .sol file (solution file)
c     SUBROUTINE BUFIN      initializes buffer for binary IO
c     SUBROUTINE BUFEX      flushes buffer at the end of binary IO
c     SUBROUTINE BUFIO2     adds integer*2 words to buffer for binary IO
c     SUBROUTINE BUFIO4     adds integer*4 words to buffer for binary IO
c     SUBROUTINE BUFIOR     adds real*4 words to buffer for binary IO (currently adds real*8)
c     SUBROUTINE BUFIOA     adds character words to buffer for binary IO
c     SUBROUTINE BUFIOC     adds complex*4 words to buffer for binary IO (currently complex*8)
c     SUBROUTINE BUFIOD     adds real*8 words to buffer for binary IO
c     SUBROUTINE BUFIOCD    adds complex*8 words to buffer for binary IO
c     SUBROUTINE BUFIOL     adds logical*4 words to buffer for binary IO
c     SUBROUTINE BUFIO      adds words to buffer for binary IO (called by other BUFIO routines)
c     SUBROUTINE XFERW2     routine swaps bytes of data between two arrays (called by BUFIO)
c     SUBROUTINE SYSIOW     performs binary read or write operation
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
c     LOGICAL FUNCTION LSEP .TRUE. if ACHAR is a separator. Converts all uppercases to lower cases
c              ENTRY LSEPUL like LSEP, but maintains upper and lower case
c     FUNCTION AFUPC        changes lower case to upper case
c     SUBROUTINE INLINE     reads ALINE and sets pointers (calls TIDY)
c     SUBROUTINE RALINE     reads string from input buffer ALINE
c     SUBROUTINE REVAR      shell for calling RVAR
c     SUBROUTINE RAVAR      return a single character from ALINE
c     CHARACTER(1) FUNCTION ACHAR    return a single character from ALINE
c     COMPLEX FUNCTION CVAR return a complex number from ALINE
c     SUBROUTINE RILINE     return a string from ALINE
c     SUBROUTINE ISETVR     shell for calling IVAR
c     SUBROUTINE SWITCH     driver for SWTCHS
c     SUBROUTINE SWTCHS     reassign logical units for input and output functions
c     SUBROUTINE LUSWAP     swap lugical units, called by SWTCHS
C
C     NOTE: grid related IO routines are in GFGRID.FOR
c           trace related IO routines are in GFTRACE.FOR
c
c
c
c --------------------
c GFGRID.FOR contains:
c --------------------
c
c     SUBROUTINE GRINIT     grid initialization routine (for contour plots)
c     SUBROUTINE GRIDIN     driver for filling grid with function values
c     SUBROUTINE GFFILL     adds function values to grid
c          ENTRY PFILL      same as GFFILL, but plots a dot map
c     SUBROUTINE MINMAX     determines minimum and maximum values in the grid
c     SUBROUTINE GRIDIO     reads or writes contents of common /GRID/ to .SOL file
c     SUBROUTINE RAIO       reads or writes the (grid) array RA to a binary .GRI file
c     SUBROUTINE GRIO       reads or writes contents of common /GRID/ to .GRI file
c     SUBROUTINE GRTEMPIO   reads data from .GRI files into dummy storage for subtracting grids
c     SUBROUTINE RATEMPIO   reads or writes the array (grid) RA2 to a file
c     SUBROUTINE DIFFGRID   reads a grid and subtracts subtracts it from the current grid
c     SUBROUTINE GRIDOUT    writes the contents of common /GRID/ formatted for debugging purposes
c
c
c ---------------------
c GFTRACE.FOR contains:
c ---------------------
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
c
c --------------------
c GFSERV.FOR contains:
c --------------------
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
c
c ---------------------
c GFEXTRA.FOR contains:
c ---------------------
c
c     BLOCK DATA PLSETD          initializes plotting parameters for single screen graphics
c     SUBROUTINE CURGET          reads characters in ALINE while in cursor mode
c          ENTRY CURGETUL        same as CURGET, but maintain upper and lower case characters in ALINE
c     SUBROUTINE CURDAT          brings up cursor for data retrieval
c     SUBROUTINE PLSETCUR        resets text cursor and blanks lines
c     SUBROUTINE SETCUR          sets cursor by use of ANSI.SYS routine
c     SUBROUTINE GRAPHEDIT       for adding contour levels to contours and text for labeling features
c     SUBROUTINE PLOTID          plots program and user identification
c     SUBROUTINE PLOTANNOTATION  adds annotations to plot from buffers when plotting on a hard copy device
c     SUBROUTINE CURSOR          brings up cross-hair for setting plot windows
c     SUBROUTINE CVRXIX          convert world coordinates to screen coordinates and vice versa
c     INTEGER(4) FUNCTION IRQLOC generate cross-hairs
c     SUBROUTINE DRAW            draw a vector on graphics device or write to surfer boundary line file
c     SUBROUTINE DRAWW           same as DRAW, but not encapsulated (see plot routines GFPLOTS and GFPLOTSS)
c     SUBROUTINE DRCIRC          draws a circle
c     SUBROUTINE GFCIRCLE        called by DRCIRC
c     SUBROUTINE SCALE           sets scale for plotting
c     SUBROUTINE DOTS            pixel translation for HP plotters
c     SUBROUTINE PLOTON          activate plotting mode
c          ENTRY PLOTOF          deactivate plotting mode
c     series or routines that form graphoria graphics encapsulations in GSS graphics routines
c     SUBROUTINE LABCOL          determines color and line type from a analytic element label
c     SUBROUTINE LAYOUT          draws layout of analytic elements
c     SUBROUTINE LOGO            selects proper logo for graphics screen
c     SUBROUTINE PLGRID          plots grid on the graphics device
c     SUBROUTINE GSSIN           activate plotting (GSS plotting calls)
c          ENTRY GSSOUT          deactivate plotting (GSS plotting calls)
c     SUBROUTINE CHANGECOLOR     color adjustments (GSS plotting calls)
c          ENTRY CHANGELINETYPE  line type adjustment (GSS plotting calls)
c          ENTRY CHANGELINEWIDTH line width adjustment (GSS plotting calls)
c          ENTRY CHANGEMARKER    select marker (GSS plotting calls)
c     SUBROUTINE GFPLOT          plot line on graphics device in pixels (encapsulated)
c          ENTRY GFPLOTT         plot line on graphics device in pixels (not encapsulated)
c     SUBROUTINE GFPLOTS         plot line on graphics device in world coordinates (encapsulated)
c          ENTRY GFPLOTSS        plot line on graphics device in world coordinates (not encapsulated)
c     SUBROUTINE INSECT          intersection of lines for graphics clipping
c     SUBROUTINE PAGE            clear graphics screen
c     SUBROUTINE TONE            sound tone
c     SUBROUTINE HIPIN           activate HP plotter (HP plotting routines)
c     SUBROUTINE HIPOUT          deactivate HP plotter (HP plotting routines)
c     SUBROUTINE HPXY            plot a line (HP plotting routines)
c     SUBROUTINE GFNEWPEN        change pen (HP plotting routines)
c     SUBROUTINE NEWLINE         select new line type (HP plotting routines)
c     SUBROUTINE NEWCOL          select new color (HP plotting routines)
c     SUBROUTINE PLOTIO          reads or writes PLOT.INC common blocks to .SOL file
c     SUBROUTINE PLOTOUT         writes contents of PLOT.INC common blocks (formatted) for debugging
c     SUBROUTINE CHECK           interactive driver to check boundary conditions
c     SUBROUTINE PIEZOMETERS     printing or plotting of piezometer data for calibration
c     SUBROUTINE PIEZLEGEND      printing of difference in modeled/observed head on screen
c     SUBROUTINE ZMARK           generate vertical elevation ticmark
c     SUBROUTINE TMARK           generate time increment ticmark
c     SUBROUTINE CLEARSCREEN     clears screen (interactive program use)
c
c
c

c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     always update MODIFICATION DATE on line 235
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c
c     Update revision data in format # 1200 (line # 642 and 654)
c     Update expiration date in line # 215
c     Update number of equations in line # 69 and 97
c     Update program version/identification in lines 36 - 64
c
C              ********************************
C              *                              *
C              *       G  F  L  O  W  1       *
C              *                              *
C              ********************************
C
C     Single aquifer Dupuit-Forhheimer model for conjunctive surface water
C     and groundwater flow in both confined and unconfined aquifers. 
C
      IMPLICIT NONE
      INTEGER NEQCOM,NEQMAX,NCHA,I,NCHAR1,NCHAR2,IDUM,
     &        NEQ,IERR,IMXSZE,IVAR,ngridsize
      CHARACTER(20) AIDENT1,AIDENT2
      CHARACTER(127) ACLINE
      CHARACTER(8) aDate
      CHARACTER(10)aTime
      COMMON /IDENT/ AIDENT1,AIDENT2,NCHAR1,NCHAR2
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'TRACOM.INC'
      INCLUDE 'MAIN.INC'
c
      NEQMAX=6000             ! maximum number of equations   (increased from 5000 to 6000 on 5/5/2017)
      ngridsize=200           ! set (maximum) grid size
c
      CALL INITS
C
C     parse the GFLOW1.EXE command line
C
      CALL GETCL (ACLINE)
      NCHA=NBLANK(ACLINE)
      DO 3 I=1,132
      ALINE(I)=' '
  3   CONTINUE      
      IF (NCHA.GT.0) THEN
        DO 5 I=1,NCHA
        ALINE(I)=ACLINE(I:I)
  5     CONTINUE
        CALL TIDY
        IDUM=IVAR(1)
        IF (LERROR.OR.LMISS) THEN  ! no number of equations found
          NEQCOM=NEQMAX
          LERROR=.FALSE.
          LMISS=.FALSE.
          CALL GETFN(1)  ! try to read filename as first parameter
        ELSE
          NEQCOM=MIN(IDUM,NEQMAX) ! apply number of equations
          CALL GETFN(2)           ! try to read filename as second parameter
        ENDIF
      ELSE  ! no command line parameters found
        NEQCOM=NEQMAX
        LMISS=.TRUE.
      ENDIF
c
c     The next statement renders the GUI estimate of the number of equations obsolete.
C     The parsing of the command line for a number of equations is maintained for reasons
c     of backward compatibility.
c     Use of the command line number of equations can be reinstated by deleting the next statement.
c
      neqcom=ngridsize
c
      AFILE2=AFILE              ! save filename for use in setup routine
      NEQ=MAX(NEQCOM,ngridsize) ! ensure minimum grid size for contouring
C
          call date_and_time(aDate,aTime)    ! get date and time to stamp .sol, .mtr, and .dec files
          aDateTime=aDate//aTime(1:8)
c
c REVISE THE LOGIC ABOVE TO IGNORE COMMANDLINE INPUT OF THE NUMBER OF EQUATIONS
c
      neq=ngridsize ! ignore all logic for getting the number of equations from the command line or setting a maximum
      CALL SETUP(NEQ,ngridsize)  ! enter the driver routine for GFLOW
      STOP
      END
C
C -----------------------------------------------------------------------------
C
      SUBROUTINE SETUP (neq,igrsize)
C
C -----------------------------------------------------------------------------
C
      IMPLICIT NONE
      LOGICAL LDATACHANGE,LBAD,LESCAP,LTWL,lDirectFromDiskTemp,
     &        LDBBOTTOMOFF,LDBBOTTOMON,LFLSLAKEITERATIONS,lfleakage,
     &        L1,L2,L3,L4,L5,L6,L7,L8,lflksolving
      INTEGER IGRSIZE,neq,IPIV,KTYPE,NCHA,IEXT,ILAB,I,J,
     &        JUMP,ITER,ITEROUT,JUMPOLD,IDUM3,IVAR,iterm1,ierr,
     &        isolvetype,ireporting,ISTATUS,nequation,matsize,
     &        iticks,iticks1,iticks2,idum,neq_initial,nstr,ndrscrsize,
     &        idim_dra
      REAL(8) DRA,DRSCR1,DRSCR2,DRFAC,RQ,RDUM,RFHEAD,
     &        RDUM1,RDUM2RQI,RFPOT,RVAR,RFHGHT,RFBASE,RFCONPT,RQI,RDUM2,
     &        RFINTERFACE
      COMPLEX(8) CZI,CALPH,CDUM,CZ,CVAR,COMEGA
      CHARACTER(1) AWORD(211),APAR(16)
      CHARACTER(256) AMESSAGE
      CHARACTER(8) aDate
      CHARACTER(10)aTime
      CHARACTER(16) aDateTimeStamp
c      CHARACTER(8) ADATE
      DIMENSION RQI(3),idim_dra(2)
      ALLOCATABLE DRA(:,:),DRSCR1(:),DRSCR2(:),DRFAC(:,:),CALPH(:),
     &            CZI(:),KTYPE(:),IPIV(:)
      INCLUDE 'COM3D.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MAIN.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'MATCH.INC'     
      INCLUDE 'DWMN.INC'
      INCLUDE 'TRACOM.INC'
      EXTERNAL RFPOT,RFPSI,RFHEAD,RFDSCH
      DATA LDATACHANGE /.FALSE./
      DATA AWORD  /
     .             '?',' ',
     .             'G','R','I','D',' ',  'S','W','I','T',' ',
     .             'A','Q','U','I',' ',  'W','E','L','L',' ',
     .             'L','I','N','E',' ',  'I','N','H','O',' ',
     .             'L','A','Y','O',' ',  'S','O','L','V',' ',
     .             'C','H','E','C',' ',  'P','L','O','T',' ',
     .             'L','O','A','D',' ',  'S','A','V','E',' ',
     .             'S','T','O','P',' ',  'D','I','S','C',' ',
     .             'T','W','E','L',' ',  'P','P','W','E',' ',
     .             'P','O','T','E',' ',  'H','E','A','D',' ',
     .             'S','D','3','D',' ',  'S','P','E','C',' ',
     .             'V','E','L','O',' ',  'T','R','A','C',' ',
     .             'P','A','G','E',' ',  'D','E','B','U',' ',
     .             'C','U','R','S',' ',  'C','O','M','M',' ',
     .             'M','A','P',' ',      'W','A','T','E',' ',
     .             'T','I','T','L',' ',  'O','M','E','G',' ',
     .             'D','A','T','A',' ',  'V','I','E','W',' ',
     .             'D','O','S',' ',      'P','O','N','D',' ',
     .             'M','O','D','E',' ',  'S','I','N','K',' ',
     .             'E','X','T','R',' ',  'L','E','A','K',' ',
     .             'T','I','M','E',' ',  'F','A','S','T',' ',
     .             'I','N','T','E',' ',  'B','F','N','A',' ',
     .             ATERM/
      DATA APAR   /'G','R','O','U',' ',
     .             'B','A','S','E',' ',
     .             'C','O','N','J',' ',
     .             ATERM/
c
      save      ! IMPORTANT TO KEEP THE ALLOCATABLE ARRAYS!!!!
c
c      NXMAX=MIN0(NXMAX,neq) ! reduce max. (contour) grid size for small # equations
c      NYMAX=MIN0(NYMAX,neq)
c
C
C     Inintial allocation or matrix arrays necessary for contouring
C
      neq_initial=neq ! initial number of equations
      nstr=neq        ! initial number of strength parameters (unknowns or degrees of freedom)
      ndrscrsize=MAX(neq,1+(nstr*nstr)/31) ! in case we will use this also for matrix solution
      ALLOCATE (DRA(NEQ,NEQ),DRSCR1(ndrscrsize),
     &          DRSCR2(NEQ),DRFAC(4,NEQ),CALPH(NEQ),CZI(NEQ),
     &          KTYPE(NEQ),IPIV(NEQ),STAT=IERR)
      if (ierr.ne.0) then
      write (*,9190)  ! error in allocating arrays
       AMESS(1)='Error in initial allocation of grid arrays.'
       AMESS(2)='No ERROR.LOG file opened yet.'
       AMESS(3)='Insufficient RAM available. Execution halted.'
       CALL HALT(3) ! stop program execution for batch version
      end if
C
      LSOL=.FALSE.
      IF (LMISS.OR.LERROR) THEN ! failed to read a filename
        NCHA=0
      ELSE                      ! get filename length for parsing
        NCHA=NBLANK(AFILE2)
      ENDIF
      LMISS=.FALSE.
      LERROR=.FALSE.
      iext=min(ncha-1,3)
      IF (NCHA.NE.0) THEN   ! determine of .dat or .sol file and act accordingly
      IF (AFILE2(NCHA-iext:NCHA).EQ.'.SOL'.OR.
     &    AFILE2(NCHA-iext:NCHA).EQ.'.sol') THEN
      ASSIGN 650 TO ILAB
      ALINE(1)='R'
      ALINE(2)='E'
      ALINE(3)='A'
      ELSE
      ASSIGN 200 TO ILAB
      ALINE(1)='S'            ! assume command line is .dat file
      ALINE(2)='W'
      ALINE(3)='I'
      ENDIF
      ALINE(4)=' '
      J=0
      DO  5 I=5,NCHA+5
      J=J+1
      ALINE(I)=AFILE2(J:J)
   5  CONTINUE
      CALL TIDY
      GOTO ILAB
      ENDIF
      write (ilume,1200)
      if (lucon) CALL INLINE
 10   IF(LERROR) WRITE(ILUER,2000) ALINE2
      IF(LMISS) WRITE(ILUER,2001) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      LCURS=.FALSE.
      IF (.NOT.LSOL) LGRCOM=.FALSE.
      if (lucon) WRITE (ILUME,1300)
      CALL INLINE
  11  CALL MATCH(AWORD,1,JUMP,LBAD)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      IF(.NOT.LBAD) GOTO 15
      GOTO (10,13), JUMP
  13  CONTINUE
      WRITE(ILUER,3000) ALINE2
      LERROR=.FALSE.
      GOTO 10
  15  GOTO (100,150,200,250,300,350,400,450,500,550,
     &      600,650,700,750,960,850,950,960,960,900,
     &      960,960,990,992,993,1100,50,1150,1160,1170,
     &      960,1180,1190,1195,800,1197,800,1000,425,875,
     &      490,960,970),JUMP
C
C   Command summary
C
  50  AFILE='GFSUM.HLP       '
c      CALL HELP                    ! NOT AVAILABLE IN BATCH MODE
      GOTO 10 
C 
C   Help 
C 
 100  AFILE='GFHLP.HLP       '
c      CALL HELP                    ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C   Grid 
C
 150  CONTINUE
      write (*,1110)
 1110 FORMAT (/,' Filling grid for contour plot.')
      idim_dra=UBOUND (dra)
c      write (iluer,5678) idim_dra(1:2),igrsize
c 5678 format ('setup5678: max. dimension DRA: DRA(',i5,',',i5,')',
c     &         ' igrsize=',i5)
      if (idim_dra(2).lt.igrsize) then ! may happen for highly asymmetric matrices
        deallocate (dra)
        ALLOCATE (dra(igrsize,igrsize),stat=ierr)
        if (ierr.ne.0) then
          write (*,9192)  ! error in allocating arrays
          AMESS(1)='Error in allocating contour grid array.'
          AMESS(2)='Insufficient RAM available. Execution halted.'
          CALL HALT(2) ! stop program execution for batch version
        end if
      end if
c      call timer(iticks1)
      CALL GRIDIN (DRA,IGRSIZE)
c      call timer(iticks2)
c      iticks=iticks2-iticks1
c      write (iluer,1234) iticks
c 1234 format (' setup1234: time for grid=',i10)
      GOTO 10
C
C   Switch module
C 
 200  CALL SWITCH
c
c
c
       write (ilume,1111)
 1111  format
     &(' GFLOW Solver, latest update January 24, 2018',//,
c
     & ' Review errors in boundary conditions and close this file.',//)
c
c
c
      IF (.not.lucon) WRITE (*,1112)
 1112 format (/,' Reading input data.')
      GOTO 10
C
C   Input given parameter
C
 250  CALL GVIN (LSOL)
      IF (loadsol.AND..NOT.LSOL) LDATACHANGE=.TRUE.
      GOTO 10
C
C   Input Wells
C
 300  CALL WLIN (LSOL,DRA,IGRSIZE,DRSCR1)
      IF (loadsol.AND..NOT.LSOL) LDATACHANGE=.TRUE.
      GOTO 10
C
C   Input line sinks
C
 350  CALL LSIN (LSOL,DRA,IGRSIZE,DRSCR1)
      IF (loadsol.AND..NOT.LSOL) LDATACHANGE=.TRUE.
      GOTO 10
C
C   Input doublets
C
 400  CALL DBIN (LSOL,DRA,IGRSIZE,DRSCR1)
      IF (loadsol.AND..NOT.LSOL) LDATACHANGE=.TRUE.
      IF (LERROR) GOTO 10
      GOTO 10
C
C   Leakage elements
C
 425  CONTINUE
      matsize=neq_initial ! actual matrix size has not yet been determined
      CALL LKIN (LSOL,matsize,
     &  dra,drscr1,drscr2,drfac,czi,calph,ktype,nstr)
      GOTO 10
      LBAD=.TRUE.
      GOTO 10
C
C      Layout
C
 450  CALL PLAPAR(LESCAP)
      IF (LESCAP) GOTO 10
C      CALL PLOTON          ! NOT AVAILABLE IN BATCH MODE
C      CALL LAYOUT          ! NOT AVAILABLE IN BATCH MODE
C      CALL PLOTOF          ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
c      ENDIF
c      GOTO 451
C
C   Solve groundwater flow or streamflow problem
C
 490  CONTINUE   !                                  FASTSOLVE command
      isolvetype=1 ! default to solution over disk
      ireporting=0 ! default to no reporting
      GOTO 501
 500  CONTINUE   !                                  SOLVE command
      isolvetype=0 ! default to non-linear equations regenerating matrix and decoposing it at each iteration
      ireporting=1 ! default to reporting at each step
      
c
c     determine matrix dimensions
c
 501  Number_equations=0
      Number_strengths=0
      call lsmatsize(Number_equations,Number_strengths)
      !write (iluer,1006) Number_equations,Number_strengths
 1006 format (' setup4: LS:  Number_equations,Number_strengths ',2(i5))
      CALL DBCOLLOCATION_PREP() ! generate proper collocation points
      call dbkfac()
      call dbmatsize(Number_equations,Number_strengths)
      !write (iluer,1007) Number_equations,Number_strengths
 1007 format (' setup4: DB: Number_equations,Number_strengths ',2(i5))
      call wlmatsize(Number_equations,Number_strengths)
      !write (iluer,1005) Number_equations,Number_strengths
 1005 format (' setup4: WL: Number_equations,Number_strengths ',2(i5))
      call w3matsize(Number_equations,Number_strengths)
      !write (iluer,1008) Number_equations,Number_strengths
 1008 format (' setup4: W3: Number_equations,Number_strengths ',2(i5))
      call lkmatsize(Number_equations,Number_strengths)
      !write (iluer,1009) Number_equations,Number_strengths
 1009 format (' setup4: LK: Number_equations,Number_strengths ',2(i5))
      call gvmatsize(Number_equations,Number_strengths)
      !write (iluer,1004) Number_equations,Number_strengths
 1004 format (' setup4: GV: Number_equations,Number_strengths ',2(i5))

c
c     check to see if currently allocated arrays are sufficiently large else
c     deallocate and reallocate
c
      NEQ=Number_equations
      NSTR=Number_strengths
      idim_dra=UBOUND (dra)
      if (neq.gt.idim_dra(1).or.nstr.gt.idim_dra(2)) THEN ! must reallocate arrays
        l1=allocated(dra)
        l2=allocated(drscr1)
        l3=allocated(drscr2)
        l4=allocated(drfac)
        l5=allocated(calph)
        l6=allocated(czi)
        l7=allocated(ktype)
        l8=allocated(ipiv)
        IF (L1) DEALLOCATE (DRA)
        IF (L2) DEALLOCATE (DRSCR1)
        IF (L3) DEALLOCATE (DRSCR2)
        IF (L4) DEALLOCATE (DRFAC)
        IF (L5) DEALLOCATE (CALPH)
        IF (L6) DEALLOCATE (CZI)
        IF (L7) DEALLOCATE (KTYPE)
        IF (L8) DEALLOCATE (IPIV)
        ndrscrsize=MAX(neq,1+(nstr*nstr)/31)
        ALLOCATE (DRA(NEQ,NSTR),DRSCR1(ndrscrsize),
     &            DRSCR2(NEQ),DRFAC(4,NEQ),CALPH(NEQ),CZI(NEQ),
     &            KTYPE(NEQ),IPIV(NEQ),STAT=IERR)
        if (ierr.ne.0) then
          write (*,9191)  ! error in allocating arrays
          AMESS(1)='Error in allocating matrix arrays.'
          AMESS(2)='Insufficient RAM available. Execution halted.'
          CALL HALT(2) ! stop program execution for batch version
        end if
      end if
      if (lflksolving()) then
c     Open the file "gfmfnonitor.dat" to be used to monitor performance of
c     the conjunctive GFLOW - MODFLOW solution process.
c     Note: This file must exist and be empty prior to running GFLOW!
c
      open (UNIT=ilutmp,FILE='gfmfmonitor.dat',RECL=132,
     &      POSITION='append',iostat=ierr)
      if (ierr.ne.0) THEN ! failed to open file, abort
       call IOSTAT_MSG (ierr,amessage)
       write (iluer,1010) ierr,amessage
 1010  format ('Failed to open -gfmfmonitor.dat- in Setup.',/,
     & 'Error number is ',i5,' Error message is:',/,
     & a256,/)
       amess(1)='Could not open -gfmfmonitor.dat-'
       amess(2)='See IO error message in -error.log-'
       amess(3)='Program execution aborted.'
       call halt (3)
      end if
      end if
c
c     now parse the Solve command
c
c
c      solve (number iterations inner loop) (number of iterations outer loop) (solution type) (reporting flag)
c
      IF (LTWL()) then
       GOTO 10 ! do not solve if transient wells are present and time>0
       write (ilume,5050)
       write (*,5050)
      endif
      CALL MATCH (APAR,2,JUMP,LBAD)  ! start of parameter parsing
      IF (LBAD) THEN    ! no stream flow
        LBAD=.FALSE.
        JUMP=1
        ITER=IVAR(2)
        ITEROUT=IVAR(3)
        if (lerror) then ! to accomodate old WhAEM GUI command without isolve and ireporting codes.
         lerror=.false.
         lmiss=.false.
         iterout=0
        else
         isolvetype=ivar(4)
         ireporting=ivar(5)
        endif
      ELSE             ! conjunctive GW and stream flow
        ITER=IVAR(3)
        ITEROUT=IVAR(4)
         isolvetype=ivar(5)
         ireporting=ivar(6)
      ENDIF
      IF (LERROR) THEN
        LMISS=.FALSE.
        LERROR=.FALSE.
        ITER=0
        ITEROUT=0
        isolvetype=0 ! default to non-linear equations regenerating matrix and decomposing it at each iteration
        ireporting=1 ! default to reporting at each step
        WRITE (ILUER,4050)
      ENDIF
C                         ! end of parameter parsing, now setup solution procedure
      NOUTERLOOP=NOUTERLOOP+ITEROUT
      lDirectFromDisk=isolvetype.EQ.1
c     recreate matrix and decompose at every step when isolvetype=0
      lGaussSeidel=isolvetype.EQ.2
c     no error reporting at BC during iterations when ireporting=0
      lErrorReport=ireporting.EQ.1
c     assign LU 10 and 11, but 11 is not used when lDirectFromDisk=.false.
      if (nsol.eq.0) then ! create matrix files and open
         lwrongfile=.TRUE. ! set, since we test on it in SOLUT
         afile=TRIM(abasename) // ".mtr"
         open (UNIT=10,IOSTAT=istatus,ERR=512,FILE=afile,            ! using LU 10 for coefficient matrix file
     &        STATUS='REPLACE',ACCESS='TRANSPARENT',FORM='BINARY')
         afile=TRIM(abasename) // ".dec"
         open (UNIT=11,IOSTAT=istatus,ERR=513,FILE=afile,
     &        STATUS='REPLACE',ACCESS='TRANSPARENT',FORM='BINARY')  ! using LU 11 for decomposed matrix file
      else                                   ! solution exisits, just open files
         afile=TRIM(abasename) // ".mtr"
         open (UNIT=10,IOSTAT=istatus,ERR=514,FILE=afile,
     &        STATUS='OLD',ACCESS='TRANSPARENT',FORM='BINARY')
         afile=TRIM(abasename) // ".dec"
         open (UNIT=11,IOSTAT=istatus,ERR=515,FILE=afile,
     &        STATUS='OLD',ACCESS='TRANSPARENT',FORM='BINARY')
         lWrongFile=.false.
         read (UNIT=10,IOSTAT=istatus,ERR=516) aDateTimeStamp,nEquation
         if (aDateTimeStamp.ne.aDateTime) lWrongfile=.true.! incompatible .mtr file matrix must be regenerated
         if (lDirectFromDisk) then
          lwrongfile=.true.
          read (UNIT=11,IOSTAT=istatus,ERR=503) aDateTimeStamp,nEquation
          lwrongfile=.false. ! decomposed matrix file used before, now see if during latest run
          if (aDateTimeStamp.ne.aDateTime) lWrongFile=.true.! incompatible .dec file matrix must be regenerated
  503     CONTINUE ! could not read datestamp and # of equations, previous run must have been without .dec
         end if
         REWIND(10)
         REWIND(11)
         if (lWrongFile) then
c
          call date_and_time(aDate,aTime)    ! get date and time to stamp new .sol, .mtr, and .dec files
          aDateTime=aDate//aTime(1:8)
c
         end if
      end if
      call gfread_convergence_file() ! read the file "converge.tab"
c                    --------------------------------------------------
      IF (JUMP.EQ.1) THEN                 ! solve GROUNDWATER FLOW ONLY
c                    --------------------------------------------------
      CALL NOCONJUNCTIVE  ! set flag for baseflow calculations .FALSE.
      if (ldbbottomOFF()) THEN ! base jump temporarily converted to k jump
      write (*,6010) iter
      if (.not.lucon) write (ilume,6010) iter
      end if
      WRITE (*,6100)
      if (.not.lucon) write (ilume,6100)
      CALL SOLUT
     &(DRA,neq,nstr,DRSCR1,DRSCR2,DRFAC,CALPH,CZI,KTYPE,IPIV,ITER)
      if (ldbbottomON ()) then
      DO i=1,3    ! add 3 iterations
      lDirectFromDiskTemp=lDirectFromDisk
        if (i.eq.1.or..not.lquit) then ! only execute if solution did not yet converge
         WRITE (*,6100)
         if (.not.lucon) write (ilume,6100)
         if (i.eq.1) lDirectFromDisk=.false.  ! force new matrix decompsition
      CALL SOLUT
     &(DRA,neq,nstr,DRSCR1,DRSCR2,DRFAC,CALPH,CZI,KTYPE,IPIV,1)
         lDirectFromDisk=lDirectFromDiskTemp ! restore solution type
        endif
      enddo
      endif
c      write (iluer,1007)
c 1007 format ('gfmod7: writing MODFLOW grid files from GW solution.')
c      call lkwrite_upperheads()
c      call lkwrite_leakage()
c      call lkwrite_resistances()
c      call lkwrite_error()
      LGRCOM=.FALSE.    ! contour grid no longer compatible
      IF (NSOL.LT.0) THEN ! error detected during solution process
        NSOL=0
        if (lDirectFromDisk) then  ! close matrix files
          CLOSE (10)
          CLOSE (11)
        end if
        GOTO 10
      ENDIF
      IF (lquit) THEN
      WRITE (ILUME,6210)
      ENDIF
      if (.not.LerrorReport) then
      CALL SOLUT
     &(DRA,neq,nstr,DRSCR1,DRSCR2,DRFAC,CALPH,CZI,KTYPE,IPIV,0)   ! no solution, but give an error report
      end if
      LSOL=.TRUE.
      ENDIF   ! end solve GW flow only
c                     --------------------------------------------
      IF (JUMP.EQ.2) THEN                 ! solve STREAM FLOW ONLY
c                     --------------------------------------------
      IF (.NOT.LSOL) THEN
      WRITE (ILUER,6000)
        if (lDirectFromDisk) then    ! close matrix files
          CLOSE (10)
          close (11)
        end if
      GOTO 10
      ENDIF
      write (*,6200)
      if (.not.lucon) WRITE (ILUME,6200)
      iterm1=iter-1
      CALL GENSTREAMFLOW (ITERM1,lErrorReport,nsol,
     &                          ilstablelength,niterarray,relaxarray)
      call genstreamflow (1,.TRUE.,nsol, ! do last iteration with an error report
     &                          ilstablelength,niterarray,relaxarray)
      CALL NOCONJUNCTIVE
      IF (ITER.GT.0) THEN
      WRITE (ILUME,6300)
      LSOL=.FALSE.
      LGRCOM=.FALSE.    ! grid no longer compatible
c      CALL TIMER(ILAYER) ! ILAYER used as a code to identify a solution
      ENDIF
      ENDIF         ! end solve streamflow only
c                   -----------------------------------------------
      IF (JUMP.EQ.3) THEN                 ! solve CONJUNCTIVE GW/SW
c                   -----------------------------------------------
      IF (ITER.LT.2) THEN
       ITER=2
       WRITE (ILUER,6005)
      ENDIF
      if (ldbbottomOFF()) THEN ! base jump temporarily converted to k jump
       write (*,6010) iter
       if (.not.lucon) write (ilume,6010) iter
      end if
      call lsread_relaxation_file  ! read the file "relax.tab" if present
      call set_lakedone()   ! initialize logical LAKEDONE
c
c     start timing the solution procedure
c
  505 call timer(iticks1)
c
      DO I=1,ITER
      WRITE (*,6100)
      if (.not.lucon) write (ilume,6100)
      CALL SOLUT
     &(DRA,neq,nstr,DRSCR1,DRSCR2,DRFAC,CALPH,CZI,KTYPE,IPIV,1)
      LGRCOM=.FALSE.    ! grid no longer compatible
      IF (NSOL.LT.0) THEN ! error detected during solution process
        NSOL=0
        if (lDirectFromDisk) then   ! close matrix files
          CLOSE (10)
          close (11)
        end if
        GOTO 10
      ENDIF
      IF (lquit.and.i.ge.minimum_iterations) THEN
        WRITE (ILUME,6210)
        GOTO 510
      ENDIF
      IF (I.NE.ITER) THEN
        WRITE (*,6200)
        if (.not.lucon) WRITE (ILUME,6200)
        CALL GENSTREAMFLOW (1,lErrorReport,nsol,
     &                          ilstablelength,niterarray,relaxarray)
      ENDIF
      enddo
c
c     end timing solution procedure
c
 510  continue
      call timer(iticks2)
      iticks=iticks2-iticks1
      write (ilume,2222) iticks
 2222 format (' TOTAL SOLUTION TIME =                             ',i10,
     & ' E-2 seconds')
c
c     add extra solutions with base jump instead of surrogate k jump
c
      if (ldbbottomON ()) then
        lDirectFromDiskTemp=lDirectFromDisk
        DO I=1,3                ! adding extra iterations
         WRITE (*,6100)
         if (.not.lucon) write (ilume,6100)
         if (i.eq.1) lDirectFromDisk=.false.  ! force new matrix decomposition
      CALL SOLUT
     &(DRA,neq,nstr,DRSCR1,DRSCR2,DRFAC,CALPH,CZI,KTYPE,IPIV,1)
         lDirectFromDisk=lDirectFromDiskTemp ! restore solution type
         LGRCOM=.FALSE.    ! contour grid no longer compatible
         IF (NSOL.LT.0) THEN ! error detected during solution process
           NSOL=0
           if (lDirectFromDisk) then  ! close matrix files
             CLOSE (10)
             close (11)
           end if
           GOTO 10
         ENDIF
        IF (lquit) THEN
          WRITE (ILUME,6210)
          GOTO 511
        ENDIF
         IF (I.NE.ITER) THEN
           WRITE (*,6200)
           if (.not.lucon) WRITE (ILUME,6200)
           CALL GENSTREAMFLOW (1,lErrorReport,nsol,
     &                          ilstablelength,niterarray,relaxarray)
         ENDIF
        END do
      endif
c
c     end of logic to add iterations for base jump
c
  511 call lslakewaterbalance(.FALSE.,lErrorReport) ! write lake water balance at end of inner loop
c
      IF (LFLSLAKEITERATIONS(NOUTERLOOP).and..not.lakedone) GOTO 505
c
c     Solution procedure is complete, wrap up calculations and error reporting
c
      CALL SOLUT
     &(DRA,neq,nstr,DRSCR1,DRSCR2,DRFAC,CALPH,CZI,KTYPE,IPIV,0)    ! with number of iterations set to zero do not
c                                  resolve, but offer an error report
      WRITE (ILUME,6205)
      CALL GENSTREAMFLOW (0,.TRUE.,nsol,  ! calculate streamflow without corrections and
     &                          ilstablelength,niterarray,relaxarray)
c                                       set lErrorReport=true to force a report
      call lslakewaterbalance(.TRUE.,lErrorReport) ! write lake water balance at end of solution
      LSOL=.TRUE.
      ENDIF
      IF (LUCON) THEN
      WRITE (ILUME,6220)
      CALL INLINE
      ENDIF  ! end of conjunctive surface water and groundwater loop.
        if (lDirectFromDisk) then ! close matrix files
          CLOSE (10)
          close (11)
        end if
      call lkwrite_upperheads()
      call lkwrite_leakage()
      call lkwrite_resistances()
      call lkwrite_error()
      call lkgfmfmonitor()
      GOTO 10
C
C     failure to create matrix files, halt program
C
 512  call iostat_msg(istatus,amessage)
      write (*,9097) amessage ! error in allocating file
      write (ilume,9097) amessage
      AMESS(1)='Error allocating matrix file on disk.'
      AMESS(2)='See IO error in Message.log file. Execution halted.'
      CALL HALT(2) ! stop program execution for batch version
 513  call iostat_msg(istatus,amessage)
      write (*,9098) amessage ! error in allocating file
      write (ilume,9098) amessage
      AMESS(1)='Error allocating decomposed matrix file on disk.'
      AMESS(2)='See IO error in Message.log file. Execution halted.'
      CALL HALT(2) ! stop program execution for batch version
 514  call iostat_msg(istatus,amessage)
      write (*,9097) amessage ! error in allocating file
      write (ilume,9097) amessage
      AMESS(1)='Error opening matrix file on disk.'
      AMESS(2)='See IO error in Message.log file. Execution halted.'
      CALL HALT(2) ! stop program execution for batch version
 515  call iostat_msg(istatus,amessage)
      write (*,9098) amessage ! error in allocating file
      write (ilume,9098) amessage
      AMESS(1)='Error opening decomposed matrix file on disk.'
      AMESS(2)='See IO error in Message.log file. Execution halted.'
      CALL HALT(2) ! stop program execution for batch version
 516  call iostat_msg(istatus,amessage)
      write (*,9099) amessage ! error in allocating file
      write (ilume,9099) amessage
      AMESS(1)='Error reading from matrix file on disk.'
      AMESS(2)='See IO error in Message.log file. Execution halted.'
      CALL HALT(2) ! stop program execution for batch version
c
C   Check module
C
 550  IF (.NOT.LSOL) THEN
      WRITE (ILUER,6003)
      GOTO 10
      ENDIF
C      CALL CHECK (DRA,IGRSIZE,DRSCR1)  ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C   Contour Head, Potential, or Stream Functions
C
 600  IF (.NOT.LSOL) GOTO 1050
C      IF (.NOT.LGRCOM) THEN          ! NOT AVAILABLE IN BATCH MODE
C      WRITE (ILUER,3050)             ! NOT AVAILABLE IN BATCH MODE
C      GOTO 10                        ! NOT AVAILABLE IN BATCH MODE
C      ENDIF                          ! NOT AVAILABLE IN BATCH MODE
C      CALL PGRPAR(LESCAP)            ! NOT AVAILABLE IN BATCH MODE
C      IF (LESCAP) GOTO 10            ! NOT AVAILABLE IN BATCH MODE
C 630  CALL CNPLOT (DRA,IGRSIZE,DRSCR1) ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     load a solution file
C 
 650  write (*,1113)
 1113 format (/,' Read a solution from the disk.')
      CALL BIO (2,'READ ')
      write (*,1214)
 1214 format (/,' Solution has been read.')
      loadsol=.TRUE.
      LINALREADY=.TRUE.
c      nsol=0               ! this reset is necessary to make a resolve work correctly
      CDUM=(1.0D21,1.0D21)      ! reset previous point memory
      CALL DBPREP(CDUM)
      RDUM=RFPOT(CDUM)
      GOTO 10
C
C     save a solution file
C
 700  write (*,1114)
 1114 format (/,' Write a solution to the disk')
      CALL BIO (2,'WRITE')
      GOTO 10
c
c     stop program execution
c
  750 STOP
c
C   Disc sink (2D)  "Pond" command
C
 800  CALL PDIN (LSOL,DRA,IGRSIZE,DRSCR1)
      IF (loadsol.AND..NOT.LSOL) LDATACHANGE=.TRUE.
      GOTO 10
C
C   Transient Wells
C
 850  CALL TWIN
      IF (loadsol.AND..NOT.LSOL) LDATACHANGE=.TRUE.
      GOTO 10
C
C     Set the time for transient wells
C
  875 RDUM=RVAR(2)
      IF (LERROR.OR.LMISS) GOTO 10
      CALL TIME (RDUM)
      WRITE (ILUME,9090) RDUM
      CALL SETLOCALTRANSMISSIVITY
      GOTO 10
C
C   SD3D (disc sink) (3D)
C
 900  CALL DISCIN (LSOL,DRA,IGRSIZE,DRSCR1)
      IF (loadsol.AND..NOT.LSOL) LDATACHANGE=.TRUE.
      GOTO 10
C
C   Partial. Penetrat. Wells (3D Function)
C
 950  CALL W3IN (LSOL,DRA,IGRSIZE,DRSCR1)
      IF (loadsol.AND..NOT.LSOL) LDATACHANGE=.TRUE.
      GOTO 10
C
C   Potential/head/discharge/spec.discharge/velocity,interface at CZ, R3DZ
C
 960  IF (.NOT.LSOL) GOTO 1050
      CZ=CVAR(2)
      IF (LERROR) GOTO 10
      RDUM=RVAR(4)
      IF (LERROR) THEN
      R3DZ=RFBASE(CZ)+RFHGHT(CZ)    ! default z-value at saturated aquifer top
      LERROR=.FALSE.
      LMISS=.FALSE.
      ELSE
      R3DZ=RDUM
      ENDIF
      JUMPOLD=JUMP
c      CALL TWTIME (CZ) ! skip for batch version
      IF (JUMPOLD.EQ. 18) THEN
C ----------------------------Potential      
      RDUM=RFPOT(CZ)
      WRITE (ILUME,7000) CZ,R3DZ,RDUM
      IF (LUOUTFILE) WRITE (ILUOUT,7000) CZ,R3DZ,RDUM
      ENDIF
      IF (JUMPOLD.EQ.31) THEN
C-----------------------------Complex potential
      CDUM=COMEGA (CZ)+RFCONPT(CZ)
      WRITE (ILUME,7001) CZ,CDUM
      IF (LUOUTFILE) WRITE (ILUOUT,7001) CZ,CDUM
      ENDIF
      IF (JUMPOLD.EQ.19) THEN
C ----------------------------Head      
      RDUM=RFHEAD(CZ)
      WRITE (ILUME,8000) CZ,R3DZ,RDUM
      IF (LUOUTFILE) WRITE (ILUOUT,8000) CZ,R3DZ,RDUM
      ENDIF
      IF (JUMPOLD.EQ.15) THEN
C ----------------------------Discharge vector      
      CALL DISCH (CZ,RQI)
      WRITE (ILUME,8060) CZ,RQI(1),RQI(2)
      IF (LUOUTFILE) WRITE (ILUOUT,8060) CZ,RQI(1),RQI(2)
      ENDIF
      IF (JUMPOLD.EQ.21) THEN
C ----------------------------Specific discharge vector      
      CALL SDISCH (CZ,RQI)
      WRITE (ILUME,8070) CZ,R3DZ,RQI
      IF (LUOUTFILE) WRITE (ILUOUT,8070) CZ,R3DZ,RQI
      ENDIF
      IF (JUMPOLD.EQ.22) THEN
C ----------------------------Velocity vector      
      CALL VELOC(CZ,RQI)
      WRITE (ILUME,9050) CZ,R3DZ,RQI
      IF (LUOUTFILE) WRITE (ILUOUT,9050) CZ,R3DZ,RQI
      ENDIF
      IF (JUMPOLD.EQ.42) THEN
C ---------------------------Interface elevation
      RDUM=RFINTERFACE (CZ)
      WRITE (ILUME,9060) CZ,RDUM
      IF (LUOUTFILE) WRITE (ILUOUT,9060) CZ,RDUM
      END IF
      GOTO 10
C
C     reading the BaseFileName
C
 970  CALL GETFN(2)
      IF (LERROR.OR.LMISS) THEN
        ABASENAME='GFLOW'  ! no basefilename found
        LERROR=.FALSE.
        LMISS=.FALSE.
        WRITE (ILUER,9095)
      ELSE
        ABASENAME=TRIM(AFILE)   ! found basefilename, now store it
c        OPEN (UNIT=2,IOSTAT=IERR,FILE=ABASENAME,STATUS='SCRATCH')
c        IF (IERR.NE.0) THEN  ! filename does not work, replace
c          WRITE (ILUER,9096) ABASENAME
c          ABASENAME='GFLOW'
c        END IF
c        CLOSE (2)         ! delete the test file
      END IF
      GOTO 10
C
C   Tracing streamlines
C
 990  IF (.NOT.LSOL) GOTO 1050
      CZ=(0.0D0,0.0D0)
c      CALL TWTIME (CZ) ! skip for batch version
      write (*,1115)
 1115 format (/,' Tracing streamlines.')
      CALL STREAM (DRA,IGRSIZE,DRSCR1)
      GOTO 10
C
C   Clear page
C
 992  CONTINUE
C      CALL CLEARSCREEN            ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C   Debug
C
 993  CALL DEBUG
      GOTO 10
C
C    Write *.XTR file
C
 1000 CONTINUE
      write (*,1116)
 1116 format (/,' Creating Extract file.')
      CALL EXTRACT(LSOL)
      GOTO 10
C
C   Cursor
C
1100  continue
c      LCURS=.TRUE.                            ! NOT AVAILABLE IN BATCH MODE
c      CALL PLAPAR (LESCAP)                    ! NOT AVAILABLE IN BATCH MODE
c      IF (LESCAP) GOTO 10                     ! NOT AVAILABLE IN BATCH MODE
c      LCONT=LGRCOM                            ! NOT AVAILABLE IN BATCH MODE
c      CALL CURDAT (DRA,IGRSIZE,DRSCR1,LCONT)    ! NOT AVAILABLE IN BATCH MODE
c      LCURS=.FALSE.                           ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     map
C
 1150 continue
c      CALL MAP                 ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     groundwatershed 
C
 1160 CONTINUE
c      CALL WSIN                ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     Title
C      
 1170 CALL GETFN (2)
      IF (LMISS) THEN
        LMISS=.FALSE.
        LERROR=.FALSE.
        WRITE (ILUME,4000) ATITLE
        GOTO 10
      ENDIF
      ATITLE=AFILE
      IF (loadsol) LDATACHANGE=.TRUE.
      GOTO 10
C
C     Write new input data file
C
 1180 CALL DATIO
      LDATACHANGE=.FALSE. 
      GOTO 10
C
C     View ASCII files
C
 1190 CALL GETFN (2)
      IF (LERROR.OR.LMISS) THEN
        LMISS=.FALSE.
        LERROR=.FALSE.
        WRITE (ILUME,1390)
        CALL INLINE
        IF (ALINE(1).EQ.'Q') GOTO 10
        CALL GETFN(1)
        IF (LERROR) GOTO 10
      ENDIF
c      CALL VIEW                            ! NOT AVAILABLE IN BATCH MODE
      GOTO 10      
C
C     DOS escape
C      
 1195 WRITE (ILUME,1400)
c      CALL SYSTEM ('COMMAND')              ! NOT AVAILABLE IN BATCH MODE
      GOTO 10
C
C     Model origin (hidden command needed for GAEP and GFLOW 2000)
C      
 1197 RDUM1=RVAR(2)
      RDUM2=RVAR(3)
      IDUM3=IVAR(4)
      IF (LERROR.OR.LMISS) GOTO 10
      RMODORIGX=RDUM1
      RMODORIGY=RDUM2
      IMODORCODE=IDUM3
      GOTO 10
C
C   Error message - No solution
C 
 1050 WRITE (ILUER,6000)
      GOTO 10
C
1200  FORMAT(' Program GFLOW1 (rev. 10/21/2014)',/
     &       ' Copyright (C) H.M. Haitjema 1991 - 2014',/
     &       ' Haitjema Software',/
     &       ' 2738 Brigs Bend',/
     &       ' Bloomington, IN 47401',/
     &       ' phone: (812) 336 2464 fax: (812) 336 2508',///
     &       ' EXTENDED VERSION 1.10',///
     &       ' The graphics cursor in this program requires a MOUSE.',/
     &       ' The graphics screens require ANSI.SYS loaded.',/
     &       ' The program is designed for file IO without path names,'/
     &       ' and using fixed filename extensions. See manual for',/
     &       ' proper program setup and disk organization.',//
     &       ' Maximum number of equations = ',I5,////,
     &       ' Press <Enter> to continue.')
1300  FORMAT (/,
     &' -------------------- M A I N  P R O G R A M -----------------'
     &                                                             ,//,
     &' --INPUT/OUTPUT--      --ANALYTIC ELEMENTS--  --MISCELLANEOUS-',
     &                                                                /,
     &' SWitch (filename)     WEll                   TItle [title]',/,
     &' SAve   (filename)     PPwell                 AQuifer',/,
     &' LOad   (filename)     TWell                  MAp',/,
     &' DAta   (filename)     SInkdisc               DOs',/,
     &' VIew   (file.ext)     SD3d                   STop',/,
     &' EXtract               LInesink               COmmand summary',/,
     &'                       INhomogeneity          <F1> = help',//,
     &' --SOLUTION--          --NUMERICAL OUTPUT--   --GRAPHICAL OUTPUT'
     &,'--',/,
     &' SOlve [GRoun] [it.]   HEad      (x,y[,z])    LAyout',/,
     &' SOlve BAseflow        POtential (x,y[,z])    GRid',/,
     &' SOlve COnjunc [it.]   OMega     (x,y)        PLot',/,
     &' CHeck                 DIscharge (x,y)        CUrsor',/,
     &'                       SPecific  (x,y[,z])    TRace',/,
     &'                       VElocity  (x,y[,z])',/,' >')
C
1390  FORMAT (' Give name of ASCII file, press <Esc> to cancel.',/,' >')     
1400  FORMAT (' You escaped to DOS, type EXIT to return to GFLOW1.')     
1600  FORMAT (' Press <F7> for hardcopy or <Esc> to exit graphics.')
2000  FORMAT (' ***ILLEGAL OR MISSING PARAMETERS:',/,' ',80A1)
2001  FORMAT (' ***MISSING PARAMETER(S):',/,' ',80A1)
3000  FORMAT (' ***ILLEGAL COMMAND:',/,' ',80A1)
3050  FORMAT (' ***ERROR: no grid or incompatible grid!')
4000  FORMAT (' title=  ',A16,' to change give TITLE (title) command.'
     &,/)
4050  FORMAT (' ***ERROR: incorrect or incomplete iteration parameters',
     &          ' in SOLVE or FASTSOLVE command.')
5000  FORMAT (' You may still type SAVE, else press <Enter>!'/)
5050  format (' Transient well(s) present, no new solution created.')
6000  FORMAT (' ***ERROR: NO GROUNDWATER SOLUTION!',
     &' GIVE <SOLVE> COMMAND.'/)
6003  FORMAT (' ***ERROR: no solution, no access to check module.',/)
6005  FORMAT (' ***WARNING: # iterations raised to 2.')
6010  format (/,' Note: base jumps converted to k jumps for first ',I3,
     &    ' iterations.',/,
     &    ' Base jumps restored during 3 extra iterations.',/)
6100  FORMAT (/' Solving groundwater flow problem .......')
6110  FORMAT (' ***ERROR solve command not available in demo version',/)
6200  FORMAT (/' Solving stream flow problem .............')     
6205  FORMAT (/' Ending with stream flow calculations without',
     &' corrrecting for negative flows.',/,
     &' You may click on line-sink vertices of gray stream sections to',
     &' verify that',/,
     &' residual negative stream flows are not too large.')
6210  FORMAT (/,' Convergence reached on groundwater solution; '
     &          'iterations aborted.')
6220  FORMAT (' You may list or plot errors in boundary conditions in',
     &' the CHECK module.',/,
     &' Press <Enter> to continue.')
6300  FORMAT (' Always end solution procedure with:',
     &' SOLVE BASE 0')
6500  FORMAT (' Input data may have changed!',/,
     &' Do you want to create a new input data file?',/,
     &' Type YES or NO >')     
6600  FORMAT (' ***ERROR in opening ',A16,' iostat=',I5) 
7000  FORMAT (' Pot. at ',2(F11.1,1X),F11.3,' = ',G14.7)
7001  FORMAT (' Omega at ',2(F11.1,1X),' = ',2(G14.7,1X))
8000  FORMAT (' Head at ',2(F11.1,1X),F11.3,' = ',G14.7)
8060  FORMAT (' Location (x,y) = ',2(G14.7,1X),/,
     &' Discharge (Qx,Qy) = ',2(G14.7,1X))
8070  FORMAT (' Location (x,y,z) = ',3(G14.7,1X),/,
     &' Specific Discharge (qx,qy,qz) = ',3(G14.7,1X))
9050  FORMAT (' Location (x,y,z) = ',3(G14.7,1X),/,
     &' Velocity (Vx,Vy,Vz) = ',3(G14.7,1X))
 9060 FORMAT (' Interface elevation at ',G14.7,1x,G14.7,' = ',G14.7)
9070  FORMAT (' ***ERROR: program expired, request update!!',/,
     &' Press <Enter> to return to DOS')
9080  FORMAT (' Command is only available in "extended version".')
 9090 FORMAT (' Time for transient solutions is set to ',E14.7)
 9095 FORMAT (' Warning: No BaseFilename found, using "GFLOW" ')
 9096 FORMAT (' ERROR in SETUP: ',a8,' is an illegal DOS filename.',/,
     &        ' It has been replaced by "GFLOW"')
 9097 format (' Error in OPEN routine for coefficient matrix IO:',/
     &        ,a80,/,' Program execution aborted.')
 9098 format (' Error in OPEN routine for decomposed matrix IO:',/
     &        ,a80,/,' program execution aborted.')
 9099 format (' Error in reading from matrix file:',/
     &        ,a80,/,' Program execution aborted.')
 9190 format (' Error in allocating matrix arrays:',/,
     &        ' Probably insufficient RAM available',/,
     &        ' Cannot allocate enough memory for a contour grid',/,
     &        ' Free up allocatable RAM, or increase RAM.',/,
     &        ' Program execution aborted.')
 9191 format (' Error in allocating matrix arrays:',/,
     &        ' Probably insufficient RAM available',/,
     &' Not enough RAM for current number of equations.',/,
     &        ' Free up allocatable RAM, or increase RAM.',/,
     &        ' Program execution aborted.')
 9192 format (' Error in allocating grid array:',/,
     &        ' Probably insufficient RAM available',/,
     &        ' Free up allocatable RAM, or increase RAM.',/,
     &        ' Program execution aborted.')
C
      END
C
C -----------------------------------------------------------------------------
C
      SUBROUTINE INITS
C
C -----------------------------------------------------------------------------
C
      IMPLICIT COMPLEX (C), LOGICAL (L)
      INCLUDE 'MAIN.INC'
      INCLUDE 'MATCH.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'COM3D.INC'
      LINALREADY=.FALSE.
      RPI=3.141592654
      RPI2=2.0*RPI
      RO2PI=1.0/RPI2
      ROPI=1.0/RPI
      CPI=CMPLX(0.0,RPI)
      C2PI=CMPLX(0.0,RPI2)
      CI=(0.0,1.0)
      CIM=(0.0,-1.0)
      LSOL=.FALSE.
      LDISPL=.FALSE.
      NSOL=0
      nouterloop=0
      IMODORIGX=0
      IMODORIGY=0
      IMODORCODE=0
      lGaussSeidel=.false.
      lDirectFromDisk=.true.
      lErrorReport=.true.
      abasename='GFLOW'
      loadsol=.false.
      aDateTime='Not_Set'
      RETURN
      END



