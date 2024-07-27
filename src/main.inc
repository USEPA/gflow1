      INTEGER(4) NSOL,nouterloop,nlsiterarray,nlsrelaxarray,
     &           niterarray,ilstablelength,
     &           Number_equations,Number_strengths
      LOGICAL LSOL,LDISPL,LINALREADY,lErrorReport,lGaussSeidel,
     &              lDirectFromDisk,loadsol,lWrongFile
      REAL(8) RPI2,RO2PI,ROPI,RERMAX,relaxarray
      COMPLEX(8) CPI,C2PI,CI,CIM
      CHARACTER(80) AHELPPATH,ACURRENTPATH
      CHARACTER(8)  abasename
      character(16) aDateTime
      CHARACTER(2)  AHELPDRIVE,ACURRENTDRIVE
      parameter (nlsiterarray=100,nlsrelaxarray=100)
      COMMON/MAINC/ RPI2,RO2PI,ROPI,CPI,C2PI,CI,CIM,RERMAX,LDISPL,
     &              LSOL,NSOL,nouterloop,lErrorReport,lGaussSeidel,
     &              lDirectFromDisk,abasename,aDateTime,loadsol,
     &              lWrongFile,niterarray(nlsiterarray),
     &              relaxarray(nlsrelaxarray),ilstablelength,
     &              Number_equations,Number_strengths
      COMMON/MAIND/ LINALREADY
      COMMON/MAINE/ AHELPPATH,ACURRENTPATH,AHELPDRIVE,ACURRENTDRIVE
C
C     RPI2              2 pi
C     RO2PI             1/(2pi)
C     ROPI              1/pi
C     CPI               (0.0,pi)
C     C2PI              (0.0,2pi)
C     CI                (0.0,1.0)
C     CIM               (0.0,-1.0)
C     LSOL              true if a solution has been generated
C     NSOL              number of solutions (iterations)
C     RERMAX            maximum error in all boundary conditions
C     LDISPL            true when errors are to be plotted, else they will be listed
C     LINALREADY        true when data has been entered
C     AHELPPATH         path to the help files
C     ACURRENTPATH      path to the current directory
C     AHELPDRIVE        drive where help files reside
C     ACURRENTDRIVE     current drive
c     abasename         basename defined by the GUI
c     aDateTime         date and time stamp used in saving .sol, .mtr and .dec files
c     loadsol           true if a solution has been loaded from disk
c     lWrongFile        true if DateTimeStamp of files do not correspond, create new matrix
c     lErrorReport      true is errors in boundary conditions are to be reported during iterations
c     lGaussSeidel      true if Gauss Seidel type iterative solution procedure is to be used
c     lDirectFromDisk   true if matrix and decomposed matrix are to be stored and Sherman-Morrison to be used
c     niterarray        array of iterations numbers for change in relaxation for streamflow solution
c     relaxarray        array of relaxation factors associated with niterarray
c     ilstablelength    length of niterarray and relaxarray
c     Number_equations  number of equations as determined by "matsize" routines called in main program.
c     Number_strengths  number of strength parameters (degrees of freedom) as determined by "matsize"
c                       routines called in main program.
C
c
      