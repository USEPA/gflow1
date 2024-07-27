      INTEGER(4) ILPNT,NERR,NPBAD,NCHAR,ilinelength
      LOGICAL LERROR,LMISS
      CHARACTER(1 )ALINE,ALINE2,ATERM
      CHARACTER(16) AFILE,AFILE2
C
      PARAMETER (ATERM='$',ilinelength=132)
C
      COMMON /MATCM/ ILPNT(40),NERR,LERROR,LMISS,NPBAD,NCHAR
      COMMON /MATCM2/ AFILE
      COMMON /MATCM1/ ALINE(ilinelength)
      COMMON /MATCM3/ ALINE2(ilinelength)
      COMMON /MATCH4/ AFILE2
C
C     ILPNT          pointer array for command words and parameters
C     NERR
C     NPBAD
C     NCHAR
C     ILINELENGTH    number of characters read from the DOS buffer (keyboard) or file
C     LERROR         true if illegal parameter
C     LMISS          true if parameter not found
C     ALINE          buffer in which line is read
C     ALINE2         duplicate buffer for line
C     ATERM          command string termination character
C     AFILE          buffer in which filename is read
C     AFILE2         duplicate filename buffer
C
C

