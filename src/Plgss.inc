      INTEGER(4) IIX,MAXCOOR,MAXPOINT
      PARAMETER (MAXCOOR=2000,MAXPOINT=MAXCOOR/2)      
      COMMON /LIBUF/ IIX(0:MAXCOOR)
C
C     MAXCOOR     maximum number of coordinates
C     MAXPOINT    maximum number of points
C     IIX         coordinate array
C      
