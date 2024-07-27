      INTEGER(4) NPZMX,NPIEZ
      REAL RHPIEZ,RMINPZ,RMAXPZ,RMINPZ0,RMAXPZ0,RMAXABS,RHMX
      COMPLEX CZPIEZ
      CHARACTER(16) APZLAB
      PARAMETER (NPZMX=1000)
      COMMON /CHCK/ NPIEZ,CZPIEZ(NPZMX),RHPIEZ(NPZMX),APZLAB(NPZMX),
     &              RMINPZ,RMAXPZ,RMINPZ0,RMAXPZ0,RMAXABS,RHMX
C
C     NPZMX       maximum number of piezometers allowed
C     NPIEZ       current number of piezometers stored
C     CZPIEZ      location (x,y) of piezometers
C     RHPIEZ      piezometric heads in the piezometers
C     APZLAB      labels of the piezometers
C     RMINPZ      minimum value of modeled - observed head to be plotted
C     RMAXPZ      maximum value of modeled - observed head to be plotted
C     RMINPZ0     minimum value of modeled - observed head found
C     RMAXPZ0     maximum value of modeled - observed head found
C     RMAXABS     absolute value of maximum difference
C     RHMX        maximum marker size
C
      