      INTEGER(4) NTW,NTWMX
      REAL(8) RTWCT,RTWRD0,RTWST0,RPI4,RTWR0,RTWST,RTWQ,RTWH0,RTWSTO,
     &        RTWTRANS
      COMPLEX(8) CTWZ
      CHARACTER(16) ATWLAB
c
      PARAMETER (NTWMX=500)   ! number of transient wells (including images)
c
      COMMON /TWELL/ NTW,RTWCT,RTWRD0,RTWST0,RPI4,
     &               CTWZ(NTWMX),RTWR0(NTWMX),RTWST(NTWMX),RTWQ(NTWMX),
     &               RTWH0(NTWMX),RTWSTO(NTWMX),RTWTRANS(NTWMX),
     &               ATWLAB(NTWMX)
C
C     NTWMX       maximum number of transient wells
C     NTW         number of transient wells
C     RTWCT       current time
C     RTWRD0      default well radius
C     CTWZ        well coordinates
C     RTWR0       well radius
C     RTWST       starting time of well discharge (pumping)
C     RTWQ        well discharge
C     RTWTRANS    local transmissivity using local k and local saturated aquifer thickness
C     RTWH0       local saturated aquifer thickness at the onset of pumping
C     RTWSTO      (local) storage coefficient
C     ATWLAB      label for well of maximum 16 characters

     
