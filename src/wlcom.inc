      INTEGER(4) NWL,NWLMX,nwldischarges
      REAL(8) RWLQ,RWLH,RWLR,RWLD0,RWLERRH,RPI2,RWLRD0,DWLPOT
      COMPLEX(8) CWLZ
      LOGICAL LWLH
      CHARACTER(16) AWLAB
c
      PARAMETER (NWLMX=2000)         ! maximum number of wells
c
      COMMON /WELL/ NWL,nwldischarges,RPI2,RWLRD0,
     &              CWLZ(NWLMX),RWLQ(NWLMX),RWLH(NWLMX),RWLR(NWLMX),
     &              RWLD0(3,NWLMX),RWLERRH(NWLMX),AWLAB(NWLMX),
     &              LWLH(NWLMX),DWLPOT(NWLMX)
C
C     NWLMX       total number of wells available
C     NWL         total number of wells entered
C     nwldischarges number of unkown wells discharges in the matrix
C     RWLRD0      default well radius
C     CWLZ        center of the well
C     RWLQ        discharge of the well (pumping rate)
C     RWLH        head in the well (0 when discharge specified)
C     RWLR        radius of the well
C     RWLD0       vector (3D) from well center to point CZ
C                 (used in WLNEAR for determining end of streamline)
C     RWLERRH     relative error in head at the well
C     AWLAB       label of the well
C     LWLH        .TRUE. when head specified well
C     DWLPOT      discharge potentials (array) at control point
C
       
