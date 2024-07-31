      INTEGER(4) NPD,NPDHD,NPDRC,IPDHD,IPDRC,NPDMX
      LOGICAL LPDGIV,LPDHD,LPDRC,LPDRCH,LPDPER
      REAL(8) RPI2,RPDC0,RPDD0,RPDS,RPDH,RPDR,RPDRES,RPDEP,RPDSB
      COMPLEX(8) CPDZ
      CHARACTER(16) APDLAB
c
      PARAMETER (NPDMX=700)            ! maximum number of 2D sink discs
c
      COMMON/POND/NPD,NPDHD,NPDRC,
     &            LPDHD,LPDRC,LPDRCH,LPDPER,
     &            RPI2,RPDC0,RPDD0,
     &            IPDHD(NPDMX),IPDRC(NPDMX),
     &            LPDGIV(NPDMX),APDLAB(NPDMX),
     &            CPDZ(NPDMX),RPDS(NPDMX),RPDH(NPDMX),RPDR(NPDMX),
     &            RPDRES(NPDMX),RPDEP(NPDMX),RPDSB(NPDMX)
C
C     NPDMX    maximum number of ponds
C     NPD      total number of ponds
C     NPDHD    number of head specified ponds
C     NPDRC    number of recharge specified ponds
C     LPDHD    true for plotting head specified ponds
C     LPDRC    true for plotting recharge specified ponds
C     LPDRCH   true for highlighting recharging head specified ponds
C     LPDPER   true for highlighting percolating head specified ponds
C     RPI2     6.2831....
C     RPDC0    default resistance of layer underneath head specified ponds
C     RPDD0    default depth of bottom of resistance layer below water table
C     IPDHD    pointer to head specified ponds
C     IPDRC    pointer to recharge specified ponds
C     LPDGIV   true when a head specified pond percolates (becomes "given" and
C              is eliminated from the matrix solution procedure)
C     APDLAB   label of pond
C     CPDZ     center of pond
C     RPDS     net strength (exfiltration rate) of pond (top+bottom)
C     RPDH     specified head at pond center
C     RPDR     radius of pond
C     RPDRES   resistance of layer underneath pond
C     RPDEP    depth of bottom of resistance layer below water table in pond
C     RPDSB    exfiltration rate at aquifer bottom (only non zero for selected
C              discharge specified ponds

