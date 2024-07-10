c
c WARNING: R3DCPZ has a fixed dimension which should not be less
c          than the number of equations in the matrix !!!!!
c
      INTEGER(4) IMODORCODE
      LOGICAL L3DREV,L3DEND,L3DPL,L3DHOR,L3DVER
      REAL(8) R3DH,R3DK,R3DPOR,R3DZ,R3DSTO,R3DCPZ,R3DX1,R3DY1,R3DZ1,
     &     R3DX2,R3DY2,R3DZ2,R3DZA,RETARDATION,RHALFLIFE,
     &     RMODORIGX,RMODORIGY
      COMMON/COM3D/ R3DH,R3DK,R3DPOR,R3DSTO,R3DZ,R3DCPZ(2000),
     &              RETARDATION,RHALFLIFE,
     &              IMODORCODE,RMODORIGX,RMODORIGY,
     &              L3DREV,L3DEND,R3DX1,R3DY1,R3DZ1,R3DX2,R3DY2,R3DZ2,
     &              L3DPL,L3DHOR,L3DVER
      COMMON/TEM3D/ R3DZA(300)
C
C     R3DH          aquifer thickness in 3D zone  (set in 3D modules, but not used)
C     R3DK          aquifer hydraulic conductivity in 3D zone (set in 3D modules, but not used)
C     R3DPOR        aquifer porosity in 3D zone  (set in 3D modules, but not used)
C     R3DZ          elevation of point at which calculations are carried out
C     R3DSTO        aquifer storage coefficient in 3D zone (not set or used)
C     R3DCPZ        array of elevations for 3D grids (for 3D contour plots)
C     R3DX1         lower left x-coordinate for 3D grid
C     R3DY1         lower left y-coordinate for 3D grid
C     R3DZ1         lower left z-coordinate for 3D grid
C     R3DX2         upper right x-coordinate for 3D grid
C     R3DY2         upper right y-coordinate for 3D grid
C     R3DZ2         upper right z-coordinate for 3D grid
C     R3DZA         elevation belonging to CZA in matrix generation routines
C     RETARDATION   retardation coefficient for equilibrium adsorption
C     RHALFLIFE     half-life for first order decay
C     RMODORIGX     x-coordinate in world coordinates for model origin (stored for GUI, not used in solver)
C     RMODORIGY     y-coordinate in world coordinates for model origin (stored for GUI, not used in solver)
C     L3DREV        true if pathline tracing occurs backward in time
C     L3DEND        true is pathline has arrived at a boundary (analytic element, aquifer boundary or window boundary)
C     L3DPL         true if in 3D plotting mode
C     L3DHOR        true if a grid is generated on a horizontal cross-section through the aquifer
C     L3DVER        true if a grid is generated on a vertical cross-section through the aquifer
C
      