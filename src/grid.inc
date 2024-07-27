      INTEGER(4) NGRWMX,NX,NY,IX1,IX2,IY1,IY2,NDOTX,NDOTY,NXMAX,NYMAX,
     &         ILAYER,ISFCNT,ISFREC,NGRWIN,NGRDUM,NLEVEL,NLEVELS,NTEXTS,
     &         IGRIDCODE
      LOGICAL LPHI,LPSI,LHEAD,LRESP,LSCRN,LPRINT,LPLOT,LHPLOT,LGRID,
     &        LCURS,LSURF,LGRCOM,LDISCH,LFLOWNET,LDOTMAP,LSUBTRACT,
     &        LSLIDES,LEAKANCE,LEAKAGE,LBASE,LTOP,LHEIGHT,LCONDUCTIVITY,
     &        LPOROSITY,LRECHARGE,LINTELV
      REAL(8) RX1,RY1,RX2,RY2,RMAX,RMIN,RSCALE,RXG,RYG,RSTEP,RVFAC,
     &  RXSCALE,RYSCALE ,RXMN,RYMN,RXMX,RYMX,RDS0,RWX1,RWY1,RWX2,RWY2,
     &     RMINLEVEL,RINCREMENT,RMINLVLPSI,NLEVELPSI

      COMPLEX(8) CZWIN,CZBL0
      CHARACTER(16) ATITLE
      PARAMETER (NGRWMX=10)      
      COMMON /GRID/ RX1,RY1,RX2,RY2,RMAX,RMIN,
     &              NX,NY,RSCALE,RXG(1000),RYG(1000),IX1,IY1,
     &              IX2,IY2,RSTEP,LPHI,LPSI,LHEAD,NDOTX,NDOTY,
     &              NXMAX,NYMAX,LRESP(25),ILAYER,LSCRN,LPRINT,
     &              LPLOT,LHPLOT,LGRID,RVFAC,RXSCALE,RYSCALE,
     &              RXMN,RYMN,RXMX,RYMX,LCURS,RDS0,
     &              LSURF,ISFCNT,ISFREC,CZBL0,
     &              RWX1,RWY1,RWX2,RWY2,LGRCOM,LDISCH
      COMMON /GRIDC/ ATITLE
      COMMON /GRID2/ NGRWIN,CZWIN(4,NGRWMX)
      COMMON /GRID3/ LFLOWNET,LDOTMAP,LSUBTRACT,LSLIDES,NLEVELS,NTEXTS,
     &               LEAKAGE,LEAKANCE,NGRDUM(3),LBASE,LTOP,LHEIGHT,
     &               LCONDUCTIVITY,LPOROSITY,LRECHARGE,LINTELV
      COMMON /GRID4/ RMINLEVEL,RINCREMENT,NLEVEL
      COMMON /GRID5/ RMINLVLPSI,NLEVELPSI
      COMMON /GRID6/ IGRIDCODE
C     
C RX1,RY1,RX2,RY2 lower left and upper right corner of window
C RMAX, RMIN      maximum and minimum values in the grid
C NX,NY           number of grid points in x- and y-direction
C RSCALE          scaling factor to translate user coordinates into
C                 graphics display area coordinates
C RXG,RYG         arrays for x- and y-coordinates in the grid
C IX1,IY1,IX2,IY2 lower left and upper right corner of graphics display area
C RSTEP           current stepsize (delta s) along streamlines (for tracing)
C ,LHEAD potential, streamfunction, head to be contoured
C NDOTX, NDOTY    maximum graphics display units available in x_ and y_direction
C NXMAX,NYMAX     maximum grid points allowable in x- and y_direction
C LRESP(i)        .TRUE. when i th command in PLOPAR has been used
C IlAYER          layer number (currently not in use)
C LSCRN,LPRINT   .TRUE. when screen or printer selected
C LPLOT,LHPLOT   .TRUE. when (GSS) plotter or HIPLOT selected
C LGRID          .TRUE. when a grid is available for contouring
C RVFAC          vertical scaling factor 
C RXSCALE        horizontal scale from user coordinates to graphics display area
C RYSCALE        vertical scale from user coordinates to graphics display area
C RXMN,RXMX      minimum and maximum x-values of analytic element coordinates
C RYMN,RYMX      minimum and maximum y-values of analytic element coordinates
C LCURS          cursor selected
C RDS0           default stepsize
C LSURF          .TRUE. when surfer file writing is selected
C ISFCNT,ISFREC  current record to be written, record to write number of points
C                in, used to create surfer boundary files
C CZBL0          previous point in boundary line file writing
C RWX1,RWY1      lower left corner of grid window (stored with grid)
C RWX2,RWY2      upper right corner of grid window (stored with grid)
C LGRCOM        .TRUE. when grid window and current window correspond
C LDISCH        absolute value of discharge vector to be contoured
C ATITLE        title of problem being solved
C NGRWIN        number of windows stored for ``select option" in window
C CZWIN(4,I)    corner coordinates of I th window
C LFLOWNET      flownet to be contoured
C LDOTMAP       dots to be plotted on the screen during grid generation
C NGRDUM        dummy array (reserved for future use)
C NLEVELS       number of contour levels stored in level buffer in GRAPHEDIT
C NTEXTS        number of text items stored in text buffer in GRAPHEDIT
C IGRIDCODE     code stored in ILAYER to identify solution
C LEAKANCE      apparent leakance factor to be contoured
C LEAKAGE       actual leakage to be contoured
C LBASE         aquifer base elevation to be contoured
C LTOP          aquifer top elevation to be contoured
C LHEIGHT       saturated aquifer height to be contoured
C LCONDUCTIVITY hydraulic conductivity to be contoured
C LPOROSITY     porosity to be contoured
C LRECHARGE     recharge to be contoured
C LINTELV       interface elevation to be contoured

