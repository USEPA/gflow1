      LOGICAL LRAINX,LRAINY,LRAINR,LINTERFACE
      REAL(8) RK,RH,RHEAD0,RBASE,RPOTC,RCON,RPI,RQ0,
     &        RAINX,RAINY,RAIN,RPOR,RSPECIFICGRAVITYFRESH,
     &        RSPECIFICGRAVITYSALT,RSEALEVEL,RGVFAC1,RGVFAC2
      COMPLEX(8) CREFZ,CUNALP
      COMMON/GIVEN/RK,RH,RHEAD0,RBASE,RPOTC,CREFZ,RCON,RPI,RQ0,CUNALP,
     &             RAINX,RAINY,RAIN,LRAINX,LRAINY,LRAINR,RPOR,
     &             RSPECIFICGRAVITYFRESH,RSPECIFICGRAVITYSALT,RSEALEVEL,
     &             LINTERFACE,RGVFAC1,RGVFAC2
c
C     RK         regional hydraulic conductivity (outside all inhomogeneity domains)
C     RH         regional aquifer thickness (outside all inhomogeneity domains)
C     RHEAD0     head at the reference point
C     RBASE      regional aquifer base elevation (outside all inhomogeneity domains)
C     RPOTC      potential at the reference point
C     CREFZ      location of the reference point
C     RCON       integration constant for the complex potential function
C     RPI        the number pi
C     RQ0        ambient uniform flow rate
C     CUNALP     i*alpha, where alpha is the direction of the uniform flow measured with respect to the x-axis
C     RAINX      areal infiltration due to rainfall flowing in x-direction (not in use)
C     RAINY      areal infiltration due to rainfall flowing in y-direction (not in use)
C     RAIN       areal infiltration due to rainfall flowing in radial direction (not in use)
C     LRAINX     true if rainfall flows in x-direction (not in use)
C     LRAINY     true if rainfall flows in y-direction (not in use)
C     LRAINR     true if rainfall flows in radial direction (not in use)
C     RPOR       regional (effective) porosity (outside all inhomogeneity domains)
C     RSPECIFICGRAVITYFRESH  specific gravity of fresh water
C     RSPECIFICGRAVITYSALT   specific gravity of salt water
C     RSEALEVEL              average local sea level (same datum as heads, e.g. msl)
C     LINTERFACE             true if interface flow to be included
C     RGVFAC1    spec. grav. salt water/spec. grav. fresh water
C     RGVFAC2    spec. grav. salt water/(spec. grav. salt water - spec. grav. fresh water)
C
     