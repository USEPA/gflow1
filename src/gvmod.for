C     Last change:  HMH  18 Jul 2002    4:36 pm
c       This file contains the following routines and functions:
c
c       GVDAT       initializes data in given module common blocks
c       GVIN        input of aquifer data, etc. from .dat file
c
c -------------------
c gvmat.for contains:
c -------------------
c
c       GVCZI     generates reference control point
c       GVMAT     generates matrix coefficients for reference point
c       GVKNO     generates known vector for reference point
C       GVSUB     substitutes solution vector (element) into integration constant
c  
c -------------------
c gvfun.for contains:
c -------------------
c
c       CFUFOM         calculates the complex potential due to uniform flow
c       GVQI           calculates the discharge vector due to uniform flow
c       RFRAINPT       real discharge potential for a recharge function (not in use)
c       RFCONPT        adds integration constant to the real discharge potential
C
c ---------------------
c gvnflow.for contains:
c ---------------------
c
c       RFNFGV   calculates flow across line from CZ1 to CZ2 due to uniform flow
C
c ------------------
c gvio.for contains:
c ------------------
c
c     GVIO            reads or writes contents of GVCOM.INC to .sol file
c     GVDATIO         writes all given (aquifer) data to a ".dat" file
c     GVOUT           writes contents of GVCOM formatted ti ILUOUT for debugging purposes
C
c ---------------------
c gvcheck.for contains:
c ---------------------
c
c     GVERROR   reports the error in head at the reference point
C
c --------------------
c gvserv.for contains:
c --------------------
c
c      RFTOP      returns elevation of the confining layer at CZ
c      RFHEDP     returns the HEAD belonging to the potential RPOT
c      GVPAR      returns given commonblock data as arguments
c      RFBASE     returns the aquifer base at CZ
c      RFPOTH     returns the POTENTIAL belonging to the head RHEDIN
c      RFPERM     returns the hydraulic conductivity at CZ
c      RFPOR      returns the porosity at CZ
c      RFHGHT     returns the SATURATED AQUIFER THICKNESS at CZ
c      RGVREFDIST returns the distance from CZ to the reference point
C
c ---------------------
c gvextra.for contains:
c ---------------------
c
c     GVCHECK
C
C
C --------------------------------------------------------------------
C
	 BLOCK DATA GVDAT
C
C --------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'GVCOM.INC'
      INCLUDE 'LUSYS.INC'
      DATA RK,RH,RHEAD0,RBASE,RPOTC,CREFZ/3*1.0D0,2*0.0D0,(0.0D0,0.0D0)/
      DATA RCON,RPI,RQ0,CUNALP /0.5D0,3.141592653589793D0,
     &                          0.0D0,(0.0D0,0.0D0)/
      DATA RAINX,RAINY,RAIN,LRAINX,LRAINY,LRAINR /3*0.0D0,3*.FALSE./
      DATA RSPECIFICGRAVITYFRESH,RSPECIFICGRAVITYSALT /2*0.0D0/
      DATA RSEALEVEL /0.0D0/
      DATA LINTERFACE /.FALSE./
      DATA RPOR /0.2D0/
      END
C
C --------------------------------------------------------------------
C
      SUBROUTINE GVIN (LSOL)
C
C --------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) JUMP
      LOGICAL LSOL,LBAD
      REAL(8) RDUM,RVAR,RDUMX,RDUMY,SQROOT,RALPH,
     &        RDUM1,RDUM2
      COMPLEX(8) CDUM,CVAR
      CHARACTER(1) AWORD(51)
      INCLUDE 'GVCOM.INC'
      INCLUDE 'LUSYS.INC'
      INCLUDE 'MATCH.INC'
      DATA AWORD /'P','E','R','M',' ',
     .            'T','H','I','C',' ',
     .            'R','E','F','E',' ',
     .            'Q','U','I','T',' ',
     .            'U','N','I','F',' ',
     .            'R','A','I','N',' ',
     .            'P','O','R','O',' ',
     .            'B','O','T','T',' ',
     .            'B','A','S','E',' ',
     .            'I','N','T','E',' ',
     .            ATERM/
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
  10  IF (LERROR.OR.LMISS) WRITE (ILUER,2000) ALINE2
      LERROR=.FALSE.
      LMISS=.FALSE.
      if (lucon)  then
      WRITE (ILUME,1000) RK
      WRITE (ILUME,1010) RPOR
      WRITE (ILUME,1001) RH
      WRITE (ILUME,1004) RBASE      
      WRITE (ILUME,1002) CREFZ,RHEAD0
      WRITE (ILUME,1003) REAL(RQ0*EXP(CUNALP)),AIMAG(RQ0*EXP(CUNALP))
      IF (LRAINX) WRITE (ILUME,1005) RAINX,RAIN
      IF (LRAINY) WRITE (ILUME,1006) RAINY,RAIN
      IF (LRAINR) WRITE (ILUME,1007) RAINX,RAINY,RAIN
      WRITE (ILUME,1009)
      endif
      CALL INLINE
  11  CALL MATCH (AWORD,1,JUMP,LBAD)
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      IF (.NOT.LBAD) GOTO 15
c      CALL CLEARSCREEN                 ! NOT AVAILABLE IN BATCH MODE
      GOTO (10,13),JUMP
  13  WRITE (ILUER,3000) ALINE2
      LERROR=.FALSE.
      GOTO 10
  15  GOTO (100,200,300,400,500,600,700,800,800,900),JUMP
C
C     HYDRAULIC CONDUCTIVITY
C
 100  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RK=RDUM
c      CALL DBKFAC       ! calculate new coefficients for inhomogeneities (not needed in batch mode)
      LSOL=.FALSE.
      GOTO 10
C
C     AQUIFER HEIGHT
C
 200  RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RH=RDUM
      LSOL=.FALSE.
      GOTO 10
C
C     REFERENCE POINT
C
 300  CDUM=CVAR(2)
      RDUM=RVAR(4)
      IF (LERROR) GOTO 10
      CREFZ=CDUM
      RHEAD0=RDUM
      LSOL=.FALSE.
      GOTO 10
C
C     RETURN
C
 400  CONTINUE
      RETURN
C
C     UNIFORM FLOW
C
 500  RDUMX=RVAR(2)
      RDUMY=RVAR(3)
      IF (LERROR) GOTO 10
      RQ0=SQROOT(RDUMX*RDUMX+RDUMY*RDUMY)
      if (rq0.lt.1.0e-25) then
      ralph=0.0
      else
      RALPH=ATAN2(RDUMY,RDUMX)
      ENDif
      CUNALP=CMPLX(0.0,RALPH)
      LSOL=.FALSE.
      GOTO 10
C
C     RAIN
C
 600  GOTO 10        !  COMMAND IS BLOCKED OFF, USE RECHARGE INHOMOGENEITY INSTEAD
C      LRAINX=.FALSE.
C      LRAINY=.FALSE.
C      LRAINR=.FALSE.
C      IF (LERROR) WRITE (ILUER,2000)
C      LERROR=.FALSE.
C      LMISS=.FALSE.
C      if (lucon) WRITE (ILUME,4000)
C      CALL INLINE
C      IF (ALINE(1).EQ.'X') LRAINX=.TRUE.
C      IF (ALINE(1).EQ.'Y') LRAINY=.TRUE.
C      IF (ALINE(1).EQ.'R') LRAINR=.TRUE.
C      IF (ALINE(1).EQ.'C') GOTO 10
C      RDUM=RVAR(2)
C      RDUM1=RVAR(3)
C      IF (LERROR) GOTO 600
C      IF (LRAINX) THEN
C      RAINX=RDUM
C      RAIN=RDUM1
C      LSOL=.FALSE.
C      GOTO 10
C      ENDIF
C      IF (LRAINY) THEN
C      RAINY=RDUM
C      RAIN=RDUM1
C      LSOL=.FALSE.
C      GOTO 10
C      ENDIF
C      IF (LRAINR) THEN
C      RDUM2=RVAR(4)
C      IF (LERROR) GOTO 600
C      RAINX=RDUM
C      RAINY=RDUM1
C      RAIN=RDUM2
C      LSOL=.FALSE.
C      GOTO 10
C      ENDIF
C      if (lucon) WRITE (ILUME,3000)
C      GOTO 600
C
C     POROSITY
C
  700 RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      IF (RDUM.LE.0.0.OR.RDUM.GT.1.0) THEN
      WRITE (ILUER,5000)
      GOTO 10
      ENDIF
      RPOR=RDUM
      GOTO 10
C
C     aquifer base
C
  800 RDUM=RVAR(2)
      IF (LERROR) GOTO 10
      RBASE=RDUM
      LSOL=.FALSE.
      GOTO 10
c
c     interface
c
  900 RDUM=RVAR(2)
      RDUM1=RVAR(3)
      RDUM2=RVAR(4)
      IF (LERROR) GOTO 10
      LERROR=.TRUE.
      IF (RDUM.LE.0.0D0) THEN ! illegal fresh water specific gravity
      WRITE (ILUER,6001) RDUM
      GOTO 10
      END IF
      IF (RDUM1.LE.0.0D0) THEN ! illegal salt water specific gravity
      WRITE (ILUER,6002) RDUM2
      GOTO 10
      END IF
      LERROR=.FALSE.
      RSPECIFICGRAVITYFRESH=RDUM
      RSPECIFICGRAVITYSALT=RDUM1
      RSEALEVEL=RDUM2
      RGVFAC1=RSPECIFICGRAVITYSALT/RSPECIFICGRAVITYFRESH
      RGVFAC2=RSPECIFICGRAVITYSALT/
     &       (RSPECIFICGRAVITYSALT-RSPECIFICGRAVITYFRESH)
      LINTERFACE=.TRUE. ! activate interface flow
      GOTO 10
C
 1000 FORMAT ('       ---------- AQUIFER module ----------- ',/,
     &        ' PERMEABILITY ',G14.7)
 1001 FORMAT (' THICKNESS ',G14.7)
 1002 FORMAT (' REFERENCE x,y ',2F11.1,' head ',G14.7)
 1003 FORMAT (' UNIFLOW Qx ',G14.7,' Qy ',G14.7)
 1004 FORMAT (' BASE ',G11.4)
 1005 FORMAT (' Rain flows in x-dir. X0=',G14.7,' N=',G14.7)
 1006 FORMAT (' Rain flows in y-dir. Y0=',G14.7,' N=',G14.7)
 1007 FORMAT (' Rain flows radially X0=',G14.7,' Y0=',G14.7,' N=',G14.7)
 1009 FORMAT (' <Esc> or QUIT',/,' >')
 1010 FORMAT (' POROSITY ',G14.7)
 2000 FORMAT(' ***ILLEGAL or MISSING PARAMETER(S) in aquifer module:',/,
     &       ' ',80A1)
 3000 FORMAT (' ***ILLEGAL COMMAND in aquifer module:',/,
     &        ' ',80A1)
 4000 FORMAT (' Enter if rain should flow off:'/
     &        ' in x-dir.  <X>  (X0)  (N)'/
     &        ' in y-dir.  <Y>  (Y0)  (N)'/
     &        ' radially <R>  (X0)  (Y0)  (N)'/
     &        ' Type CLEAR to eliminate rain'/)
 5000 FORMAT (' ***ERROR: porosity must be between 0 and 1, try again.')
 6001 FORMAT (' ***ERROR: illegal fresh water specific gravity: ',D14.7)
 6002 FORMAT (' ***ERROR: illegal salt water specific gravity: ',D14.7)
C
      END
