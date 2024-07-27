C     Last change:  HMH   8 Aug 2005    6:47 pm
c     This file contains the following routines and functions:
c
c     GVIO            reads or writes contents of GVCOM.INC to .sol file
c     GVDATIO         writes all given (aquifer) data to a ".dat" file
c     GVOUT           writes contents of GVCOM formatted ti ILUOUT for debugging purposes
C
C-------------------------------------------------------------------------------------------------
C
      SUBROUTINE GVIO (ICODE,ILU,RVERSION,ierr)
C
C-------------------------------------------------------------------------------------------------
C
C     Routine reads or writes contents of GVCOM.INC common blocks to
C     an external file.
C     ICODE=41 write
C     ICODE=73 read
C
      IMPLICIT NONE
      INTEGER(4) ICODE,ILU,IERR
      REAL(8) RVERSION
      INCLUDE 'GVCOM.INC'
C
      if (ierr.ne.0) return
      CALL BUFIOR (RK,5,ILU,ICODE,IERR)
      CALL BUFIOC (CREFZ,1,ILU,ICODE,IERR)
      CALL BUFIOR (RCON,3,ILU,ICODE,IERR)
      CALL BUFIOC (CUNALP,1,ILU,ICODE,IERR)
      CALL BUFIOR (RAINX,3,ILU,ICODE,IERR)
      CALL BUFIOL (LRAINX,3,ILU,ICODE,IERR)
      CALL BUFIOR (RPOR,1,ILU,ICODE,IERR)
      IF (RVERSION.LE.6.0) RETURN
      CALL BUFIOR (RSPECIFICGRAVITYFRESH,3,ILU,ICODE,IERR)
      CALL BUFIOL (LINTERFACE,1,ILU,ICODE,IERR)
      CALL BUFIOR (RGVFAC1,2,ILU,ICODE,IERR)
C
      RETURN
      END
C
C-------------------------------------------------------------------------------------------------
C
      SUBROUTINE GVDATIO (ILU)
C
C-------------------------------------------------------------------------------------------------
C
C     Routine writes all given (aquifer) data to a ".dat" file.
C     Routine is called in DATIO.
C
      IMPLICIT NONE
      INTEGER(4) ILU
      COMPLEX CFLOW,CEXP
      INCLUDE 'GVCOM.INC'
      INCLUDE 'LUSYS.INC'
C
      WRITE (ILU,1000)
      WRITE (ILU,2000) RK
      WRITE (ILU,3000) RH
      WRITE (ILU,4000) RBASE
      WRITE (ILU,5000) RPOR
      IF (LINTERFACE) THEN
      WRITE (ILU,5050) RSPECIFICGRAVITYFRESH,RSPECIFICGRAVITYSALT,
     &                 RSEALEVEL
      END IF
      WRITE (ILU,6000) CREFZ,RHEAD0
      CFLOW=RQ0*EXP(CUNALP)
      WRITE (ILU,7000) CFLOW
      WRITE (ILU,8000)
      RETURN
C     
 1000 FORMAT (' aquifer')
 2000 FORMAT (' permeability ',G14.7)
 3000 FORMAT (' thickness    ',G14.7)
 4000 FORMAT (' base         ',G14.7)
 5000 FORMAT (' porosity     ',G14.7)
 5050 FORMAT (' interface    ',3G14.7)
 6000 FORMAT (' reference    ',2G14.7,2X,G14.7)
 7000 FORMAT (' uniflow      ',G14.7,2X,G14.7)
 8000 FORMAT (' quit')
      END
C
C-------------------------------------------------------------------------------------------------
C
       SUBROUTINE GVOUT (ICALL)
C
C-------------------------------------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER(4) ICALL
      CHARACTER(1) ADUM
      INCLUDE  'GVCOM.INC'
      INCLUDE  'LUSYS.INC'
      SAVE
      WRITE  (ILUOUT,1000)  RK,RH,RHEAD0,RBASE
      WRITE  (ILUOUT,1002)  RPOTC,CREFZ
      WRITE  (ILUOUT,1004)  RCON,RPI,RQ0
      WRITE  (ILUOUT,1005)  CUNALP
      WRITE  (ILUOUT,1006)  RAINX,RAINY,RAIN
      WRITE  (ILUOUT,1008)  LRAINX,LRAINY,LRAINR
      WRITE  (ILUOUT,1010)  RPOR
      WRITE  (ILUOUT,1020)  LINTERFACE
      WRITE  (ILUOUT,1030)  RSPECIFICGRAVITYFRESH,RSPECIFICGRAVITYSALT,
     &                      RSEALEVEL
      WRITE  (ILUOUT,2000)  ICALL
      IF (LUOUTFILE) WRITE (ILUME,2000) ICALL
      READ (ILUIN,3000)  ADUM
      IF (ADUM.EQ.'Q'.OR.ADUM.EQ.'q') RETURN
      RETURN
 1000 FORMAT  (' GVOUT: RK,RH,RHEAD0,RBASE ', 4(G14.7))
 1002 FORMAT  (' GVOUT: RPOTC,CREFZ  ',3(G14.7))
 1004 FORMAT  (' GVOUT: RCON,RPI,RQ0  ',3(G14.7))
 1005 FORMAT  (' GVOUT: CUNALP  ',G14.7,G14.7)
 1006 FORMAT  (' GVOUT: RAINX,RAINY,RAIN  ',3(G14.7))
 1008 FORMAT  (' GVOUT: LRAINX,LRAINY,LRAINR  ',3L3)
 1010 FORMAT  (' GVOUT: RPOR  ',G14.7)
 1020 FORMAT  (' GVOUT: LINTERFACE  ',L3)
 1030 FORMAT  (' GVOUT: RSPECIFICGRAVITYFRESH,RSPECIFICGRAVITYSALT, ',
     &         'RSEALEVEL ',/,3(G14.7))
 2000 FORMAT  (' GVOUT: ICALL=  ',I3,' press <Enter> to continue.')
 3000 FORMAT  (A1)
      END
